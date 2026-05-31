"""
InterPro REST API client — uses the /entry/all/ endpoint to fetch ALL
signatures (InterPro entries + Pfam + PANTHER + CATH-Gene3D + SSF +
PIRSF + NCBIfam + SMART + PROSITE + CDD) in a single call per protein.

This gives much richer coverage than the previous /entry/InterPro/ +
/entry/unintegrated/ split:
  - Each member-database signature comes with its own e-value, score,
    model, and protein locations
  - InterPro integrated entries (IPRxxxxx) come with member_databases
    dict showing which signatures they aggregate
  - Hierarchy, GO terms, cross_references are preserved

The raw JSON is saved as outputs/caches/interpro_raw/{acc}.json — single source
of truth. The consolidate.py parser reads directly from these files.

Endpoint: https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}/
"""
from __future__ import annotations

import json
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

from ..schema import EnrichmentRecord, EnrichmentSource, EnrichmentType
from ..utils.checkpoint import save as cp_save, load as cp_load

log = logging.getLogger(__name__)

# Single endpoint that returns ALL signatures (InterPro + member databases)
# covering a given UniProt protein. page_size=200 covers every realistic case.
INTERPRO_ALL = (
    "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{acc}/"
    "?page_size=200"
)

_HEADERS = {
    "User-Agent": "bioinformatics_pipeline/1.0 (academic project)",
    "Accept": "application/json",
}


def _fetch_url(url: str, max_retries: int = 5, backoff_base: float = 2.0,
               timeout: int = 60) -> dict:
    """
    Fetch a single InterPro URL with proper status code handling:
      - 200 OK              → return parsed JSON
      - 204 No Content      → return empty results (no match for this protein)
      - 404 Not Found       → return empty results
      - 429 Rate Limited    → retry with exponential backoff
      - 5xx Server Error    → retry with exponential backoff
      - other 4xx           → raise exception
    """
    for attempt in range(max_retries):
        try:
            resp = requests.get(url, headers=_HEADERS, timeout=timeout)
        except requests.exceptions.Timeout:
            wait = backoff_base ** attempt
            log.warning("Timeout on %s — waiting %.1fs", url, wait)
            time.sleep(wait)
            continue
        except requests.exceptions.ConnectionError as e:
            wait = backoff_base ** attempt
            log.warning("Connection error on %s: %s — waiting %.1fs", url, e, wait)
            time.sleep(wait)
            continue

        if resp.status_code == 200:
            return resp.json()
        if resp.status_code in (204, 404):
            return {"results": [], "count": 0}
        if resp.status_code == 429:
            wait = backoff_base ** attempt
            log.warning("Rate limited (429) on %s — waiting %.1fs", url, wait)
            time.sleep(wait)
            continue
        if resp.status_code >= 500:
            wait = backoff_base ** attempt
            log.warning("Server error %d on %s — waiting %.1fs", resp.status_code, url, wait)
            time.sleep(wait)
            continue
        resp.raise_for_status()

    raise RuntimeError(f"Failed to GET {url} after {max_retries} attempts")


def _fetch_one(acc: str, raw_dir: Path) -> dict | None:
    """
    Fetch the complete InterPro coverage for one UniProt accession.
    Saves raw JSON to raw_dir/{acc}.json.

    The saved JSON wraps the API response with the accession so consolidate.py
    has the protein ID directly available:
        {
            "accession": "Q46505",
            "data": <raw API response from /entry/all/protein/uniprot/Q46505/>
        }
    """
    raw_file = raw_dir / f"{acc}.json"

    # If already cached locally, read from disk
    if raw_file.exists():
        try:
            with open(raw_file, encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass  # corrupt file — re-fetch

    try:
        data = _fetch_url(INTERPRO_ALL.format(acc=acc))
    except Exception as e:
        log.debug("InterPro /entry/all/ for %s failed: %s", acc, e)
        return None

    wrapped = {
        "accession": acc,
        "data": data,
    }

    raw_file.parent.mkdir(parents=True, exist_ok=True)
    with open(raw_file, "w", encoding="utf-8") as f:
        json.dump(wrapped, f, ensure_ascii=False)

    return wrapped


def fetch_proteins(
    accessions: list[str],
    batch_size: int = 1,
    delay: float = 0.0,
    max_workers: int = 10,
    checkpoint_path: Path | None = None,
    raw_dir: Path | None = None,
) -> list[EnrichmentRecord]:
    """
    Fetch complete InterPro /entry/all/ coverage for each protein.
    Returns a minimal list of EnrichmentRecords for enrich.py summary stats.
    Full extraction is done by interpro_extract.py reading raw JSON directly.
    """
    if raw_dir is None:
        raw_dir = Path("outputs/caches/interpro_raw")
    raw_dir.mkdir(parents=True, exist_ok=True)

    done_accs: set[str] = set()
    records: list[EnrichmentRecord] = []

    if checkpoint_path:
        cp = cp_load(checkpoint_path)
        if cp:
            done_accs = set(cp.get("done", []))
            records = [EnrichmentRecord(**r) for r in cp.get("records", [])]
            log.info("InterPro: resuming — %d / %d already done",
                     len(done_accs), len(accessions))

    remaining = [a for a in accessions if a not in done_accs]
    total = len(accessions)

    if not remaining:
        log.info("InterPro: all %d accessions already cached", total)
        return records

    completed = len(done_accs)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_fetch_one, acc, raw_dir): acc
            for acc in remaining
        }
        for future in as_completed(futures):
            acc = futures[future]
            wrapped = future.result()
            if wrapped:
                new_records = _parse_entry_minimal(wrapped)
                records.extend(new_records)
            done_accs.add(acc)
            completed += 1

            if completed % 100 == 0 or completed == total:
                log.info("InterPro: %d / %d proteins fetched", completed, total)
                if checkpoint_path:
                    cp_save(checkpoint_path, {
                        "done": list(done_accs),
                        "records": [r.model_dump() for r in records],
                    })

    log.info("InterPro: done — %d records from %d accessions", len(records), total)
    return records


def _parse_entry_minimal(wrapped: dict) -> list[EnrichmentRecord]:
    """
    Minimal parser — extracts just enough for enrich.py summary stats.
    Full extraction is in interpro_extract.py used by consolidate.py.

    Produces one record per signature found, distinguishing:
      - INTERPRO_ENTRY for entries with source_database == "interpro" (IPR...)
      - INTERPRO_UNINTEGRATED_SIGNATURE for member-db signatures (Pfam, etc.)
    """
    records: list[EnrichmentRecord] = []
    acc = wrapped.get("accession", "")
    if not acc:
        return records

    data = wrapped.get("data") or {}
    results = data.get("results") or []

    for result in results:
        meta = result.get("metadata") or {}
        sig_acc = meta.get("accession", "")
        if not sig_acc:
            continue

        source_db = (meta.get("source_database") or "").lower()
        is_interpro = (source_db == "interpro")

        etype = (EnrichmentType.INTERPRO_ENTRY if is_interpro
                 else EnrichmentType.INTERPRO_UNINTEGRATED_SIGNATURE)

        records.append(EnrichmentRecord(
            uniprot_accession=acc,
            source=EnrichmentSource.INTERPRO,
            enrichment_type=etype,
            value=sig_acc,
            label=meta.get("name") or None,
            extras={
                "type":              meta.get("type")           or None,
                "source_database":   meta.get("source_database") or None,
                "integrated_to":     meta.get("integrated")     or None,
            },
        ))

    return records