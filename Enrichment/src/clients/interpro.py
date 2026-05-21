"""
InterPro REST API client.

Fetches the complete InterPro entry coverage for each protein and saves the
raw JSON response to outputs/interpro_raw/{acc}.json — nothing is filtered
or lost.

The raw JSON is the single source of truth. The consolidate.py script reads
directly from these raw files and extracts whatever fields are needed.

Uses ThreadPoolExecutor for parallel fetching (default 10 workers) and supports
checkpointing every 100 proteins for resumable runs.

API documentation: https://interpro-documentation.readthedocs.io/en/latest/
Endpoint used: https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/UniProt/{acc}/
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

# Returns ALL InterPro entries covering a given UniProt protein.
INTERPRO_ENTRY = (
    "https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/UniProt/{acc}/"
    "?page_size=200"
)

# Unintegrated member-database signatures.
INTERPRO_UNINTEGRATED = (
    "https://www.ebi.ac.uk/interpro/api/entry/unintegrated/protein/UniProt/{acc}/"
    "?page_size=200"
)

_HEADERS = {
    "User-Agent": "bioinformatics_pipeline/1.0 (academic project)",
    "Accept": "application/json",
}


def _fetch_url(url: str, max_retries: int = 5, backoff_base: float = 2.0,
               timeout: int = 60) -> dict:
    """
    Fetch a single InterPro URL with proper handling of:
      - 200 OK → return parsed JSON
      - 204 No Content → return empty results (no match for this protein)
      - 404 Not Found → return empty results
      - 429 Rate Limited → retry with exponential backoff
      - 5xx Server Error → retry with exponential backoff
      - other 4xx → raise exception
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
            # No InterPro matches for this protein — legitimate.
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
        # 4xx other than 429 / 404 — no point retrying
        resp.raise_for_status()

    raise RuntimeError(f"Failed to GET {url} after {max_retries} attempts")


def _fetch_one(acc: str, raw_dir: Path) -> dict | None:
    """
    Fetch a single InterPro coverage entry for a UniProt accession and save
    raw JSON to raw_dir/{acc}.json. Returns the parsed dict or None on error.

    The saved JSON wraps both the integrated entries and the unintegrated
    signatures so consolidate.py has everything in one place.
    """
    raw_file = raw_dir / f"{acc}.json"

    # If already cached locally, read from disk
    if raw_file.exists():
        try:
            with open(raw_file, encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            pass  # corrupt file — re-fetch

    integrated:   dict | None = None
    unintegrated: dict | None = None

    try:
        integrated = _fetch_url(INTERPRO_ENTRY.format(acc=acc))
    except Exception as e:
        log.debug("InterPro integrated: %s failed: %s", acc, e)

    try:
        unintegrated = _fetch_url(INTERPRO_UNINTEGRATED.format(acc=acc))
    except Exception as e:
        log.debug("InterPro unintegrated: %s failed: %s", acc, e)

    if integrated is None and unintegrated is None:
        return None

    data = {
        "accession": acc,
        "integrated":   integrated   if integrated   is not None else {"results": [], "count": 0},
        "unintegrated": unintegrated if unintegrated is not None else {"results": [], "count": 0},
    }

    raw_file.parent.mkdir(parents=True, exist_ok=True)
    with open(raw_file, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False)

    return data


def fetch_proteins(
    accessions: list[str],
    batch_size: int = 1,
    delay: float = 0.0,
    max_workers: int = 10,
    checkpoint_path: Path | None = None,
    raw_dir: Path | None = None,
) -> list[EnrichmentRecord]:
    """
    Fetch complete InterPro coverage for each protein, saving raw JSON.
    Returns a minimal list of EnrichmentRecords for enrich.py summary stats.
    """
    if raw_dir is None:
        raw_dir = Path("outputs/interpro_raw")
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
            data = future.result()
            if data:
                new_records = _parse_entry_minimal(data)
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


def _parse_entry_minimal(data: dict) -> list[EnrichmentRecord]:
    """Minimal parser — extracts just enough for enrich.py summary stats."""
    records: list[EnrichmentRecord] = []
    acc = data.get("accession", "")
    if not acc:
        return records

    def _emit(etype: EnrichmentType, value: str, label: str | None = None,
              extras: dict | None = None) -> None:
        if not value:
            return
        records.append(EnrichmentRecord(
            uniprot_accession=acc,
            source=EnrichmentSource.INTERPRO,
            enrichment_type=etype,
            value=value,
            label=label,
            extras=extras or {},
        ))

    # ── Integrated InterPro entries ─────────────────────────────────────
    integrated = data.get("integrated", {}) or {}
    for result in integrated.get("results", []):
        meta = result.get("metadata", {}) or {}
        ipr_id = meta.get("accession", "")
        if not ipr_id:
            continue
        _emit(
            EnrichmentType.INTERPRO_ENTRY,
            ipr_id,
            label=meta.get("name") or None,
            extras={
                "type":              meta.get("type")             or None,
                "source_database":   meta.get("source_database")  or None,
                "integrated_to":     meta.get("integrated")       or None,
                "member_databases":  list((meta.get("member_databases") or {}).keys()),
                "go_terms_count":    len(meta.get("go_terms") or []),
            },
        )

    # ── Unintegrated member-database signatures ─────────────────────────
    unintegrated = data.get("unintegrated", {}) or {}
    for result in unintegrated.get("results", []):
        meta = result.get("metadata", {}) or {}
        sig_id = meta.get("accession", "")
        if not sig_id:
            continue
        _emit(
            EnrichmentType.INTERPRO_UNINTEGRATED_SIGNATURE,
            sig_id,
            label=meta.get("name") or None,
            extras={
                "type":            meta.get("type")            or None,
                "source_database": meta.get("source_database") or None,
            },
        )

    return records