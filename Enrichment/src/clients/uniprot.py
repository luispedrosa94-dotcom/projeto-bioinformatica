"""
UniProt REST API client.

Fetches the complete UniProt entry for each protein and saves the raw JSON
response to outputs/caches/uniprot_raw/{acc}.json — nothing is filtered or lost.

The raw JSON is the single source of truth. The consolidate.py script reads
directly from these raw files and extracts whatever fields are needed.

Uses ThreadPoolExecutor for parallel fetching (default 10 workers) and supports
checkpointing every 100 proteins for resumable runs.
"""
from __future__ import annotations

import json
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from ..schema import EnrichmentRecord, EnrichmentSource, EnrichmentType
from ..utils.rate_limit import get_with_retry
from ..utils.checkpoint import save as cp_save, load as cp_load

log = logging.getLogger(__name__)

UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}"

_HEADERS = {
    "User-Agent": "bioinformatics_pipeline/1.0 (academic project)",
    "Accept": "application/json",
}


def _fetch_one(acc: str, raw_dir: Path) -> dict | None:
    """
    Fetch a single UniProt entry and save raw JSON to raw_dir/{acc}.json.
    Returns the parsed dict or None on error.
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
        resp = get_with_retry(UNIPROT_ENTRY.format(acc=acc), headers=_HEADERS)
        data = resp.json()
        raw_file.parent.mkdir(parents=True, exist_ok=True)
        with open(raw_file, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)
        return data
    except Exception as e:
        log.debug("UniProt: %s failed: %s", acc, e)
        return None


def fetch_proteins(
    accessions: list[str],
    batch_size: int = 1,
    delay: float = 0.0,
    max_workers: int = 10,
    checkpoint_path: Path | None = None,
    raw_dir: Path | None = None,
) -> list[EnrichmentRecord]:
    """
    Fetch complete UniProt entries, saving raw JSON for each protein.
    Also returns EnrichmentRecords for backward compatibility with enrich.py.
    """
    if raw_dir is None:
        raw_dir = Path("outputs/caches/uniprot_raw")
    raw_dir.mkdir(parents=True, exist_ok=True)

    done_accs: set[str] = set()
    records: list[EnrichmentRecord] = []

    if checkpoint_path:
        cp = cp_load(checkpoint_path)
        if cp:
            done_accs = set(cp.get("done", []))
            records = [EnrichmentRecord(**r) for r in cp.get("records", [])]
            log.info("UniProt: resuming — %d / %d already done",
                     len(done_accs), len(accessions))

    remaining = [a for a in accessions if a not in done_accs]
    total = len(accessions)

    if not remaining:
        log.info("UniProt: all %d accessions already cached", total)
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
                log.info("UniProt: %d / %d proteins fetched", completed, total)
                if checkpoint_path:
                    cp_save(checkpoint_path, {
                        "done": list(done_accs),
                        "records": [r.model_dump() for r in records],
                    })

    log.info("UniProt: done — %d records from %d accessions", len(records), total)
    return records


def _parse_entry_minimal(entry: dict) -> list[EnrichmentRecord]:
    """
    Minimal parser — extracts just enough for enrich.py summary stats
    and GO aspect resolution. Full extraction is done by consolidate.py
    reading the raw JSON directly.
    """
    records: list[EnrichmentRecord] = []
    acc = entry.get("primaryAccession", "")
    if not acc:
        return records

    def _emit(etype, value, label=None, extras=None):
        if not value:
            return
        records.append(EnrichmentRecord(
            uniprot_accession=acc,
            source=EnrichmentSource.UNIPROT,
            enrichment_type=etype,
            value=value,
            label=label,
            extras=extras or {},
        ))

    # Reviewed status
    entry_type = entry.get("entryType", "").lower()
    is_reviewed = "reviewed" in entry_type and "unreviewed" not in entry_type
    _emit(EnrichmentType.REVIEWED_STATUS, "reviewed" if is_reviewed else "unreviewed")

    # Annotation score
    ann_score = entry.get("annotationScore")
    if ann_score is not None:
        _emit(EnrichmentType.ANNOTATION_SCORE, str(ann_score))

    # Protein existence
    _emit(EnrichmentType.PROTEIN_EXISTENCE, entry.get("proteinExistence", ""))

    # Organism
    organism = entry.get("organism", {})
    org_name = organism.get("scientificName", "")
    if org_name:
        _emit(EnrichmentType.ORGANISM, org_name, extras={
            "common_name": organism.get("commonName", "") or None,
            "taxon_id": str(organism.get("taxonId", "")) or None,
            "lineage": organism.get("lineage", []),
        })

    # Protein name
    try:
        pname = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        try:
            pname = entry["proteinDescription"]["submittedNames"][0]["fullName"]["value"]
        except (KeyError, TypeError, IndexError):
            pname = None
    if pname:
        _emit(EnrichmentType.PROTEIN_NAME, pname)

    # Gene name
    try:
        genes = entry.get("genes", [])
        if genes:
            gene_name = genes[0]["geneName"]["value"]
            _emit(EnrichmentType.GENE_NAME, gene_name)
    except (KeyError, TypeError, IndexError):
        pass

    # GO terms for aspect resolution
    aspect_map = {
        "F:": EnrichmentType.GO_MF,
        "P:": EnrichmentType.GO_BP,
        "C:": EnrichmentType.GO_CC,
    }
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        go_id   = xref.get("id", "")
        props   = {p["key"]: p["value"] for p in xref.get("properties", [])}
        go_name = props.get("GoTerm", "")
        evidence = props.get("GoEvidenceType", "")
        prefix  = go_name[:2] if len(go_name) >= 2 else ""
        etype   = aspect_map.get(prefix)
        if etype and go_id:
            ev_parts = evidence.split(":", 1)
            _emit(etype, go_id, label=go_name[2:].strip(), extras={
                "evidence_code": ev_parts[0] if ev_parts else None,
                "evidence_source": ev_parts[1] if len(ev_parts) > 1 else None,
            })

    # Keywords
    for kw in entry.get("keywords", []):
        kw_name = kw.get("name", "")
        if kw_name:
            _emit(EnrichmentType.KEYWORD, kw_name)

    # Subcellular location
    for comment in entry.get("comments", []):
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for loc in comment.get("subcellularLocations", []):
            loc_val = loc.get("location", {}).get("value", "")
            if loc_val:
                _emit(EnrichmentType.SUBCELLULAR_LOCATION, loc_val)

    # Function description
    for comment in entry.get("comments", []):
        if comment.get("commentType") != "FUNCTION":
            continue
        for text in comment.get("texts", []):
            val = text.get("value", "")
            if val:
                _emit(EnrichmentType.FUNCTION_DESCRIPTION, val)

    return records
