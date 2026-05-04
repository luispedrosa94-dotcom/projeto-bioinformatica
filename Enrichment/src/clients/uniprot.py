"""
UniProt REST API client.

Uses ThreadPoolExecutor to fetch multiple proteins in parallel (default: 10
workers), making it ~10x faster than sequential queries while staying within
UniProt's rate limits.

Includes checkpoint support: progress is saved every 100 proteins so the
pipeline can resume after interruption without redoing completed work.

GO aspect resolution is performed locally using the GO term → aspect map
extracted directly from UniProt records (F:/P:/C: prefixes), without any
additional API call.
"""
from __future__ import annotations

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


def _fetch_one(acc: str) -> list[EnrichmentRecord]:
    """Fetch and parse a single UniProt entry. Returns [] on error."""
    try:
        resp = get_with_retry(UNIPROT_ENTRY.format(acc=acc), headers=_HEADERS)
        return _parse_entry(resp.json())
    except Exception as e:
        log.debug("UniProt: %s failed: %s", acc, e)
        return []


def fetch_proteins(
    accessions: list[str],
    batch_size: int = 1,      # unused — kept for API compatibility
    delay: float = 0.0,       # unused — threading handles concurrency
    max_workers: int = 10,
    checkpoint_path: Path | None = None,
) -> list[EnrichmentRecord]:
    """
    Fetch protein annotations from UniProt using parallel threads.

    - max_workers=10 → ~10 simultaneous requests → ~10x faster than sequential
    - Saves a checkpoint every 100 proteins
    - On restart, skips already-fetched proteins automatically
    """
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
        log.info("UniProt: all %d accessions already in checkpoint", total)
        return records

    completed = len(done_accs)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_fetch_one, acc): acc for acc in remaining}
        for future in as_completed(futures):
            acc = futures[future]
            new_records = future.result()
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


def _parse_entry(entry: dict) -> list[EnrichmentRecord]:
    """Parse a single UniProt JSON entry into EnrichmentRecords."""
    records: list[EnrichmentRecord] = []
    acc = entry.get("primaryAccession", "")
    if not acc:
        return records

    # Reviewed status
    entry_type = entry.get("entryType", "").lower()
    is_reviewed = "reviewed" in entry_type and "unreviewed" not in entry_type
    records.append(EnrichmentRecord(
        uniprot_accession=acc,
        source=EnrichmentSource.UNIPROT,
        enrichment_type=EnrichmentType.REVIEWED_STATUS,
        value="reviewed" if is_reviewed else "unreviewed",
    ))

    # Protein name (reviewed → recommendedName, unreviewed → submittedName)
    try:
        pname = entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        try:
            pname = entry["proteinDescription"]["submittedNames"][0]["fullName"]["value"]
        except (KeyError, TypeError, IndexError):
            pname = None
    if pname:
        records.append(EnrichmentRecord(
            uniprot_accession=acc,
            source=EnrichmentSource.UNIPROT,
            enrichment_type=EnrichmentType.PROTEIN_NAME,
            value=pname,
        ))

    # Gene name
    try:
        genes = entry.get("genes", [])
        if genes:
            gene_name = genes[0]["geneName"]["value"]
            records.append(EnrichmentRecord(
                uniprot_accession=acc,
                source=EnrichmentSource.UNIPROT,
                enrichment_type=EnrichmentType.GENE_NAME,
                value=gene_name,
            ))
    except (KeyError, TypeError, IndexError):
        pass

    # GO terms — aspect encoded as prefix F: / P: / C:
    aspect_map = {
        "F:": EnrichmentType.GO_MF,
        "P:": EnrichmentType.GO_BP,
        "C:": EnrichmentType.GO_CC,
    }
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        go_id    = xref.get("id", "")
        props    = {p["key"]: p["value"] for p in xref.get("properties", [])}
        go_name  = props.get("GoTerm", "")
        evidence = props.get("GoEvidenceType", "")
        prefix   = go_name[:2] if len(go_name) >= 2 else ""
        etype    = aspect_map.get(prefix)
        if etype and go_id:
            records.append(EnrichmentRecord(
                uniprot_accession=acc,
                source=EnrichmentSource.UNIPROT,
                enrichment_type=etype,
                value=go_id,
                label=go_name[2:].strip(),
                extras={"evidence": evidence},
            ))

    # Subcellular location
    try:
        for comment in entry.get("comments", []):
            if comment.get("commentType") != "SUBCELLULAR LOCATION":
                continue
            for loc in comment.get("subcellularLocations", []):
                loc_val = loc.get("location", {}).get("value", "")
                if loc_val:
                    records.append(EnrichmentRecord(
                        uniprot_accession=acc,
                        source=EnrichmentSource.UNIPROT,
                        enrichment_type=EnrichmentType.SUBCELLULAR_LOCATION,
                        value=loc_val,
                    ))
    except (KeyError, TypeError):
        pass

    # Function description (first only)
    try:
        for comment in entry.get("comments", []):
            if comment.get("commentType") != "FUNCTION":
                continue
            for text in comment.get("texts", []):
                val = text.get("value", "")
                if val:
                    records.append(EnrichmentRecord(
                        uniprot_accession=acc,
                        source=EnrichmentSource.UNIPROT,
                        enrichment_type=EnrichmentType.FUNCTION_DESCRIPTION,
                        value=val,
                    ))
                    break
    except (KeyError, TypeError):
        pass

    # Keywords
    for kw in entry.get("keywords", []):
        kw_name = kw.get("name", "")
        if kw_name:
            records.append(EnrichmentRecord(
                uniprot_accession=acc,
                source=EnrichmentSource.UNIPROT,
                enrichment_type=EnrichmentType.KEYWORD,
                value=kw_name,
            ))

    return records
