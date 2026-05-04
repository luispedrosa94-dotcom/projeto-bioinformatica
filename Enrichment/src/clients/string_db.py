"""
STRING API client.

STRING provides protein-protein association data.

STRING API docs: https://string-db.org/help/api/
"""
from __future__ import annotations

import logging
import time

from ..schema import EnrichmentRecord, EnrichmentSource, EnrichmentType
from ..utils.rate_limit import post_with_retry

log = logging.getLogger(__name__)

STRING_API = "https://string-db.org/api"

# Minimum combined score to keep an interaction (0-1000)
# 400 = medium confidence, 700 = high confidence
DEFAULT_MIN_SCORE = 400

_HEADERS = {
    "User-Agent": "bioinformatics_pipeline/1.0 (academic project)",
    "Content-Type": "application/x-www-form-urlencoded",
}


def map_identifiers(
    accessions: list[str],
    species: int = 0,
    batch_size: int = 500,
    delay: float = 1.0,
) -> dict[str, str]:
    """
    Map UniProt accessions to STRING protein identifiers.

    Returns {uniprot_accession: string_id}.
    species=0 lets STRING auto-detect species from the identifier.
    """
    mapping: dict[str, str] = {}
    total = len(accessions)

    for i in range(0, total, batch_size):
        batch = accessions[i: i + batch_size]
        log.info("STRING: mapping identifiers batch %d-%d of %d",
                 i + 1, min(i + batch_size, total), total)

        # STRING expects identifiers newline-separated
        data = {
            "identifiers": "\n".join(batch),
            "limit": 1,
            "echo_query": 1,
            "caller_identity": "bioinformatics_project",
        }
        if species:
            data["species"] = species

        try:
            resp = post_with_retry(
                f"{STRING_API}/json/get_string_ids",
                data=data,
                headers=_HEADERS,
            )
            items = resp.json()
        except Exception as e:
            log.error("STRING mapping batch %d failed: %s", i, e)
            continue

        for item in items:
            query_acc = item.get("queryItem", "")
            string_id = item.get("stringId", "")
            if query_acc and string_id:
                mapping[query_acc] = string_id

        if i + batch_size < total:
            time.sleep(delay)

    log.info("STRING: mapped %d / %d accessions to STRING IDs", len(mapping), total)
    return mapping


def fetch_interactions(
    string_ids: list[str],
    uniprot_map: dict[str, str],
    min_score: int = DEFAULT_MIN_SCORE,
    batch_size: int = 500,
    delay: float = 1.0,
) -> list[EnrichmentRecord]:
    """
    Fetch protein-protein interactions from STRING.
    uniprot_map: {uniprot_accession: string_id}
    """
    records: list[EnrichmentRecord] = []
    total = len(string_ids)
    # Reverse: string_id → uniprot_accession
    rev_map = {v: k for k, v in uniprot_map.items()}

    for i in range(0, total, batch_size):
        batch = string_ids[i: i + batch_size]
        log.info("STRING: fetching interactions batch %d-%d of %d",
                 i + 1, min(i + batch_size, total), total)

        data = {
            "identifiers": "\n".join(batch),
            "required_score": min_score,
            "caller_identity": "bioinformatics_project",
        }

        try:
            resp = post_with_retry(
                f"{STRING_API}/json/network",
                data=data,
                headers=_HEADERS,
            )
            edges = resp.json()
        except Exception as e:
            log.error("STRING interactions batch %d failed: %s", i, e)
            continue

        for edge in edges:
            sid_a = edge.get("stringId_A", "")
            sid_b = edge.get("stringId_B", "")
            score = edge.get("score", 0)
            pname_b = edge.get("preferredName_B", sid_b)

            acc_a = rev_map.get(sid_a, sid_a)
            if acc_a and sid_b:
                records.append(EnrichmentRecord(
                    uniprot_accession=acc_a,
                    source=EnrichmentSource.STRING,
                    enrichment_type=EnrichmentType.INTERACTION_PARTNER,
                    value=sid_b,
                    label=pname_b,
                    score=round(score / 1000.0, 4),
                    extras={
                        "string_id_a": sid_a,
                        "string_id_b": sid_b,
                        "combined_score": score,
                    },
                ))

        if i + batch_size < total:
            time.sleep(delay)

    log.info("STRING: fetched %d interaction records", len(records))
    return records


def fetch_enrichment(
    string_ids: list[str],
    species: int = 0,
    delay: float = 1.0,
) -> list[EnrichmentRecord]:
    """
    Fetch functional enrichment from STRING for the full protein set.
    Results are set-level (tagged with uniprot_accession='SET_LEVEL').
    """
    records: list[EnrichmentRecord] = []
    log.info("STRING: fetching functional enrichment for %d proteins", len(string_ids))

    data = {
        "identifiers": "\n".join(string_ids),
        "caller_identity": "bioinformatics_project",
    }
    if species:
        data["species"] = species

    try:
        resp = post_with_retry(
            f"{STRING_API}/json/enrichment",
            data=data,
            headers=_HEADERS,
        )
        items = resp.json()
    except Exception as e:
        log.error("STRING enrichment failed: %s", e)
        return records

    for item in items:
        term_id = item.get("term", "")
        if not term_id:
            continue
        records.append(EnrichmentRecord(
            uniprot_accession="SET_LEVEL",
            source=EnrichmentSource.STRING,
            enrichment_type=EnrichmentType.FUNCTIONAL_ENRICHMENT,
            value=term_id,
            label=item.get("description", ""),
            score=item.get("fdr"),
            extras={
                "category": item.get("category", ""),
                "p_value": item.get("p_value"),
                "fdr": item.get("fdr"),
                "number_of_genes": item.get("number_of_genes"),
            },
        ))

    log.info("STRING: fetched %d enrichment terms", len(records))
    return records
