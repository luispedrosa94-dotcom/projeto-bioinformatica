"""
Helper functions for extracting InterPro information from the raw JSON
files saved by the interpro client.

Used by consolidate.py to build the 'interpro' section of each
protein_profile entry. Preserves all available information from the
InterPro7 API response.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
import json


def _safe_get(d: dict | None, *path, default=None):
    """Safe nested dictionary access — never crashes if a key is missing."""
    cur = d if d is not None else {}
    for k in path:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(k)
        if cur is None:
            return default
    return cur


def _extract_entry(result: dict) -> dict:
    """
    Convert one raw InterPro result block into a canonical dict that
    preserves all fields the API exposes.

    Each `result` from /entry/InterPro/protein/UniProt/{acc}/ has:
      - metadata        : info about the InterPro entry itself
      - proteins        : a list with the location info for this protein
    """
    meta = result.get("metadata", {}) or {}
    proteins = result.get("proteins", []) or []

    # Location info (positions in the protein) lives inside proteins[*]
    protein_blocks = []
    for prot in proteins:
        protein_blocks.append({
            "accession":          prot.get("accession") or None,
            "protein_length":     prot.get("protein_length"),
            "source_database":    prot.get("source_database") or None,
            "organism":           prot.get("organism") or None,
            "entry_protein_locations": prot.get("entry_protein_locations") or [],
        })

    # GO terms attached to this InterPro entry, by category
    go_terms = meta.get("go_terms") or []
    go_by_category = {
        "molecular_function": [],
        "biological_process": [],
        "cellular_component": [],
    }
    category_map = {
        "molecular_function": "molecular_function",
        "biological_process": "biological_process",
        "cellular_component": "cellular_component",
        # InterPro API also uses short codes
        "F": "molecular_function",
        "P": "biological_process",
        "C": "cellular_component",
    }
    for g in go_terms:
        cat_raw = g.get("category", {}) or {}
        cat_name = ""
        if isinstance(cat_raw, dict):
            cat_name = cat_raw.get("name") or cat_raw.get("code") or ""
        else:
            cat_name = str(cat_raw)
        cat = category_map.get(cat_name, "molecular_function")
        go_by_category[cat].append({
            "go_id":       g.get("identifier") or g.get("id") or None,
            "name":        g.get("name") or None,
            "category":    cat,
            "category_raw": cat_name or None,
        })

    # Member databases (Pfam, SMART, PROSITE, etc. signatures inside this entry)
    member_dbs = meta.get("member_databases") or {}
    member_dbs_clean = {}
    for db_name, sigs in member_dbs.items():
        if isinstance(sigs, dict):
            # Format: {"PF00001": "Name of family", ...}
            member_dbs_clean[db_name] = [
                {"accession": acc, "name": name}
                for acc, name in sigs.items()
            ]
        else:
            member_dbs_clean[db_name] = sigs

    # Hierarchy (parent / children InterPro entries)
    hierarchy = meta.get("hierarchy") or {}

    # Cross-references (PDB, Reactome, MetaCyc, etc.)
    cross_refs = meta.get("cross_references") or {}

    # Literature references
    literature = meta.get("literature") or {}

    return {
        "accession":         meta.get("accession") or None,
        "name":              meta.get("name") or None,
        "short_name":        meta.get("name_short") or meta.get("short_name") or None,
        "description":       meta.get("description") or [],
        "type":              meta.get("type") or None,
        "source_database":   meta.get("source_database") or None,
        "integrated":        meta.get("integrated") or None,
        "member_databases":  member_dbs_clean,
        "go_terms":          go_by_category,
        "hierarchy":         hierarchy,
        "cross_references":  cross_refs,
        "literature":        literature,
        "entry_id":          meta.get("entry_id") or None,
        "entry_date":        meta.get("entry_date") or None,
        "proteins":          protein_blocks,
    }


def _extract_unintegrated_signature(result: dict) -> dict:
    """
    Convert one raw unintegrated-signature result block into a canonical dict.

    These are member-database signatures (Pfam, SMART, etc.) that are NOT
    yet integrated into an InterPro entry. They still convey useful domain
    information for the protein.
    """
    meta = result.get("metadata", {}) or {}
    proteins = result.get("proteins", []) or []

    protein_blocks = []
    for prot in proteins:
        protein_blocks.append({
            "accession":          prot.get("accession") or None,
            "protein_length":     prot.get("protein_length"),
            "source_database":    prot.get("source_database") or None,
            "organism":           prot.get("organism") or None,
            "entry_protein_locations": prot.get("entry_protein_locations") or [],
        })

    return {
        "accession":         meta.get("accession") or None,
        "name":              meta.get("name") or None,
        "short_name":        meta.get("name_short") or meta.get("short_name") or None,
        "type":              meta.get("type") or None,
        "source_database":   meta.get("source_database") or None,
        "go_terms":          meta.get("go_terms") or [],
        "cross_references":  meta.get("cross_references") or {},
        "literature":        meta.get("literature") or {},
        "entry_id":          meta.get("entry_id") or None,
        "entry_date":        meta.get("entry_date") or None,
        "proteins":          protein_blocks,
    }


def extract_interpro_section(raw_data: dict) -> dict:
    """
    Build the 'interpro' section of one protein profile from the raw JSON
    saved by the InterPro client.

    Output schema:
    {
      "entries": [ { ...integrated InterPro entries... } ],
      "unintegrated_signatures": [ { ...unintegrated signatures... } ],
      "summary": {
        "total_integrated_entries":         int,
        "total_unintegrated_signatures":    int,
        "by_type":                          { "Domain": n, "Family": n, ... },
        "member_databases_used":            [ "pfam", "smart", ... ],
        "go_terms_count":                   { "BP": n, "MF": n, "CC": n },
      }
    }
    """
    integrated_raw   = raw_data.get("integrated",   {}) or {}
    unintegrated_raw = raw_data.get("unintegrated", {}) or {}

    entries = [
        _extract_entry(r) for r in integrated_raw.get("results", []) or []
    ]
    unintegrated = [
        _extract_unintegrated_signature(r) for r in unintegrated_raw.get("results", []) or []
    ]

    # ── Summary statistics ──────────────────────────────────────────────
    type_counter = Counter()
    for e in entries:
        if e.get("type"):
            type_counter[e["type"]] += 1

    member_dbs_used: set[str] = set()
    for e in entries:
        for db in (e.get("member_databases") or {}).keys():
            member_dbs_used.add(db)

    go_counts = {"molecular_function": 0, "biological_process": 0, "cellular_component": 0}
    for e in entries:
        for cat, terms in (e.get("go_terms") or {}).items():
            go_counts[cat] = go_counts.get(cat, 0) + len(terms)

    summary = {
        "total_integrated_entries":      len(entries),
        "total_unintegrated_signatures": len(unintegrated),
        "by_type":                       dict(type_counter),
        "member_databases_used":         sorted(member_dbs_used),
        "go_terms_count":                go_counts,
    }

    return {
        "entries": entries,
        "unintegrated_signatures": unintegrated,
        "summary": summary,
    }


def load_interpro_raw(raw_dir: Path, accession: str) -> dict | None:
    """Load and return the raw InterPro JSON file for one accession, or None."""
    raw_file = raw_dir / f"{accession}.json"
    if not raw_file.exists():
        return None
    try:
        with open(raw_file, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return None