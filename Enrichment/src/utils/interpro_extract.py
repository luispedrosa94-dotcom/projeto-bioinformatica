"""
Parser for InterPro /entry/all/ raw JSON.

Reads outputs/caches/interpro_raw/{acc}.json (saved by interpro.py client) and
builds the 'interpro' section of one protein profile.

Preserves ALL fields the API exposes:
  - Per entry: accession, name, source_database, type, integrated,
    member_databases, go_terms, hierarchy, cross_references, literature,
    description, entry_id, entry_date
  - Per protein location: fragments (start, end, dc-status),
    representative, model, score, subfamily, protein_length, organism,
    in_alphafold, in_bfvd
  - Summary statistics aggregating across all entries
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


def _extract_protein_locations(proteins_list: list) -> list[dict]:
    """
    Extract the per-protein location info from the 'proteins' array.
    Each entry from /entry/all/ contains a proteins[] list — usually with
    one element (the protein we queried). We preserve all useful fields.
    """
    blocks = []
    for prot in proteins_list or []:
        blocks.append({
            "accession":             prot.get("accession") or None,
            "protein_length":        prot.get("protein_length"),
            "source_database":       prot.get("source_database") or None,
            "organism":              prot.get("organism") or None,
            "in_alphafold":          prot.get("in_alphafold"),
            "in_bfvd":               prot.get("in_bfvd"),
            "entry_protein_locations": prot.get("entry_protein_locations") or [],
        })
    return blocks


def _normalise_go_terms(go_terms_raw) -> dict[str, list[dict]]:
    """
    Convert the raw 'go_terms' list into a dict grouped by GO category:
        {
          "molecular_function": [{"go_id": "GO:...", "name": "...", ...}, ...],
          "biological_process": [...],
          "cellular_component": [...],
        }

    The API uses either short codes (F/P/C) or full names. We handle both.
    """
    by_category = {
        "molecular_function": [],
        "biological_process": [],
        "cellular_component": [],
    }
    if not go_terms_raw:
        return by_category

    category_map = {
        "molecular_function": "molecular_function",
        "biological_process": "biological_process",
        "cellular_component": "cellular_component",
        "F": "molecular_function",
        "P": "biological_process",
        "C": "cellular_component",
    }

    for g in go_terms_raw:
        cat_raw = g.get("category", {}) or {}
        if isinstance(cat_raw, dict):
            cat_key = cat_raw.get("code") or cat_raw.get("name") or ""
        else:
            cat_key = str(cat_raw)
        cat = category_map.get(cat_key, "molecular_function")
        by_category[cat].append({
            "go_id":        g.get("identifier") or g.get("id") or None,
            "name":         g.get("name") or None,
            "category":     cat,
            "category_raw": cat_key or None,
        })

    return by_category


def _flatten_member_databases(member_dbs: dict) -> dict[str, list[dict]]:
    """
    Convert the raw 'member_databases' nested dict into a cleaner structure:
        {
            "pfam":       [{"accession": "PF00890", "name": "FAD binding domain"}, ...],
            "ncbifam":    [{"accession": "TIGR02061", "name": "..."}, ...],
        }

    Empty/None member_databases become an empty dict.
    """
    if not member_dbs:
        return {}
    clean = {}
    for db_name, sigs in member_dbs.items():
        if isinstance(sigs, dict):
            clean[db_name] = [
                {"accession": acc, "name": name}
                for acc, name in sigs.items()
            ]
        else:
            clean[db_name] = sigs
    return clean


def _extract_entry(result: dict) -> dict:
    """
    Convert one raw entry block from the /entry/all/ response into a
    canonical dict preserving all useful fields.
    """
    meta = result.get("metadata") or {}
    proteins = result.get("proteins") or []

    return {
        # Identity
        "accession":         meta.get("accession") or None,
        "name":              meta.get("name") or None,
        "short_name":        meta.get("name_short") or meta.get("short_name") or None,
        "description":       meta.get("description") or [],
        "type":              meta.get("type") or None,
        "source_database":   meta.get("source_database") or None,

        # Integration relationship
        "integrated":        meta.get("integrated") or None,         # signature → IPR
        "member_databases":  _flatten_member_databases(meta.get("member_databases")),  # IPR → signatures

        # Biology
        "go_terms":          _normalise_go_terms(meta.get("go_terms")),

        # Cross-references / hierarchy / literature
        "hierarchy":         meta.get("hierarchy") or {},
        "cross_references":  meta.get("cross_references") or {},
        "literature":        meta.get("literature") or {},

        # Versioning
        "entry_id":          meta.get("entry_id") or None,
        "entry_date":        meta.get("entry_date") or None,

        # Per-protein evidence (positions, e-values, models, subfamily)
        "proteins":          _extract_protein_locations(proteins),
    }


def _build_summary(entries: list[dict]) -> dict:
    """
    Build summary statistics across all entries for this protein.
    Useful for the Streamlit app and for quick filtering.
    """
    by_source_db = Counter()
    by_type      = Counter()

    interpro_integrated = []   # IPRxxxxx entries (source_database == "interpro")
    unintegrated_sigs   = []   # signatures with integrated == None

    # GO term collection — unique by go_id across all entries
    go_unique = {}             # go_id → {go_id, name, category, sources: [...]}
    go_counts_by_cat = {
        "molecular_function": 0,
        "biological_process": 0,
        "cellular_component": 0,
    }

    member_dbs_used = set()

    for e in entries:
        src_db = e.get("source_database") or "unknown"
        by_source_db[src_db] += 1

        etype = e.get("type") or "unknown"
        by_type[etype] += 1

        if src_db == "interpro":
            interpro_integrated.append(e.get("accession"))
        elif e.get("integrated") is None:
            unintegrated_sigs.append(e.get("accession"))

        # Member databases referenced by this entry
        for db in (e.get("member_databases") or {}).keys():
            member_dbs_used.add(db)

        # GO terms — collect unique, track which entries cite each
        for cat, terms in (e.get("go_terms") or {}).items():
            for term in terms:
                gid = term.get("go_id")
                if not gid:
                    continue
                if gid not in go_unique:
                    go_unique[gid] = {
                        "go_id":    gid,
                        "name":     term.get("name"),
                        "category": cat,
                        "sources":  [],
                    }
                    go_counts_by_cat[cat] = go_counts_by_cat.get(cat, 0) + 1
                src_entry = e.get("accession")
                if src_entry and src_entry not in go_unique[gid]["sources"]:
                    go_unique[gid]["sources"].append(src_entry)

    return {
        "total_entries":             len(entries),
        "total_interpro_integrated": len(interpro_integrated),
        "total_unintegrated":        len(unintegrated_sigs),
        "by_source_database":        dict(by_source_db),
        "by_type":                   dict(by_type),
        "member_databases_used":     sorted(member_dbs_used),
        "go_terms_count":            go_counts_by_cat,
        "go_terms_list":             list(go_unique.values()),
        "interpro_integrated_ids":   sorted(set(filter(None, interpro_integrated))),
        "unintegrated_ids":          sorted(set(filter(None, unintegrated_sigs))),
    }


def extract_interpro_section(raw_data: dict) -> dict:
    """
    Build the 'interpro' section of one protein profile from the raw JSON
    saved by the interpro.py client.

    Input raw_data shape (from /entry/all/protein/uniprot/{acc}/):
        {
            "accession": "Q46505",
            "data": {
                "count": 17,
                "results": [
                    {"metadata": {...}, "proteins": [...]},
                    ...
                ]
            }
        }

    Output:
        {
            "entries": [ { ...all 17 entries, fully expanded... } ],
            "summary": { ...aggregated statistics... }
        }
    """
    data = raw_data.get("data") or {}
    results = data.get("results") or []

    entries = [_extract_entry(r) for r in results]
    summary = _build_summary(entries)

    return {
        "entries": entries,
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