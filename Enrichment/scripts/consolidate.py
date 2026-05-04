"""
Stage 2b — Consolidation

Aggregates all annotations (Stage 1 + Stage 2) into a single wide-format
table with one row per protein. This is the direct input to the LLM
harmonization stage (Stage 3).

For each protein, produces:
  - Identity: accession, name, gene, reviewed status, organism hint
  - Functional: GO terms (BP/MF/CC), EC numbers, KEGG pathways/KOs
  - Structural: domains (COG, Pfam, TIGR, SMART, KOG)
  - Contextual: subcellular location, function description, keywords
  - Evidence: which tools annotated it, confidence distribution
  - STRING: set-level functional enrichment terms

Usage:
    python scripts/consolidate.py --config configs/default.yaml
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def resolve_path(base: Path, p: str) -> Path:
    pp = Path(p)
    return pp if pp.is_absolute() else (base / pp).resolve()


def _dedupe(items: list) -> list:
    """Deduplicate preserving order."""
    seen = set()
    out = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _best_value(records: list[dict], annotation_type: str) -> str | None:
    """Return the value of the first HIGH-confidence record of given type, or any."""
    matches = [r for r in records if r.get("annotation_type") == annotation_type]
    if not matches:
        return None
    high = [r for r in matches if r.get("confidence_level") == "high"]
    return (high or matches)[0]["value"]


def _all_values(records: list[dict], annotation_type: str,
                min_confidence: str | None = None) -> list[str]:
    """Return all unique values for a given annotation type."""
    order = {"high": 0, "medium": 1, "low": 2, "unknown": 3}
    matches = [r for r in records if r.get("annotation_type") == annotation_type]
    if min_confidence:
        cutoff = order[min_confidence]
        matches = [r for r in matches if order.get(r.get("confidence_level", "unknown"), 3) <= cutoff]
    matches.sort(key=lambda r: order.get(r.get("confidence_level", "unknown"), 3))
    return _dedupe([r["value"] for r in matches])


def _all_values_with_sources(records: list[dict], annotation_type: str) -> list[dict]:
    """Return values with their sources and confidence for a given type."""
    seen: dict[str, dict] = {}
    order = {"high": 0, "medium": 1, "low": 2, "unknown": 3}
    for r in records:
        if r.get("annotation_type") != annotation_type:
            continue
        val = r["value"]
        tool = r.get("source_tool", "").replace("SourceTool.", "")
        conf = r.get("confidence_level", "unknown").replace("ConfidenceLevel.", "")
        if val not in seen:
            seen[val] = {"value": val, "sources": [], "confidence": conf, "label": r.get("label")}
        if tool not in seen[val]["sources"]:
            seen[val]["sources"].append(tool)
        # Keep best confidence
        if order.get(conf, 3) < order.get(seen[val]["confidence"], 3):
            seen[val]["confidence"] = conf
    return sorted(seen.values(), key=lambda x: order.get(x["confidence"], 3))


def _overall_confidence(records: list[dict]) -> str:
    """Derive overall annotation confidence from the distribution."""
    counts = defaultdict(int)
    for r in records:
        conf = r.get("confidence_level", "unknown").replace("ConfidenceLevel.", "")
        counts[conf] += 1
    total = sum(counts.values())
    if total == 0:
        return "unknown"
    high_frac = counts["high"] / total
    medium_frac = counts["medium"] / total
    if high_frac >= 0.5:
        return "high"
    if high_frac + medium_frac >= 0.5:
        return "medium"
    return "low"


def consolidate(
    proteins_path: Path,
    annotations_path: Path,
    uniprot_path: Path,
    string_path: Path | None,
    output_path: Path,
) -> None:
    log = logging.getLogger("consolidate")

    log.info("Loading data...")
    with open(proteins_path) as f:
        proteins = json.load(f)
    with open(annotations_path) as f:
        annotations = json.load(f)
    with open(uniprot_path) as f:
        uniprot_records = json.load(f)

    string_enrichment: list[dict] = []
    if string_path and string_path.exists():
        with open(string_path) as f:
            string_enrichment = json.load(f)

    log.info("Loaded %d proteins, %d annotations, %d UniProt records",
             len(proteins), len(annotations), len(uniprot_records))

    # Index annotations and UniProt records by accession
    ann_by_acc: dict[str, list[dict]] = defaultdict(list)
    for r in annotations:
        ann_by_acc[r["uniprot_accession"]].append(r)

    uni_by_acc: dict[str, list[dict]] = defaultdict(list)
    for r in uniprot_records:
        uni_by_acc[r["uniprot_accession"]].append(r)

    # STRING enrichment is set-level — same for all proteins
    string_terms = [
        {"term": r["value"], "label": r.get("label", ""), "fdr": r.get("score")}
        for r in string_enrichment
        if r.get("enrichment_type") == "functional_enrichment"
    ]

    consolidated: list[dict] = []

    for prot in proteins:
        acc = prot["uniprot_accession"]
        ann = ann_by_acc.get(acc, [])
        uni = uni_by_acc.get(acc, [])

        # ── Identity ────────────────────────────────────────────────────
        reviewed = next(
            (r["value"] for r in uni if r.get("enrichment_type") == "reviewed_status"),
            "unreviewed"
        )
        protein_name = next(
            (r["value"] for r in uni if r.get("enrichment_type") == "protein_name"),
            _best_value(ann, "protein_description")
        )
        gene_name = next(
            (r["value"] for r in uni if r.get("enrichment_type") == "gene_name"),
            _best_value(ann, "preferred_gene_name")
        )
        subcellular_location = _dedupe([
            r["value"] for r in uni if r.get("enrichment_type") == "subcellular_location"
        ])
        function_description = next(
            (r["value"] for r in uni if r.get("enrichment_type") == "function_description"),
            _best_value(ann, "function_cc")
        )
        keywords = _dedupe([
            r["value"] for r in uni if r.get("enrichment_type") == "keyword"
        ])

        # ── GO terms ─────────────────────────────────────────────────────
        # Merge GO terms from Stage 1 and UniProt, deduplicating
        go_bp = _all_values_with_sources(ann, "GO_BP")
        go_mf = _all_values_with_sources(ann, "GO_MF")
        go_cc = _all_values_with_sources(ann, "GO_CC")

        # Add UniProt GO terms (already have aspect)
        uni_go_by_type = {"GO_BP": [], "GO_MF": [], "GO_CC": []}
        for r in uni:
            if r.get("enrichment_type") in uni_go_by_type:
                uni_go_by_type[r["enrichment_type"]].append({
                    "value": r["value"],
                    "sources": ["uniprot"],
                    "confidence": "high",  # UniProt GO terms are curated
                    "label": r.get("label"),
                })

        def _merge_go(stage1: list[dict], uniprot_go: list[dict]) -> list[dict]:
            seen_vals = {r["value"] for r in stage1}
            merged = stage1.copy()
            for r in uniprot_go:
                if r["value"] in seen_vals:
                    # Add UniProt as additional source to existing record
                    existing = next(x for x in merged if x["value"] == r["value"])
                    if "uniprot" not in existing["sources"]:
                        existing["sources"].append("uniprot")
                        existing["confidence"] = "high"  # curated evidence upgrades confidence
                else:
                    merged.append(r)
                    seen_vals.add(r["value"])
            return merged

        go_bp = _merge_go(go_bp, uni_go_by_type["GO_BP"])
        go_mf = _merge_go(go_mf, uni_go_by_type["GO_MF"])
        go_cc = _merge_go(go_cc, uni_go_by_type["GO_CC"])

        # ── EC numbers ───────────────────────────────────────────────────
        ec_numbers = _all_values_with_sources(ann, "EC")

        # ── KEGG ─────────────────────────────────────────────────────────
        kegg_ko       = _all_values(ann, "KEGG_ko")
        kegg_pathways = _all_values(ann, "KEGG_pathway")
        kegg_modules  = _all_values(ann, "KEGG_module")

        # ── Domains ──────────────────────────────────────────────────────
        domains = {
            "cog":    _all_values_with_sources(ann, "cog"),
            "kog":    _all_values_with_sources(ann, "kog"),
            "pfam":   _all_values_with_sources(ann, "pfam"),
            "tigrfam":_all_values_with_sources(ann, "tigrfam"),
            "smart":  _all_values_with_sources(ann, "smart"),
        }
        cog_categories = _all_values(ann, "cog_category")

        # ── Evidence summary ─────────────────────────────────────────────
        tools = _dedupe([
            r["source_tool"].replace("SourceTool.", "")
            for r in ann
        ])
        conf_counts = defaultdict(int)
        for r in ann:
            c = r.get("confidence_level", "unknown").replace("ConfidenceLevel.", "")
            conf_counts[c] += 1
        overall_conf = _overall_confidence(ann)

        # ── Assemble record ──────────────────────────────────────────────
        consolidated.append({
            # Identity
            "uniprot_accession": acc,
            "original_id": prot["original_id"],
            "db_source": prot["db_source"],
            "entry_name": prot.get("entry_name"),
            "in_poorly_annotated_subset": prot.get("in_poorly_annotated_subset", False),

            # Curated metadata (UniProt)
            "reviewed": reviewed,
            "protein_name": protein_name,
            "gene_name": gene_name,
            "subcellular_location": subcellular_location,
            "function_description": function_description,
            "keywords": keywords,

            # GO terms (merged Stage 1 + UniProt)
            "go_bp": go_bp,
            "go_mf": go_mf,
            "go_cc": go_cc,

            # Functional identifiers
            "ec_numbers": ec_numbers,
            "kegg_ko": kegg_ko,
            "kegg_pathways": kegg_pathways,
            "kegg_modules": kegg_modules,

            # Domain annotations
            "domains": domains,
            "cog_categories": cog_categories,

            # Evidence summary
            "tools": tools,
            "annotation_count": len(ann),
            "confidence_distribution": dict(conf_counts),
            "overall_confidence": overall_conf,

            # STRING set-level enrichment (same for all proteins in the set)
            "string_set_enrichment": string_terms,
        })

    log.info("Consolidated %d protein records", len(consolidated))

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(consolidated, f, indent=2, ensure_ascii=False)
    log.info("Written → %s", output_path)

    # Summary
    reviewed_count = sum(1 for r in consolidated if r["reviewed"] == "reviewed")
    with_name      = sum(1 for r in consolidated if r["protein_name"])
    with_go        = sum(1 for r in consolidated if r["go_bp"] or r["go_mf"] or r["go_cc"])
    with_ec        = sum(1 for r in consolidated if r["ec_numbers"])
    with_func      = sum(1 for r in consolidated if r["function_description"])
    conf_dist      = defaultdict(int)
    for r in consolidated:
        conf_dist[r["overall_confidence"]] += 1

    print("\n=== Consolidation summary ===")
    print(f"Proteins:                  {len(consolidated)}")
    print(f"  Reviewed (SwissProt):    {reviewed_count}")
    print(f"  Unreviewed (TrEMBL):     {len(consolidated) - reviewed_count}")
    print(f"With protein name:         {with_name}")
    print(f"With GO terms:             {with_go}")
    print(f"With EC number:            {with_ec}")
    print(f"With function description: {with_func}")
    print(f"\nOverall confidence:")
    for c in ["high", "medium", "low", "unknown"]:
        print(f"  {c}: {conf_dist[c]}")
    print(f"\nSTRING set-level enrichment: {len(string_terms)} terms")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/default.yaml")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    cfg_path = Path(args.config)
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    base = cfg_path.parent

    output_root = resolve_path(base, cfg["output_root"])

    consolidate(
        proteins_path   = resolve_path(base, cfg["proteins_path"]),
        annotations_path= resolve_path(base, cfg["annotations_path"]),
        uniprot_path    = output_root / "uniprot_enrichment.json",
        string_path     = output_root / "string_enrichment.json",
        output_path     = output_root / "protein_profiles.json",
    )


if __name__ == "__main__":
    main()
