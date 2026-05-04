"""
Stage 2 — Enrichment

Usage:
    python scripts/enrich.py --config configs/default.yaml
    python scripts/enrich.py --config configs/default.yaml --scope poorly_annotated
    python scripts/enrich.py --config configs/default.yaml --skip-string
    python scripts/enrich.py --config configs/default.yaml --workers 20
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.clients.uniprot import fetch_proteins
from src.clients.string_db import map_identifiers, fetch_interactions, fetch_enrichment
from src.schema import EnrichmentRecord, EnrichmentType


def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def resolve_path(base: Path, p: str) -> Path:
    pp = Path(p)
    return pp if pp.is_absolute() else (base / pp).resolve()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/default.yaml")
    parser.add_argument("--scope", choices=["all", "poorly_annotated"], default=None)
    parser.add_argument("--skip-string", action="store_true")
    parser.add_argument("--workers", type=int, default=10,
                        help="Parallel workers for UniProt (default: 10)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
    log = logging.getLogger("enrich")

    cfg_path = Path(args.config)
    cfg = load_config(cfg_path)
    base = cfg_path.parent

    annotations_path = resolve_path(base, cfg["annotations_path"])
    proteins_path    = resolve_path(base, cfg["proteins_path"])
    output_root      = resolve_path(base, cfg["output_root"])
    output_root.mkdir(parents=True, exist_ok=True)

    # Checkpoint directory
    cp_dir = output_root / "checkpoints"
    cp_dir.mkdir(exist_ok=True)

    api_cfg = cfg.get("api", {})
    scope   = args.scope or cfg.get("scope", "all")

    # ── Load Stage 1 outputs ──────────────────────────────────────────────
    log.info("Loading Stage 1 outputs...")
    with open(annotations_path) as f:
        annotations = json.load(f)
    with open(proteins_path) as f:
        proteins = json.load(f)
    log.info("Loaded %d annotations, %d proteins", len(annotations), len(proteins))

    # ── Scope filter ──────────────────────────────────────────────────────
    if scope == "poorly_annotated":
        scope_accs = {p["uniprot_accession"] for p in proteins if p["in_poorly_annotated_subset"]}
        log.info("Scope: poorly_annotated — %d proteins", len(scope_accs))
    else:
        scope_accs = {p["uniprot_accession"] for p in proteins}
        log.info("Scope: all — %d proteins", len(scope_accs))

    accessions = sorted(scope_accs)

    # ── Step 1: UniProt ───────────────────────────────────────────────────
    log.info("=== Step 1: UniProt API (%d parallel workers) ===", args.workers)
    uniprot_records = fetch_proteins(
        accessions,
        max_workers=args.workers,
        checkpoint_path=cp_dir / "uniprot.json",
    )

    # ── Step 2: Resolve GO_unknown aspects ───────────────────────────────
    log.info("=== Step 2: GO aspect resolution (from UniProt data) ===")

    # UniProt already returns GO terms with aspect prefix (F:/P:/C:)
    # Extract a GO term → aspect map from the UniProt records we already fetched
    go_aspect_map: dict[str, str] = {}
    prefix_to_aspect = {
        EnrichmentType.GO_MF: "MF",
        EnrichmentType.GO_BP: "BP",
        EnrichmentType.GO_CC: "CC",
    }
    for rec in uniprot_records:
        if rec.enrichment_type in prefix_to_aspect:
            go_aspect_map[rec.value] = prefix_to_aspect[rec.enrichment_type]

    log.info("GO aspect map built from UniProt: %d terms", len(go_aspect_map))

    # Apply resolution to GO_unknown annotations
    resolved_count = 0
    still_unknown  = 0
    aspect_to_type = {"BP": "GO_BP", "MF": "GO_MF", "CC": "GO_CC"}
    for record in annotations:
        if record["annotation_type"] == "GO_unknown" and record["uniprot_accession"] in scope_accs:
            aspect = go_aspect_map.get(record["value"])
            if aspect:
                record["annotation_type"] = aspect_to_type[aspect]
                resolved_count += 1
            else:
                still_unknown += 1

    log.info("GO resolution: %d resolved, %d still unknown", resolved_count, still_unknown)

    # ── Step 3: STRING ────────────────────────────────────────────────────
    string_records: list[EnrichmentRecord] = []

    if not args.skip_string:
        log.info("=== Step 3: STRING API ===")
        string_id_map = map_identifiers(
            accessions,
            batch_size=api_cfg.get("string_batch_size", 500),
            delay=api_cfg.get("string_delay", 1.0),
        )
        log.info("STRING: mapped %d / %d accessions", len(string_id_map), len(accessions))

        if string_id_map:
            string_ids = list(string_id_map.values())
            interaction_records = fetch_interactions(
                string_ids,
                uniprot_map=string_id_map,
                min_score=api_cfg.get("string_min_score", 400),
                batch_size=api_cfg.get("string_batch_size", 500),
                delay=api_cfg.get("string_delay", 1.0),
            )
            string_records.extend(interaction_records)

            enrichment_records = fetch_enrichment(
                string_ids,
                delay=api_cfg.get("string_delay", 1.0),
            )
            string_records.extend(enrichment_records)
    else:
        log.info("Skipping STRING (--skip-string)")

    # ── Write outputs ─────────────────────────────────────────────────────
    # 1. Updated annotations (GO_unknown resolved)
    with open(output_root / "annotations.json", "w", encoding="utf-8") as f:
        json.dump(annotations, f, indent=2, ensure_ascii=False)
    log.info("Updated annotations → %s", output_root / "annotations.json")

    # 2. UniProt enrichment records
    uniprot_out = [r.model_dump() for r in uniprot_records]
    with open(output_root / "uniprot_enrichment.json", "w", encoding="utf-8") as f:
        json.dump(uniprot_out, f, indent=2, ensure_ascii=False)
    log.info("UniProt enrichment → %d records", len(uniprot_out))

    # 3. STRING records
    if string_records:
        string_out = [r.model_dump() for r in string_records]
        with open(output_root / "string_enrichment.json", "w", encoding="utf-8") as f:
            json.dump(string_out, f, indent=2, ensure_ascii=False)
        log.info("STRING enrichment → %d records", len(string_out))

    # 4. GO aspect map (audit trail)
    with open(output_root / "go_aspect_map.json", "w", encoding="utf-8") as f:
        json.dump(go_aspect_map, f, indent=2, ensure_ascii=False)
    log.info("GO aspect map → %d terms", len(go_aspect_map))

    # ── Summary ───────────────────────────────────────────────────────────
    from collections import Counter
    print("\n=== Enrichment summary ===")
    print(f"Proteins queried:        {len(accessions)}")
    print(f"UniProt records:         {len(uniprot_records)}")
    print(f"GO_unknown resolved:     {resolved_count} ({still_unknown} still unknown)")
    print(f"STRING interactions:     {sum(1 for r in string_records if r.enrichment_type == EnrichmentType.INTERACTION_PARTNER)}")
    print(f"STRING enrichment terms: {sum(1 for r in string_records if r.enrichment_type == EnrichmentType.FUNCTIONAL_ENRICHMENT)}")
    print("\nUniProt records by type:")
    for t, n in sorted(Counter(r.enrichment_type.value for r in uniprot_records).items(), key=lambda x: -x[1]):
        print(f"  {t}: {n}")


if __name__ == "__main__":
    main()
