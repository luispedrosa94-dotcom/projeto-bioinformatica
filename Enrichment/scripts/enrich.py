"""
Stage 2 — Enrichment

Fetches the complete UniProt entry for each protein and saves the raw JSON
to outputs/uniprot_raw/{acc}.json. Also fetches InterPro coverage and saves
to outputs/interpro_raw/{acc}.json. Resolves GO_unknown aspects using the
GO term → aspect map extracted from UniProt records.

Usage:
    python scripts/enrich.py --config configs/default.yaml
    python scripts/enrich.py --config configs/default.yaml --workers 20
    python scripts/enrich.py --config configs/default.yaml --scope poorly_annotated
    python scripts/enrich.py --config configs/default.yaml --skip-interpro
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.clients.uniprot import fetch_proteins as fetch_uniprot
from src.clients.interpro import fetch_proteins as fetch_interpro
from src.schema import EnrichmentType


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
    parser.add_argument("--workers", type=int, default=10,
                        help="Parallel workers for UniProt/InterPro (default: 10)")
    parser.add_argument("--skip-interpro", action="store_true",
                        help="Skip InterPro fetching (faster for testing)")
    parser.add_argument("--skip-uniprot", action="store_true",
                        help="Skip UniProt fetching (useful when only adding InterPro)")
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

    cp_dir = output_root / "checkpoints"
    cp_dir.mkdir(exist_ok=True)

    scope = args.scope or cfg.get("scope", "all")

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
    uniprot_records = []
    if args.skip_uniprot:
        log.info("=== Step 1: UniProt API (SKIPPED) ===")
    else:
        log.info("=== Step 1: UniProt API (%d parallel workers) ===", args.workers)
        uniprot_records = fetch_uniprot(
            accessions,
            max_workers=args.workers,
            checkpoint_path=cp_dir / "uniprot.json",
            raw_dir=output_root / "uniprot_raw",
        )

    # ── Step 2: Resolve GO_unknown aspects ───────────────────────────────
    log.info("=== Step 2: GO aspect resolution (from UniProt data) ===")

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

    # ── Step 3: InterPro ──────────────────────────────────────────────────
    interpro_records = []
    if args.skip_interpro:
        log.info("=== Step 3: InterPro API (SKIPPED) ===")
    else:
        log.info("=== Step 3: InterPro API (%d parallel workers) ===", args.workers)
        interpro_records = fetch_interpro(
            accessions,
            max_workers=args.workers,
            checkpoint_path=cp_dir / "interpro.json",
            raw_dir=output_root / "interpro_raw",
        )

    # ── Write outputs ─────────────────────────────────────────────────────
    with open(output_root / "annotations.json", "w", encoding="utf-8") as f:
        json.dump(annotations, f, indent=2, ensure_ascii=False)
    log.info("Updated annotations → %s", output_root / "annotations.json")

    uniprot_out = [r.model_dump() for r in uniprot_records]
    with open(output_root / "uniprot_enrichment.json", "w", encoding="utf-8") as f:
        json.dump(uniprot_out, f, indent=2, ensure_ascii=False)
    log.info("UniProt enrichment → %d records", len(uniprot_out))

    with open(output_root / "go_aspect_map.json", "w", encoding="utf-8") as f:
        json.dump(go_aspect_map, f, indent=2, ensure_ascii=False)
    log.info("GO aspect map → %d terms", len(go_aspect_map))

    if not args.skip_interpro:
        interpro_out = [r.model_dump() for r in interpro_records]
        with open(output_root / "interpro_enrichment.json", "w", encoding="utf-8") as f:
            json.dump(interpro_out, f, indent=2, ensure_ascii=False)
        log.info("InterPro enrichment → %d records", len(interpro_out))

    # ── Summary ───────────────────────────────────────────────────────────
    from collections import Counter
    print("\n=== Enrichment summary ===")
    print(f"Proteins queried:         {len(accessions)}")
    print(f"UniProt records:          {len(uniprot_records)}")
    print(f"InterPro records:         {len(interpro_records)}")
    print(f"GO_unknown resolved:      {resolved_count} ({still_unknown} still unknown)")
    print(f"Raw UniProt files:        {len(list((output_root / 'uniprot_raw').glob('*.json')))}")
    print(f"Raw InterPro files:       {len(list((output_root / 'interpro_raw').glob('*.json')))}")
    if uniprot_records:
        print("\nUniProt records by type:")
        for t, n in sorted(Counter(r.enrichment_type.value for r in uniprot_records).items(), key=lambda x: -x[1]):
            print(f"  {t}: {n}")
    if interpro_records:
        print("\nInterPro records by type:")
        for t, n in sorted(Counter(r.enrichment_type.value for r in interpro_records).items(), key=lambda x: -x[1]):
            print(f"  {t}: {n}")


if __name__ == "__main__":
    main()