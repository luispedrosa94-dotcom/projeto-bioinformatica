"""
Run normalization: read raw outputs, apply scope filter, write
annotations.parquet and proteins.parquet.

Usage:
    python -m scripts.normalize --config configs/default.yaml
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import replace
from pathlib import Path

import pandas as pd

# Make src/ importable when run as a script
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.config import load_config, load_scope_accessions
from src.parsers.base import filter_to_scope
from src.parsers.clean import CleanParser
from src.parsers.deepfri import DeepFRIParser
from src.parsers.deepgo2 import DeepGO2Parser
from src.parsers.eggnog import EggnogParser
from src.parsers.recognizer import RecognizerParser
from src.parsers.upimapi import UpimapiParser
from src.schema import ProteinRecord
from src.utils.ids import parse_protein_id


# Registry of available parsers. Add new parsers here as they're implemented.
# NOTE: Foldseek and ColabFold are out of scope per the supervisor's decision —
# only homology- and ML-based tools are processed in this iteration.
PARSER_REGISTRY = {
    "upimapi": UpimapiParser,
    "eggnog": EggnogParser,
    "clean": CleanParser,
    "deepgo2": DeepGO2Parser,
    "deepfri": DeepFRIParser,
    "recognizer": RecognizerParser,
}


def _build_parser(tool_name: str, cfg) -> object:
    """Instantiate a parser, passing tool-specific config."""
    if tool_name == "upimapi":
        return UpimapiParser(raw_data_root=cfg.raw_data_root, thresholds=cfg.upimapi_thresholds)
    if tool_name == "eggnog":
        return EggnogParser(raw_data_root=cfg.raw_data_root, thresholds=cfg.eggnog_thresholds)
    if tool_name == "clean":
        return CleanParser(raw_data_root=cfg.raw_data_root, thresholds=cfg.clean_thresholds)
    if tool_name == "deepgo2":
        return DeepGO2Parser(
            raw_data_root=cfg.raw_data_root,
            thresholds=cfg.deepgo2_thresholds,
            min_score=cfg.deepgo2_min_score,
        )
    if tool_name == "deepfri":
        return DeepFRIParser(
            raw_data_root=cfg.raw_data_root,
            thresholds=cfg.deepfri_thresholds,
            min_score=cfg.deepfri_min_score,
        )
    if tool_name == "recognizer":
        return RecognizerParser(
            raw_data_root=cfg.raw_data_root,
            thresholds=cfg.recognizer_thresholds,
            top_n_per_db=cfg.recognizer_top_n_per_db,
        )
    # Future parsers receive their own thresholds from cfg here.
    raise KeyError(tool_name)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/default.yaml")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
    log = logging.getLogger("normalize")

    cfg = load_config(args.config)
    log.info("Config loaded: scope=%s, output=%s", cfg.scope_mode, cfg.output_root)

    scope = load_scope_accessions(cfg)
    if scope is None:
        log.info("Scope: ALL proteins (no filter)")
    else:
        log.info("Scope: %d accessions (mode=%s)", len(scope), cfg.scope_mode)

    cfg.output_root.mkdir(parents=True, exist_ok=True)
    (cfg.output_root / "01_normalization").mkdir(exist_ok=True)

    # Run each enabled parser
    all_records: list = []
    for tool_name, enabled in cfg.tools_enabled.items():
        if not enabled:
            log.info("Skipping %s (disabled)", tool_name)
            continue
        if tool_name not in PARSER_REGISTRY:
            log.warning("No parser implemented for %s — skipping", tool_name)
            continue
        instance = _build_parser(tool_name, cfg)
        recs = instance.parse()
        recs_in_scope = filter_to_scope(recs, scope)
        log.info(
            "%s: %d records parsed, %d in scope",
            tool_name, len(recs), len(recs_in_scope),
        )
        all_records.extend(recs_in_scope)

    if not all_records:
        log.warning("No records produced — exiting.")
        return

    # Build annotations DataFrame
    annotations_df = pd.DataFrame([r.model_dump() for r in all_records])
    # Convert enum values to strings for serialization
    annotations_df["source_tool"] = annotations_df["source_tool"].astype(str)
    annotations_df["annotation_type"] = annotations_df["annotation_type"].astype(str)
    annotations_df["score_type"] = annotations_df["score_type"].astype(str)
    annotations_df["confidence_level"] = annotations_df["confidence_level"].astype(str)
    # raw_extras stays as native dict — JSON serializes it as a nested object,
    # which is more readable than a stringified JSON inside a string field.

    annotations_path = cfg.output_root / "01_normalization/annotations.json"
    annotations_df.to_json(
        annotations_path,
        orient="records",
        indent=2,
        force_ascii=False,
    )
    log.info("Wrote %d annotation rows → %s", len(annotations_df), annotations_path)

    # Build proteins table (one row per unique accession).
    # Compute the poorly-annotated set independently of the active scope
    # so we can flag proteins regardless of how scope was set.
    poorly_set: set[str] = set()
    if cfg.poorly_annotated_file and cfg.poorly_annotated_file.exists():
        poorly_cfg = replace(cfg, scope_mode="poorly_annotated")
        poorly_set = load_scope_accessions(poorly_cfg) or set()

    proteins: dict[str, ProteinRecord] = {}
    for r in all_records:
        if r.uniprot_accession in proteins:
            continue
        parsed = parse_protein_id(r.original_id)
        proteins[r.uniprot_accession] = ProteinRecord(
            uniprot_accession=r.uniprot_accession,
            original_id=r.original_id,
            db_source=parsed.db_source,
            entry_name=parsed.entry_name,
            in_poorly_annotated_subset=r.uniprot_accession in poorly_set,
        )

    proteins_df = pd.DataFrame([p.model_dump() for p in proteins.values()])
    proteins_path = cfg.output_root / "01_normalization/proteins.json"
    proteins_df.to_json(
        proteins_path,
        orient="records",
        indent=2,
        force_ascii=False,
    )
    log.info("Wrote %d protein rows → %s", len(proteins_df), proteins_path)

    # Quick QC summary
    print("\n=== Normalization summary ===")
    print(f"Proteins:     {len(proteins_df)}")
    print(f"Annotations:  {len(annotations_df)}")
    print("\nAnnotations by tool:")
    print(annotations_df.groupby("source_tool").size().to_string())
    print("\nAnnotations by type:")
    print(annotations_df.groupby("annotation_type").size().sort_values(ascending=False).to_string())
    print("\nAnnotations by confidence level:")
    print(annotations_df.groupby("confidence_level").size().to_string())
    print("\nConfidence × tool:")
    print(
        annotations_df.groupby(["source_tool", "confidence_level"]).size()
        .unstack(fill_value=0)
        .to_string()
    )


if __name__ == "__main__":
    main()
