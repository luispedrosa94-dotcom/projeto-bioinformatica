"""
Configuration loading and scope filtering.

The scope determines which subset of proteins enters the pipeline. Supported
modes:

    all                — every protein found in any input file
    poorly_annotated   — the curated subset of poorly-annotated proteins
    custom             — a user-provided list of UniProt accessions

Per-tool confidence thresholds for the high/medium/low classification are
also loaded here and passed to parsers via PipelineConfig.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from .utils.confidence import HomologyThresholds, MLThresholds, DEFAULT_CLEAN
from .utils.ids import _ACCESSION_RE


@dataclass
class PipelineConfig:
    # Where the raw tool outputs live
    raw_data_root: Path
    # Where normalized parquet files will be written
    output_root: Path

    # Scope filter
    scope_mode: str = "all"
    poorly_annotated_file: Optional[Path] = None
    custom_accessions: list[str] = field(default_factory=list)

    # Per-tool parser switches
    tools_enabled: dict[str, bool] = field(default_factory=dict)

    # Parser-level knobs
    recognizer_top_n_per_db: int = 5
    deepfri_min_score: float = 0.0
    deepgo2_min_score: float = 0.0

    # Confidence-level thresholds, per tool
    upimapi_thresholds: HomologyThresholds = field(default_factory=HomologyThresholds)
    eggnog_thresholds: HomologyThresholds = field(default_factory=HomologyThresholds)
    recognizer_thresholds: HomologyThresholds = field(default_factory=HomologyThresholds)
    deepfri_thresholds: MLThresholds = field(default_factory=MLThresholds)
    deepgo2_thresholds: MLThresholds = field(default_factory=MLThresholds)
    clean_thresholds: MLThresholds = field(default_factory=lambda: DEFAULT_CLEAN)


def _homology_thresholds_from(d: Optional[dict]) -> HomologyThresholds:
    if not d:
        return HomologyThresholds()
    return HomologyThresholds(
        high=float(d.get("high", 1e-50)),
        medium=float(d.get("medium", 1e-10)),
    )


def _ml_thresholds_from(d: Optional[dict], default: MLThresholds) -> MLThresholds:
    if not d:
        return default
    return MLThresholds(
        high=float(d.get("high", default.high)),
        medium=float(d.get("medium", default.medium)),
    )


def load_config(path: str | Path) -> PipelineConfig:
    """Load a YAML config and resolve relative paths against the config file's directory."""
    path = Path(path)
    with path.open("r") as f:
        raw = yaml.safe_load(f)

    base = path.parent

    def _resolve(p: Optional[str]) -> Optional[Path]:
        if p is None:
            return None
        pp = Path(p)
        return pp if pp.is_absolute() else (base / pp).resolve()

    confidence = raw.get("confidence_thresholds", {}) or {}

    return PipelineConfig(
        raw_data_root=_resolve(raw["raw_data_root"]),
        output_root=_resolve(raw["output_root"]),
        scope_mode=raw.get("scope", {}).get("mode", "all"),
        poorly_annotated_file=_resolve(raw.get("scope", {}).get("poorly_annotated_file")),
        custom_accessions=raw.get("scope", {}).get("custom_accessions", []) or [],
        tools_enabled=raw.get("tools_enabled", {}),
        recognizer_top_n_per_db=raw.get("parsers", {}).get("recognizer_top_n_per_db", 5),
        deepfri_min_score=raw.get("parsers", {}).get("deepfri_min_score", 0.0),
        deepgo2_min_score=raw.get("parsers", {}).get("deepgo2_min_score", 0.0),
        upimapi_thresholds=_homology_thresholds_from(confidence.get("upimapi")),
        eggnog_thresholds=_homology_thresholds_from(confidence.get("eggnog")),
        recognizer_thresholds=_homology_thresholds_from(confidence.get("recognizer")),
        deepfri_thresholds=_ml_thresholds_from(confidence.get("deepfri"), MLThresholds()),
        deepgo2_thresholds=_ml_thresholds_from(confidence.get("deepgo2"), MLThresholds()),
        clean_thresholds=_ml_thresholds_from(confidence.get("clean"), DEFAULT_CLEAN),
    )


def load_scope_accessions(cfg: PipelineConfig) -> Optional[set[str]]:
    """Return the set of UniProt accessions to keep, or None if scope='all'."""
    mode = cfg.scope_mode.lower()
    if mode == "all":
        return None

    if mode == "custom":
        return {a.strip().upper() for a in cfg.custom_accessions if a.strip()}

    if mode == "poorly_annotated":
        if not cfg.poorly_annotated_file or not cfg.poorly_annotated_file.exists():
            raise FileNotFoundError(
                f"scope=poorly_annotated requires a valid poorly_annotated_file; got {cfg.poorly_annotated_file}"
            )
        accessions: set[str] = set()
        with cfg.poorly_annotated_file.open("r") as f:
            for line in f:
                for m in _ACCESSION_RE.finditer(line.upper()):
                    accessions.add(m.group(1))
        return accessions

    raise ValueError(f"Unknown scope mode: {cfg.scope_mode!r}")
