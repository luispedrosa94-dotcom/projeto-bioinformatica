"""
Confidence-level mapping rules.

Each tool reports confidence on its own scale. To support the supervisor's
requirement of categorical levels (`high` / `medium` / `low`), every parser
must call one of these helpers. Thresholds are configurable from YAML so
they can be re-tuned without code changes.

Default thresholds are the starting points discussed at the supervision
meeting and are NOT meant to be final. They should be refined empirically
(e.g. by inspecting score distributions per tool on the working dataset)
and justified in the final article's methods section.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from ..schema import ConfidenceLevel


@dataclass(frozen=True)
class HomologyThresholds:
    """For e-value-based tools (UPIMAPI, eggNOG-mapper, reCOGnizer).

    Smaller e-value = stronger evidence. `high` is the upper bound of the
    HIGH band (inclusive); `medium` is the upper bound of the MEDIUM band.
    Anything above `medium` is LOW.
    """
    high: float = 1e-50
    medium: float = 1e-10


@dataclass(frozen=True)
class MLThresholds:
    """For ML-based tools (DeepFRI, DeepGO2, CLEAN).

    Larger score = stronger evidence. `high` is the lower bound of the HIGH
    band (inclusive); `medium` is the lower bound of the MEDIUM band.
    Anything below `medium` is LOW.
    """
    high: float = 0.7
    medium: float = 0.3


@dataclass(frozen=True)
class StructuralThresholds:
    """For structure-search (Foldseek) and structure-prediction (ColabFold).

    Currently OUT OF SCOPE per supervisor decision — kept here for future
    re-inclusion without schema migration.
    """
    foldseek_high_bitscore: float = 100.0
    foldseek_medium_bitscore: float = 50.0
    foldseek_high_evalue: float = 1e-10
    plddt_high: float = 90.0
    plddt_medium: float = 70.0


# ---- Per-tool defaults registry ----
# Per-tool overrides are applied on top of these by reading the YAML config.

DEFAULT_HOMOLOGY = HomologyThresholds()
DEFAULT_ML = MLThresholds()

# Tool-specific overrides discussed at the meeting:
DEFAULT_CLEAN = MLThresholds(high=0.5, medium=0.1)


def from_evalue(
    evalue: Optional[float],
    thresholds: HomologyThresholds = DEFAULT_HOMOLOGY,
) -> ConfidenceLevel:
    """Map a homology-search e-value to a categorical confidence level."""
    if evalue is None:
        return ConfidenceLevel.UNKNOWN
    try:
        e = float(evalue)
    except (TypeError, ValueError):
        return ConfidenceLevel.UNKNOWN
    if e <= thresholds.high:
        return ConfidenceLevel.HIGH
    if e <= thresholds.medium:
        return ConfidenceLevel.MEDIUM
    return ConfidenceLevel.LOW


def from_ml_score(
    score: Optional[float],
    thresholds: MLThresholds = DEFAULT_ML,
) -> ConfidenceLevel:
    """Map an ML confidence score (0..1, larger = better) to a categorical level."""
    if score is None:
        return ConfidenceLevel.UNKNOWN
    try:
        s = float(score)
    except (TypeError, ValueError):
        return ConfidenceLevel.UNKNOWN
    if s >= thresholds.high:
        return ConfidenceLevel.HIGH
    if s >= thresholds.medium:
        return ConfidenceLevel.MEDIUM
    return ConfidenceLevel.LOW
