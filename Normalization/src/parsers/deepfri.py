"""
DeepFRI parser.

DeepFRI is a deep-learning function predictor that takes a protein sequence
and predicts Gene Ontology terms — specifically, Molecular Function (MF)
terms only. Unlike DeepGO2 which spreads its predictions across three
files (BP/CC/MF), DeepFRI focuses on MF and produces a single CSV.

The CSV has a header and four columns:
    Protein, GO_term/EC_number, Score, GO_term/EC_number name

Despite the column name suggesting both GO and EC, in practice all entries
are GO terms (because DeepFRI's MF model only predicts GOs). The 'name'
column gives human-readable descriptions ("lyase activity", "hydrolase
activity, acting on ester bonds", etc.) which we use as the `label` field.

Note: there is an auxiliary `DeepFRI_MF_pred_scores.json` file containing
all class scores (not just top hits). We ignore it because the predictions
CSV already contains the filtered, ranked top hits — which is what we want.

Confidence: ML score in [0, 1]. We use the same default ML thresholds
as DeepGO2 (high>=0.7, medium>=0.3) which can be tuned per-tool in the YAML.

evidence_rank: when a protein has multiple predictions, sort by score
descending and assign rank 1 to the most confident.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from ..schema import AnnotationRecord, AnnotationType, ScoreType, SourceTool
from ..utils.confidence import MLThresholds, from_ml_score
from ..utils.ids import parse_protein_id
from .base import BaseParser


class DeepFRIParser(BaseParser):
    tool_name = "deepfri"

    rel_path = Path("deepfri") / "DeepFRI_MF_predictions.csv"

    def __init__(
        self,
        raw_data_root: Path,
        thresholds: MLThresholds | None = None,
        min_score: float = 0.0,
    ):
        super().__init__(raw_data_root)
        self.thresholds = thresholds or MLThresholds()
        # min_score: predictions strictly below this are dropped at parse time.
        # Default 0.0 keeps everything (filtering deferred to later stages).
        self.min_score = float(min_score)

    def parse(self) -> list[AnnotationRecord]:
        path = self.raw_data_root / self.rel_path
        if not path.exists():
            self.log.warning("DeepFRI file not found: %s — skipping", path)
            return []

        # Read all rows, grouping by protein for ranking.
        # The CSV uses standard quoting so labels with commas
        # (e.g. "hydrolase activity, acting on ester bonds") are handled
        # automatically by the csv module.
        by_protein: dict[str, list[tuple[str, float, str]]] = defaultdict(list)
        rows_processed = 0
        skipped_low_score = 0

        with path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows_processed += 1
                raw_id = (row.get("Protein") or "").strip()
                go_term = (row.get("GO_term/EC_number") or "").strip()
                score_str = (row.get("Score") or "").strip()
                label = (row.get("GO_term/EC_number name") or "").strip()

                if not raw_id or not go_term:
                    continue
                try:
                    score = float(score_str)
                except ValueError:
                    continue

                if score < self.min_score:
                    skipped_low_score += 1
                    continue

                # Defensive: skip non-GO entries if any ever appear
                # (the column header allows EC numbers in principle).
                if not go_term.startswith("GO:"):
                    continue

                by_protein[raw_id].append((go_term, score, label or None))  # type: ignore[arg-type]

        records: list[AnnotationRecord] = []

        # For each protein, sort by score descending and assign rank
        for raw_id, preds in by_protein.items():
            parsed = parse_protein_id(raw_id)
            # Sort by score descending; ties broken by GO term for determinism
            preds.sort(key=lambda x: (-x[1], x[0]))

            for rank, (go_term, score, label) in enumerate(preds, start=1):
                confidence = from_ml_score(score, self.thresholds)
                records.append(
                    AnnotationRecord(
                        uniprot_accession=parsed.accession,
                        original_id=parsed.original,
                        source_tool=SourceTool.DEEPFRI,
                        annotation_type=AnnotationType.GO_MF,
                        value=go_term,
                        label=label,
                        score=score,
                        score_type=ScoreType.CONFIDENCE,
                        confidence_level=confidence,
                        evidence_rank=rank,
                        raw_extras={},
                    )
                )

        self.log.info(
            "DeepFRI: parsed %d records from %d rows (%d below min_score=%.2f)",
            len(records), rows_processed, skipped_low_score, self.min_score,
        )
        return records
