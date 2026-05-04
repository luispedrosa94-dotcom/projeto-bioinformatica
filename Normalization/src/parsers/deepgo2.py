"""
DeepGO2 parser.

DeepGO2 is a deep-learning Gene Ontology predictor. It produces three separate
TSV files, one per GO aspect:

    subset_f_preds_bp.tsv   →  Biological Process predictions
    subset_f_preds_cc.tsv   →  Cellular Component predictions
    subset_f_preds_mf.tsv   →  Molecular Function predictions

Each file is headerless with three tab-separated columns:
    <protein_id>  <GO_term>  <score>

Unlike UPIMAPI and eggNOG (where GO aspect must be resolved later via UniProt),
DeepGO2 already separates aspects by file. We use this directly: each row
becomes an annotation tagged GO_BP, GO_CC, or GO_MF based on its source file.

Confidence: ML score in [0, 1]. We use CLEAN-independent ML thresholds
(default: high>=0.7, medium>=0.3) which can be tuned per-tool in the YAML.

evidence_rank: when a protein has multiple predictions for the same aspect,
we sort by score descending and assign rank 1 to the most confident.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from ..schema import AnnotationRecord, AnnotationType, ScoreType, SourceTool
from ..utils.confidence import MLThresholds, from_ml_score
from ..utils.ids import parse_protein_id
from .base import BaseParser


# Mapping from filename suffix to canonical AnnotationType.
# Order matters: it determines the order in which files are processed.
_ASPECT_FILES = [
    ("subset_f_preds_bp.tsv", AnnotationType.GO_BP),
    ("subset_f_preds_cc.tsv", AnnotationType.GO_CC),
    ("subset_f_preds_mf.tsv", AnnotationType.GO_MF),
]


class DeepGO2Parser(BaseParser):
    tool_name = "deepgo2"

    rel_path = Path("deepgo2")  # directory containing the three TSV files

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
        dir_path = self.raw_data_root / self.rel_path
        if not dir_path.is_dir():
            self.log.warning("DeepGO2 directory not found: %s — skipping", dir_path)
            return []

        records: list[AnnotationRecord] = []
        total_rows = 0
        skipped_low_score = 0

        for filename, ann_type in _ASPECT_FILES:
            file_path = dir_path / filename
            if not file_path.exists():
                self.log.warning("DeepGO2 file missing: %s — skipping aspect", file_path)
                continue

            # Read all rows for this aspect, grouping by protein for ranking
            by_protein: dict[str, list[tuple[str, float]]] = defaultdict(list)
            file_rows = 0

            with file_path.open("r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if len(row) < 3:
                        continue
                    raw_id, go_term, score_str = row[0].strip(), row[1].strip(), row[2].strip()
                    if not raw_id or not go_term:
                        continue
                    try:
                        score = float(score_str)
                    except ValueError:
                        continue
                    file_rows += 1

                    if score < self.min_score:
                        skipped_low_score += 1
                        continue

                    by_protein[raw_id].append((go_term, score))

            total_rows += file_rows

            # For each protein, sort by score descending and assign rank
            for raw_id, preds in by_protein.items():
                parsed = parse_protein_id(raw_id)
                # Sort by score descending; ties broken by GO term for determinism
                preds.sort(key=lambda x: (-x[1], x[0]))

                for rank, (go_term, score) in enumerate(preds, start=1):
                    confidence = from_ml_score(score, self.thresholds)
                    records.append(
                        AnnotationRecord(
                            uniprot_accession=parsed.accession,
                            original_id=parsed.original,
                            source_tool=SourceTool.DEEPGO2,
                            annotation_type=ann_type,
                            value=go_term,
                            label=None,  # DeepGO2 does not give human-readable GO names
                            score=score,
                            score_type=ScoreType.CONFIDENCE,
                            confidence_level=confidence,
                            evidence_rank=rank,
                            raw_extras={},
                        )
                    )

        self.log.info(
            "DeepGO2: parsed %d records from %d rows (%d below min_score=%.2f)",
            len(records), total_rows, skipped_low_score, self.min_score,
        )
        return records
