"""
CLEAN parser.

CLEAN is a deep-learning enzyme function predictor: given a protein sequence,
it predicts EC numbers with associated confidence scores. The output format is
a headerless CSV with a variable number of columns per row:

    tr|A0B5W5|A0B5W5_METTP,EC:1.2.7.4/0.9367
    tr|A0B6R9|A0B6R9_METTP,EC:2.7.13.3/0.0005,EC:4.2.3.55/0.0004

Each row starts with the protein ID, followed by 1+ predictions in the form
'EC:<number>/<score>'. We emit one AnnotationRecord per prediction.

Confidence handling (per supervisor meeting and Decision 1 of the design):
  - We KEEP all predictions, including score=0.0000.
  - The categorical confidence_level reflects the score using CLEAN-specific
    thresholds (high>=0.5, medium>=0.1) which are more lenient than the
    generic ML thresholds because CLEAN scores tend to be lower overall.
  - Filtering of low-confidence predictions is deferred to later stages.

evidence_rank preserves the order in which CLEAN reports predictions for a
given protein (1 = most confident according to CLEAN's own ordering).
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator

from ..schema import AnnotationRecord, AnnotationType, ScoreType, SourceTool
from ..utils.confidence import MLThresholds, DEFAULT_CLEAN, from_ml_score
from ..utils.ids import parse_protein_id
from .base import BaseParser


def _split_ec_field(field: str) -> tuple[str | None, float | None]:
    """Parse one CLEAN prediction field 'EC:<number>/<score>'.

    Returns (ec_number, score) where ec_number has the 'EC:' prefix stripped,
    or (None, None) if the field is malformed.

    Examples:
        'EC:1.2.7.4/0.9367'  -> ('1.2.7.4', 0.9367)
        'EC:3.6.1.-/0.5'     -> ('3.6.1.-', 0.5)
        ''                   -> (None, None)
        'malformed'          -> (None, None)
    """
    field = field.strip()
    if not field:
        return None, None

    # Strip optional 'EC:' prefix
    if field.startswith("EC:"):
        field = field[3:]

    # Split on the LAST '/' to be robust if EC numbers ever contain '/'
    # (they shouldn't, but cheap safety).
    if "/" not in field:
        return None, None

    ec_part, score_part = field.rsplit("/", 1)
    ec_part = ec_part.strip()
    score_part = score_part.strip()

    if not ec_part:
        return None, None

    try:
        score = float(score_part)
    except ValueError:
        return None, None

    return ec_part, score


class CleanParser(BaseParser):
    tool_name = "clean"

    rel_path = Path("clean") / "subset_f_maxsep.csv"

    def __init__(
        self,
        raw_data_root: Path,
        thresholds: MLThresholds | None = None,
    ):
        super().__init__(raw_data_root)
        self.thresholds = thresholds or DEFAULT_CLEAN

    def parse(self) -> list[AnnotationRecord]:
        path = self.raw_data_root / self.rel_path
        if not path.exists():
            self.log.warning("CLEAN file not found: %s — skipping", path)
            return []

        records: list[AnnotationRecord] = []
        rows_processed = 0
        skipped_malformed = 0

        # Use the csv module directly because the file has variable columns
        # per row (pandas would struggle with that without extra work).
        with path.open("r", newline="", encoding="utf-8") as f:
            reader = csv.reader(f)
            for row in reader:
                if not row:
                    continue
                rows_processed += 1

                raw_id = row[0].strip()
                if not raw_id:
                    continue
                parsed = parse_protein_id(raw_id)

                # Each subsequent column is one EC prediction
                for rank, field in enumerate(row[1:], start=1):
                    ec, score = _split_ec_field(field)
                    if ec is None:
                        skipped_malformed += 1
                        continue

                    confidence = from_ml_score(score, self.thresholds)

                    records.append(
                        AnnotationRecord(
                            uniprot_accession=parsed.accession,
                            original_id=parsed.original,
                            source_tool=SourceTool.CLEAN,
                            annotation_type=AnnotationType.EC,
                            value=ec,
                            label=None,  # CLEAN does not provide a human-readable label
                            score=score,
                            score_type=ScoreType.CONFIDENCE,
                            confidence_level=confidence,
                            evidence_rank=rank,
                            raw_extras={},
                        )
                    )

        self.log.info(
            "CLEAN: parsed %d EC predictions from %d rows (%d malformed fields skipped)",
            len(records), rows_processed, skipped_malformed,
        )
        return records
