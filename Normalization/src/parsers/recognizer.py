"""
reCOGnizer parser.

reCOGnizer runs DIAMOND against multiple CDD sub-databases (COG, KOG, Pfam,
TIGRFAM, SMART, PRK, cd, NF, and others) and merges all hits into a single
TSV: reCOGnizer_results.tsv. Each row is one alignment hit of one protein
against one database entry.

Columns of interest:
    qseqid             — protein ID (pipe-form: sp|ACC|ENTRY)
    DB ID              — hit identifier, prefix determines the database
    Protein description — human-readable label for the hit (from CDD)
    EC number          — EC annotation transferred from the hit (COG/TIGR only)
    evalue             — alignment e-value (smaller = better; NaN for 'cl' rows)
    bitscore           — alignment bit-score
    Functional category — COG/KOG functional category string
    KO                 — KEGG Orthology term (only present for COG hits)

DB ID prefix → AnnotationType mapping
--------------------------------------
    COG    → COG
    KOG    → KOG
    pfam   → PFAM
    TIGR   → TIGRFAM
    smart  → SMART
    PRK    → PRK
    cd     → CDD          (CDD specific domains, no description in practice)
    cl     → skip         (CDD superfamily clusters; NaN evalue/description)
    NF, PLN, PTZ, CHL, PHA, MTH → NCBI_CURATED

Additional records emitted per qualifying row (when field is non-empty):
    • EC number        → AnnotationType.EC
    • KO               → AnnotationType.KEGG_KO
    • Functional category → AnnotationType.COG_CATEGORY
    • Protein description → AnnotationType.PROTEIN_DESCRIPTION
      (only emitted when the description is a meaningful free-text string,
       not when it is identical to the DB ID — reCOGnizer sometimes
       copies the ID as placeholder)

top_n_per_db: for each (protein, DB prefix) pair, only the N lowest-evalue
hits are kept. This limits noise from databases like KOG where a protein
can have dozens of weak hits.

Confidence: evalue-based (HomologyThresholds), same as UPIMAPI and eggNOG.
'cl' rows have NaN evalue and are skipped entirely — they carry no annotation
payload and would all land in UNKNOWN confidence.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Iterator

import pandas as pd

from ..schema import AnnotationRecord, AnnotationType, ConfidenceLevel, ScoreType, SourceTool
from ..utils.confidence import HomologyThresholds, from_evalue
from ..utils.ids import parse_protein_id
from .base import BaseParser

log = logging.getLogger(__name__)

# ---- DB prefix → AnnotationType ----

_DB_PREFIX_MAP: dict[str, AnnotationType] = {
    "COG": AnnotationType.COG,
    "KOG": AnnotationType.KOG,
    "pfam": AnnotationType.PFAM,
    "TIGR": AnnotationType.TIGR,
    "smart": AnnotationType.SMART,
    "PRK": AnnotationType.PRK,
    "cd": AnnotationType.CDD,
    # NCBI CDD sub-databases
    "NF": AnnotationType.NCBI_CURATED,
    "PLN": AnnotationType.NCBI_CURATED,
    "PTZ": AnnotationType.NCBI_CURATED,
    "CHL": AnnotationType.NCBI_CURATED,
    "PHA": AnnotationType.NCBI_CURATED,
    "MTH": AnnotationType.NCBI_CURATED,
}

# Prefixes to skip entirely (no annotation payload).
_SKIP_PREFIXES: frozenset[str] = frozenset({"cl"})


def _safe_str(x) -> str:
    """Coerce pandas value to clean string; NaN / '-' → empty string."""
    if x is None:
        return ""
    if isinstance(x, float) and pd.isna(x):
        return ""
    s = str(x).strip()
    if s in ("-", "nan"):
        return ""
    return s


def _safe_float(x) -> float | None:
    """Coerce to float; return None on failure."""
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _db_prefix(db_id: str) -> str:
    """Extract the alphabetic prefix of a DB ID (e.g. 'COG1155' → 'COG')."""
    for i, ch in enumerate(db_id):
        if ch.isdigit():
            return db_id[:i]
    return db_id  # all-alpha (shouldn't happen in practice)


def _split_ec(raw: str) -> Iterator[str]:
    """Yield cleaned EC numbers from reCOGnizer's EC field.

    The field contains a single EC number (e.g. '3.6.3.14') or is empty.
    Partial EC numbers with dashes are valid and preserved.
    """
    s = raw.strip()
    if not s or s == "-":
        return
    yield s


class RecognizerParser(BaseParser):
    tool_name = "recognizer"

    rel_path = Path("recognizer_results") / "reCOGnizer_results.tsv"

    def __init__(
        self,
        raw_data_root: Path,
        thresholds: HomologyThresholds | None = None,
        top_n_per_db: int = 5,
    ):
        super().__init__(raw_data_root)
        self.thresholds = thresholds or HomologyThresholds()
        # top_n_per_db: keep only the N best (lowest evalue) hits per
        # (protein, DB prefix) to avoid noise from multi-hit databases.
        self.top_n_per_db = int(top_n_per_db)

    def parse(self) -> list[AnnotationRecord]:
        path = self.raw_data_root / self.rel_path
        if not path.exists():
            self.log.warning("reCOGnizer file not found: %s — skipping", path)
            return []

        df = pd.read_csv(path, sep="\t", low_memory=False)

        # Normalise column names (strip leading/trailing whitespace)
        df.columns = [c.strip() for c in df.columns]

        rows_total = len(df)
        rows_skipped_cl = 0
        rows_skipped_no_type = 0

        records: list[AnnotationRecord] = []

        # Group by (protein, db_prefix) for top-N filtering.
        # We collect rows into a dict keyed by (raw_id, prefix), then
        # sort by evalue and keep top_n_per_db before emitting records.
        groups: dict[tuple[str, str], list[pd.Series]] = defaultdict(list)

        for _, row in df.iterrows():
            raw_id = _safe_str(row.get("qseqid", ""))
            db_id = _safe_str(row.get("DB ID", ""))

            if not raw_id or not db_id:
                continue

            prefix = _db_prefix(db_id)

            if prefix in _SKIP_PREFIXES:
                rows_skipped_cl += 1
                continue

            if prefix not in _DB_PREFIX_MAP:
                rows_skipped_no_type += 1
                continue

            groups[(raw_id, prefix)].append(row)

        # For each (protein, prefix) group: sort by evalue asc, keep top_n
        for (raw_id, prefix), rows in groups.items():
            # Sort: rows with evalue=0 first (perfect hits), then ascending.
            # NaN evalue (shouldn't happen after cl-skip, but guard anyway).
            def _evalue_key(r):
                e = _safe_float(r.get("evalue"))
                return (1, 1e300) if e is None else (0, e)

            rows_sorted = sorted(rows, key=_evalue_key)
            rows_top = rows_sorted[: self.top_n_per_db]

            parsed = parse_protein_id(raw_id)
            annotation_type = _DB_PREFIX_MAP[prefix]

            for rank, row in enumerate(rows_top, start=1):
                db_id = _safe_str(row.get("DB ID", ""))
                description = _safe_str(row.get("Protein description", ""))
                ec_raw = _safe_str(row.get("EC number", ""))
                ko_raw = _safe_str(row.get("KO", ""))
                func_cat = _safe_str(row.get("Functional category", ""))
                evalue = _safe_float(row.get("evalue"))
                bitscore = _safe_float(row.get("bitscore"))

                confidence = from_evalue(evalue, self.thresholds)

                # Primary record: the DB domain / family hit
                records.append(
                    AnnotationRecord(
                        uniprot_accession=parsed.accession,
                        original_id=parsed.original,
                        source_tool=SourceTool.RECOGNIZER,
                        annotation_type=annotation_type,
                        value=db_id,
                        label=description or None,
                        score=evalue,
                        score_type=ScoreType.EVALUE,
                        confidence_level=confidence,
                        evidence_rank=rank,
                        raw_extras={"bitscore": bitscore},
                    )
                )

                # Secondary records — additional payloads carried by the row.
                # These inherit the same evalue/confidence as the primary hit.

                # EC number (COG, TIGR rows)
                for ec in _split_ec(ec_raw):
                    records.append(
                        AnnotationRecord(
                            uniprot_accession=parsed.accession,
                            original_id=parsed.original,
                            source_tool=SourceTool.RECOGNIZER,
                            annotation_type=AnnotationType.EC,
                            value=ec,
                            label=None,
                            score=evalue,
                            score_type=ScoreType.EVALUE,
                            confidence_level=confidence,
                            evidence_rank=rank,
                            raw_extras={"source_db": db_id},
                        )
                    )

                # KEGG Orthology (COG rows only)
                if ko_raw:
                    # Strip 'ko:' prefix if present
                    ko_val = ko_raw.removeprefix("ko:").strip()
                    if ko_val:
                        records.append(
                            AnnotationRecord(
                                uniprot_accession=parsed.accession,
                                original_id=parsed.original,
                                source_tool=SourceTool.RECOGNIZER,
                                annotation_type=AnnotationType.KEGG_KO,
                                value=ko_val,
                                label=None,
                                score=evalue,
                                score_type=ScoreType.EVALUE,
                                confidence_level=confidence,
                                evidence_rank=rank,
                                raw_extras={"source_db": db_id},
                            )
                        )

                # COG functional category
                if func_cat:
                    records.append(
                        AnnotationRecord(
                            uniprot_accession=parsed.accession,
                            original_id=parsed.original,
                            source_tool=SourceTool.RECOGNIZER,
                            annotation_type=AnnotationType.COG_CATEGORY,
                            value=func_cat.strip(),
                            label=None,
                            score=evalue,
                            score_type=ScoreType.EVALUE,
                            confidence_level=confidence,
                            evidence_rank=rank,
                            raw_extras={"source_db": db_id},
                        )
                    )

                # Protein description — emit only when non-trivial
                # (not empty, not identical to the DB ID)
                if description and description != db_id:
                    records.append(
                        AnnotationRecord(
                            uniprot_accession=parsed.accession,
                            original_id=parsed.original,
                            source_tool=SourceTool.RECOGNIZER,
                            annotation_type=AnnotationType.PROTEIN_DESCRIPTION,
                            value=description,
                            label=None,
                            score=evalue,
                            score_type=ScoreType.EVALUE,
                            confidence_level=confidence,
                            evidence_rank=rank,
                            raw_extras={"source_db": db_id},
                        )
                    )

        self.log.info(
            "reCOGnizer: parsed %d records from %d rows "
            "(%d cl-rows skipped, %d unknown-prefix skipped)",
            len(records),
            rows_total,
            rows_skipped_cl,
            rows_skipped_no_type,
        )
        return records
