"""
UPIMAPI parser.

UPIMAPI runs DIAMOND against UniProt and joins each best hit with UniProt
metadata, producing one row per query protein. From a single row we extract
several distinct annotations:

    Protein names                 → PROTEIN_DESCRIPTION (one record)
    Function [CC]                 → FUNCTION_CC          (one record)
    Gene Ontology (GO)            → GO_UNKNOWN           (one per GO term)
    EC number                     → EC                   (one per EC code)
    Pfam                          → PFAM                 (one per Pfam family)
    Protein families              → PROTEIN_FAMILY       (one record)

Note on confidence: UPIMAPI annotations are inherited from a homologous
UniProt entry, so the BLAST e-value/bitscore of the match describes the
TRANSFER confidence — not the annotation's intrinsic confidence. We attach
the e-value as the score for every annotation derived from the row.

GO aspect (BP / MF / CC) is NOT directly available in UPIMAPI output —
the field mixes all three. We tag them as GO_UNKNOWN here and resolve the
aspect later during the UniProt API enrichment step.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Iterator

import pandas as pd

from ..schema import AnnotationRecord, AnnotationType, ScoreType, SourceTool
from ..utils.confidence import HomologyThresholds, from_evalue
from ..utils.ids import parse_protein_id
from .base import BaseParser


# ---- Field-level parsers (each handles one UPIMAPI column) ----

_GO_TERM_RE = re.compile(r"\[GO:(\d{7})\]")
# EC pattern. We use lookarounds instead of \b because \b doesn't work
# correctly around the dash (-) used for unspecified EC sublevels.
_EC_RE = re.compile(
    r"(?<![\w.])(\d+\.\d+\.\d+\.\d+|\d+\.\d+\.\d+\.-|\d+\.\d+\.-\.-|\d+\.-\.-\.-)(?![\w.])"
)


def _split_go_terms(go_field: str) -> Iterator[tuple[str, str]]:
    """UPIMAPI GO format example:
        'cell outer membrane [GO:0009279]; iron ion transport [GO:0006826]'
    Yields (go_id_full, label) pairs.
    """
    if not go_field or not isinstance(go_field, str):
        return
    for chunk in go_field.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        m = _GO_TERM_RE.search(chunk)
        if not m:
            continue
        go_id = f"GO:{m.group(1)}"
        # Label is everything before the [GO:...] bracket
        label = _GO_TERM_RE.sub("", chunk).strip().rstrip(",;").strip()
        yield go_id, label or None  # type: ignore[misc]


def _split_ec_numbers(ec_field: str) -> Iterator[str]:
    """UPIMAPI EC field can be empty, a single EC, or multiple separated by ';' or ','."""
    if not ec_field or not isinstance(ec_field, str):
        return
    seen: set[str] = set()
    for m in _EC_RE.finditer(ec_field):
        ec = m.group(1)
        if ec not in seen:
            seen.add(ec)
            yield ec


def _split_pfam(pfam_field: str) -> Iterator[str]:
    """Pfam field example: 'PF14322;PF07980;'"""
    if not pfam_field or not isinstance(pfam_field, str):
        return
    for chunk in pfam_field.split(";"):
        chunk = chunk.strip()
        if chunk and chunk.startswith("PF"):
            yield chunk


def _safe_float(x) -> float | None:
    """Coerce to float, returning None for blanks / non-numeric / NaN."""
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    if pd.isna(f):
        return None
    return f


def _safe_str(x) -> str:
    """Coerce to a clean string, treating NaN and None as empty."""
    if x is None:
        return ""
    if isinstance(x, float) and pd.isna(x):
        return ""
    return str(x).strip()


# ---- Parser ----

class UpimapiParser(BaseParser):
    tool_name = "upimapi"

    # Path relative to raw_data_root
    rel_path = Path("upimapi") / "UPIMAPI_results.tsv"

    EXPECTED_COLS = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "Entry", "Entry Name", "Protein names", "Gene Ontology (GO)",
        "Function [CC]", "EC number", "Annotation", "Protein families", "Pfam",
    ]

    def __init__(
        self,
        raw_data_root: Path,
        thresholds: HomologyThresholds | None = None,
    ):
        super().__init__(raw_data_root)
        self.thresholds = thresholds or HomologyThresholds()

    def parse(self) -> list[AnnotationRecord]:
        path = self.raw_data_root / self.rel_path
        if not path.exists():
            self.log.warning("UPIMAPI file not found: %s — skipping", path)
            return []

        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

        # Sanity: schema check
        missing = [c for c in self.EXPECTED_COLS if c not in df.columns]
        if missing:
            self.log.warning("UPIMAPI missing expected columns: %s", missing)

        records: list[AnnotationRecord] = []

        for row_idx, row in df.iterrows():
            raw_id = _safe_str(row.get("qseqid"))
            if not raw_id:
                continue
            parsed = parse_protein_id(raw_id)

            evalue = _safe_float(row.get("evalue"))
            bitscore = _safe_float(row.get("bitscore"))
            pident = _safe_float(row.get("pident"))

            # Confidence level shared by every annotation derived from this
            # row, since they all inherit from the same homology hit.
            row_confidence = from_evalue(evalue, self.thresholds)

            # Carry homology-transfer metadata in raw_extras for every record
            transfer_meta = {
                "subject_entry": _safe_str(row.get("Entry")),
                "subject_entry_name": _safe_str(row.get("Entry Name")),
                "pident": pident,
                "bitscore": bitscore,
                "evalue": evalue,
            }

            def _emit(ann_type: AnnotationType, value: str, label: str | None = None) -> None:
                if not value:
                    return
                records.append(
                    AnnotationRecord(
                        uniprot_accession=parsed.accession,
                        original_id=parsed.original,
                        source_tool=SourceTool.UPIMAPI,
                        annotation_type=ann_type,
                        value=value,
                        label=label,
                        score=evalue,
                        score_type=ScoreType.EVALUE if evalue is not None else ScoreType.NONE,
                        confidence_level=row_confidence,
                        evidence_rank=1,  # UPIMAPI gives one best hit per protein
                        raw_extras=transfer_meta.copy(),
                    )
                )

            # 1. Protein names → description
            protein_names = _safe_str(row.get("Protein names"))
            if protein_names:
                _emit(AnnotationType.PROTEIN_DESCRIPTION, protein_names)

            # 2. Function [CC] (UniProt curated function note)
            function_cc = _safe_str(row.get("Function [CC]"))
            if function_cc:
                _emit(AnnotationType.FUNCTION_CC, function_cc)

            # 3. Protein families
            families = _safe_str(row.get("Protein families"))
            if families:
                _emit(AnnotationType.PROTEIN_FAMILY, families)

            # 4. GO terms (aspect unknown)
            go_field = _safe_str(row.get("Gene Ontology (GO)"))
            for go_id, go_label in _split_go_terms(go_field):
                _emit(AnnotationType.GO_UNKNOWN, go_id, label=go_label)

            # 5. EC numbers
            ec_field = _safe_str(row.get("EC number"))
            for ec in _split_ec_numbers(ec_field):
                _emit(AnnotationType.EC, ec)

            # 6. Pfam IDs
            pfam_field = _safe_str(row.get("Pfam"))
            for pfam in _split_pfam(pfam_field):
                _emit(AnnotationType.PFAM, pfam)

        self.log.info("UPIMAPI: parsed %d annotation records from %d rows", len(records), len(df))
        return records
