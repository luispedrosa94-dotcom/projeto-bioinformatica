"""
eggNOG-mapper parser.

eggNOG-mapper assigns each query protein to an orthologous group, then
transfers functional annotations from the group's consensus. From a single
row we extract many canonical records, each tagged with the same e-value
based confidence (since they all inherit from the same orthology hit).

Fields retained (13 of 21 — see Decision 1 in design notes):
  COG_category, Description, Preferred_name, GOs, EC,
  KEGG_ko, KEGG_Pathway, KEGG_Module, KEGG_Reaction, BRITE,
  CAZy, PFAMs, eggNOG_OGs

Fields dropped (low coverage or redundant):
  seed_ortholog, max_annot_lvl, KEGG_rclass, KEGG_TC, BiGG_Reaction

Like UPIMAPI, GO terms come without aspect (BP/MF/CC) so we tag them as
GO_unknown and resolve aspect later via the UniProt API (Stage 2 — Enrichment).
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


# ---- Field-level helpers ----

# eggNOG_OGs format: "COG1152@1|root,arCOG02428@2157|Archaea,2XTBY@28890|Euryarchaeota"
# We extract just the OG IDs (before @), discarding the taxonomic context (Decision 2).
_OG_RE = re.compile(r"([^,@]+)@\d+\|[^,]+")

# EC pattern (same as UPIMAPI parser; reused conceptually but redefined locally
# to keep parsers independent)
_EC_RE = re.compile(
    r"(?<![\w.])(\d+\.\d+\.\d+\.\d+|\d+\.\d+\.\d+\.-|\d+\.\d+\.-\.-|\d+\.-\.-\.-)(?![\w.])"
)


def _safe_str(x) -> str:
    """Coerce to a clean string. eggNOG uses '-' for missing values, treat as empty."""
    if x is None:
        return ""
    if isinstance(x, float) and pd.isna(x):
        return ""
    s = str(x).strip()
    if s == "-":
        return ""
    return s


def _safe_float(x) -> float | None:
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    if pd.isna(f):
        return None
    return f


def _split_csv(field: str) -> Iterator[str]:
    """Split a comma-separated field, dropping empties and stripping whitespace."""
    if not field:
        return
    seen: set[str] = set()
    for chunk in field.split(","):
        c = chunk.strip()
        if c and c != "-" and c not in seen:
            seen.add(c)
            yield c


def _split_eggnog_ogs(field: str) -> Iterator[str]:
    """Extract OG IDs (e.g. 'COG1152', 'arCOG02428') from the OG_ID@taxlevel|name format."""
    if not field:
        return
    seen: set[str] = set()
    for m in _OG_RE.finditer(field):
        og = m.group(1).strip()
        if og and og not in seen:
            seen.add(og)
            yield og


def _split_kegg_pathway(field: str) -> Iterator[str]:
    """KEGG_Pathway field can contain both 'ko*****' and 'map*****' for the same
    pathway (one is the reference, one is the map). We keep both — they are
    semantically distinct identifiers — but deduplicate exact repeats.
    """
    yield from _split_csv(field)


def _split_brite(field: str) -> Iterator[str]:
    """BRITE: comma-separated KEGG hierarchy IDs (ko00000, ko00001, ...)."""
    yield from _split_csv(field)


def _split_pfams(field: str) -> Iterator[str]:
    """eggNOG PFAMs come as names ('Fer4_7,Fer4_9,Prismane'), not IDs.
    Different from UPIMAPI which gives IDs (PF14322). We keep names as values
    and let the UniProt enrichment stage map them to IDs if needed.
    """
    yield from _split_csv(field)


def _split_ec(field: str) -> Iterator[str]:
    """eggNOG EC field can be a single EC or several separated by commas."""
    if not field:
        return
    seen: set[str] = set()
    for m in _EC_RE.finditer(field):
        ec = m.group(1)
        if ec not in seen:
            seen.add(ec)
            yield ec


# ---- Parser ----

class EggnogParser(BaseParser):
    tool_name = "eggnog"

    rel_path = Path("eggnogmapper_results") / "eggnog_mapper_results.emapper.annotations"

    # Columns we care about (must match the eggNOG-mapper v2 header exactly)
    EXPECTED_COLS = [
        "#query", "seed_ortholog", "evalue", "score",
        "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
        "Preferred_name", "GOs", "EC",
        "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
        "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
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
            self.log.warning("eggNOG file not found: %s — skipping", path)
            return []

        # eggNOG annotations file has metadata lines starting with '##' followed
        # by the header line starting with '#query'. We use comment='#' filter
        # but then we'd lose the header — so we read manually.
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            comment="#",       # skip comment lines
            header=None,
            names=self.EXPECTED_COLS,
        )

        # The file always ends with a few comment lines (## summary stats) which
        # comment='#' already drops. Sanity-check that we have data.
        if df.empty:
            self.log.warning("eggNOG file produced 0 data rows — check format")
            return []

        records: list[AnnotationRecord] = []

        for _, row in df.iterrows():
            raw_id = _safe_str(row.get("#query"))
            if not raw_id:
                continue
            parsed = parse_protein_id(raw_id)

            evalue = _safe_float(row.get("evalue"))
            bitscore = _safe_float(row.get("score"))   # eggNOG calls bitscore 'score'

            row_confidence = from_evalue(evalue, self.thresholds)

            # Shared metadata for every record from this row
            transfer_meta = {
                "seed_ortholog": _safe_str(row.get("seed_ortholog")),
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
                        source_tool=SourceTool.EGGNOG,
                        annotation_type=ann_type,
                        value=value,
                        label=label,
                        score=evalue,
                        score_type=ScoreType.EVALUE if evalue is not None else ScoreType.NONE,
                        confidence_level=row_confidence,
                        evidence_rank=1,
                        raw_extras=transfer_meta.copy(),
                    )
                )

            # --- Free-text identity fields ---
            description = _safe_str(row.get("Description"))
            if description:
                _emit(AnnotationType.PROTEIN_DESCRIPTION, description)

            preferred_name = _safe_str(row.get("Preferred_name"))
            if preferred_name:
                _emit(AnnotationType.PREFERRED_GENE_NAME, preferred_name)

            # --- COG category (single letter, e.g. 'C' for energy production) ---
            cog_cat = _safe_str(row.get("COG_category"))
            if cog_cat:
                # Some rows have multiple letters (e.g. 'CR') for proteins that
                # span multiple categories. Emit each separately.
                for letter in cog_cat:
                    if letter.isalpha():
                        _emit(AnnotationType.COG_CATEGORY, letter)

            # --- Orthologous groups (Decision 2: just the IDs) ---
            ogs = _safe_str(row.get("eggNOG_OGs"))
            for og in _split_eggnog_ogs(ogs):
                _emit(AnnotationType.EGGNOG_OG, og)

            # --- GO terms (aspect unknown, Decision 5) ---
            gos = _safe_str(row.get("GOs"))
            for go_id in _split_csv(gos):
                if go_id.startswith("GO:"):
                    _emit(AnnotationType.GO_UNKNOWN, go_id)

            # --- EC numbers ---
            ec_field = _safe_str(row.get("EC"))
            for ec in _split_ec(ec_field):
                _emit(AnnotationType.EC, ec)

            # --- KEGG identifiers (4 fields) ---
            for ko in _split_csv(_safe_str(row.get("KEGG_ko"))):
                _emit(AnnotationType.KEGG_KO, ko)
            for path_id in _split_kegg_pathway(_safe_str(row.get("KEGG_Pathway"))):
                _emit(AnnotationType.KEGG_PATHWAY, path_id)
            for mod in _split_csv(_safe_str(row.get("KEGG_Module"))):
                _emit(AnnotationType.KEGG_MODULE, mod)
            for rxn in _split_csv(_safe_str(row.get("KEGG_Reaction"))):
                _emit(AnnotationType.KEGG_REACTION, rxn)

            # --- BRITE classifications ---
            for brite in _split_brite(_safe_str(row.get("BRITE"))):
                _emit(AnnotationType.BRITE, brite)

            # --- CAZy families (rare but valuable when present) ---
            for cazy in _split_csv(_safe_str(row.get("CAZy"))):
                _emit(AnnotationType.CAZY, cazy)

            # --- Pfam domains (names, not IDs — see _split_pfams docstring) ---
            for pfam in _split_pfams(_safe_str(row.get("PFAMs"))):
                _emit(AnnotationType.PFAM, pfam)

        self.log.info(
            "eggNOG: parsed %d annotation records from %d rows",
            len(records), len(df),
        )
        return records
