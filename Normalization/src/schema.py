"""
Canonical schema for normalized protein annotations.

Every parser in src/parsers/ must produce records conforming to AnnotationRecord.
This is the single source of truth for the normalization stage of the pipeline.
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field, field_validator


# ---- Controlled vocabularies ----
# These enums constrain what each parser can emit. Adding a new tool or
# annotation type requires extending these enums explicitly — this prevents
# silent schema drift.

class SourceTool(str, Enum):
    UPIMAPI = "upimapi"
    RECOGNIZER = "recognizer"
    EGGNOG = "eggnog"
    DEEPFRI = "deepfri"
    DEEPGO2 = "deepgo2"
    CLEAN = "clean"
    FOLDSEEK = "foldseek"
    COLABFOLD = "colabfold"


class AnnotationType(str, Enum):
    # Gene Ontology (aspect known)
    GO_BP = "GO_BP"
    GO_MF = "GO_MF"
    GO_CC = "GO_CC"
    # Gene Ontology (aspect not known at parse time — resolved later)
    GO_UNKNOWN = "GO_unknown"

    # Functional identifiers
    EC = "EC"
    KEGG_KO = "KEGG_ko"
    KEGG_PATHWAY = "KEGG_pathway"
    KEGG_MODULE = "KEGG_module"
    KEGG_REACTION = "KEGG_reaction"
    BRITE = "brite"

    # Domains / families
    PFAM = "pfam"
    SMART = "smart"
    TIGR = "tigrfam"
    COG = "cog"
    KOG = "kog"
    NCBI_CURATED = "ncbi_curated"
    PRK = "prk"
    CDD = "cdd"
    CAZY = "cazy"

    # Categories / classifications
    COG_CATEGORY = "cog_category"
    EGGNOG_OG = "eggnog_og"

    # Free-text descriptions
    PROTEIN_DESCRIPTION = "protein_description"
    PROTEIN_FAMILY = "protein_family"
    FUNCTION_CC = "function_cc"
    PREFERRED_GENE_NAME = "preferred_gene_name"

    # Structure-related
    STRUCTURAL_HIT = "structural_hit"
    PLDDT_MEAN = "plddt_mean"
    PLDDT_MAX = "plddt_max"
    PLDDT_MIN = "plddt_min"
    PROTEIN_LENGTH = "protein_length"


class ScoreType(str, Enum):
    EVALUE = "evalue"
    BITSCORE = "bitscore"
    PIDENT = "pident"           # percent identity
    CONFIDENCE = "confidence"    # 0..1 ML confidence
    PLDDT = "plddt"              # AlphaFold per-residue confidence
    NONE = "none"


class ConfidenceLevel(str, Enum):
    """Categorical confidence abstraction over heterogeneous score scales.

    Different tools express confidence on incompatible scales (small e-value =
    good, large ML score = good — opposite directions). The categorical level
    gives downstream stages — particularly the LLM harmonizer — a uniform
    ordinal signal without losing the raw `score` for auditing.

    Per-tool mapping rules live in src/utils/confidence.py.
    """
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    UNKNOWN = "unknown"


# ---- Records ----

class AnnotationRecord(BaseModel):
    """One piece of evidence about one protein from one tool.

    Long-format design choice: each annotation (one GO term, one EC number,
    one Pfam domain) is a separate row. This makes downstream merging,
    enrichment, and LLM context assembly trivial.
    """
    # Identity
    uniprot_accession: str = Field(..., description="Canonical UniProt accession, e.g. A0B9K2")
    original_id: str = Field(..., description="ID as it appeared in the source file")
    source_tool: SourceTool

    # Annotation payload
    annotation_type: AnnotationType
    value: str = Field(..., description="The annotation itself (term, ID, description, ...)")

    # Optional human-readable label (e.g. GO term name, EC enzyme name)
    label: Optional[str] = None

    # Confidence
    score: Optional[float] = None
    score_type: ScoreType = ScoreType.NONE
    confidence_level: ConfidenceLevel = ConfidenceLevel.UNKNOWN

    # When a tool emits multiple ranked hits per protein, this preserves order
    evidence_rank: Optional[int] = None

    # Anything tool-specific that doesn't fit the canonical fields
    raw_extras: dict[str, Any] = Field(default_factory=dict)

    @field_validator("uniprot_accession")
    @classmethod
    def _validate_accession(cls, v: str) -> str:
        # UniProt accessions: 6 or 10 alphanumeric chars, starting with [O,P,Q] or [A-N,R-Z]
        # We allow lowercase too and uppercase it.
        v = v.strip().upper()
        if not v or len(v) not in (6, 10):
            # Don't hard-fail — emit but warn via downstream QC.
            # Some tools may have non-standard IDs we want to preserve.
            pass
        return v


class ProteinRecord(BaseModel):
    """One row per unique protein in the working dataset."""
    uniprot_accession: str
    original_id: str
    db_source: str               # 'sp' (SwissProt) or 'tr' (TrEMBL) — or 'unknown'
    entry_name: Optional[str] = None  # e.g. AATA_METTP
    in_poorly_annotated_subset: bool = False
