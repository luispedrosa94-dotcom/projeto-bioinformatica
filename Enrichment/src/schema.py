"""
Canonical schema for the Enrichment stage.

Each record represents one piece of information retrieved from an external
API (UniProt or STRING) for a given protein. The format mirrors the
AnnotationRecord from Stage 1 so both can be merged downstream.
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field


class EnrichmentSource(str, Enum):
    UNIPROT = "uniprot"
    STRING = "string"


class EnrichmentType(str, Enum):
    # UniProt fields
    REVIEWED_STATUS = "reviewed_status"         # 'reviewed' or 'unreviewed'
    PROTEIN_NAME = "protein_name"               # recommended full name
    GENE_NAME = "gene_name"                     # primary gene name
    GO_BP = "GO_BP"
    GO_MF = "GO_MF"
    GO_CC = "GO_CC"
    EC = "EC"
    SUBCELLULAR_LOCATION = "subcellular_location"
    KEYWORD = "keyword"
    FUNCTION_DESCRIPTION = "function_description"

    # GO aspect resolution (GO_unknown → GO_BP/MF/CC)
    GO_ASPECT = "go_aspect"

    # STRING fields
    INTERACTION_PARTNER = "interaction_partner"
    INTERACTION_SCORE = "interaction_score"
    FUNCTIONAL_ENRICHMENT = "functional_enrichment"


class EnrichmentRecord(BaseModel):
    """One piece of API-retrieved information for one protein."""

    uniprot_accession: str = Field(..., description="Canonical UniProt accession")
    source: EnrichmentSource
    enrichment_type: EnrichmentType
    value: str = Field(..., description="The enrichment value")
    label: Optional[str] = None
    score: Optional[float] = None              # e.g. STRING combined score
    extras: dict[str, Any] = Field(default_factory=dict)


class GOAspectRecord(BaseModel):
    """Maps a GO term to its aspect (BP / MF / CC)."""
    go_term: str
    aspect: str    # 'biological_process' | 'molecular_function' | 'cellular_component'
    aspect_short: str  # 'BP' | 'MF' | 'CC'
