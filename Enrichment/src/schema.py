"""
Canonical schema for the Enrichment stage.

Each record represents one piece of information retrieved from an external
API (UniProt, InterPro) for a given protein. The format mirrors the
AnnotationRecord from Stage 1 so both can be merged downstream.
"""
from __future__ import annotations

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field


class EnrichmentSource(str, Enum):
    UNIPROT = "uniprot"
    INTERPRO = "interpro"


class EnrichmentType(str, Enum):
    # UniProt fields — identity
    REVIEWED_STATUS = "reviewed_status"       # 'reviewed' | 'unreviewed'
    ANNOTATION_SCORE = "annotation_score"     # 1.0–5.0 quality score
    PROTEIN_EXISTENCE = "protein_existence"   # evidence level for protein existence
    PROTEIN_NAME = "protein_name"             # recommended/submitted full name
    ALTERNATIVE_NAME = "alternative_name"     # alternative protein names
    GENE_NAME = "gene_name"                   # primary gene name
    GENE_NAME_SYNONYM = "gene_name_synonym"   # synonymous gene names
    ORGANISM = "organism"                     # scientific name + taxon + lineage
    ENTRY_NAME = "entry_name"                 # UniProtKB entry name (e.g. ATPA_HUMAN)

    # UniProt fields — function
    FUNCTION_DESCRIPTION = "function_description"  # free-text functional description
    CATALYTIC_ACTIVITY = "catalytic_activity"       # EC + reaction + Rhea/ChEBI IDs
    PATHWAY = "pathway"                             # metabolic pathway description
    SUBUNIT = "subunit"                             # quaternary structure
    SIMILARITY = "similarity"                       # protein family membership
    GO_BP = "GO_BP"
    GO_MF = "GO_MF"
    GO_CC = "GO_CC"
    EC = "EC"                                       # EC number from cross-references

    # UniProt fields — localisation & structure
    SUBCELLULAR_LOCATION = "subcellular_location"
    FEATURE_DOMAIN = "feature_domain"               # domain with sequence positions
    FEATURE_ACTIVE_SITE = "feature_active_site"     # active site with position
    FEATURE_BINDING_SITE = "feature_binding_site"   # binding site with ligand + position

    # UniProt fields — metadata
    KEYWORD = "keyword"
    SEQUENCE = "sequence"                           # amino acid sequence + length + mass
    REFERENCE = "reference"                         # literature reference / PMID

    # GO aspect resolution (GO_unknown → GO_BP/MF/CC)
    GO_ASPECT = "go_aspect"

    # InterPro fields
    INTERPRO_ENTRY = "interpro_entry"                            # integrated InterPro entry (IPRxxxxxx)
    INTERPRO_UNINTEGRATED_SIGNATURE = "interpro_unintegrated_signature"  # signature not yet integrated




class EnrichmentRecord(BaseModel):
    """One piece of API-retrieved information for one protein."""

    uniprot_accession: str = Field(..., description="Canonical UniProt accession")
    source: EnrichmentSource
    enrichment_type: EnrichmentType
    value: str = Field(..., description="The enrichment value")
    label: Optional[str] = None
    score: Optional[float] = None              # e.g. ML prediction score
    extras: dict[str, Any] = Field(default_factory=dict)


class GOAspectRecord(BaseModel):
    """Maps a GO term to its aspect (BP / MF / CC)."""
    go_term: str
    aspect: str    # 'biological_process' | 'molecular_function' | 'cellular_component'
    aspect_short: str  # 'BP' | 'MF' | 'CC'