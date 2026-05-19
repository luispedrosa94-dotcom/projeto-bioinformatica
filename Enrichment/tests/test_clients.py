"""
Tests for the Enrichment stage clients.

These tests use mocked HTTP responses so they run offline without hitting
the real APIs. Integration tests (requiring network) are marked with
@pytest.mark.integration and skipped by default.
"""
from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.schema import EnrichmentSource, EnrichmentType, GOAspectRecord


# ── UniProt parser tests (no network needed) ──────────────────────────────

class TestUniProtEntryParser:
    """Test _parse_entry with synthetic UniProt JSON payloads."""

    def _make_entry(self, acc="A0B9K2", reviewed=True):
        """Minimal synthetic UniProt entry."""
        return {
            "primaryAccession": acc,
            "entryType": "UniProtKB reviewed (Swiss-Prot)" if reviewed else "UniProtKB unreviewed (TrEMBL)",
            "proteinDescription": {
                "recommendedName": {
                    "fullName": {"value": "Archaeal ATPase subunit A"}
                }
            },
            "genes": [{"geneName": {"value": "atpA"}}],
            "uniProtKBCrossReferences": [
                {
                    "database": "GO",
                    "id": "GO:0005524",
                    "properties": [
                        {"key": "GoTerm", "value": "F:ATP binding"},
                        {"key": "GoEvidenceType", "value": "IEA"},
                    ],
                },
                {
                    "database": "GO",
                    "id": "GO:0006754",
                    "properties": [
                        {"key": "GoTerm", "value": "P:ATP biosynthetic process"},
                        {"key": "GoEvidenceType", "value": "IMP"},
                    ],
                },
                {
                    "database": "GO",
                    "id": "GO:0005886",
                    "properties": [
                        {"key": "GoTerm", "value": "C:plasma membrane"},
                        {"key": "GoEvidenceType", "value": "IDA"},
                    ],
                },
            ],
            "comments": [
                {
                    "commentType": "SUBCELLULAR LOCATION",
                    "subcellularLocations": [
                        {"location": {"value": "Cell membrane"}}
                    ],
                },
                {
                    "commentType": "FUNCTION",
                    "texts": [{"value": "Catalyzes the hydrolysis of ATP."}],
                },
            ],
            "keywords": [
                {"name": "ATP-binding"},
                {"name": "Hydrolase"},
            ],
        }

    def test_reviewed_status(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry(reviewed=True))
        status = next(r for r in records if r.enrichment_type == EnrichmentType.REVIEWED_STATUS)
        assert status.value == "reviewed"

    def test_unreviewed_status(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry(reviewed=False))
        status = next(r for r in records if r.enrichment_type == EnrichmentType.REVIEWED_STATUS)
        assert status.value == "unreviewed"

    def test_protein_name(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry())
        name = next(r for r in records if r.enrichment_type == EnrichmentType.PROTEIN_NAME)
        assert name.value == "Archaeal ATPase subunit A"

    def test_gene_name(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry())
        gene = next(r for r in records if r.enrichment_type == EnrichmentType.GENE_NAME)
        assert gene.value == "atpA"

    def test_go_aspect_mf(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry())
        mf = [r for r in records if r.enrichment_type == EnrichmentType.GO_MF]
        assert any(r.value == "GO:0005524" for r in mf)
        go_rec = next(r for r in mf if r.value == "GO:0005524")
        assert go_rec.label == "ATP binding"

    def test_go_aspect_bp(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry())
        bp = [r for r in records if r.enrichment_type == EnrichmentType.GO_BP]
        assert any(r.value == "GO:0006754" for r in bp)

    def test_go_aspect_cc(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry())
        cc = [r for r in records if r.enrichment_type == EnrichmentType.GO_CC]
        assert any(r.value == "GO:0005886" for r in cc)

    def test_subcellular_location(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry())
        loc = next(r for r in records if r.enrichment_type == EnrichmentType.SUBCELLULAR_LOCATION)
        assert loc.value == "Cell membrane"

    def test_function_description(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry())
        func = next(r for r in records if r.enrichment_type == EnrichmentType.FUNCTION_DESCRIPTION)
        assert "ATP" in func.value

    def test_keywords(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry(self._make_entry())
        kws = [r for r in records if r.enrichment_type == EnrichmentType.KEYWORD]
        kw_values = {r.value for r in kws}
        assert "ATP-binding" in kw_values
        assert "Hydrolase" in kw_values

    def test_missing_entry_graceful(self):
        """An empty entry should not crash — just produce a reviewed status record."""
        from src.clients.uniprot import _parse_entry
        records = _parse_entry({"primaryAccession": "X00000", "entryType": "UniProtKB unreviewed (TrEMBL)"})
        assert any(r.enrichment_type == EnrichmentType.REVIEWED_STATUS for r in records)

    def test_empty_accession_skipped(self):
        from src.clients.uniprot import _parse_entry
        records = _parse_entry({})
        assert records == []
