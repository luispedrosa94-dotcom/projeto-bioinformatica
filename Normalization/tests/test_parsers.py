"""Tests covering ID parsing and the UPIMAPI parser against tiny real-data samples."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pytest

from src.utils.ids import parse_protein_id
from src.utils.confidence import (
    HomologyThresholds,
    MLThresholds,
    DEFAULT_CLEAN,
    from_evalue,
    from_ml_score,
)
from src.schema import ConfidenceLevel
from src.parsers.upimapi import _split_go_terms, _split_ec_numbers, _split_pfam


class TestIDParser:
    def test_pipe_form_swissprot(self):
        p = parse_protein_id("sp|A0B9K2|AATA_METTP")
        assert p.accession == "A0B9K2"
        assert p.db_source == "sp"
        assert p.entry_name == "AATA_METTP"

    def test_pipe_form_trembl(self):
        p = parse_protein_id("tr|X5DZ82|X5DZ82_9BACT")
        assert p.accession == "X5DZ82"
        assert p.db_source == "tr"
        assert p.entry_name == "X5DZ82_9BACT"

    def test_underscore_form_simple(self):
        p = parse_protein_id("sp_A0B9K2_AATA_METTP")
        assert p.accession == "A0B9K2"
        assert p.db_source == "sp"

    def test_underscore_form_with_alphafold_suffix(self):
        p = parse_protein_id(
            "sp_A0B9K2_AATA_METTP_unrelaxed_rank_001_alphafold2_model_1_seed_000"
        )
        assert p.accession == "A0B9K2"

    def test_bare_accession(self):
        p = parse_protein_id("A0B9K2")
        assert p.accession == "A0B9K2"

    def test_long_accession(self):
        p = parse_protein_id("sp|A0A5N4D8M3|TRYP_PIG")
        assert p.accession == "A0A5N4D8M3"


class TestUpimapiFieldParsers:
    def test_go_terms_multiple(self):
        result = list(_split_go_terms(
            "cell outer membrane [GO:0009279]; iron ion transport [GO:0006826]"
        ))
        assert len(result) == 2
        assert result[0] == ("GO:0009279", "cell outer membrane")
        assert result[1] == ("GO:0006826", "iron ion transport")

    def test_go_terms_empty(self):
        assert list(_split_go_terms("")) == []
        assert list(_split_go_terms(None)) == []  # type: ignore[arg-type]

    def test_ec_numbers(self):
        ecs = list(_split_ec_numbers("EC:7.1.2.2"))
        assert ecs == ["7.1.2.2"]

    def test_ec_numbers_partial(self):
        ecs = list(_split_ec_numbers("3.6.1.- ; 1.2.-.-"))
        assert "3.6.1.-" in ecs
        assert "1.2.-.-" in ecs

    def test_ec_dedup(self):
        ecs = list(_split_ec_numbers("7.1.2.2; 7.1.2.2; 1.1.1.1"))
        assert ecs == ["7.1.2.2", "1.1.1.1"]

    def test_pfam(self):
        pfams = list(_split_pfam("PF14322;PF07980;"))
        assert pfams == ["PF14322", "PF07980"]


class TestConfidenceMapping:
    """Confidence mapping rules — these are starting thresholds that may be
    refined later, so the tests intentionally pin only the ordinal contract."""

    def test_evalue_high(self):
        assert from_evalue(1e-100) == ConfidenceLevel.HIGH
        assert from_evalue(0.0) == ConfidenceLevel.HIGH
        assert from_evalue(1e-50) == ConfidenceLevel.HIGH  # boundary inclusive

    def test_evalue_medium(self):
        assert from_evalue(1e-30) == ConfidenceLevel.MEDIUM
        assert from_evalue(1e-10) == ConfidenceLevel.MEDIUM  # boundary inclusive

    def test_evalue_low(self):
        assert from_evalue(1e-5) == ConfidenceLevel.LOW
        assert from_evalue(1.0) == ConfidenceLevel.LOW

    def test_evalue_unknown(self):
        assert from_evalue(None) == ConfidenceLevel.UNKNOWN
        assert from_evalue("not a number") == ConfidenceLevel.UNKNOWN  # type: ignore[arg-type]

    def test_evalue_custom_thresholds(self):
        strict = HomologyThresholds(high=1e-100, medium=1e-50)
        assert from_evalue(1e-60, strict) == ConfidenceLevel.MEDIUM
        assert from_evalue(1e-30, strict) == ConfidenceLevel.LOW

    def test_ml_score_high(self):
        assert from_ml_score(0.95) == ConfidenceLevel.HIGH
        assert from_ml_score(0.7) == ConfidenceLevel.HIGH  # boundary inclusive

    def test_ml_score_medium(self):
        assert from_ml_score(0.5) == ConfidenceLevel.MEDIUM
        assert from_ml_score(0.3) == ConfidenceLevel.MEDIUM  # boundary inclusive

    def test_ml_score_low(self):
        assert from_ml_score(0.1) == ConfidenceLevel.LOW
        assert from_ml_score(0.0) == ConfidenceLevel.LOW

    def test_ml_score_clean_thresholds(self):
        # CLEAN-specific thresholds: high=0.5, medium=0.1
        assert from_ml_score(0.6, DEFAULT_CLEAN) == ConfidenceLevel.HIGH
        assert from_ml_score(0.3, DEFAULT_CLEAN) == ConfidenceLevel.MEDIUM
        assert from_ml_score(0.05, DEFAULT_CLEAN) == ConfidenceLevel.LOW

    def test_ordinal_consistency(self):
        # As e-value increases (worse), confidence should be monotone non-increasing
        levels = [from_evalue(e) for e in [1e-100, 1e-30, 1e-5, 1.0]]
        order = {ConfidenceLevel.HIGH: 3, ConfidenceLevel.MEDIUM: 2, ConfidenceLevel.LOW: 1}
        ranks = [order[l] for l in levels]
        assert ranks == sorted(ranks, reverse=True)


class TestEggnogFieldParsers:
    """Tests for the eggNOG-specific helpers."""

    def test_split_csv_basic(self):
        from src.parsers.eggnog import _split_csv
        assert list(_split_csv("a,b,c")) == ["a", "b", "c"]

    def test_split_csv_strips_whitespace(self):
        from src.parsers.eggnog import _split_csv
        assert list(_split_csv("a, b , c")) == ["a", "b", "c"]

    def test_split_csv_dedupes(self):
        from src.parsers.eggnog import _split_csv
        assert list(_split_csv("a,b,a,c,b")) == ["a", "b", "c"]

    def test_split_csv_handles_dash(self):
        from src.parsers.eggnog import _split_csv
        # eggNOG uses '-' for missing — treated as empty field
        assert list(_split_csv("")) == []

    def test_split_eggnog_ogs_simple(self):
        from src.parsers.eggnog import _split_eggnog_ogs
        result = list(_split_eggnog_ogs("arCOG01151@1|root,arCOG01151@2157|Archaea"))
        # Same OG appears twice at different tax levels — should dedupe
        assert result == ["arCOG01151"]

    def test_split_eggnog_ogs_multiple(self):
        from src.parsers.eggnog import _split_eggnog_ogs
        ogs = list(_split_eggnog_ogs(
            "COG1152@1|root,arCOG02428@2157|Archaea,2XTBY@28890|Euryarchaeota,2N9CY@224756|Methanomicrobia"
        ))
        assert ogs == ["COG1152", "arCOG02428", "2XTBY", "2N9CY"]

    def test_split_ec_eggnog(self):
        from src.parsers.eggnog import _split_ec
        assert list(_split_ec("1.2.7.4")) == ["1.2.7.4"]

    def test_split_ec_multiple(self):
        from src.parsers.eggnog import _split_ec
        ecs = list(_split_ec("1.2.7.4,3.6.1.-"))
        assert "1.2.7.4" in ecs
        assert "3.6.1.-" in ecs

    def test_safe_str_handles_dash(self):
        from src.parsers.eggnog import _safe_str
        # eggNOG uses '-' for missing values
        assert _safe_str("-") == ""
        assert _safe_str("real_value") == "real_value"
        assert _safe_str("  spaced  ") == "spaced"


class TestCleanFieldParser:
    """Tests for the CLEAN-specific helpers."""

    def test_split_ec_basic(self):
        from src.parsers.clean import _split_ec_field
        ec, score = _split_ec_field("EC:1.2.7.4/0.9367")
        assert ec == "1.2.7.4"
        assert score == 0.9367

    def test_split_ec_zero_score(self):
        from src.parsers.clean import _split_ec_field
        ec, score = _split_ec_field("EC:3.6.1.54/0.0000")
        assert ec == "3.6.1.54"
        assert score == 0.0

    def test_split_ec_partial_number(self):
        from src.parsers.clean import _split_ec_field
        # Partial EC numbers (e.g. 3.6.1.-) still parse correctly
        ec, score = _split_ec_field("EC:3.6.1.-/0.5")
        assert ec == "3.6.1.-"
        assert score == 0.5

    def test_split_ec_strips_prefix(self):
        from src.parsers.clean import _split_ec_field
        # Both prefixed and non-prefixed forms work
        ec_with, _ = _split_ec_field("EC:1.1.1.1/0.5")
        ec_without, _ = _split_ec_field("1.1.1.1/0.5")
        assert ec_with == ec_without == "1.1.1.1"

    def test_split_ec_empty(self):
        from src.parsers.clean import _split_ec_field
        ec, score = _split_ec_field("")
        assert ec is None and score is None

    def test_split_ec_malformed(self):
        from src.parsers.clean import _split_ec_field
        # No slash separator — malformed
        ec, score = _split_ec_field("EC:1.1.1.1")
        assert ec is None and score is None

    def test_split_ec_non_numeric_score(self):
        from src.parsers.clean import _split_ec_field
        ec, score = _split_ec_field("EC:1.1.1.1/notanumber")
        assert ec is None and score is None

    def test_split_ec_strips_whitespace(self):
        from src.parsers.clean import _split_ec_field
        ec, score = _split_ec_field("  EC:1.1.1.1/0.5  ")
        assert ec == "1.1.1.1"
        assert score == 0.5


class TestDeepGO2Parser:
    """Tests for DeepGO2 parser logic — uses tmp files to test end-to-end."""

    def _make_files(self, tmp_path):
        """Create a minimal DeepGO2 directory layout in a temp folder."""
        d = tmp_path / "deepgo2"
        d.mkdir()
        (d / "subset_f_preds_bp.tsv").write_text(
            "tr|A0B9K2|AATA_METTP\tGO:0006754\t0.85\n"
            "tr|A0B9K2|AATA_METTP\tGO:0009987\t0.32\n"
            "tr|A0B9K2|AATA_METTP\tGO:0050896\t0.12\n"
        )
        (d / "subset_f_preds_cc.tsv").write_text(
            "tr|A0B9K2|AATA_METTP\tGO:0005886\t0.71\n"
        )
        (d / "subset_f_preds_mf.tsv").write_text(
            "tr|A0B9K2|AATA_METTP\tGO:0005524\t0.95\n"
            "tr|A0B9K2|AATA_METTP\tGO:0016887\t0.55\n"
        )
        return tmp_path

    def test_aspect_assignment(self, tmp_path):
        """Each file's records should be tagged with the correct GO aspect."""
        from src.parsers.deepgo2 import DeepGO2Parser
        from src.schema import AnnotationType
        self._make_files(tmp_path)
        parser = DeepGO2Parser(raw_data_root=tmp_path)
        records = parser.parse()
        bp = [r for r in records if r.annotation_type == AnnotationType.GO_BP]
        cc = [r for r in records if r.annotation_type == AnnotationType.GO_CC]
        mf = [r for r in records if r.annotation_type == AnnotationType.GO_MF]
        assert len(bp) == 3
        assert len(cc) == 1
        assert len(mf) == 2

    def test_evidence_rank_sorted_by_score(self, tmp_path):
        """Within a (protein, aspect), rank should reflect score order."""
        from src.parsers.deepgo2 import DeepGO2Parser
        from src.schema import AnnotationType
        self._make_files(tmp_path)
        parser = DeepGO2Parser(raw_data_root=tmp_path)
        records = parser.parse()
        bp = [r for r in records if r.annotation_type == AnnotationType.GO_BP]
        bp.sort(key=lambda r: r.evidence_rank)
        # Highest score (0.85, GO:0006754) should be rank 1
        assert bp[0].evidence_rank == 1
        assert bp[0].value == "GO:0006754"
        assert bp[0].score == 0.85
        # Lowest (0.12) should be rank 3
        assert bp[2].evidence_rank == 3
        assert bp[2].score == 0.12

    def test_min_score_filtering(self, tmp_path):
        """Predictions below min_score should be dropped."""
        from src.parsers.deepgo2 import DeepGO2Parser
        from src.schema import AnnotationType
        self._make_files(tmp_path)
        parser = DeepGO2Parser(raw_data_root=tmp_path, min_score=0.5)
        records = parser.parse()
        bp = [r for r in records if r.annotation_type == AnnotationType.GO_BP]
        # Only one BP prediction has score >= 0.5 (the 0.85 one)
        assert len(bp) == 1
        assert bp[0].value == "GO:0006754"

    def test_confidence_levels(self, tmp_path):
        """Confidence categories should follow ML thresholds."""
        from src.parsers.deepgo2 import DeepGO2Parser
        from src.schema import ConfidenceLevel
        self._make_files(tmp_path)
        parser = DeepGO2Parser(raw_data_root=tmp_path)
        records = parser.parse()
        # Score 0.95 → high; 0.32 → medium; 0.12 → low (using default ML thresholds)
        scores_to_levels = {r.score: r.confidence_level for r in records}
        assert scores_to_levels[0.95] == ConfidenceLevel.HIGH
        assert scores_to_levels[0.32] == ConfidenceLevel.MEDIUM
        assert scores_to_levels[0.12] == ConfidenceLevel.LOW

    def test_missing_directory(self, tmp_path):
        """Should not crash when the deepgo2 directory does not exist."""
        from src.parsers.deepgo2 import DeepGO2Parser
        parser = DeepGO2Parser(raw_data_root=tmp_path)
        records = parser.parse()
        assert records == []


class TestDeepFRIParser:
    """Tests for the DeepFRI parser."""

    def _make_csv(self, tmp_path):
        """Create a minimal DeepFRI CSV in a temp folder."""
        d = tmp_path / "deepfri"
        d.mkdir()
        # Note the quoted label with comma — tests CSV quoting
        (d / "DeepFRI_MF_predictions.csv").write_text(
            "Protein,GO_term/EC_number,Score,GO_term/EC_number name\n"
            "tr|A0B5T7|A0B5T7_METTP,GO:0016829,0.85,lyase activity\n"
            "tr|A0B5T7|A0B5T7_METTP,GO:0016788,\"0.45\",\"hydrolase activity, acting on ester bonds\"\n"
            "tr|A0B5T7|A0B5T7_METTP,GO:0016810,0.10,\"hydrolase activity, acting on carbon-nitrogen bonds\"\n"
            "tr|X5DZ82|X5DZ82_9BACT,GO:0005524,0.95,ATP binding\n"
        )
        return tmp_path

    def test_basic_parse(self, tmp_path):
        from src.parsers.deepfri import DeepFRIParser
        from src.schema import AnnotationType, SourceTool
        self._make_csv(tmp_path)
        parser = DeepFRIParser(raw_data_root=tmp_path)
        records = parser.parse()
        assert len(records) == 4
        # All records should be GO_MF
        assert all(r.annotation_type == AnnotationType.GO_MF for r in records)
        # All records should be tagged with the deepfri tool
        assert all(r.source_tool == SourceTool.DEEPFRI for r in records)

    def test_label_preservation(self, tmp_path):
        """Labels with commas should be preserved correctly."""
        from src.parsers.deepfri import DeepFRIParser
        self._make_csv(tmp_path)
        parser = DeepFRIParser(raw_data_root=tmp_path)
        records = parser.parse()
        labels_by_go = {r.value: r.label for r in records}
        assert labels_by_go["GO:0016829"] == "lyase activity"
        # Quoted label with comma must round-trip intact
        assert labels_by_go["GO:0016788"] == "hydrolase activity, acting on ester bonds"
        assert labels_by_go["GO:0005524"] == "ATP binding"

    def test_evidence_rank_sorted_by_score(self, tmp_path):
        from src.parsers.deepfri import DeepFRIParser
        self._make_csv(tmp_path)
        parser = DeepFRIParser(raw_data_root=tmp_path)
        records = parser.parse()
        a0_records = [r for r in records if r.uniprot_accession == "A0B5T7"]
        a0_records.sort(key=lambda r: r.evidence_rank)
        # Highest score (0.85) should be rank 1
        assert a0_records[0].evidence_rank == 1
        assert a0_records[0].value == "GO:0016829"
        # Lowest (0.10) should be rank 3
        assert a0_records[2].evidence_rank == 3
        assert a0_records[2].value == "GO:0016810"

    def test_min_score_filtering(self, tmp_path):
        from src.parsers.deepfri import DeepFRIParser
        self._make_csv(tmp_path)
        parser = DeepFRIParser(raw_data_root=tmp_path, min_score=0.5)
        records = parser.parse()
        # Only score >= 0.5 should pass: 0.85 and 0.95
        assert len(records) == 2
        scores = sorted(r.score for r in records)
        assert scores == [0.85, 0.95]

    def test_confidence_levels(self, tmp_path):
        from src.parsers.deepfri import DeepFRIParser
        from src.schema import ConfidenceLevel
        self._make_csv(tmp_path)
        parser = DeepFRIParser(raw_data_root=tmp_path)
        records = parser.parse()
        scores_to_levels = {r.score: r.confidence_level for r in records}
        # Default thresholds: high>=0.7, medium>=0.3
        assert scores_to_levels[0.95] == ConfidenceLevel.HIGH
        assert scores_to_levels[0.85] == ConfidenceLevel.HIGH
        assert scores_to_levels[0.45] == ConfidenceLevel.MEDIUM
        assert scores_to_levels[0.10] == ConfidenceLevel.LOW

    def test_missing_file(self, tmp_path):
        from src.parsers.deepfri import DeepFRIParser
        parser = DeepFRIParser(raw_data_root=tmp_path)
        records = parser.parse()
        assert records == []


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


class TestRecognizerParser:
    """Tests for the reCOGnizer parser."""

    def _make_tsv(self, tmp_path):
        """Minimal reCOGnizer TSV with representative DB types."""
        d = tmp_path / "recognizer_results"
        d.mkdir()
        lines = [
            "qseqid\tDB ID\tProtein description\tEC number\tCDD ID\ttaxonomic_range_name\ttaxonomic_range\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tGeneral functional category\tFunctional category\tKO",
            "sp|A0B9K2|AATA_METTP\tCOG1155\tATPase subunit A\t3.6.3.14\tCDD:440769\t\t\t62.0\t582.0\t213.0\t3.0\t1.0\t575.0\t2.0\t583.0\t0.0\t1080.0\tMETABOLISM\tEnergy production and conversion\tK02117",
            "sp|A0B9K2|AATA_METTP\tKOG1350\tF0F1-type ATP synthase beta subunit\t\tCDD:229292\t\t\t31.0\t299.0\t183.0\t8.0\t205.0\t497.0\t175.0\t458.0\t1.32e-28\t117.0\tMETABOLISM\tEnergy production and conversion\t",
            "sp|A0B9K2|AATA_METTP\tpfam00006\tATP synthase alpha/beta nucleotide binding\t\tCDD:999\t\t\t55.0\t300.0\t50.0\t1.0\t1.0\t300.0\t1.0\t300.0\t9.43e-107\t317.0\t\t\t",
            "sp|X5DZ82|X5DZ82_9BACT\tTIGR01043\tATP synthase subunit A\t7.1.2.2\tCDD:888\t\t\t60.0\t580.0\t200.0\t2.0\t1.0\t580.0\t1.0\t580.0\t0.0\t962.0\t\tEnergy production and conversion\t",
            "sp|A0B9K2|AATA_METTP\tcl38909\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t",
            "sp|A0B5T7|A0B5T7_METTP\tCOG0148\tEnolase\t4.2.1.11\tCDD:777\t\t\t90.0\t430.0\t10.0\t0.0\t1.0\t430.0\t1.0\t430.0\t1e-5\t200.0\tMETABOLISM\tCarbohydrate transport and metabolism\t",
        ]
        (d / "reCOGnizer_results.tsv").write_text("\n".join(lines))
        return tmp_path

    def test_basic_parse(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        assert len(records) > 0

    def test_cl_rows_skipped(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        assert "cl38909" not in [r.value for r in records]

    def test_db_type_mapping(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        from src.schema import AnnotationType
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        types_by_value = {r.value: r.annotation_type for r in records}
        assert types_by_value.get("COG1155") == AnnotationType.COG
        assert types_by_value.get("KOG1350") == AnnotationType.KOG
        assert types_by_value.get("pfam00006") == AnnotationType.PFAM
        assert types_by_value.get("TIGR01043") == AnnotationType.TIGR

    def test_secondary_records_ec(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        from src.schema import AnnotationType
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        ec_values = {r.value for r in records if r.annotation_type == AnnotationType.EC}
        assert "3.6.3.14" in ec_values
        assert "7.1.2.2" in ec_values

    def test_secondary_records_ko(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        from src.schema import AnnotationType
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        ko_records = [r for r in records if r.annotation_type == AnnotationType.KEGG_KO]
        assert any(r.value == "K02117" for r in ko_records)

    def test_secondary_records_cog_category(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        from src.schema import AnnotationType
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        cat_values = {r.value for r in records if r.annotation_type == AnnotationType.COG_CATEGORY}
        assert "Energy production and conversion" in cat_values

    def test_confidence_levels(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        from src.schema import AnnotationType, ConfidenceLevel
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        cog_rec = next(r for r in records if r.value == "COG1155")
        assert cog_rec.confidence_level == ConfidenceLevel.HIGH
        low_rec = next(r for r in records if r.value == "COG0148")
        assert low_rec.confidence_level == ConfidenceLevel.LOW

    def test_top_n_per_db_filtering(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        from src.schema import AnnotationType
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path, top_n_per_db=1)
        records = parser.parse()
        kog_recs = [r for r in records
                    if r.annotation_type == AnnotationType.KOG
                    and r.uniprot_accession == "A0B9K2"]
        assert len(kog_recs) == 1

    def test_evidence_rank_assigned(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        from src.schema import AnnotationType
        self._make_tsv(tmp_path)
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        cog_recs = [r for r in records
                    if r.annotation_type == AnnotationType.COG
                    and r.uniprot_accession == "A0B9K2"]
        assert cog_recs[0].evidence_rank == 1

    def test_missing_file(self, tmp_path):
        from src.parsers.recognizer import RecognizerParser
        parser = RecognizerParser(raw_data_root=tmp_path)
        records = parser.parse()
        assert records == []
