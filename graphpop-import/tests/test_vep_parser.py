"""Tests for vep_parser.py — VEP/SnpEff annotation parsing."""

import csv
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from graphpop_import.vep_parser import (
    GENE_HEADER,
    HAS_CONSEQUENCE_HEADER,
    VEPParser,
    _PRED_SCORE_RE,
    _load_variant_id_set,
    _parse_pred_score,
)


# ---------------------------------------------------------------------------
# _parse_pred_score tests
# ---------------------------------------------------------------------------


class TestParsePredScore:
    def test_deleterious_with_score(self):
        """Standard VEP format: 'deleterious(0.02)' -> ('deleterious', '0.02')."""
        pred, score = _parse_pred_score("deleterious(0.02)")
        assert pred == "deleterious"
        assert score == "0.02"

    def test_tolerated_with_score(self):
        pred, score = _parse_pred_score("tolerated(0.45)")
        assert pred == "tolerated"
        assert score == "0.45"

    def test_benign_polyphen(self):
        pred, score = _parse_pred_score("benign(0.001)")
        assert pred == "benign"
        assert score == "0.001"

    def test_probably_damaging(self):
        pred, score = _parse_pred_score("probably_damaging(0.998)")
        assert pred == "probably_damaging"
        assert score == "0.998"

    def test_empty_string(self):
        """Empty string returns ('', '')."""
        pred, score = _parse_pred_score("")
        assert pred == ""
        assert score == ""

    def test_none_like_empty(self):
        """Falsy empty value returns ('', '')."""
        pred, score = _parse_pred_score("")
        assert pred == ""
        assert score == ""

    def test_prediction_only_no_parentheses(self):
        """Value without parentheses returns (value, '')."""
        pred, score = _parse_pred_score("deleterious")
        assert pred == "deleterious"
        assert score == ""

    def test_score_zero(self):
        pred, score = _parse_pred_score("tolerated(0)")
        assert pred == "tolerated"
        assert score == "0"

    def test_score_one(self):
        pred, score = _parse_pred_score("damaging(1.0)")
        assert pred == "damaging"
        assert score == "1.0"


# ---------------------------------------------------------------------------
# _load_variant_id_set tests
# ---------------------------------------------------------------------------


class TestLoadVariantIdSet:
    def test_loads_ids_from_csv(self, tmp_path):
        """Reads first column of variant_nodes.csv (skips header)."""
        csv_path = tmp_path / "variant_nodes.csv"
        csv_path.write_text(
            "variantId:ID(Variant),:LABEL,chr,pos:long\n"
            "chr22:100:A:T,Variant,chr22,100\n"
            "chr22:200:G:C,Variant,chr22,200\n"
            "chr22:300:C:A,Variant,chr22,300\n"
        )
        ids = _load_variant_id_set(csv_path)
        assert ids == {"chr22:100:A:T", "chr22:200:G:C", "chr22:300:C:A"}

    def test_empty_csv(self, tmp_path):
        """CSV with only header returns empty set."""
        csv_path = tmp_path / "variant_nodes.csv"
        csv_path.write_text("variantId:ID(Variant),:LABEL\n")
        ids = _load_variant_id_set(csv_path)
        assert ids == set()


# ---------------------------------------------------------------------------
# VEPParser._resolve_variant_id tests (static method)
# ---------------------------------------------------------------------------


class TestResolveVariantId:
    def test_snp_direct_match(self):
        """SNP allele matches directly."""
        alt_to_vid = {"T": "chr1:100:A:T"}
        matched_alts = {"T"}
        result = VEPParser._resolve_variant_id("T", "A", ["T"], alt_to_vid, matched_alts)
        assert result == "chr1:100:A:T"

    def test_deletion_dash(self):
        """VEP uses '-' for deletions."""
        alt_to_vid = {"A": "chr1:100:AT:A"}
        matched_alts = {"A"}
        result = VEPParser._resolve_variant_id("-", "AT", ["A"], alt_to_vid, matched_alts)
        assert result == "chr1:100:AT:A"

    def test_insertion_bases(self):
        """VEP uses inserted bases (no ref padding) for insertions."""
        alt_to_vid = {"AGC": "chr1:100:A:AGC"}
        matched_alts = {"AGC"}
        result = VEPParser._resolve_variant_id("GC", "A", ["AGC"], alt_to_vid, matched_alts)
        assert result == "chr1:100:A:AGC"

    def test_no_match_returns_none(self):
        """Unmatched allele returns None."""
        alt_to_vid = {"T": "chr1:100:A:T"}
        matched_alts = {"T"}
        result = VEPParser._resolve_variant_id("G", "A", ["T"], alt_to_vid, matched_alts)
        assert result is None

    def test_not_in_matched_alts(self):
        """Direct match exists but allele is not in matched_alts."""
        alt_to_vid = {"T": "chr1:100:A:T", "G": "chr1:100:A:G"}
        matched_alts = {"G"}  # only G is matched
        result = VEPParser._resolve_variant_id("T", "A", ["T", "G"], alt_to_vid, matched_alts)
        assert result is None


# ---------------------------------------------------------------------------
# VEPParser._resolve_ann_variant_id tests (static method for SnpEff)
# ---------------------------------------------------------------------------


class TestResolveAnnVariantId:
    def test_direct_match(self):
        """SnpEff ANN allele matches VCF ALT directly."""
        alt_to_vid = {"T": "chr1:100:A:T"}
        matched_alts = {"T"}
        result = VEPParser._resolve_ann_variant_id("T", alt_to_vid, matched_alts)
        assert result == "chr1:100:A:T"

    def test_not_matched(self):
        """Allele not in matched_alts returns None."""
        alt_to_vid = {"T": "chr1:100:A:T"}
        matched_alts = set()
        result = VEPParser._resolve_ann_variant_id("T", alt_to_vid, matched_alts)
        assert result is None

    def test_unknown_allele(self):
        alt_to_vid = {"T": "chr1:100:A:T"}
        matched_alts = {"T"}
        result = VEPParser._resolve_ann_variant_id("G", alt_to_vid, matched_alts)
        assert result is None


# ---------------------------------------------------------------------------
# VEPParser._parse_csq_format_from_desc tests
# ---------------------------------------------------------------------------


class TestParseCSQFormat:
    def test_standard_csq_description(self):
        """Standard VEP CSQ header description is parsed correctly."""
        desc = (
            'Consequence annotations from Ensembl VEP. '
            'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|'
            'BIOTYPE|SIFT|PolyPhen|CADD_PHRED|REVEL'
        )
        field_idx = VEPParser._parse_csq_format_from_desc(desc)
        assert field_idx["allele"] == 0
        assert field_idx["consequence"] == 1
        assert field_idx["impact"] == 2
        assert field_idx["symbol"] == 3
        assert field_idx["gene"] == 4
        assert field_idx["sift"] == 8
        assert field_idx["polyphen"] == 9
        assert field_idx["cadd"] == 10
        assert field_idx["revel"] == 11

    def test_missing_required_field_raises(self):
        """Missing required field (Gene) raises ValueError."""
        desc = "Format: Allele|Consequence|IMPACT|SYMBOL"
        with pytest.raises(ValueError, match="missing required field"):
            VEPParser._parse_csq_format_from_desc(desc)

    def test_no_format_prefix_raises(self):
        """Description without 'Format: ' raises ValueError."""
        desc = "Some description without format info"
        with pytest.raises(ValueError, match="no 'Format: '"):
            VEPParser._parse_csq_format_from_desc(desc)


# ---------------------------------------------------------------------------
# VEPParser._parse_ann_format tests
# ---------------------------------------------------------------------------


class TestParseANNFormat:
    def test_ann_fixed_layout(self):
        """ANN format returns fixed index map regardless of description."""
        desc = "Functional annotations: 'Allele | Annotation | ...'"
        field_idx = VEPParser._parse_ann_format(desc)
        assert field_idx["allele"] == 0
        assert field_idx["consequence"] == 1
        assert field_idx["impact"] == 2
        assert field_idx["symbol"] == 3
        assert field_idx["gene"] == 4
        assert field_idx["feature_type"] == 5
        assert field_idx["feature"] == 6
        assert field_idx["biotype"] == 7
        # SIFT/PolyPhen should NOT be present for SnpEff
        assert "sift" not in field_idx
        assert "polyphen" not in field_idx


# ---------------------------------------------------------------------------
# _PRED_SCORE_RE regex tests
# ---------------------------------------------------------------------------


class TestPredScoreRegex:
    def test_matches_standard_format(self):
        m = _PRED_SCORE_RE.match("deleterious(0.02)")
        assert m is not None
        assert m.group(1) == "deleterious"
        assert m.group(2) == "0.02"

    def test_no_match_for_bare_word(self):
        m = _PRED_SCORE_RE.match("deleterious")
        assert m is None

    def test_no_match_for_empty(self):
        m = _PRED_SCORE_RE.match("")
        assert m is None

    def test_underscore_in_prediction(self):
        m = _PRED_SCORE_RE.match("probably_damaging(0.999)")
        assert m is not None
        assert m.group(1) == "probably_damaging"
        assert m.group(2) == "0.999"
