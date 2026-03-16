"""Tests for ploidy detection, ploidy-aware allele counting, and panel sex loading."""

import sys
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from graphpop_import.vcf_parser import (
    _build_ploidy_packed,
    _detect_ploidy,
    _load_panel,
    _vectorized_gt_pack,
)


# ---------------------------------------------------------------------------
# Helpers: mock cyvcf2 Variant for _detect_ploidy
# ---------------------------------------------------------------------------


def _mock_variant(genotypes_list):
    """Create a mock cyvcf2 Variant with given genotypes list.

    Each entry is [a0, a1, is_phased] (diploid) or [a0, is_phased] (haploid).
    """
    v = MagicMock()
    v.genotypes = genotypes_list
    return v


# ---------------------------------------------------------------------------
# _detect_ploidy tests
# ---------------------------------------------------------------------------


class TestDetectPloidy:
    def test_all_diploid(self):
        """All samples have 3-element genotypes → all_diploid, no haploid flags."""
        v = _mock_variant([
            [0, 0, True],   # HOM_REF diploid
            [0, 1, True],   # HET diploid
            [1, 1, True],   # HOM_ALT diploid
        ])
        mode, flags = _detect_ploidy(v)
        assert mode == "all_diploid"
        assert flags.sum() == 0
        assert len(flags) == 3

    def test_all_haploid(self):
        """All samples have 2-element genotypes → all_haploid."""
        v = _mock_variant([
            [0, True],   # REF haploid
            [1, True],   # ALT haploid
            [0, False],  # REF haploid unphased
            [1, True],   # ALT haploid
        ])
        mode, flags = _detect_ploidy(v)
        assert mode == "all_haploid"
        assert flags.sum() == 4
        assert all(flags)

    def test_mixed(self):
        """Mix of diploid and haploid → mixed mode with correct flags."""
        v = _mock_variant([
            [0, 0, True],   # diploid female
            [0, 1, True],   # diploid female (HET)
            [0, True],      # haploid male
            [1, True],      # haploid male
            [0, True],      # haploid male
        ])
        mode, flags = _detect_ploidy(v)
        assert mode == "mixed"
        expected = [False, False, True, True, True]
        np.testing.assert_array_equal(flags, expected)


# ---------------------------------------------------------------------------
# Ploidy-aware allele counting (inline tests via _stream logic)
# ---------------------------------------------------------------------------


class TestAlleleCounting:
    def test_haploid_counts(self):
        """All-haploid: 3 REF + 2 ALT → ac=2, an=5."""
        # cyvcf2 gt_types for haploid: REF→0 (HOM_REF), ALT→3 (HOM_ALT)
        gt_types = np.array([0, 0, 0, 3, 3], dtype=np.int8)

        # Haploid counting: ac = n_hom_alt, an = n_called
        n_miss = int(np.sum(gt_types == 2))
        n_hom_alt = int(np.sum(gt_types == 3))
        n_called = len(gt_types) - n_miss

        ac = n_hom_alt          # 2
        an = n_called           # 5
        af = ac / an if an > 0 else 0.0

        assert ac == 2
        assert an == 5
        assert abs(af - 0.4) < 1e-10

    def test_mixed_chrX_counts(self):
        """Mixed chrX: 2 diploid females + 3 haploid males.

        Females: 1 HET (gt=1), 1 HOM_REF (gt=0)
        Males: 1 ALT (gt=3), 2 REF (gt=0)
        Expected: ac = (1 + 0) + 1 = 2 diploid + 1 haploid = 2
                  Actually: ac_dip = n_het + 2*n_hom_alt_dip = 1 + 0 = 1
                            ac_hap = n_hom_alt_hap = 1
                            ac = 1 + 1 = 2
                  an = 2*2 + 3 = 7 (2 diploid females × 2 alleles + 3 haploid males × 1)
        """
        gt_types = np.array([1, 0, 3, 0, 0], dtype=np.int8)
        haploid_flags = np.array([False, False, True, True, True])

        idx = np.arange(5)
        hap_k = haploid_flags[idx]

        # Diploid subset
        dip_mask = ~hap_k
        gt_dip = gt_types[dip_mask]
        n_miss_dip = int(np.sum(gt_dip == 2))
        n_het_dip = int(np.sum(gt_dip == 1))
        n_hom_alt_dip = int(np.sum(gt_dip == 3))
        n_called_dip = len(gt_dip) - n_miss_dip
        ac_dip = n_het_dip + 2 * n_hom_alt_dip
        an_dip = 2 * n_called_dip

        # Haploid subset
        gt_hap = gt_types[hap_k]
        n_miss_hap = int(np.sum(gt_hap == 2))
        n_hom_alt_hap = int(np.sum(gt_hap == 3))
        n_called_hap = len(gt_hap) - n_miss_hap
        ac_hap = n_hom_alt_hap
        an_hap = n_called_hap

        ac = ac_dip + ac_hap
        an = an_dip + an_hap
        het_count = n_het_dip

        assert ac_dip == 1, f"ac_dip={ac_dip}"
        assert an_dip == 4, f"an_dip={an_dip}"
        assert ac_hap == 1, f"ac_hap={ac_hap}"
        assert an_hap == 3, f"an_hap={an_hap}"
        assert ac == 2
        assert an == 7
        assert het_count == 1

    def test_mixed_with_missing(self):
        """Mixed ploidy with missing samples excluded correctly."""
        # 2 diploid (1 HET, 1 MISSING) + 2 haploid (1 ALT, 1 MISSING)
        gt_types = np.array([1, 2, 3, 2], dtype=np.int8)
        haploid_flags = np.array([False, False, True, True])

        idx = np.arange(4)
        hap_k = haploid_flags[idx]

        # Diploid
        gt_dip = gt_types[~hap_k]
        n_miss_dip = int(np.sum(gt_dip == 2))
        n_het_dip = int(np.sum(gt_dip == 1))
        n_hom_alt_dip = int(np.sum(gt_dip == 3))
        n_called_dip = len(gt_dip) - n_miss_dip
        ac_dip = n_het_dip + 2 * n_hom_alt_dip
        an_dip = 2 * n_called_dip

        # Haploid
        gt_hap = gt_types[hap_k]
        n_miss_hap = int(np.sum(gt_hap == 2))
        n_hom_alt_hap = int(np.sum(gt_hap == 3))
        n_called_hap = len(gt_hap) - n_miss_hap
        ac_hap = n_hom_alt_hap
        an_hap = n_called_hap

        ac = ac_dip + ac_hap
        an = an_dip + an_hap

        # 1 called diploid (HET) → ac_dip=1, an_dip=2
        # 1 called haploid (ALT) → ac_hap=1, an_hap=1
        assert ac == 2
        assert an == 3


# ---------------------------------------------------------------------------
# _load_panel sex tests
# ---------------------------------------------------------------------------


class TestLoadPanelSex:
    def test_panel_with_gender(self, tmp_path):
        """1000G panel format with 'gender' column → sex dict populated."""
        panel = tmp_path / "panel.tsv"
        panel.write_text(
            "sample\tpop\tsuper_pop\tgender\n"
            "S1\tGBR\tEUR\tmale\n"
            "S2\tGBR\tEUR\tfemale\n"
            "S3\tYRI\tAFR\tmale\n"
        )
        sample_to_pop, sample_to_sex = _load_panel(panel, "superpopulation")
        assert sample_to_pop == {"S1": "EUR", "S2": "EUR", "S3": "AFR"}
        assert sample_to_sex == {"S1": 1, "S2": 2, "S3": 1}

    def test_ped_with_sex_numeric(self, tmp_path):
        """PED format with numeric 'Sex' column (1=male, 2=female)."""
        ped = tmp_path / "panel.ped"
        ped.write_text(
            "SampleID\tPopulation\tSuperpopulation\tSex\n"
            "NA001\tGBR\tEUR\t1\n"
            "NA002\tGBR\tEUR\t2\n"
            "NA003\tCHS\tEAS\t1\n"
        )
        sample_to_pop, sample_to_sex = _load_panel(ped, "superpopulation")
        assert sample_to_pop == {"NA001": "EUR", "NA002": "EUR", "NA003": "EAS"}
        assert sample_to_sex == {"NA001": 1, "NA002": 2, "NA003": 1}

    def test_panel_without_sex(self, tmp_path):
        """Panel without sex/gender column → empty sex dict."""
        panel = tmp_path / "panel.tsv"
        panel.write_text(
            "sample\tpop\tsuper_pop\n"
            "R1\tJAP\tEAS\n"
            "R2\tIND\tSAS\n"
        )
        sample_to_pop, sample_to_sex = _load_panel(panel, "superpopulation")
        assert sample_to_pop == {"R1": "EAS", "R2": "SAS"}
        assert sample_to_sex == {}


# ---------------------------------------------------------------------------
# _vectorized_gt_pack and _build_ploidy_packed
# ---------------------------------------------------------------------------


class TestPacking:
    def test_vectorized_gt_pack_roundtrip(self):
        """Vectorized packing matches manual unpacking for all genotype values."""
        gt_types = np.array([0, 1, 3, 2, 0, 3, 1, 2, 0], dtype=np.int8)
        packed = _vectorized_gt_pack(gt_types)

        # Remap: cyvcf2 [0,1,2,3] → packed [0,1,3,2]
        remap = {0: 0, 1: 1, 2: 3, 3: 2}

        for i, gt in enumerate(gt_types):
            byte_idx = i >> 2
            bit_shift = (i & 3) << 1
            extracted = (packed[byte_idx] >> bit_shift) & 0x03
            assert extracted == remap[gt], f"sample {i}: got {extracted}, expected {remap[gt]}"

    def test_ploidy_packed_roundtrip(self):
        """Ploidy packing produces correct bit layout."""
        flags = np.array([True, False, True, False, False, True, False, False, True])
        packed = _build_ploidy_packed(flags)

        for i, expected in enumerate(flags):
            byte_idx = i >> 3
            bit_idx = i & 7
            bit = (packed[byte_idx] >> bit_idx) & 1
            assert bit == int(expected), f"sample {i}: got {bit}, expected {int(expected)}"

    def test_ploidy_packed_all_diploid(self):
        """All-diploid flags produce all-zero bytes."""
        flags = np.zeros(16, dtype=bool)
        packed = _build_ploidy_packed(flags)
        assert all(b == 0 for b in packed)
        assert len(packed) == 2  # 16 samples → 2 bytes
