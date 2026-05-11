"""Unit tests for the 1000G subset-extraction helper (R0.2)."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from graphpop_bench.paper2_drivers.subset_1000g import (
    DEFAULT_PANEL_PATH,
    DEFAULT_VCF_DIR,
    PanelRow,
    SUPER_POPS,
    extract_subset_vcf,
    load_panel,
    samples_for_sub_pop,
    samples_for_super_pop,
    samples_for_super_pop_set,
    vcf_path_for_chr,
)


# ---------------------------------------------------------------------------
# Panel parsing (gated on the real panel file being on disk)
# ---------------------------------------------------------------------------

def _panel_available() -> bool:
    return DEFAULT_PANEL_PATH.exists()


@pytest.mark.skipif(
    not _panel_available(),
    reason="1000G panel TSV not present on this host",
)
def test_load_panel_returns_2504_rows():
    rows = load_panel()
    assert len(rows) == 2504


@pytest.mark.skipif(
    not _panel_available(),
    reason="1000G panel TSV not present on this host",
)
def test_load_panel_rows_have_expected_fields():
    rows = load_panel()
    first = rows[0]
    assert isinstance(first, PanelRow)
    assert first.sample.startswith("HG") or first.sample.startswith("NA")
    assert first.super_pop in set(SUPER_POPS)


@pytest.mark.skipif(
    not _panel_available(),
    reason="1000G panel TSV not present on this host",
)
def test_samples_for_super_pop_AFR_count_661():
    rows = load_panel()
    assert len(samples_for_super_pop(rows, "AFR")) == 661


@pytest.mark.skipif(
    not _panel_available(),
    reason="1000G panel TSV not present on this host",
)
def test_samples_for_super_pop_EUR_count_503():
    rows = load_panel()
    assert len(samples_for_super_pop(rows, "EUR")) == 503


@pytest.mark.skipif(
    not _panel_available(),
    reason="1000G panel TSV not present on this host",
)
def test_samples_for_sub_pop_YRI_count_108():
    rows = load_panel()
    assert len(samples_for_sub_pop(rows, "YRI")) == 108


@pytest.mark.skipif(
    not _panel_available(),
    reason="1000G panel TSV not present on this host",
)
def test_samples_for_super_pop_set_AFR_EUR_returns_1164():
    rows = load_panel()
    assert len(samples_for_super_pop_set(rows, ["AFR", "EUR"])) == 661 + 503


# ---------------------------------------------------------------------------
# Synthetic panel — runs without the real 1000G file on disk
# ---------------------------------------------------------------------------

def test_load_panel_synthetic(tmp_path):
    p = tmp_path / "panel.txt"
    p.write_text(
        "sample\tpop\tsuper_pop\tgender\n"
        "S0\tYRI\tAFR\tmale\n"
        "S1\tGBR\tEUR\tfemale\n"
        "S2\tJPT\tEAS\tfemale\n"
    )
    rows = load_panel(p)
    assert len(rows) == 3
    assert samples_for_super_pop(rows, "AFR") == ["S0"]
    assert samples_for_sub_pop(rows, "GBR") == ["S1"]
    assert samples_for_super_pop_set(
        rows, ["AFR", "EAS"]) == ["S0", "S2"]


def test_load_panel_rejects_bad_header(tmp_path):
    p = tmp_path / "bad.txt"
    p.write_text("sample\tpop\tNOPE\tgender\nS0\tYRI\tAFR\tmale\n")
    with pytest.raises(ValueError, match="header"):
        load_panel(p)


def test_load_panel_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="panel file"):
        load_panel(tmp_path / "nope.tsv")


def test_load_panel_skips_blank_rows(tmp_path):
    p = tmp_path / "p.tsv"
    p.write_text(
        "sample\tpop\tsuper_pop\tgender\n"
        "S0\tYRI\tAFR\tmale\n"
        "\n"
        "S1\tGBR\tEUR\tfemale\n"
    )
    rows = load_panel(p)
    assert len(rows) == 2


# ---------------------------------------------------------------------------
# vcf_path_for_chr
# ---------------------------------------------------------------------------

def test_vcf_path_for_chr_basic():
    p = vcf_path_for_chr(22)
    assert "chr22" in str(p)
    assert str(p).endswith(".vcf.gz")


def test_vcf_path_for_chr_string_input():
    p1 = vcf_path_for_chr("22")
    p2 = vcf_path_for_chr("chr22")
    assert p1 == p2


# ---------------------------------------------------------------------------
# Subset extraction (gated on bcftools + real VCF on disk)
# ---------------------------------------------------------------------------

def _bcftools_available() -> bool:
    return shutil.which("bcftools") is not None


def _chr22_vcf_available() -> bool:
    return vcf_path_for_chr(22).exists()


@pytest.mark.skipif(
    not (_bcftools_available()
         and _chr22_vcf_available()
         and _panel_available()),
    reason="bcftools + chr22 VCF + panel needed for extraction smoke",
)
def test_extract_subset_vcf_micro(tmp_path):
    """Tiny 22:50000000-50100000 region × 10-sample subset."""
    panel = load_panel()
    eur_samples = samples_for_super_pop(panel, "EUR")[:10]
    out = tmp_path / "chr22_EUR_micro.vcf.gz"
    result = extract_subset_vcf(
        input_vcf=vcf_path_for_chr(22),
        output_vcf=out,
        sample_ids=eur_samples,
        region="chr22:50000000-50100000",
    )
    assert result.output_vcf.exists()
    assert result.output_vcf.with_suffix(".gz.tbi").exists()
    assert result.n_samples == 10
    assert result.n_variants >= 0  # may be zero if region is sparse
    # Receipt JSON parses + has expected keys.
    receipt = json.loads(result.receipt_path.read_text())
    assert receipt["n_samples_requested"] == 10
    assert receipt["region"] == "chr22:50000000-50100000"


def test_extract_subset_vcf_empty_samples_raises(tmp_path):
    """API guard: empty sample list should not silently call bcftools."""
    fake_vcf = tmp_path / "x.vcf.gz"
    fake_vcf.write_text("")
    with pytest.raises(ValueError, match="non-empty"):
        extract_subset_vcf(
            input_vcf=fake_vcf,
            output_vcf=tmp_path / "out.vcf.gz",
            sample_ids=[],
        )


def test_extract_subset_vcf_missing_input_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="input VCF"):
        extract_subset_vcf(
            input_vcf=tmp_path / "nope.vcf.gz",
            output_vcf=tmp_path / "out.vcf.gz",
            sample_ids=["S0"],
        )
