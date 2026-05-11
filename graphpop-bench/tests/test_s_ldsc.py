"""Tests for the s-LDSC wrapper (Finucane et al. 2015).

Parser + TSV-emit + binary-detection tests run without LDSC;
the integration test skips when ``ldsc.py`` is missing on PATH.
"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import pytest

from graphpop_bench.competitors import (
    CategoryRow,
    SLdscResult,
    SLdscRunner,
    parse_log,
    parse_results,
)
from graphpop_bench.competitors.s_ldsc import (
    _build_partition_tsv,
    _build_total_tsv,
    _build_ldsc_cmd,
    _find_ldsc_binary,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_RESULTS_HEADER = (
    "Category Prop._SNPs Prop._h2 Prop._h2_std_error "
    "Enrichment Enrichment_std_error Enrichment_p "
    "Coefficient Coefficient_std_error Coefficient_z-score\n"
)

_RESULTS_BODY = (
    "baseL2_0   0.234  0.512  0.041  2.190  0.180  3.4e-12  "
    "1.2e-08  3.0e-09  4.0\n"
    "Coding_UCSC.bedL2_0   0.012  0.075  0.020  6.200  1.500  "
    "2.1e-05  4.4e-08  1.0e-08  4.4\n"
    "Conserved_LindbladTohL2_0   0.027  0.103  0.015  3.810  0.450  "
    "1.5e-16  9.0e-09  2.0e-09  4.5\n"
)

_RESULTS_FIXTURE = _RESULTS_HEADER + _RESULTS_BODY

_LOG_FIXTURE = textwrap.dedent("""\
    Reading summary statistics from sumstats.gz ...
    After merging with reference panel LD, 1185001 SNPs remain.
    After merging with regression SNP LD, 1183421 SNPs remain.

    Total Observed scale h2: 0.234 (0.012)
    Lambda GC: 1.421
    Mean Chi^2: 1.612
    Intercept: 1.011 (0.0093)
    Ratio: 0.018 (0.015)

    Reading annotations ...
""")


def _write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


# ---------------------------------------------------------------------------
# parse_results — parser unit tests
# ---------------------------------------------------------------------------

def test_parse_results_canonical_format(tmp_path):
    p = _write(tmp_path, "ldsc.results", _RESULTS_FIXTURE)
    rows = parse_results(p)
    assert len(rows) == 3
    base = rows[0]
    assert base.category == "baseL2_0"
    assert base.prop_snps == pytest.approx(0.234)
    assert base.prop_h2 == pytest.approx(0.512)
    assert base.prop_h2_se == pytest.approx(0.041)
    assert base.enrichment == pytest.approx(2.190)
    assert base.enrichment_se == pytest.approx(0.180)
    assert base.enrichment_p == pytest.approx(3.4e-12)


def test_parse_results_header_only_returns_empty(tmp_path):
    p = _write(tmp_path, "ldsc.results", _RESULTS_HEADER)
    assert parse_results(p) == []


def test_parse_results_empty_file_returns_empty(tmp_path):
    p = _write(tmp_path, "ldsc.results", "")
    assert parse_results(p) == []


def test_parse_results_skips_malformed_numeric_rows(tmp_path):
    bad = _RESULTS_HEADER + (
        "good_cat 0.1 0.2 0.01 1.5 0.1 0.001 0 0 0\n"
        "bad_cat  NA  NA  NA  NA  NA  NA     0 0 0\n"
    )
    p = _write(tmp_path, "ldsc.results", bad)
    rows = parse_results(p)
    assert len(rows) == 1
    assert rows[0].category == "good_cat"


def test_parse_results_missing_required_column_raises(tmp_path):
    """Strip Enrichment from header → parser refuses."""
    bad_header = (
        "Category Prop._SNPs Prop._h2 Prop._h2_std_error "
        "Enrichment_std_error Enrichment_p Coefficient\n"
    )
    p = _write(tmp_path, "ldsc.results", bad_header)
    with pytest.raises(ValueError, match="missing required"):
        parse_results(p)


def test_parse_results_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match=".results file"):
        parse_results(tmp_path / "nope.results")


# ---------------------------------------------------------------------------
# parse_log — log file extraction
# ---------------------------------------------------------------------------

def test_parse_log_extracts_all_metrics(tmp_path):
    p = _write(tmp_path, "ldsc.log", _LOG_FIXTURE)
    h2, h2_se, intercept, intercept_se, lambda_gc, n_snps = parse_log(p)
    assert h2 == pytest.approx(0.234)
    assert h2_se == pytest.approx(0.012)
    assert intercept == pytest.approx(1.011)
    assert intercept_se == pytest.approx(0.0093)
    assert lambda_gc == pytest.approx(1.421)
    assert n_snps == 1183421


def test_parse_log_handles_missing_lambda(tmp_path):
    text = (
        "After merging with regression SNP LD, 500 SNPs remain.\n"
        "Total Observed scale h2: 0.1 (0.01)\n"
        "Intercept: 1.0 (0.01)\n"
    )
    p = _write(tmp_path, "ldsc.log", text)
    h2, _, intercept, _, lambda_gc, n_snps = parse_log(p)
    assert h2 == pytest.approx(0.1)
    assert intercept == pytest.approx(1.0)
    assert lambda_gc is None
    assert n_snps == 500


def test_parse_log_missing_file_returns_all_none(tmp_path):
    result = parse_log(tmp_path / "missing.log")
    assert result == (None, None, None, None, None, None)


def test_parse_log_empty_file_returns_all_none(tmp_path):
    p = _write(tmp_path, "ldsc.log", "")
    assert parse_log(p) == (None, None, None, None, None, None)


# ---------------------------------------------------------------------------
# TSV emission
# ---------------------------------------------------------------------------

def test_build_partition_tsv_emits_header_and_rows(tmp_path):
    rows = [
        CategoryRow("baseL2_0", 0.234, 0.512, 0.041,
                    2.19, 0.18, 3.4e-12),
        CategoryRow("CodingL2_0", 0.012, 0.075, 0.02,
                    6.2, 1.5, 2.1e-05),
    ]
    p = tmp_path / "partition_h2.tsv"
    _build_partition_tsv(rows, p)
    lines = p.read_text().strip().splitlines()
    assert lines[0] == (
        "category\tprop_snps\tprop_h2\tprop_h2_se\t"
        "enrichment\tenrichment_se\tenrichment_p"
    )
    assert lines[1].startswith("baseL2_0\t0.234\t")
    assert lines[2].startswith("CodingL2_0\t")


def test_build_partition_tsv_empty_input(tmp_path):
    p = tmp_path / "partition_h2.tsv"
    _build_partition_tsv([], p)
    lines = p.read_text().strip().splitlines()
    assert len(lines) == 1
    assert lines[0].startswith("category\t")


def test_build_total_tsv_with_all_values(tmp_path):
    p = tmp_path / "total_h2.tsv"
    _build_total_tsv(
        total_h2=0.234, total_h2_se=0.012,
        intercept=1.011, intercept_se=0.0093,
        lambda_gc=1.421, n_snps=1183421, path=p)
    lines = p.read_text().strip().splitlines()
    assert lines[0] == (
        "total_h2\ttotal_h2_se\tintercept\tintercept_se\t"
        "lambda_gc\tn_snps"
    )
    body = lines[1].split("\t")
    assert float(body[0]) == pytest.approx(0.234)
    assert int(body[5]) == 1183421


def test_build_total_tsv_with_missing_values(tmp_path):
    p = tmp_path / "total_h2.tsv"
    _build_total_tsv(
        total_h2=0.1, total_h2_se=0.01,
        intercept=None, intercept_se=None,
        lambda_gc=None, n_snps=None, path=p)
    # Use rstrip("\n") rather than full strip(): the latter would
    # eat trailing tabs and elide the empty-trailing-column cells.
    lines = p.read_text().rstrip("\n").split("\n")
    body = lines[1].split("\t")
    # 6 columns; the last 4 are empty strings.
    assert len(body) == 6
    assert body[0] == f"{0.1:.8g}"
    assert body[1] == f"{0.01:.8g}"
    assert body[2] == "" and body[3] == "" and body[4] == ""
    assert body[5] == ""


# ---------------------------------------------------------------------------
# Binary detection + cmd construction
# ---------------------------------------------------------------------------

def test_is_available_returns_bool():
    assert isinstance(SLdscRunner.is_available(), bool)


def test_detect_version_returns_none_for_nonexistent_binary():
    assert SLdscRunner.detect_version("/no/such/binary") is None


def test_build_ldsc_cmd_minimal():
    cmd = _build_ldsc_cmd(
        "ldsc.py", Path("sumstats.gz"),
        "ref/", "wld/", None, True,
        Path("/tmp/out/ldsc"), [])
    assert cmd[0] == "ldsc.py"
    assert "--h2" in cmd and "sumstats.gz" in cmd
    assert "--ref-ld-chr" in cmd and "ref/" in cmd
    assert "--w-ld-chr" in cmd and "wld/" in cmd
    assert "--overlap-annot" in cmd
    assert "--frqfile-chr" not in cmd


def test_build_ldsc_cmd_with_frq_and_extras():
    cmd = _build_ldsc_cmd(
        "ldsc.py", Path("sumstats.gz"),
        "ref/", "wld/", "frq/", False,
        Path("/tmp/out/ldsc"),
        ["--print-coefficients"])
    assert "--frqfile-chr" in cmd and "frq/" in cmd
    assert "--overlap-annot" not in cmd
    assert "--print-coefficients" in cmd


# ---------------------------------------------------------------------------
# Runner — missing-binary + arg-validation paths
# ---------------------------------------------------------------------------

def test_run_raises_clean_error_when_ldsc_missing(tmp_path):
    runner = SLdscRunner(ldsc_binary=None)
    if runner.ldsc_binary is None:
        with pytest.raises(FileNotFoundError, match="ldsc.py"):
            runner.run(tmp_path / "sumstats.gz", tmp_path / "out",
                       ref_ld_chr="ref/", w_ld_chr="wld/")


# ---------------------------------------------------------------------------
# Integration: only when ldsc.py is on PATH (skipped otherwise)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    _find_ldsc_binary() is None,
    reason="ldsc.py not on PATH; skipping integration smoke test",
)
def test_runner_constructs_when_binary_present():
    runner = SLdscRunner()
    assert runner.ldsc_binary is not None
    # Sanity: we can probe version without crashing.
    SLdscRunner.detect_version(runner.ldsc_binary)
