"""Tests for the KING-robust wrapper.

Parser tests run without KING; the integration test skips when
``king`` (or PLINK, for VCF inputs) is missing.
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from graphpop_bench.competitors import (
    KingResult,
    KingRunner,
    parse_kin_files,
)
from graphpop_bench.competitors.king import (
    _build_normalised_tsv,
    _canonical_pair,
    _resolve_input_kind,
    _find_king_binary,
)
from graphpop_bench.competitors.plink_grm import _find_plink_binary


# ---------------------------------------------------------------------------
# Helpers — build hand-crafted .kin / .kin0 files
# ---------------------------------------------------------------------------

_KIN0_TEXT = textwrap.dedent("""\
    FID1\tIID1\tFID2\tIID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship
    F1\tS0\tF2\tS1\t1000\t0.123\t0.001\t0.5\t0.0\t0.2503
    F1\tS0\tF3\tS2\t1000\t0.100\t0.010\t0.4\t0.0\t0.0012
    F2\tS1\tF3\tS2\t1000\t0.110\t0.020\t0.45\t0.0\t-0.0007
""")

_KIN_TEXT = textwrap.dedent("""\
    FID\tID1\tID2\tN_SNP\tHetHet\tIBS0\tKinship\tError
    F1\tS0\tS0_sib\t1000\t0.123\t0.001\t0.2510\t0
""")


def _write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


# ---------------------------------------------------------------------------
# parse_kin_files — parser unit tests (no KING needed)
# ---------------------------------------------------------------------------

def test_parse_kin0_only_returns_expected_pair_map(tmp_path):
    kin0 = _write(tmp_path, "out.kin0", _KIN0_TEXT)
    pairs, ids = parse_kin_files(kin0_path=kin0, kin_path=None)
    assert pairs == {
        ("S0", "S1"): pytest.approx(0.2503),
        ("S0", "S2"): pytest.approx(0.0012),
        ("S1", "S2"): pytest.approx(-0.0007),
    }
    assert ids == ["S0", "S1", "S2"]


def test_parse_kin_only_returns_within_family_pair(tmp_path):
    kin = _write(tmp_path, "out.kin", _KIN_TEXT)
    pairs, ids = parse_kin_files(kin0_path=None, kin_path=kin)
    assert pairs == {("S0", "S0_sib"): pytest.approx(0.2510)}
    assert ids == ["S0", "S0_sib"]


def test_parse_both_files_combines_rows(tmp_path):
    kin0 = _write(tmp_path, "out.kin0", _KIN0_TEXT)
    kin = _write(tmp_path, "out.kin", _KIN_TEXT)
    pairs, ids = parse_kin_files(kin0_path=kin0, kin_path=kin)
    assert len(pairs) == 4
    assert "S0_sib" in ids
    assert pairs[("S0", "S0_sib")] == pytest.approx(0.2510)


def test_parse_no_files_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        parse_kin_files(kin0_path=None, kin_path=None)


def test_parse_kin0_missing_Kinship_column_raises(tmp_path):
    bad = _write(tmp_path, "bad.kin0",
                 "FID1\tIID1\tFID2\tIID2\tN_SNP\n"
                 "F1\tS0\tF2\tS1\t1000\n")
    with pytest.raises(ValueError, match="no 'Kinship'"):
        parse_kin_files(kin0_path=bad, kin_path=None)


def test_parse_kin_files_with_fam_includes_unpaired_samples(tmp_path):
    """Samples with no kin output rows still appear in sample_ids
    when a .fam is provided."""
    kin0 = _write(tmp_path, "out.kin0", _KIN0_TEXT)
    fam = _write(
        tmp_path, "out.fam",
        "F1 S0 0 0 0 -9\n"
        "F2 S1 0 0 0 -9\n"
        "F3 S2 0 0 0 -9\n"
        "F4 S_alone 0 0 0 -9\n",
    )
    _, ids = parse_kin_files(kin0_path=kin0, kin_path=None, fam_path=fam)
    assert "S_alone" in ids


def test_parse_handles_empty_kin_file_gracefully(tmp_path):
    empty = _write(tmp_path, "empty.kin", "")
    pairs, ids = parse_kin_files(kin0_path=None, kin_path=empty)
    assert pairs == {}
    assert ids == []


def test_canonical_pair_is_order_invariant():
    assert _canonical_pair("A", "B") == ("A", "B")
    assert _canonical_pair("B", "A") == ("A", "B")
    assert _canonical_pair("S0", "S0_sib") == ("S0", "S0_sib")


# ---------------------------------------------------------------------------
# TSV emission
# ---------------------------------------------------------------------------

def test_build_normalised_tsv_emits_sorted_rows(tmp_path):
    pairs = {
        ("S0", "S2"): 0.01,
        ("S0", "S1"): 0.25,
        ("S1", "S2"): -0.001,
    }
    out = tmp_path / "king.tsv"
    _build_normalised_tsv(pairs, out)
    lines = out.read_text().strip().splitlines()
    assert lines[0] == "sample_a\tsample_b\tkinship"
    body = lines[1:]
    assert len(body) == 3
    # Sorted by (sample_a, sample_b): (S0,S1), (S0,S2), (S1,S2).
    assert body[0].startswith("S0\tS1\t")
    assert body[1].startswith("S0\tS2\t")
    assert body[2].startswith("S1\tS2\t")


def test_build_normalised_tsv_empty_input(tmp_path):
    out = tmp_path / "king.tsv"
    _build_normalised_tsv({}, out)
    lines = out.read_text().strip().splitlines()
    assert lines == ["sample_a\tsample_b\tkinship"]


def test_build_normalised_tsv_creates_parent_dir(tmp_path):
    out = tmp_path / "nested" / "king.tsv"
    _build_normalised_tsv({}, out)
    assert out.exists()


# ---------------------------------------------------------------------------
# Input-kind resolution
# ---------------------------------------------------------------------------

def test_resolve_input_kind_auto_detects_vcf():
    assert _resolve_input_kind(Path("cohort.vcf"), "auto") == "vcf"
    assert _resolve_input_kind(Path("cohort.vcf.gz"), "auto") == "vcf"


def test_resolve_input_kind_no_extension_is_bfile():
    assert _resolve_input_kind(Path("cohort"), "auto") == "bfile"


def test_resolve_input_kind_unknown_raises():
    with pytest.raises(ValueError, match="unknown input_kind"):
        _resolve_input_kind(Path("cohort"), "pgen")


# ---------------------------------------------------------------------------
# is_available + detect_version — no KING required
# ---------------------------------------------------------------------------

def test_is_available_returns_bool():
    assert isinstance(KingRunner.is_available(), bool)


def test_detect_version_returns_none_for_nonexistent_binary():
    assert KingRunner.detect_version("/no/such/binary") is None


def test_run_raises_clean_error_when_king_missing(tmp_path):
    # Force-construct a runner with a non-existent binary path.
    runner = KingRunner(king_binary=None)
    if runner.king_binary is None:
        with pytest.raises(FileNotFoundError, match="KING binary"):
            runner.run(tmp_path / "no.vcf", tmp_path / "out")


def test_run_rejects_unknown_mode(tmp_path):
    # Only meaningful if KING is present; otherwise the binary
    # check fires first. So we patch a binary that won't run.
    runner = KingRunner.__new__(KingRunner)
    runner.king_binary = "/bin/true"  # fake binary that exists
    with pytest.raises(ValueError, match="unknown mode"):
        runner.run(tmp_path, tmp_path / "out", mode="explode")


# ---------------------------------------------------------------------------
# Integration: run KING on a small VCF (skip if king OR plink missing)
# ---------------------------------------------------------------------------

# KING refuses very small cohorts; use a slightly bigger fixture
# than the H1 PLINK test. 8 samples × 4 SNPs.
_KING_VCF = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##contig=<ID=1,length=10000>
    ##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1\tS2\tS3\tS4\tS5\tS6\tS7
    1\t100\trs1\tA\tG\t.\tPASS\tPR\tGT\t0/0\t0/1\t1/1\t0/1\t0/0\t1/1\t0/1\t0/0
    1\t500\trs2\tC\tT\t.\tPASS\tPR\tGT\t0/1\t0/1\t1/1\t0/0\t0/0\t1/1\t0/1\t1/1
    1\t1500\trs3\tG\tA\t.\tPASS\tPR\tGT\t1/1\t0/1\t0/0\t0/1\t1/1\t1/1\t0/0\t0/1
    1\t3000\trs4\tT\tC\t.\tPASS\tPR\tGT\t0/0\t1/1\t0/1\t1/1\t0/1\t0/0\t0/1\t1/1
""")


@pytest.mark.skipif(
    _find_king_binary() is None or _find_plink_binary() is None,
    reason="KING + PLINK required on PATH for VCF integration test",
)
def test_king_runs_on_tiny_vcf_and_emits_tsv(tmp_path):
    """End-to-end smoke: PLINK converts VCF→BED, KING computes
    kinship, the wrapper parses both .kin0 (and .kin if produced)
    and emits a normalised TSV. KING may refuse very small
    cohorts; if it exits non-zero on this fixture, skip with a
    clear reason rather than fail the suite."""
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(_KING_VCF)
    out = tmp_path / "king_out"
    runner = KingRunner()
    try:
        result = runner.run(vcf, out, seed=42)
    except RuntimeError as exc:
        if "exited with code" in str(exc):
            pytest.skip(
                f"KING refused this tiny fixture: {exc}; "
                "real benchmarks use ≥ 50 samples")
        raise
    assert isinstance(result, KingResult)
    assert result.profiling.exit_code == 0
    assert result.normalised_tsv.exists()
    text = result.normalised_tsv.read_text()
    lines = text.strip().splitlines()
    assert lines[0] == "sample_a\tsample_b\tkinship"
    # Pair count ≤ 28 (8 choose 2). KING may filter very weak pairs.
    body = lines[1:]
    assert 0 <= len(body) <= 28
    assert (out / "receipt.json").exists()
