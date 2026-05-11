"""Tests for the PLINK 2.0 GRM wrapper.

Parser + TSV-emission tests run without PLINK. The integration
test (which actually invokes PLINK) is skipped when no PLINK
binary is on PATH — matches the CI-friendly pattern used in
graphpop-sim's SLiM tests and graphpop-pedigree.
"""
from __future__ import annotations

import gzip
import struct
import textwrap
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.competitors import (
    PlinkGrmResult,
    PlinkGrmRunner,
    parse_grm_bin,
)
from graphpop_bench.competitors.plink_grm import (
    _build_normalised_tsv,
    _resolve_input_kind,
)


# ---------------------------------------------------------------------------
# Helpers — build a hand-crafted .grm.bin + .grm.id without PLINK
# ---------------------------------------------------------------------------

def _write_grm_bin(grm: np.ndarray, path: Path) -> None:
    """Write the lower-triangle of ``grm`` as packed float32."""
    n = grm.shape[0]
    flat = []
    for i in range(n):
        flat.extend(float(grm[i, j]) for j in range(i + 1))
    arr = np.asarray(flat, dtype=np.float32)
    arr.tofile(str(path))


def _write_grm_id(sample_ids, path: Path) -> None:
    lines = [f"{sid}\t{sid}" for sid in sample_ids]
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Parser tests (no PLINK required)
# ---------------------------------------------------------------------------

def test_parse_grm_bin_round_trips_a_hand_crafted_matrix(tmp_path):
    grm_true = np.array([
        [1.00, 0.10, 0.05, 0.02, 0.01],
        [0.10, 1.00, 0.15, 0.04, 0.03],
        [0.05, 0.15, 1.00, 0.20, 0.09],
        [0.02, 0.04, 0.20, 1.00, 0.30],
        [0.01, 0.03, 0.09, 0.30, 1.00],
    ], dtype=np.float64)
    sample_ids = ["S0", "S1", "S2", "S3", "S4"]

    _write_grm_bin(grm_true, tmp_path / "out.grm.bin")
    _write_grm_id(sample_ids, tmp_path / "out.grm.id")

    grm, ids = parse_grm_bin(
        grm_bin=tmp_path / "out.grm.bin",
        grm_id=tmp_path / "out.grm.id",
    )
    assert ids == sample_ids
    np.testing.assert_allclose(grm, grm_true, atol=1e-6)


def test_parse_grm_bin_missing_id_file_raises(tmp_path):
    _write_grm_bin(np.eye(3, dtype=np.float64), tmp_path / "out.grm.bin")
    with pytest.raises(FileNotFoundError):
        parse_grm_bin(
            grm_bin=tmp_path / "out.grm.bin",
            grm_id=tmp_path / "nope.grm.id",
        )


def test_parse_grm_bin_size_mismatch_raises(tmp_path):
    _write_grm_id(["S0", "S1", "S2"], tmp_path / "out.grm.id")
    # Wrong number of entries: 3 samples expect 3*4/2 = 6 float32s.
    np.asarray([1.0, 0.1, 1.0], dtype=np.float32).tofile(
        str(tmp_path / "out.grm.bin"))
    with pytest.raises(ValueError, match="size mismatch"):
        parse_grm_bin(
            grm_bin=tmp_path / "out.grm.bin",
            grm_id=tmp_path / "out.grm.id",
        )


def test_parse_grm_bin_single_column_id_file_falls_back(tmp_path):
    # PLINK normally writes "FID IID"; this verifies graceful handling
    # of an unusual one-column file.
    (tmp_path / "out.grm.id").write_text("S0\nS1\n")
    _write_grm_bin(np.array([[1.0, 0.0], [0.0, 1.0]]),
                   tmp_path / "out.grm.bin")
    _, ids = parse_grm_bin(
        grm_bin=tmp_path / "out.grm.bin",
        grm_id=tmp_path / "out.grm.id",
    )
    assert ids == ["S0", "S1"]


def test_parse_grm_bin_empty_id_file_raises(tmp_path):
    (tmp_path / "out.grm.id").write_text("")
    (tmp_path / "out.grm.bin").write_bytes(b"")
    with pytest.raises(ValueError, match="no sample IDs"):
        parse_grm_bin(
            grm_bin=tmp_path / "out.grm.bin",
            grm_id=tmp_path / "out.grm.id",
        )


# ---------------------------------------------------------------------------
# TSV emission tests
# ---------------------------------------------------------------------------

def test_build_normalised_tsv_emits_one_row_per_unordered_pair(tmp_path):
    grm = np.array([
        [1.00, 0.12, 0.05],
        [0.12, 1.00, 0.20],
        [0.05, 0.20, 1.00],
    ], dtype=np.float64)
    out = tmp_path / "plink_grm.tsv"
    _build_normalised_tsv(grm, ["A", "B", "C"], out)
    text = out.read_text()
    lines = text.strip().splitlines()
    assert lines[0] == "sample_a\tsample_b\tkinship"
    # Body rows: AA, AB, AC, BB, BC, CC = 6 rows (n=3 → n*(n+1)/2)
    body = lines[1:]
    assert len(body) == 6
    assert body[0].startswith("A\tA\t")
    assert body[1].startswith("A\tB\t")
    assert body[-1].startswith("C\tC\t")
    # Check one off-diagonal value carries through
    ab_value = float(body[1].split("\t")[2])
    assert ab_value == pytest.approx(0.12)


def test_build_normalised_tsv_creates_parent_dir(tmp_path):
    out = tmp_path / "nested" / "deeper" / "tsv.tsv"
    _build_normalised_tsv(np.array([[1.0]]), ["S0"], out)
    assert out.exists()


# ---------------------------------------------------------------------------
# Input-kind resolution
# ---------------------------------------------------------------------------

def test_resolve_input_kind_auto_detects_vcf():
    assert _resolve_input_kind(Path("cohort.vcf"), "auto") == "vcf"
    assert _resolve_input_kind(Path("cohort.vcf.gz"), "auto") == "vcf"


def test_resolve_input_kind_auto_treats_no_extension_as_bfile():
    assert _resolve_input_kind(Path("cohort"), "auto") == "bfile"


def test_resolve_input_kind_honours_explicit_choice():
    assert _resolve_input_kind(Path("cohort.vcf"), "bfile") == "bfile"
    assert _resolve_input_kind(Path("cohort"), "vcf") == "vcf"


def test_resolve_input_kind_rejects_unknown_kind():
    with pytest.raises(ValueError, match="unknown input_kind"):
        _resolve_input_kind(Path("cohort"), "pgen")


# ---------------------------------------------------------------------------
# is_available + detect_version (no PLINK required)
# ---------------------------------------------------------------------------

def test_is_available_returns_bool():
    # Whether plink is installed or not, must return a bool without throwing.
    assert isinstance(PlinkGrmRunner.is_available(), bool)


def test_detect_version_returns_none_for_nonexistent_binary():
    assert PlinkGrmRunner.detect_version("/no/such/binary") is None


# ---------------------------------------------------------------------------
# Integration: run PLINK on a small VCF (skip if plink is missing)
# ---------------------------------------------------------------------------

# A tiny 5-sample, 4-SNP VCF — small enough that hand-eyeballing
# the resulting GRM is meaningful.
_MIN_VCF = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##contig=<ID=1,length=10000>
    ##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS0\tS1\tS2\tS3\tS4
    1\t100\trs1\tA\tG\t.\tPASS\tPR\tGT\t0/0\t0/1\t1/1\t0/1\t0/0
    1\t500\trs2\tC\tT\t.\tPASS\tPR\tGT\t0/1\t0/1\t1/1\t0/0\t0/0
    1\t1500\trs3\tG\tA\t.\tPASS\tPR\tGT\t1/1\t0/1\t0/0\t0/1\t1/1
    1\t3000\trs4\tT\tC\t.\tPASS\tPR\tGT\t0/0\t1/1\t0/1\t1/1\t0/1
""")


@pytest.mark.skipif(
    not PlinkGrmRunner.is_available(),
    reason="plink2/plink not on PATH; integration test requires PLINK",
)
def test_plink_run_on_tiny_vcf_emits_symmetric_grm_and_tsv(tmp_path):
    """PLINK 2.0 refuses to compute a GRM when there are < 50
    samples and no precomputed allele freqs — too noisy. Pass
    `--bad-freqs` to bypass that guard rail for this smoke test.
    """
    vcf_path = tmp_path / "tiny.vcf"
    vcf_path.write_text(_MIN_VCF)
    out_dir = tmp_path / "plink_out"
    runner = PlinkGrmRunner()
    result = runner.run(
        vcf_path, out_dir, seed=42, extra_args=["--bad-freqs"])
    assert isinstance(result, PlinkGrmResult)
    assert result.profiling.exit_code == 0
    # 5 samples in the VCF.
    assert len(result.sample_ids) == 5
    # Symmetric GRM with positive diagonal.
    assert result.grm.shape == (5, 5)
    np.testing.assert_allclose(result.grm, result.grm.T, atol=1e-6)
    assert np.all(np.diag(result.grm) > 0)
    # Normalised TSV emitted with the expected number of pair rows.
    text = result.normalised_tsv.read_text()
    body = text.strip().splitlines()[1:]  # drop header
    assert len(body) == 15  # 5 * 6 / 2
    # Receipt JSON written.
    assert (out_dir / "receipt.json").exists()
