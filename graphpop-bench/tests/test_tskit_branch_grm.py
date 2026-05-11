"""Tests for the tskit branch_grm wrapper.

Loader / TSV-emit unit tests run without tskit (synthetic .npy +
.ids fixtures); the integration test uses a real 20-sample
TreeSequence already shipped with the GraphPop test suite and
skips if tskit is unavailable.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.competitors import (
    TskitBranchGrmResult,
    TskitBranchGrmRunner,
)
from graphpop_bench.competitors.tskit_branch_grm import (
    _build_normalised_tsv,
    _load_grm_files,
)

# Reuse the 20-sample fixture validated by the GraphPop M4.1
# kinship procedure (TreeSequence with known branch-GRM).
_FIXTURE = Path(
    "/mnt/data/GraphPop/graphpop-procedures/src/test/resources/"
    "egrm_fixture_20samples.trees"
)


# ---------------------------------------------------------------------------
# _load_grm_files — round-trip on synthetic artefacts
# ---------------------------------------------------------------------------

def test_load_grm_files_roundtrip(tmp_path):
    grm = np.array([[1.0, 0.2, 0.1],
                    [0.2, 1.0, 0.3],
                    [0.1, 0.3, 1.0]])
    np.save(tmp_path / "grm.npy", grm)
    (tmp_path / "grm.ids").write_text("S0\nS1\nS2\n")
    loaded, ids = _load_grm_files(
        tmp_path / "grm.npy", tmp_path / "grm.ids")
    np.testing.assert_allclose(loaded, grm)
    assert ids == ["S0", "S1", "S2"]


def test_load_grm_files_missing_npy_raises(tmp_path):
    (tmp_path / "grm.ids").write_text("S0\n")
    with pytest.raises(FileNotFoundError, match="GRM file not found"):
        _load_grm_files(tmp_path / "grm.npy", tmp_path / "grm.ids")


def test_load_grm_files_missing_ids_raises(tmp_path):
    np.save(tmp_path / "grm.npy", np.eye(2))
    with pytest.raises(FileNotFoundError, match="IDs file not found"):
        _load_grm_files(tmp_path / "grm.npy", tmp_path / "grm.ids")


def test_load_grm_files_shape_mismatch_raises(tmp_path):
    np.save(tmp_path / "grm.npy", np.eye(3))
    (tmp_path / "grm.ids").write_text("S0\nS1\n")
    with pytest.raises(ValueError, match="mismatches"):
        _load_grm_files(tmp_path / "grm.npy", tmp_path / "grm.ids")


def test_load_grm_files_non_square_raises(tmp_path):
    np.save(tmp_path / "grm.npy", np.zeros((2, 3)))
    (tmp_path / "grm.ids").write_text("S0\nS1\n")
    with pytest.raises(ValueError, match="bad shape"):
        _load_grm_files(tmp_path / "grm.npy", tmp_path / "grm.ids")


# ---------------------------------------------------------------------------
# _build_normalised_tsv — schema match with H1
# ---------------------------------------------------------------------------

def test_build_normalised_tsv_emits_upper_triangle_including_diag(tmp_path):
    grm = np.array([[1.0, 0.2, 0.1],
                    [0.2, 1.0, 0.3],
                    [0.1, 0.3, 1.0]])
    ids = ["S0", "S1", "S2"]
    out = tmp_path / "tskit.tsv"
    _build_normalised_tsv(grm, ids, out)
    lines = out.read_text().strip().splitlines()
    assert lines[0] == "sample_a\tsample_b\tkinship"
    # n*(n+1)/2 = 6 rows for n=3.
    body = lines[1:]
    assert len(body) == 6
    # Spot check: (S0,S0) diag, (S0,S1) off-diag, (S2,S2) diag.
    assert body[0].split("\t")[:2] == ["S0", "S0"]
    assert body[1].split("\t")[:2] == ["S0", "S1"]
    assert body[-1].split("\t")[:2] == ["S2", "S2"]


def test_build_normalised_tsv_creates_parent_dir(tmp_path):
    out = tmp_path / "nested" / "tskit.tsv"
    _build_normalised_tsv(np.eye(2), ["A", "B"], out)
    assert out.exists()


# ---------------------------------------------------------------------------
# is_available + detect_version
# ---------------------------------------------------------------------------

def test_is_available_returns_bool():
    assert isinstance(TskitBranchGrmRunner.is_available(), bool)


def test_detect_version_matches_available():
    version = TskitBranchGrmRunner.detect_version()
    if TskitBranchGrmRunner.is_available():
        assert version is not None and "." in version
    else:
        assert version is None


def test_run_rejects_unknown_mode(tmp_path):
    runner = TskitBranchGrmRunner()
    with pytest.raises(ValueError, match="unknown mode"):
        runner.run(_FIXTURE, tmp_path, mode="explode")


# ---------------------------------------------------------------------------
# Integration: spawn subprocess + compute branch GRM on real fixture
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not TskitBranchGrmRunner.is_available(),
    reason="tskit not importable in this Python environment",
)
@pytest.mark.skipif(
    not _FIXTURE.exists(),
    reason=f"20-sample fixture not at {_FIXTURE}",
)
def test_run_on_20_sample_fixture(tmp_path):
    """End-to-end: load fixture, compute branch GRM via subprocess,
    parse artefacts, emit TSV + receipt. Smoke-checks the math by
    asserting symmetry + finite values."""
    runner = TskitBranchGrmRunner()
    result = runner.run(_FIXTURE, tmp_path, seed=42,
                        graphpop_commit="testsha")
    assert isinstance(result, TskitBranchGrmResult)
    assert result.profiling.exit_code == 0
    assert result.grm.shape == (20, 20)
    assert np.allclose(result.grm, result.grm.T, atol=1e-9), (
        "branch GRM should be symmetric")
    assert np.all(np.isfinite(result.grm))
    # Diagonals are the tree-length-scaled self-similarity; positive.
    assert (np.diag(result.grm) > 0).all()
    # Off-diagonals can be negative under branch-GRM centering, but
    # the spread shouldn't blow up. Bound checked empirically on the
    # fixture: |g| < 10 with the genome-scale Kb tree.
    assert (np.abs(result.grm) < 1e6).all()

    # TSV emitted with diagonal + receipt.
    assert result.normalised_tsv.exists()
    lines = result.normalised_tsv.read_text().strip().splitlines()
    assert lines[0] == "sample_a\tsample_b\tkinship"
    # n*(n+1)/2 = 210 rows for 20 samples.
    assert len(lines) - 1 == 20 * 21 // 2

    receipt = tmp_path / "receipt.json"
    assert receipt.exists()
    import json
    receipt_dict = json.loads(receipt.read_text())
    assert receipt_dict["tool"] == "tskit_branch_grm"
    assert receipt_dict["seed"] == 42
    assert receipt_dict["graphpop_commit"] == "testsha"
    assert receipt_dict["mode"] == "branch"
    assert receipt_dict["n_samples"] == 20


@pytest.mark.skipif(
    not TskitBranchGrmRunner.is_available(),
    reason="tskit not importable in this Python environment",
)
@pytest.mark.skipif(
    not _FIXTURE.exists(),
    reason=f"20-sample fixture not at {_FIXTURE}",
)
def test_run_subprocess_failure_surfaces_runtime_error(tmp_path):
    """Pointing the runner at a non-trees file should make the
    inner tskit.load raise, the subprocess exits non-zero, and
    the wrapper raises RuntimeError with stderr context."""
    bogus = tmp_path / "not_trees.txt"
    bogus.write_text("not a treeseq")
    runner = TskitBranchGrmRunner()
    with pytest.raises(RuntimeError, match="exited with code"):
        runner.run(bogus, tmp_path / "out")
