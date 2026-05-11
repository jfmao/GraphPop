"""Tests for the egrm reference wrapper (Fan et al. 2022).

Loader / TSV-emit / argument-validation unit tests run without
egrm (synthetic .npy + .ids + meta.json fixtures); the
integration test uses the same 20-sample TreeSequence as the
H3 test and skips if egrm is unavailable.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.competitors import (
    EgrmResult,
    EgrmRunner,
)
from graphpop_bench.competitors.egrm import (
    _build_normalised_tsv,
    _load_egrm_files,
    _parse_float_with_inf,
)

_FIXTURE = Path(
    "/mnt/data/GraphPop/graphpop-procedures/src/test/resources/"
    "egrm_fixture_20samples.trees"
)


# ---------------------------------------------------------------------------
# _load_egrm_files — round-trip on synthetic artefacts
# ---------------------------------------------------------------------------

def _write_synthetic(tmp_path: Path, *,
                     egrm: np.ndarray,
                     ids: list[str],
                     vargrm: np.ndarray | None = None,
                     total_mu: float = 12.34) -> None:
    np.save(tmp_path / "egrm.npy", egrm)
    if vargrm is not None:
        np.save(tmp_path / "vargrm.npy", vargrm)
    (tmp_path / "egrm.ids").write_text("\n".join(ids) + "\n")
    (tmp_path / "meta.json").write_text(json.dumps(
        {"total_mu": total_mu, "compute_var": vargrm is not None}))


def test_load_egrm_files_roundtrip_with_var(tmp_path):
    egrm = np.array([[0.5, 0.1, -0.1],
                     [0.1, 0.5, 0.2],
                     [-0.1, 0.2, 0.5]])
    vargrm = np.eye(3) * 0.01
    _write_synthetic(tmp_path, egrm=egrm, vargrm=vargrm,
                     ids=["S0", "S1", "S2"], total_mu=42.5)
    e, v, mu, ids = _load_egrm_files(
        tmp_path / "egrm.npy", tmp_path / "vargrm.npy",
        tmp_path / "egrm.ids", tmp_path / "meta.json",
        compute_var=True)
    np.testing.assert_allclose(e, egrm)
    np.testing.assert_allclose(v, vargrm)
    assert mu == pytest.approx(42.5)
    assert ids == ["S0", "S1", "S2"]


def test_load_egrm_files_roundtrip_without_var(tmp_path):
    egrm = np.eye(2)
    _write_synthetic(tmp_path, egrm=egrm, ids=["A", "B"],
                     total_mu=1.0)
    e, v, _, ids = _load_egrm_files(
        tmp_path / "egrm.npy", tmp_path / "vargrm.npy",
        tmp_path / "egrm.ids", tmp_path / "meta.json",
        compute_var=False)
    np.testing.assert_allclose(e, egrm)
    assert v is None
    assert ids == ["A", "B"]


def test_load_egrm_files_missing_egrm_raises(tmp_path):
    (tmp_path / "egrm.ids").write_text("S0\n")
    (tmp_path / "meta.json").write_text('{"total_mu":0,"compute_var":false}')
    with pytest.raises(FileNotFoundError, match="eGRM file not found"):
        _load_egrm_files(tmp_path / "egrm.npy", tmp_path / "vargrm.npy",
                         tmp_path / "egrm.ids", tmp_path / "meta.json",
                         compute_var=False)


def test_load_egrm_files_missing_ids_raises(tmp_path):
    np.save(tmp_path / "egrm.npy", np.eye(2))
    (tmp_path / "meta.json").write_text('{"total_mu":0,"compute_var":false}')
    with pytest.raises(FileNotFoundError, match="IDs file not found"):
        _load_egrm_files(tmp_path / "egrm.npy", tmp_path / "vargrm.npy",
                         tmp_path / "egrm.ids", tmp_path / "meta.json",
                         compute_var=False)


def test_load_egrm_files_missing_meta_raises(tmp_path):
    np.save(tmp_path / "egrm.npy", np.eye(2))
    (tmp_path / "egrm.ids").write_text("S0\nS1\n")
    with pytest.raises(FileNotFoundError, match="meta file not found"):
        _load_egrm_files(tmp_path / "egrm.npy", tmp_path / "vargrm.npy",
                         tmp_path / "egrm.ids", tmp_path / "meta.json",
                         compute_var=False)


def test_load_egrm_files_var_expected_but_missing(tmp_path):
    _write_synthetic(tmp_path, egrm=np.eye(2), ids=["A", "B"])
    # No vargrm.npy on disk; ask for it.
    with pytest.raises(FileNotFoundError, match="varGRM file expected"):
        _load_egrm_files(tmp_path / "egrm.npy", tmp_path / "vargrm.npy",
                         tmp_path / "egrm.ids", tmp_path / "meta.json",
                         compute_var=True)


def test_load_egrm_files_shape_mismatch(tmp_path):
    _write_synthetic(tmp_path, egrm=np.eye(3),
                     ids=["S0", "S1"], total_mu=0.0)
    with pytest.raises(ValueError, match="mismatches"):
        _load_egrm_files(tmp_path / "egrm.npy", tmp_path / "vargrm.npy",
                         tmp_path / "egrm.ids", tmp_path / "meta.json",
                         compute_var=False)


# ---------------------------------------------------------------------------
# _build_normalised_tsv — schema match with H1 / H3
# ---------------------------------------------------------------------------

def test_build_normalised_tsv_emits_upper_triangle_with_diag(tmp_path):
    egrm = np.array([[0.5, 0.1, -0.1],
                     [0.1, 0.5, 0.2],
                     [-0.1, 0.2, 0.5]])
    out = tmp_path / "egrm.tsv"
    _build_normalised_tsv(egrm, ["S0", "S1", "S2"], out)
    lines = out.read_text().strip().splitlines()
    assert lines[0] == "sample_a\tsample_b\tkinship"
    # n*(n+1)/2 = 6 rows.
    body = lines[1:]
    assert len(body) == 6
    assert body[0].split("\t")[:2] == ["S0", "S0"]
    assert body[1].split("\t")[:2] == ["S0", "S1"]
    assert body[-1].split("\t")[:2] == ["S2", "S2"]


def test_build_normalised_tsv_creates_parent_dir(tmp_path):
    out = tmp_path / "nested" / "egrm.tsv"
    _build_normalised_tsv(np.eye(2), ["A", "B"], out)
    assert out.exists()


# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------

def test_run_rejects_negative_rlim(tmp_path):
    runner = EgrmRunner()
    with pytest.raises(ValueError, match="rlim"):
        runner.run(_FIXTURE, tmp_path, rlim=-1.0)


def test_run_rejects_alim_le_rlim(tmp_path):
    runner = EgrmRunner()
    with pytest.raises(ValueError, match="alim"):
        runner.run(_FIXTURE, tmp_path, rlim=10.0, alim=5.0)


def test_run_rejects_negative_left(tmp_path):
    runner = EgrmRunner()
    with pytest.raises(ValueError, match="left"):
        runner.run(_FIXTURE, tmp_path, left=-1.0)


def test_run_rejects_right_le_left(tmp_path):
    runner = EgrmRunner()
    with pytest.raises(ValueError, match="right"):
        runner.run(_FIXTURE, tmp_path, left=100.0, right=50.0)


def test_parse_float_with_inf():
    import math
    assert _parse_float_with_inf("inf") == math.inf
    assert _parse_float_with_inf("INF") == math.inf
    assert _parse_float_with_inf("infinity") == math.inf
    assert _parse_float_with_inf("3.14") == pytest.approx(3.14)
    assert _parse_float_with_inf("0") == 0.0


# ---------------------------------------------------------------------------
# is_available + detect_version
# ---------------------------------------------------------------------------

def test_is_available_returns_bool():
    assert isinstance(EgrmRunner.is_available(), bool)


def test_detect_version_matches_available():
    version = EgrmRunner.detect_version()
    if EgrmRunner.is_available():
        # egrm 0.1 has no __version__ attribute → may return None.
        # We only require: not crash + consistent type.
        assert version is None or isinstance(version, str)
    else:
        assert version is None


# ---------------------------------------------------------------------------
# Integration: spawn subprocess + compute eGRM on real fixture
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not EgrmRunner.is_available(),
    reason="egrm not importable in this Python environment",
)
@pytest.mark.skipif(
    not _FIXTURE.exists(),
    reason=f"20-sample fixture not at {_FIXTURE}",
)
def test_run_on_20_sample_fixture_with_var(tmp_path):
    """End-to-end with varGRM enabled: assert shape, symmetry,
    double-centring, finite values, TSV emission, receipt fields."""
    runner = EgrmRunner()
    result = runner.run(_FIXTURE, tmp_path, seed=42,
                        graphpop_commit="testsha")
    assert isinstance(result, EgrmResult)
    assert result.profiling.exit_code == 0
    assert result.egrm.shape == (20, 20)
    assert np.allclose(result.egrm, result.egrm.T, atol=1e-9), (
        "eGRM should be symmetric")
    assert np.all(np.isfinite(result.egrm))
    # Double-centring: row sums ≈ 0 (the egrm package subtracts both
    # row mean and column mean before returning).
    np.testing.assert_allclose(result.egrm.sum(axis=0), 0.0, atol=1e-8)
    np.testing.assert_allclose(result.egrm.sum(axis=1), 0.0, atol=1e-8)
    # varGRM: same shape, finite, symmetric.
    assert result.vargrm is not None
    assert result.vargrm.shape == (20, 20)
    assert np.all(np.isfinite(result.vargrm))
    assert np.allclose(result.vargrm, result.vargrm.T, atol=1e-9)
    assert result.total_mu > 0
    # TSV with diagonal.
    assert result.normalised_tsv.exists()
    lines = result.normalised_tsv.read_text().strip().splitlines()
    assert lines[0] == "sample_a\tsample_b\tkinship"
    assert len(lines) - 1 == 20 * 21 // 2
    # Receipt fields.
    receipt = json.loads((tmp_path / "receipt.json").read_text())
    assert receipt["tool"] == "egrm"
    assert receipt["seed"] == 42
    assert receipt["graphpop_commit"] == "testsha"
    assert receipt["compute_var"] is True
    assert receipt["n_samples"] == 20
    assert receipt["total_mu"] == pytest.approx(result.total_mu)


@pytest.mark.skipif(
    not EgrmRunner.is_available(),
    reason="egrm not importable in this Python environment",
)
@pytest.mark.skipif(
    not _FIXTURE.exists(),
    reason=f"20-sample fixture not at {_FIXTURE}",
)
def test_run_compute_var_false(tmp_path):
    """With compute_var=False, vargrm is None but eGRM + TSV still
    appear, and `vargrm.npy` is NOT on disk."""
    runner = EgrmRunner()
    result = runner.run(_FIXTURE, tmp_path, compute_var=False)
    assert result.profiling.exit_code == 0
    assert result.egrm.shape == (20, 20)
    assert result.vargrm is None
    assert not (tmp_path / "vargrm.npy").exists()
    assert result.normalised_tsv.exists()


@pytest.mark.skipif(
    not EgrmRunner.is_available(),
    reason="egrm not importable in this Python environment",
)
def test_run_subprocess_failure_surfaces_runtime_error(tmp_path):
    """Pointing at a non-trees file makes the subprocess exit non-zero."""
    bogus = tmp_path / "not_trees.txt"
    bogus.write_text("not a treeseq")
    runner = EgrmRunner()
    with pytest.raises(RuntimeError, match="exited with code"):
        runner.run(bogus, tmp_path / "out")
