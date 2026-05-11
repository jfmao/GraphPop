"""Unit tests for Paper 2 Fig 3e driver helpers."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.paper2_drivers import fig3e_panels as f3e


# ---------------------------------------------------------------------------
# collapse_haploid_to_diploid_grm
# ---------------------------------------------------------------------------

def test_collapse_haploid_to_diploid_grm_basic():
    """4 haploid → 2 diploid; each diploid pair is sum of 4 hap entries."""
    hap = np.arange(16).reshape(4, 4).astype(float)
    dipl = f3e.collapse_haploid_to_diploid_grm(hap)
    assert dipl.shape == (2, 2)
    # Diploid (0, 0) = sum of hap[0:2, 0:2] = 0+1+4+5 = 10.
    assert dipl[0, 0] == 10.0
    # Diploid (0, 1) = sum of hap[0:2, 2:4] = 2+3+6+7 = 18.
    assert dipl[0, 1] == 18.0
    # Diploid (1, 1) = sum of hap[2:4, 2:4] = 10+11+14+15 = 50.
    assert dipl[1, 1] == 50.0


def test_collapse_haploid_to_diploid_grm_odd_haploids_raise():
    hap = np.zeros((3, 3))
    with pytest.raises(ValueError, match="n_haploid"):
        f3e.collapse_haploid_to_diploid_grm(hap)


# ---------------------------------------------------------------------------
# density_nnz + correlate_off_diag
# ---------------------------------------------------------------------------

def test_density_nnz_all_nonzero():
    grm = np.ones((4, 4))
    np.fill_diagonal(grm, 0)
    # 6 off-diag pairs, all > 1e-6.
    assert f3e.density_nnz(grm, threshold=1e-6) == 1.0


def test_density_nnz_all_zero():
    grm = np.zeros((4, 4))
    assert f3e.density_nnz(grm, threshold=1e-6) == 0.0


def test_density_nnz_partial():
    """Half the off-diagonal entries above threshold."""
    grm = np.zeros((4, 4))
    grm[0, 1] = grm[1, 0] = 0.5
    grm[2, 3] = grm[3, 2] = 0.7
    # 6 off-diag pairs (strict upper), 2 above threshold → 2/6.
    assert f3e.density_nnz(grm, threshold=1e-6) == pytest.approx(2 / 6)


def test_density_nnz_rejects_nonsquare():
    with pytest.raises(ValueError, match="square"):
        f3e.density_nnz(np.zeros((3, 4)))


def test_correlate_off_diag_identical_returns_one():
    rng = np.random.default_rng(0)
    grm = rng.normal(size=(5, 5))
    grm = (grm + grm.T) / 2
    pearson, spearman = f3e.correlate_off_diag(grm, grm)
    assert pearson == pytest.approx(1.0, rel=1e-9)
    assert spearman == pytest.approx(1.0, rel=1e-9)


def test_correlate_off_diag_uncorrelated():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(20, 20))
    a = (a + a.T) / 2
    b = rng.normal(size=(20, 20))
    b = (b + b.T) / 2
    pearson, spearman = f3e.correlate_off_diag(a, b)
    # 190 off-diag pairs; |r| should be small.
    assert abs(pearson) < 0.4


def test_correlate_off_diag_constant_returns_nan():
    a = np.zeros((4, 4))
    a[0, 1] = 1
    a[1, 0] = 1
    b = np.zeros((4, 4))
    b[2, 3] = 2
    b[3, 2] = 2
    # After dropping double-zeros we have 2 non-zero entries and
    # the others were eliminated → too few pairs.
    pearson, spearman = f3e.correlate_off_diag(a, b)
    assert pearson != pearson or pearson is float("nan") or np.isnan(pearson)


# ---------------------------------------------------------------------------
# Smoke
# ---------------------------------------------------------------------------

def _msprime_egrm_plink_available() -> bool:
    try:
        import msprime  # noqa: F401
        from egrm import varGRM_C  # noqa: F401
    except ImportError:
        return False
    from graphpop_bench.competitors import PlinkGrmRunner
    return PlinkGrmRunner.is_available()


@pytest.mark.skipif(
    not _msprime_egrm_plink_available(),
    reason="msprime + egrm + PLINK required for fig3e smoke",
)
def test_compute_branch_pathway_grm_shape(tmp_path):
    """5+5 diploid sim; branch GRM is symmetric square."""
    params = f3e.Fig3ECohortParams(
        afr_diploid=5, eur_diploid=5,
        sequence_length=2_000)
    ts = f3e.simulate_fig3e_cohort(params, seed=42)
    grm = f3e.compute_branch_pathway_grm(ts, list(range(min(5, ts.num_mutations))))
    n_hap = 2 * 10
    assert grm.shape == (n_hap, n_hap)
    np.testing.assert_allclose(grm, grm.T, atol=1e-9)


@pytest.mark.skipif(
    not _msprime_egrm_plink_available(),
    reason="msprime + egrm + PLINK required for fig3e smoke",
)
def test_compute_branch_pathway_grm_empty_pathway(tmp_path):
    """Empty pathway → zero matrix of correct shape."""
    params = f3e.Fig3ECohortParams(
        afr_diploid=5, eur_diploid=5,
        sequence_length=2_000)
    ts = f3e.simulate_fig3e_cohort(params, seed=42)
    grm = f3e.compute_branch_pathway_grm(ts, [])
    assert grm.shape == (20, 20)
    np.testing.assert_array_equal(grm, np.zeros((20, 20)))


@pytest.mark.skipif(
    not _msprime_egrm_plink_available(),
    reason="msprime + egrm + PLINK required for fig3e smoke",
)
def test_run_fig3e_micro(tmp_path):
    """5+5 diploid / 3 kb / 2 pathway sizes / 1 replicate micro smoke."""
    params = f3e.Fig3ECohortParams(
        afr_diploid=5, eur_diploid=5,
        sequence_length=3_000)
    result = f3e.run_fig3e(
        params=params,
        pathway_sizes=[3, 10],
        n_replicates=1,
        threshold=1e-6,
        seed=2026,
        output_dir=tmp_path,
    )
    assert result.panel_csv.exists()
    assert result.summary_path.exists()
    with open(result.panel_csv) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    # branch_nnz_frac should be ≥ plink_nnz_frac at small M.
    rows_by_M = {int(r["pathway_size"]): r for r in rows}
    if 3 in rows_by_M:
        r = rows_by_M[3]
        assert float(r["branch_nnz_frac"]) >= float(r["plink_nnz_frac"]) - 1e-6
