"""Unit tests for Paper 2 Fig 2 driver pure-logic helpers.

Covers Welford aggregation, per-entry CI / coverage, MSE, and
panel-CSV emission. No msprime / egrm / mvn dependency — the
simulate_draw integration is gated by package availability and
skipped otherwise.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.paper2_drivers import fig2_panels as panels


# ---------------------------------------------------------------------------
# Welford
# ---------------------------------------------------------------------------

def test_welford_aggregate_mean_matches_numpy():
    rng = np.random.default_rng(0)
    draws = [rng.normal(size=(4, 4)) for _ in range(50)]
    mean, var, n = panels.welford_aggregate(draws)
    np.testing.assert_allclose(mean, np.mean(draws, axis=0),
                               rtol=1e-12, atol=1e-12)
    assert n == 50


def test_welford_aggregate_variance_matches_ddof1():
    rng = np.random.default_rng(1)
    draws = [rng.normal(size=(3, 3)) for _ in range(20)]
    _, var, _ = panels.welford_aggregate(draws)
    np.testing.assert_allclose(var, np.var(draws, axis=0, ddof=1),
                               rtol=1e-11, atol=1e-11)


def test_welford_aggregate_single_draw_variance_is_zero():
    x = np.eye(3)
    mean, var, n = panels.welford_aggregate([x])
    np.testing.assert_array_equal(mean, x)
    np.testing.assert_array_equal(var, np.zeros_like(x))
    assert n == 1


def test_welford_aggregate_empty_raises():
    with pytest.raises(ValueError, match="at least one draw"):
        panels.welford_aggregate([])


def test_welford_aggregate_shape_mismatch_raises():
    draws = [np.eye(3), np.eye(2)]
    with pytest.raises(ValueError, match="shape"):
        panels.welford_aggregate(draws)


# ---------------------------------------------------------------------------
# MSE
# ---------------------------------------------------------------------------

def test_mse_zero_on_identity():
    x = np.eye(5)
    assert panels.mse(x, x) == 0.0


def test_mse_known_value():
    a = np.array([[1.0, 2.0], [3.0, 4.0]])
    b = np.array([[0.0, 0.0], [0.0, 0.0]])
    # mean of squares: (1+4+9+16)/4 = 7.5
    assert panels.mse(a, b) == pytest.approx(7.5)


def test_mse_shape_mismatch_raises():
    with pytest.raises(ValueError, match="shape"):
        panels.mse(np.eye(3), np.eye(2))


# ---------------------------------------------------------------------------
# coverage_at_alpha
# ---------------------------------------------------------------------------

def test_coverage_at_alpha_full_when_intervals_wide():
    """All-truth-in-intervals → coverage = 1."""
    rng = np.random.default_rng(2)
    truth = np.zeros((3, 3))
    # Draws are wide-noise around 0 → CI definitely covers 0.
    draws = [rng.normal(scale=5.0, size=(3, 3)) for _ in range(200)]
    cov = panels.coverage_at_alpha(draws, truth, alpha=0.05)
    assert cov == 1.0


def test_coverage_at_alpha_zero_when_intervals_too_narrow():
    """Posterior far from truth + narrow → coverage ≈ 0."""
    rng = np.random.default_rng(3)
    truth = np.full((3, 3), 100.0)
    # Posterior centred at 0 with tiny SD; CI nowhere near 100.
    draws = [rng.normal(loc=0, scale=0.01, size=(3, 3))
             for _ in range(100)]
    cov = panels.coverage_at_alpha(draws, truth, alpha=0.05)
    assert cov == 0.0


def test_coverage_at_alpha_calibrated_near_nominal():
    """Posterior IS the truth's distribution → coverage ≈ 0.95."""
    rng = np.random.default_rng(4)
    # truth ~ N(0,1) but fixed for the run; posterior draws are
    # independent N(truth, 1) — so the 95% interval covers truth
    # ≈ 95% of the time for the entries that are sampled
    # representatively. The off-diagonal entries are independent
    # under iid sampling, so we can estimate empirical coverage.
    n = 4
    n_runs = 200
    truth = rng.normal(size=(n, n))
    draws = [rng.normal(loc=truth, scale=1.0, size=(n, n))
             for _ in range(n_runs)]
    cov = panels.coverage_at_alpha(draws, truth, alpha=0.05)
    # Tolerate a coarse-grained tolerance; expected ≈ 0.95.
    assert 0.88 <= cov <= 1.0


def test_coverage_at_alpha_invalid_alpha_raises():
    with pytest.raises(ValueError, match="alpha"):
        panels.coverage_at_alpha([np.eye(2), np.eye(2)],
                                 np.eye(2), alpha=1.5)


def test_coverage_at_alpha_too_few_draws_raises():
    with pytest.raises(ValueError, match="≥ 2"):
        panels.coverage_at_alpha([np.eye(2)], np.eye(2))


# ---------------------------------------------------------------------------
# bootstrap_ci + coverage_from_ci
# ---------------------------------------------------------------------------

def test_bootstrap_ci_shape():
    rng = np.random.default_rng(5)
    map_est = rng.normal(size=(4, 4))
    lo, hi = panels.bootstrap_ci(map_est, n_bootstrap=50, rng=rng)
    assert lo.shape == map_est.shape
    assert hi.shape == map_est.shape
    assert (hi >= lo).all()


def test_coverage_from_ci_basic():
    lo = np.zeros((3, 3))
    hi = np.full((3, 3), 1.0)
    truth = np.full((3, 3), 0.5)
    assert panels.coverage_from_ci(lo, hi, truth) == 1.0
    truth_outside = np.full((3, 3), 2.0)
    assert panels.coverage_from_ci(lo, hi, truth_outside) == 0.0


# ---------------------------------------------------------------------------
# Pair selection
# ---------------------------------------------------------------------------

def test_select_fig2a_pairs_returns_k_distinct():
    rng = np.random.default_rng(6)
    truth = rng.normal(size=(8, 8))
    pairs = panels._select_fig2a_pairs(truth, k=5)
    assert len(pairs) == 5
    assert len(set(pairs)) == 5
    for (a, b) in pairs:
        assert a < b


def test_select_fig2a_pairs_too_small_returns_empty():
    truth = np.array([[1.0]])
    assert panels._select_fig2a_pairs(truth, k=5) == []


# ---------------------------------------------------------------------------
# Panel CSV emitters
# ---------------------------------------------------------------------------

def test_write_fig2a_csv_schema(tmp_path):
    rng = np.random.default_rng(7)
    draws = [rng.normal(size=(3, 3)) for _ in range(4)]
    pairs = [(0, 1), (1, 2)]
    p = tmp_path / "fig2a.csv"
    panels.write_fig2a_csv(draws, pairs, p)
    with open(p) as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
    # 2 pairs × 4 draws = 8 rows.
    assert len(rows) == 8
    assert set(rows[0].keys()) == {
        "pair_label", "sample_a", "sample_b", "draw_idx", "value"}


def test_write_fig2b_csv_schema(tmp_path):
    rows_in = [
        {"n": 5, "mse_posterior": 0.1,
         "mse_map": 0.5, "mse_bootstrap": 0.5},
        {"n": 50, "mse_posterior": 0.01,
         "mse_map": 0.5, "mse_bootstrap": 0.5},
    ]
    p = tmp_path / "fig2b.csv"
    panels.write_fig2b_csv(rows_in, p)
    with open(p) as fh:
        rows = list(csv.DictReader(fh))
    assert [int(r["n"]) for r in rows] == [5, 50]
    assert float(rows[1]["mse_posterior"]) == pytest.approx(0.01)


def test_write_fig2c_csv_schema(tmp_path):
    rows_in = [
        {"n": 5, "coverage_posterior": 0.80,
         "coverage_bootstrap": 0.70},
        {"n": 50, "coverage_posterior": 0.95,
         "coverage_bootstrap": 0.75},
    ]
    p = tmp_path / "fig2c.csv"
    panels.write_fig2c_csv(rows_in, p)
    with open(p) as fh:
        rows = list(csv.DictReader(fh))
    assert float(rows[1]["coverage_posterior"]) == pytest.approx(0.95)


# ---------------------------------------------------------------------------
# simulate_draw — skip-if-missing integration
# ---------------------------------------------------------------------------

def _msprime_egrm_available() -> bool:
    try:
        import msprime  # noqa: F401
        from egrm import varGRM_C  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(
    not _msprime_egrm_available(),
    reason="msprime + egrm required for simulate_draw",
)
def test_simulate_draw_shape_and_symmetry():
    cohort = panels.CohortParams(n_diploid=10, sequence_length=5_000)
    m = panels.simulate_draw(cohort, seed=42)
    assert m.shape == (cohort.n_haploid, cohort.n_haploid)
    np.testing.assert_allclose(m, m.T, atol=1e-12)


# ---------------------------------------------------------------------------
# run_fig2 end-to-end (tiny config)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not _msprime_egrm_available(),
    reason="msprime + egrm required for run_fig2",
)
def test_run_fig2_emits_all_panels(tmp_path):
    """Minimal end-to-end: 5-diploid cohort, 8 ground-truth draws,
    sweep N in {2, 4, 8}. ~seconds of compute. Verifies that the
    panels emit valid CSVs with the right schemas."""
    cohort = panels.CohortParams(
        n_diploid=5, sequence_length=2_000)
    result = panels.run_fig2(
        cohort=cohort,
        output_dir=tmp_path,
        n_ground_truth=8,
        sweep_n_values=[2, 4, 8],
        n_bootstrap=20,
        seed=2026,
    )
    assert result.fig2a_csv.exists()
    assert result.fig2b_csv.exists()
    assert result.fig2c_csv.exists()
    assert result.metadata_path.exists()
    # MSE should drop monotonically (or at least at the largest N
    # be ≤ N=2 case).
    with open(result.fig2b_csv) as fh:
        rows = list(csv.DictReader(fh))
    mses = [float(r["mse_posterior"]) for r in rows]
    assert mses[-1] <= mses[0]
