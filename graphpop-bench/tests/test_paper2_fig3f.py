"""Unit tests for Paper 2 Fig 3f driver helpers."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.paper2_drivers import fig3f_panels as f3f


# ---------------------------------------------------------------------------
# cross_pop_within_tolerance
# ---------------------------------------------------------------------------

def test_cross_pop_within_tolerance_within():
    h2 = {"AFR": [0.50, 0.51], "EUR": [0.48, 0.52], "EAS": [0.49]}
    within, mean, per_pop = f3f.cross_pop_within_tolerance(h2, 0.05)
    assert within is True
    # cross_pop_mean = mean of per-pop means: (0.505 + 0.50 + 0.49) / 3.
    assert mean == pytest.approx((0.505 + 0.50 + 0.49) / 3)
    assert per_pop["EAS"] == pytest.approx(0.49)


def test_cross_pop_within_tolerance_violation():
    h2 = {"AFR": [0.50, 0.51], "EUR": [0.40, 0.41], "EAS": [0.49]}
    within, _, _ = f3f.cross_pop_within_tolerance(h2, 0.05)
    assert within is False


def test_cross_pop_within_tolerance_empty():
    within, mean, per_pop = f3f.cross_pop_within_tolerance({}, 0.05)
    assert within is False
    assert np.isnan(mean)
    assert per_pop == {}


def test_cross_pop_within_tolerance_only_empty_lists():
    h2 = {"AFR": [], "EUR": []}
    within, mean, per_pop = f3f.cross_pop_within_tolerance(h2, 0.05)
    # All per_pop means are NaN → not within tolerance.
    assert within is False


# ---------------------------------------------------------------------------
# Demography + simplify (gated on msprime)
# ---------------------------------------------------------------------------

def _msprime_available() -> bool:
    try:
        import msprime  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(not _msprime_available(), reason="msprime required")
def test_build_three_pop_demography_has_three_pops():
    params = f3f.ThreePopParams()
    demo = f3f.build_three_pop_demography(params)
    names = {p.name for p in demo.populations}
    assert names == {"AFR", "EUR", "EAS"}


@pytest.mark.skipif(not _msprime_available(), reason="msprime required")
def test_simulate_three_pop_cohort_pop_label_counts():
    params = f3f.ThreePopParams(
        n_diploid_per_pop={"AFR": 5, "EUR": 5, "EAS": 5},
        sequence_length=2_000)
    ts, pop_labels = f3f.simulate_three_pop_cohort(params, seed=42)
    assert len(pop_labels) == ts.num_samples
    # 2x diploid per pop.
    counts = {pop: pop_labels.count(pop)
              for pop in {"AFR", "EUR", "EAS"}}
    assert counts == {"AFR": 10, "EUR": 10, "EAS": 10}


@pytest.mark.skipif(not _msprime_available(), reason="msprime required")
def test_simplify_per_pop_returns_three_subsets():
    params = f3f.ThreePopParams(
        n_diploid_per_pop={"AFR": 5, "EUR": 5, "EAS": 5},
        sequence_length=2_000)
    ts, pop_labels = f3f.simulate_three_pop_cohort(params, seed=42)
    per_pop = f3f.simplify_per_pop(ts, pop_labels)
    assert set(per_pop.keys()) == {"AFR", "EUR", "EAS"}
    for pop, (ts_pop, _, _) in per_pop.items():
        assert ts_pop.num_samples == 10


@pytest.mark.skipif(not _msprime_available(), reason="msprime required")
def test_remap_child_set_drops_unmapped_nodes():
    """Synthetic mapping with -1 = dropped + valid id translation."""
    import tskit
    # node_map maps original node id → new node id (or tskit.NULL = -1).
    node_map = np.array([0, 1, tskit.NULL, 2, tskit.NULL, 3])
    children = {0, 2, 3, 4, 5}
    out = f3f.remap_child_set(children, node_map)
    # 0 → 0; 2 → -1 (dropped); 3 → 2; 4 → -1 (dropped); 5 → 3.
    assert out == {0, 2, 3}


# ---------------------------------------------------------------------------
# Smoke (gated)
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
    reason="msprime + egrm required for fig3f smoke",
)
def test_run_fig3f_micro(tmp_path):
    """3+3+3 diploid / 3 kb / 1 ARG / 2 pheno replicates micro."""
    params = f3f.ThreePopParams(
        n_diploid_per_pop={"AFR": 3, "EUR": 3, "EAS": 3},
        sequence_length=3_000)
    result = f3f.run_fig3f(
        params=params,
        n_arg_draws=1,
        n_pheno_replicates=2,
        m_pathway=5,
        true_h2=0.5,
        tolerance=0.05,
        seed=2026,
        output_dir=tmp_path,
    )
    assert result.panel_csv.exists()
    with open(result.panel_csv) as fh:
        rows = list(csv.DictReader(fh))
    # 3 pops × 1 ARG × 2 phenos = 6 rows.
    assert len(rows) == 6
    pops = {r["population"] for r in rows}
    assert pops == {"AFR", "EUR", "EAS"}
