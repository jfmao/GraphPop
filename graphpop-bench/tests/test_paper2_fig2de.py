"""Unit tests for Paper 2 Fig 2d/2e driver helpers."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.paper2_drivers import fig2de_panels as f2de


# ---------------------------------------------------------------------------
# Demography construction
# ---------------------------------------------------------------------------

def _msprime_available() -> bool:
    try:
        import msprime  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(
    not _msprime_available(),
    reason="msprime required",
)
def test_build_two_pop_demography_has_two_pops():
    params = f2de.TwoPopParams()
    demo = f2de.build_two_pop_demography(params)
    names = {p.name for p in demo.populations}
    assert names == {"AFR", "EUR"}


# ---------------------------------------------------------------------------
# stratify_pairs_by_pop
# ---------------------------------------------------------------------------

def test_stratify_pairs_by_pop_three_buckets():
    labels = ["A", "A", "B", "B"]
    buckets = f2de.stratify_pairs_by_pop(labels)
    # Off-diagonal pairs: AA (one), AB (four), BB (one).
    assert sorted(buckets.keys()) == ["A-A", "A-B", "B-B"]
    assert len(buckets["A-A"]) == 1
    assert len(buckets["A-B"]) == 4
    assert len(buckets["B-B"]) == 1


def test_stratify_pairs_by_pop_diagonal_excluded_by_default():
    """No (i, i) pair appears in any bucket."""
    labels = ["X", "Y", "Z"]
    buckets = f2de.stratify_pairs_by_pop(labels)
    for pairs in buckets.values():
        for i, j in pairs:
            assert i < j


def test_stratify_pairs_by_pop_include_diagonal_flag():
    labels = ["X", "Y"]
    buckets = f2de.stratify_pairs_by_pop(
        labels, include_diagonal=True)
    flat = [p for ps in buckets.values() for p in ps]
    # i ≤ j convention with diagonal → 3 pairs.
    assert len(flat) == 3


def test_stratify_pairs_by_pop_canonical_label_order():
    """Pop-pair label is alphabetically sorted irrespective of (i, j) order."""
    labels = ["Z", "A"]
    buckets = f2de.stratify_pairs_by_pop(labels)
    assert "A-Z" in buckets and "Z-A" not in buckets


# ---------------------------------------------------------------------------
# Integration smoke (gated by msprime + egrm)
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
    reason="msprime + egrm required for fig2de smoke",
)
def test_simulate_draw_two_pop_micro():
    """Small 5+5 diploid micro draw returns symmetric square GRM."""
    params = f2de.TwoPopParams(
        n_diploid_per_pop={"AFR": 5, "EUR": 5},
        sequence_length=2_000)
    grm, pop_labels = f2de.simulate_draw_two_pop(params, seed=42)
    n_hap = 2 * (5 + 5)
    assert grm.shape == (n_hap, n_hap)
    np.testing.assert_allclose(grm, grm.T, atol=1e-9)
    # 10 AFR-haploid + 10 EUR-haploid.
    assert pop_labels.count("AFR") == 10
    assert pop_labels.count("EUR") == 10


@pytest.mark.skipif(
    not _msprime_egrm_available(),
    reason="msprime + egrm required for run_fig2de smoke",
)
def test_run_fig2de_micro(tmp_path):
    """Tiny end-to-end: 3+3 diploid / 2 kb / N=2 posterior, N=3 truth."""
    params = f2de.TwoPopParams(
        n_diploid_per_pop={"AFR": 3, "EUR": 3},
        sequence_length=2_000)
    result = f2de.run_fig2de(
        params=params,
        n_ground_truth=3,
        n_posterior=2,
        seed=2026,
        output_dir=tmp_path,
    )
    assert result.fig2d_csv.exists()
    assert result.fig2e_csv.exists()
    assert result.metadata_path.exists()
    # Fig 2d should have one row per upper-triangle pair (incl diag).
    with open(result.fig2d_csv) as fh:
        rows_d = list(csv.DictReader(fh))
    n_hap = 2 * (3 + 3)
    assert len(rows_d) == n_hap * (n_hap + 1) // 2
    # Pop labels are AFR/EUR.
    assert {r["pop_a"] for r in rows_d} <= {"AFR", "EUR"}
    # Fig 2e summary has 3 buckets (AFR-AFR, AFR-EUR, EUR-EUR).
    with open(result.fig2e_csv) as fh:
        rows_e = list(csv.DictReader(fh))
    labels = {r["pop_pair_label"] for r in rows_e}
    assert labels == {"AFR-AFR", "AFR-EUR", "EUR-EUR"}
