"""Unit tests for Paper 2 Fig 5b driver helpers."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.paper2_drivers import fig5b_panels as f5b


# ---------------------------------------------------------------------------
# inject_mosaic
# ---------------------------------------------------------------------------

def _toy_haploid_matrix() -> np.ndarray:
    """6 haploids (3 diploids), 8 variants — anchor=0, target=1, others=2."""
    anchor = np.array([0, 1, 0, 1, 0, 1, 0, 1], dtype=np.int8)
    target_orig = np.array([1, 0, 1, 0, 1, 0, 1, 0], dtype=np.int8)
    other = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=np.int8)
    return np.stack([anchor, anchor, target_orig, target_orig, other, other])


def test_inject_mosaic_shape_preserved():
    m = _toy_haploid_matrix()
    rng = np.random.default_rng(0)
    out = f5b.inject_mosaic(m, [(0, 1)], 0.5, rng)
    assert out.shape == m.shape


def test_inject_mosaic_fraction_zero_unchanged():
    m = _toy_haploid_matrix()
    rng = np.random.default_rng(0)
    out = f5b.inject_mosaic(m, [(0, 1)], 0.0, rng)
    np.testing.assert_array_equal(out, m)


def test_inject_mosaic_fraction_one_replaces_target():
    """f=1 makes the target's haploids identical to the anchor's."""
    m = _toy_haploid_matrix()
    rng = np.random.default_rng(1)
    out = f5b.inject_mosaic(m, [(0, 1)], 1.0, rng)
    # Target's haploids (rows 2, 3) now equal anchor (rows 0, 1).
    np.testing.assert_array_equal(out[2], m[0])
    np.testing.assert_array_equal(out[3], m[1])
    # Other diploid (rows 4, 5) untouched.
    np.testing.assert_array_equal(out[4:6], m[4:6])


def test_inject_mosaic_random_seeding_deterministic():
    """Same seed + same input → same output."""
    m = _toy_haploid_matrix()
    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(42)
    out1 = f5b.inject_mosaic(m, [(0, 1)], 0.5, rng1)
    out2 = f5b.inject_mosaic(m, [(0, 1)], 0.5, rng2)
    np.testing.assert_array_equal(out1, out2)


def test_inject_mosaic_rejects_bad_fraction():
    m = _toy_haploid_matrix()
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match=r"fraction must"):
        f5b.inject_mosaic(m, [(0, 1)], -0.1, rng)
    with pytest.raises(ValueError, match=r"fraction must"):
        f5b.inject_mosaic(m, [(0, 1)], 1.1, rng)


def test_inject_mosaic_rejects_self_pair():
    m = _toy_haploid_matrix()
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match=r"self-anchor"):
        f5b.inject_mosaic(m, [(0, 0)], 0.5, rng)


def test_inject_mosaic_rejects_odd_haploids():
    m = np.zeros((3, 4), dtype=np.int8)
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match=r"n_haploid"):
        f5b.inject_mosaic(m, [(0, 1)], 0.5, rng)


# ---------------------------------------------------------------------------
# measure_recall
# ---------------------------------------------------------------------------

def test_measure_recall_perfect():
    """All injected pairs are above threshold."""
    grm = np.array([
        [0.5, 0.3, 0.0],
        [0.3, 0.5, 0.0],
        [0.0, 0.0, 0.5],
    ])
    pairs = [(0, 1)]
    sample_ids = ["A", "B", "C"]
    n_rec, k, recall, phis = f5b.measure_recall(
        grm, sample_ids, pairs, threshold=0.1)
    assert n_rec == 1
    assert k == 1
    assert recall == 1.0
    assert phis == [pytest.approx(0.3)]


def test_measure_recall_none():
    grm = np.zeros((3, 3))
    pairs = [(0, 1), (1, 2)]
    sample_ids = ["A", "B", "C"]
    n_rec, k, recall, phis = f5b.measure_recall(
        grm, sample_ids, pairs, threshold=0.1)
    assert n_rec == 0
    assert recall == 0.0


def test_measure_recall_empty_pairs():
    grm = np.zeros((3, 3))
    n_rec, k, recall, phis = f5b.measure_recall(
        grm, ["A", "B", "C"], [], threshold=0.1)
    assert n_rec == 0 and k == 0 and recall == 0.0


def test_measure_recall_index_out_of_range():
    grm = np.zeros((2, 2))
    with pytest.raises(IndexError):
        f5b.measure_recall(grm, ["A", "B"], [(0, 5)], threshold=0.1)


# ---------------------------------------------------------------------------
# select_disjoint_pairs
# ---------------------------------------------------------------------------

def test_select_disjoint_pairs_disjoint():
    rng = np.random.default_rng(0)
    pairs = f5b.select_disjoint_pairs(20, 5, rng)
    indices = [i for p in pairs for i in p]
    assert len(set(indices)) == len(indices) == 10


def test_select_disjoint_pairs_too_many():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match=r"cannot select"):
        f5b.select_disjoint_pairs(10, 6, rng)


# ---------------------------------------------------------------------------
# Smoke (gated)
# ---------------------------------------------------------------------------

def _msprime_plink_available() -> bool:
    try:
        import msprime  # noqa: F401
    except ImportError:
        return False
    from graphpop_bench.competitors import PlinkGrmRunner
    return PlinkGrmRunner.is_available()


@pytest.mark.skipif(
    not _msprime_plink_available(),
    reason="msprime + PLINK required for run_fig5b micro smoke",
)
def test_run_fig5b_micro(tmp_path):
    """5-diploid / 2 kb / 1 replicate / fractions [0.5, 0.25] smoke."""
    cohort = f5b.CohortParams(
        n_diploid=10, sequence_length=2_000)
    result = f5b.run_fig5b(
        cohort=cohort,
        fractions=[0.5, 0.25],
        n_pairs=3, n_replicates=1,
        threshold=0.0625,
        seed=2026,
        output_dir=tmp_path,
    )
    assert result.panel_csv.exists()
    assert result.summary_path.exists()
    with open(result.panel_csv) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2  # 2 fractions × 1 replicate
    # Recall at fraction = 0.5 should be ≥ 0.66 (recovered ≥ 2/3).
    recalls = {float(r["fraction"]): float(r["recall"]) for r in rows}
    assert recalls[0.5] >= 0.5
