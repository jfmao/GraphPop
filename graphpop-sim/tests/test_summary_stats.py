"""Tests for the summary-stats module."""
from __future__ import annotations

import math

import pytest

from graphpop_sim import MsprimeConfig, MsprimeRunner, compute_summary_stats

pytest.importorskip("msprime")


def test_summary_stats_basic_fields_populate():
    cfg = MsprimeConfig(
        n_diploid=8, sequence_length=20_000,
        recombination_rate=1e-7, mutation_rate=1e-6, seed=11)
    ts = MsprimeRunner(cfg).simulate()
    stats = compute_summary_stats(ts)
    assert stats.n_samples == 16
    assert stats.sequence_length == 20_000
    assert stats.n_segregating_sites > 0
    assert stats.pi > 0
    assert stats.theta_w > 0
    # Tajima's D is finite for non-degenerate sims.
    assert math.isfinite(stats.tajimas_d)


def test_summary_stats_no_mutations_yields_zero_theta_and_nan_d(tmp_path):
    cfg = MsprimeConfig(
        n_diploid=8, sequence_length=10_000,
        recombination_rate=1e-7, mutation_rate=0.0, seed=2)
    ts = MsprimeRunner(cfg).simulate()
    stats = compute_summary_stats(ts)
    assert stats.n_segregating_sites == 0
    assert stats.theta_w == 0.0
    assert math.isnan(stats.tajimas_d)


def test_as_vector_returns_named_subset():
    cfg = MsprimeConfig(
        n_diploid=8, sequence_length=10_000,
        recombination_rate=1e-7, mutation_rate=1e-6, seed=2)
    ts = MsprimeRunner(cfg).simulate()
    stats = compute_summary_stats(ts)
    v = stats.as_vector(["pi", "theta_w"])
    assert v.shape == (2,)


def test_mean_fst_emerges_when_two_populations_provided():
    cfg = MsprimeConfig(
        n_diploid=10, sequence_length=20_000,
        recombination_rate=1e-7, mutation_rate=1e-6, seed=42)
    ts = MsprimeRunner(cfg).simulate()
    samples = list(ts.samples())
    half = len(samples) // 2
    stats = compute_summary_stats(
        ts, populations=[samples[:half], samples[half:]])
    # Mean F_ST should be a finite number (could be 0 for panmixia).
    assert stats.mean_fst is not None
    assert math.isfinite(stats.mean_fst)


def test_mean_fst_is_none_with_singleton_population():
    cfg = MsprimeConfig(
        n_diploid=4, sequence_length=5_000,
        recombination_rate=1e-7, mutation_rate=1e-6, seed=3)
    ts = MsprimeRunner(cfg).simulate()
    samples = list(ts.samples())
    stats = compute_summary_stats(
        ts, populations=[[samples[0]], samples[1:]])
    assert stats.mean_fst is None
