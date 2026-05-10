"""Tests for the rejection-ABC framework."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from graphpop_sim import (
    AbcResult,
    MsprimeConfig,
    MsprimeRunner,
    PriorSpec,
    SummaryStats,
    compute_summary_stats,
    load_priors,
    posterior_to_tsv,
    run_abc,
)

pytest.importorskip("msprime")


def test_load_priors_from_yaml(tmp_path):
    p = tmp_path / "priors.yaml"
    p.write_text(
        "mutation_rate:\n"
        "  distribution: log_uniform\n"
        "  low: -9\n"
        "  high: -6\n"
        "recombination_rate:\n"
        "  distribution: uniform\n"
        "  low: 1e-9\n"
        "  high: 1e-7\n"
    )
    priors = load_priors(p)
    assert len(priors) == 2
    names = {p.name for p in priors}
    assert names == {"mutation_rate", "recombination_rate"}


def test_prior_sample_log_uniform_returns_in_range():
    import numpy as np

    rng = np.random.default_rng(42)
    p = PriorSpec("mu", "log_uniform", low=-8, high=-6)
    for _ in range(20):
        s = p.sample(rng)
        assert 1e-8 <= s <= 1e-6


def test_prior_unknown_distribution_raises():
    import numpy as np

    rng = np.random.default_rng(0)
    p = PriorSpec("x", "trapezoidal", 0, 1)
    with pytest.raises(ValueError):
        p.sample(rng)


def _fake_simulator_factory(target_pi: float):
    """Return a simulator that produces stats trending toward target_pi
    based on a single param 'tune' in the config-extra space.

    Wired via the seed escape: we treat 'seed' as the tunable.
    """

    def sim(cfg: MsprimeConfig):
        ts = MsprimeRunner(cfg).simulate()
        return ts

    return sim


def test_run_abc_accepts_exactly_epsilon_fraction():
    obs = SummaryStats(pi=1e-3, theta_w=1e-3, tajimas_d=0.0,
                        mean_fst=None, n_samples=16,
                        n_segregating_sites=10,
                        sequence_length=10_000)
    priors = [
        PriorSpec("mutation_rate", "log_uniform", -8, -6),
    ]
    base = MsprimeConfig(n_diploid=4, sequence_length=5_000,
                          recombination_rate=1e-7, seed=42)
    result = run_abc(
        priors=priors, observed=obs, n_sim=20,
        epsilon=0.25, base_config=base, seed=42,
    )
    assert result.n_total == 20
    # 25 % of 20 = 5 accepted draws.
    assert result.n_accepted == 5
    assert len(result.accepted) == 5
    assert len(result.distances) == 5
    # Distances must be monotonically non-decreasing among accepted.
    for a, b in zip(result.distances, result.distances[1:]):
        assert a <= b


def test_posterior_to_tsv_round_trip(tmp_path):
    result = AbcResult(
        accepted=[
            {"mu": 1e-7, "seed": 1},
            {"mu": 2e-7, "seed": 2},
        ],
        distances=[0.1, 0.2],
        n_total=10,
        n_accepted=2,
        epsilon=0.2,
        summary_keys=["pi"],
    )
    path = tmp_path / "posterior.tsv"
    n = posterior_to_tsv(result, path)
    assert n == 2
    lines = path.read_text().splitlines()
    assert lines[0] == "distance\tmu\tseed"
    assert lines[1].startswith("0.1\t")


def test_run_abc_with_at_least_one_accepted_even_at_tiny_epsilon():
    obs = SummaryStats(pi=1e-3, theta_w=1e-3, tajimas_d=0.0,
                        mean_fst=None, n_samples=8,
                        n_segregating_sites=10,
                        sequence_length=2_000)
    priors = [
        PriorSpec("mutation_rate", "log_uniform", -8, -6),
    ]
    base = MsprimeConfig(n_diploid=2, sequence_length=2_000,
                          recombination_rate=1e-7, seed=42)
    result = run_abc(
        priors=priors, observed=obs, n_sim=5,
        epsilon=0.001,
        base_config=base, seed=42,
    )
    # Even with epsilon=0.001 we floor at 1 accepted.
    assert result.n_accepted >= 1
