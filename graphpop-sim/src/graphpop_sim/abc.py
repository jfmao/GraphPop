"""Rejection-ABC framework (Beaumont, Zhang & Balding 2002).

Run N simulations under a prior, compute summary statistics from
each, accept the top-ε by Euclidean distance to observed stats,
return the posterior sample. Pure-Python, deterministic given a
seed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import yaml

from .msprime_runner import MsprimeConfig, MsprimeRunner
from .summary_stats import SummaryStats, compute_summary_stats


@dataclass
class PriorSpec:
    """Per-parameter prior distribution.

    distribution : "uniform" or "log_uniform"
    low / high   : bounds (linear scale for uniform; log10 scale for
                    log_uniform).
    """

    name: str
    distribution: str
    low: float
    high: float

    def sample(self, rng: np.random.Generator) -> float:
        if self.distribution == "uniform":
            return float(rng.uniform(self.low, self.high))
        if self.distribution == "log_uniform":
            return float(10.0 ** rng.uniform(self.low, self.high))
        raise ValueError(f"unknown distribution: {self.distribution}")


@dataclass
class AbcResult:
    """ABC posterior sample."""

    accepted: list[dict[str, Any]] = field(default_factory=list)
    distances: list[float] = field(default_factory=list)
    n_total: int = 0
    n_accepted: int = 0
    epsilon: float = 0.05
    summary_keys: list[str] = field(default_factory=list)


def load_priors(path: str | Path) -> list[PriorSpec]:
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}
    return [
        PriorSpec(
            name=k,
            distribution=str(v.get("distribution", "uniform")),
            low=float(v["low"]),
            high=float(v["high"]),
        )
        for k, v in raw.items()
    ]


def _build_config(
    params: dict[str, Any],
    base: MsprimeConfig | None = None,
) -> MsprimeConfig:
    """Build an MsprimeConfig by overlaying ABC-sampled params."""
    cfg = base or MsprimeConfig()
    # Top-level scalar overrides take precedence; demography is the
    # one non-scalar field that callers can pass as a list-of-dicts.
    cfg = MsprimeConfig(
        n_diploid=int(params.get("n_diploid", cfg.n_diploid)),
        sequence_length=int(params.get(
            "sequence_length", cfg.sequence_length)),
        recombination_rate=float(params.get(
            "recombination_rate", cfg.recombination_rate)),
        mutation_rate=float(params.get(
            "mutation_rate", cfg.mutation_rate)),
        demography=params.get("demography", cfg.demography),
        seed=int(params.get("seed", cfg.seed)),
    )
    return cfg


def run_abc(
    priors: list[PriorSpec],
    observed: SummaryStats,
    n_sim: int,
    epsilon: float = 0.05,
    summary_keys: list[str] | None = None,
    base_config: MsprimeConfig | None = None,
    simulator: Callable[[MsprimeConfig], "Any"] | None = None,
    seed: int = 42,
) -> AbcResult:
    """Rejection-ABC.

    priors        : list of PriorSpec — each defines a sampled parameter.
    observed      : SummaryStats from the observed dataset.
    n_sim         : number of prior draws.
    epsilon       : fraction of closest simulations to accept (0 < ε ≤ 1).
    summary_keys  : subset of stats to use as the distance vector
                     (default ["pi", "theta_w", "tajimas_d"]).
    base_config   : msprime defaults to overlay sampled params on.
    simulator     : callable(MsprimeConfig) → tskit.TreeSequence;
                     defaults to MsprimeRunner.simulate. Override for
                     deterministic tests.
    seed          : RNG seed for prior sampling.
    """
    summary_keys = summary_keys or ["pi", "theta_w", "tajimas_d"]
    rng = np.random.default_rng(seed)
    obs_vec = observed.as_vector(summary_keys)

    if simulator is None:
        def simulator(cfg: MsprimeConfig):  # noqa: WPS430
            return MsprimeRunner(cfg).simulate()

    draws: list[dict[str, Any]] = []
    distances: list[float] = []
    for i in range(n_sim):
        params = {p.name: p.sample(rng) for p in priors}
        # Stamp a unique per-draw seed so msprime is deterministic
        # but draws differ.
        params.setdefault("seed", int(rng.integers(0, 2**31 - 1)))
        cfg = _build_config(params, base_config)
        ts = simulator(cfg)
        stats = compute_summary_stats(ts)
        sim_vec = stats.as_vector(summary_keys)
        d = float(np.linalg.norm(sim_vec - obs_vec))
        draws.append(params)
        distances.append(d)

    order = np.argsort(distances)
    n_accepted = max(1, int(round(n_sim * epsilon)))
    accepted_idx = order[:n_accepted]
    return AbcResult(
        accepted=[draws[i] for i in accepted_idx],
        distances=[distances[i] for i in accepted_idx],
        n_total=n_sim,
        n_accepted=n_accepted,
        epsilon=epsilon,
        summary_keys=summary_keys,
    )


def posterior_to_tsv(result: AbcResult, path: str | Path) -> int:
    """Write the accepted posterior sample to a TSV file. Returns
    the number of rows written.
    """
    if not result.accepted:
        return 0
    keys = sorted(result.accepted[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(["distance", *keys]) + "\n")
        for d, params in zip(result.distances, result.accepted):
            row = [f"{d:.6g}"] + [str(params.get(k, "")) for k in keys]
            fh.write("\t".join(row) + "\n")
    return len(result.accepted)
