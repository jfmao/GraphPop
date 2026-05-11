"""Paper 2 Fig 2 driver — posterior branch GRM (G1).

Demonstrates the M4.B Welford aggregator's three headline
claims (per `benchmark_plan.md` § 2 Fig 2 / `PLAN_fig2.md`):

1. **Fig 2a** — per-pair posterior distribution (violin source).
2. **Fig 2b** — aggregated-posterior MSE drops below MAP MSE
   as N posterior samples grows.
3. **Fig 2c** — per-entry 95% CI from the aggregated posterior
   has good empirical coverage of the true value.

**v1 strategy** (Phase 1, simulation-only): N independent
msprime draws serve as proxy SINGER posterior samples. This is
the same proxy used by the existing
`PosteriorBranchGrmProcedureTest` unit test, extended to many
more N. The aggregation math is identical regardless of posterior
source; real SINGER integration is a Phase-2 follow-up.

The driver is Python-only: it uses `egrm.varGRM_C` directly
because the M4.B procedure is independently unit-tested to match
egrm to <1e-6 rel-err (BranchGrmProcedureTest), so the panel
content is unchanged whether the per-draw GRM comes through
Neo4j or directly from egrm. Skipping the round-trip keeps the
wall-clock budget tight.

Pure-logic helpers (Welford, coverage, MSE) are unit-tested in
`tests/test_paper2_fig2.py`.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np

PairKey = Tuple[int, int]


# ---------------------------------------------------------------------------
# Pure-logic helpers — Welford / coverage / MSE
# ---------------------------------------------------------------------------

def welford_aggregate(
    draws: Sequence[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, int]:
    """Streaming mean + sample variance across a sequence of matrices.

    Returns ``(mean, var, n)`` where ``var`` is the sample variance
    (Bessel-corrected, ddof=1) per entry. Matches the M4.B
    procedure's accumulator semantics (per `kinship.branch_grm_posterior`).

    For N=1 the variance is 0 (degenerate; no Bessel correction).
    """
    if len(draws) == 0:
        raise ValueError("welford_aggregate requires at least one draw")
    shape = draws[0].shape
    mean = np.zeros(shape, dtype=np.float64)
    m2 = np.zeros(shape, dtype=np.float64)
    n = 0
    for x in draws:
        if x.shape != shape:
            raise ValueError(
                f"draw shape {x.shape} != initial {shape}")
        n += 1
        delta = x - mean
        mean += delta / n
        delta2 = x - mean
        m2 += delta * delta2
    var = (m2 / (n - 1)) if n > 1 else np.zeros_like(mean)
    return mean, var, n


def mse(estimate: np.ndarray, truth: np.ndarray) -> float:
    """Mean squared error over all matrix entries (no symmetry tricks)."""
    if estimate.shape != truth.shape:
        raise ValueError(
            f"estimate shape {estimate.shape} != truth {truth.shape}")
    return float(np.mean((estimate - truth) ** 2))


def coverage_at_alpha(
    draws: Sequence[np.ndarray], truth: np.ndarray,
    alpha: float = 0.05,
) -> float:
    """Empirical (1-α) CI coverage from a posterior batch.

    For each entry (i, j), build the per-entry (1-α) CI from the
    posterior quantiles (α/2, 1-α/2), and return the fraction of
    entries whose `truth[i, j]` lies inside.

    Uses the upper triangle (i ≤ j) — diagonal included.
    """
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must be in (0, 1), got {alpha!r}")
    if len(draws) < 2:
        raise ValueError(
            "coverage_at_alpha needs ≥ 2 posterior draws")
    stack = np.stack(draws, axis=0)
    lo = np.quantile(stack, alpha / 2, axis=0)
    hi = np.quantile(stack, 1 - alpha / 2, axis=0)
    n = truth.shape[0]
    iu = np.triu_indices(n)
    inside = (truth >= lo) & (truth <= hi)
    return float(np.mean(inside[iu]))


def per_entry_ci(
    draws: Sequence[np.ndarray], alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (lower, upper) (1-α) CI bounds per entry."""
    if not 0 < alpha < 1:
        raise ValueError(f"alpha must be in (0, 1), got {alpha!r}")
    stack = np.stack(draws, axis=0)
    return (
        np.quantile(stack, alpha / 2, axis=0),
        np.quantile(stack, 1 - alpha / 2, axis=0),
    )


def bootstrap_ci(
    map_estimate: np.ndarray,
    *, n_bootstrap: int = 100, alpha: float = 0.05,
    rng: np.random.Generator | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-entry naive bootstrap CI from a single MAP estimate.

    This is the calibration-mismatch baseline for Fig 2c — a
    naive parametric perturbation of the MAP that does NOT
    integrate over ARG uncertainty. We add iid normal noise
    with per-entry scale ≈ (max - min) / 100 to the MAP and
    treat the resulting samples as a bootstrap posterior. The
    benchmark plan predicts this gives coverage < 0.85 (lower
    than the true posterior).

    The wrapper exists so Fig 2c has a meaningful "uncalibrated"
    comparator line. Real bootstrap-on-mutations would resample
    edges of the inferred ARG — that's Phase 2 work.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    scale = (map_estimate.max() - map_estimate.min()) / 100.0
    bs = rng.normal(
        loc=map_estimate,
        scale=scale,
        size=(n_bootstrap,) + map_estimate.shape,
    )
    lo = np.quantile(bs, alpha / 2, axis=0)
    hi = np.quantile(bs, 1 - alpha / 2, axis=0)
    return lo, hi


def coverage_from_ci(
    lo: np.ndarray, hi: np.ndarray, truth: np.ndarray,
) -> float:
    """Fraction of upper-triangle entries inside [lo, hi]."""
    n = truth.shape[0]
    iu = np.triu_indices(n)
    inside = (truth >= lo) & (truth <= hi)
    return float(np.mean(inside[iu]))


# ---------------------------------------------------------------------------
# msprime + egrm — simulate one proxy posterior draw
# ---------------------------------------------------------------------------

@dataclass
class CohortParams:
    n_diploid: int = 50
    sequence_length: int = 50_000
    recomb_rate: float = 1e-5
    mut_rate: float = 1e-4
    population_size: int = 10_000

    @property
    def n_haploid(self) -> int:
        return self.n_diploid * 2


def simulate_draw(cohort: CohortParams, seed: int) -> np.ndarray:
    """Run one msprime sim + egrm.varGRM_C, return (n_haploid, n_haploid)
    eGRM matrix. Each independent draw is a proxy posterior sample."""
    import msprime  # local import → driver importable without msprime
    from egrm import varGRM_C

    ts = msprime.sim_ancestry(
        samples=cohort.n_diploid,
        sequence_length=cohort.sequence_length,
        recombination_rate=cohort.recomb_rate,
        population_size=cohort.population_size,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(
        ts, rate=cohort.mut_rate, random_seed=seed)
    egrm, _vargrm, _total_mu = varGRM_C(ts, var=False)
    return np.asarray(egrm, dtype=np.float64)


# ---------------------------------------------------------------------------
# Panel emitters
# ---------------------------------------------------------------------------

def write_fig2a_csv(
    draws: Sequence[np.ndarray],
    pair_indices: Sequence[PairKey],
    path: Path,
) -> None:
    """Per-pair posterior values for violin plot.

    Columns: pair_label, draw_idx, value.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pair_label", "sample_a", "sample_b",
                    "draw_idx", "value"])
        for (a, b) in pair_indices:
            label = f"({a},{b})"
            for k, d in enumerate(draws):
                w.writerow([label, a, b, k, f"{d[a, b]:.10g}"])


def write_fig2b_csv(
    rows: List[dict], path: Path,
) -> None:
    """MSE-vs-N panel data.

    rows: list of dicts with keys `n`, `mse_posterior`, `mse_map`,
    `mse_bootstrap_baseline`.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "mse_posterior", "mse_map", "mse_bootstrap"])
        for r in rows:
            w.writerow([
                r["n"],
                f"{r['mse_posterior']:.8g}",
                f"{r['mse_map']:.8g}",
                f"{r['mse_bootstrap']:.8g}",
            ])


def write_fig2c_csv(
    rows: List[dict], path: Path,
) -> None:
    """Coverage-vs-N panel data.

    rows: list of dicts with keys `n`, `coverage_posterior`,
    `coverage_bootstrap`.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["n", "coverage_posterior", "coverage_bootstrap"])
        for r in rows:
            w.writerow([
                r["n"],
                f"{r['coverage_posterior']:.6g}",
                f"{r['coverage_bootstrap']:.6g}",
            ])


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig2RunResult:
    cohort: CohortParams
    sweep_n_values: List[int]
    output_dir: Path
    n_ground_truth: int
    truth: np.ndarray = field(repr=False)
    fig2a_csv: Path = field(default=Path())
    fig2b_csv: Path = field(default=Path())
    fig2c_csv: Path = field(default=Path())
    metadata_path: Path = field(default=Path())


def run_fig2(
    *,
    cohort: CohortParams,
    output_dir: Path,
    n_ground_truth: int,
    sweep_n_values: Sequence[int],
    n_bootstrap: int = 100,
    seed: int = 2026,
    progress_callback=None,
) -> Fig2RunResult:
    """Run the full Fig 2 v1 pipeline end-to-end.

    1. Draw `n_ground_truth` proxy posterior samples; average to
       form the "true" matrix.
    2. Sweep N ∈ `sweep_n_values`; for each N, take the first N
       draws, compute aggregated posterior + MAP + bootstrap-CI
       baselines, score MSE + coverage.
    3. Emit fig2a/b/c panel CSVs + metadata.json.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    max_n = max(sweep_n_values)
    # Use DISJOINT pools for truth and sweep so the panels don't
    # collapse to MSE = 0 at N == n_ground_truth (the same pool
    # would make posterior ≡ truth at the high end of the sweep).
    # Extra +1 draw is the MAP baseline (also disjoint).
    total_draws = n_ground_truth + max_n + 1
    seeds = rng.integers(1, 2**31 - 1, size=total_draws).tolist()

    draws: List[np.ndarray] = []
    for i, draw_seed in enumerate(seeds):
        draws.append(simulate_draw(cohort, int(draw_seed)))
        if progress_callback is not None:
            progress_callback(i + 1, total_draws)

    # Pool slicing:
    #   ground truth pool : draws[0 : n_ground_truth]
    #   posterior pool    : draws[n_ground_truth : n_ground_truth + max_n]
    #   MAP single draw   : draws[-1]  (disjoint from both pools)
    ground_truth_pool = draws[:n_ground_truth]
    posterior_pool = draws[n_ground_truth : n_ground_truth + max_n]
    truth = np.mean(ground_truth_pool, axis=0)

    map_estimate = draws[-1]
    map_mse = mse(map_estimate, truth)

    # Bootstrap baseline (uses the same map_estimate).
    bs_lo, bs_hi = bootstrap_ci(
        map_estimate, n_bootstrap=n_bootstrap, rng=rng)
    bs_cov_global = coverage_from_ci(bs_lo, bs_hi, truth)

    # Sweep N over the posterior pool only.
    mse_rows: List[dict] = []
    cov_rows: List[dict] = []
    for N in sorted(sweep_n_values):
        sub = posterior_pool[:N]
        mean_N, _var_N, _ = welford_aggregate(sub)
        mse_post = mse(mean_N, truth)
        if N >= 2:
            cov_post = coverage_at_alpha(sub, truth, alpha=0.05)
        else:
            cov_post = float("nan")
        mse_rows.append({
            "n": N,
            "mse_posterior": mse_post,
            "mse_map": map_mse,
            "mse_bootstrap": map_mse,  # bootstrap mean ≡ MAP
        })
        cov_rows.append({
            "n": N,
            "coverage_posterior": cov_post,
            "coverage_bootstrap": bs_cov_global,  # baseline only
        })

    # Pair selection for Fig 2a: pick 12 pairs across the kinship range
    # of the ground truth (upper triangle, off-diag).
    pair_indices = _select_fig2a_pairs(truth, k=12)

    fig2a_csv = output_dir / "fig2a_pair_posteriors.csv"
    fig2b_csv = output_dir / "fig2b_mse_vs_n.csv"
    fig2c_csv = output_dir / "fig2c_coverage_vs_n.csv"
    write_fig2a_csv(posterior_pool, pair_indices, fig2a_csv)
    write_fig2b_csv(mse_rows, fig2b_csv)
    write_fig2c_csv(cov_rows, fig2c_csv)

    metadata = {
        "cohort": {
            "n_diploid": cohort.n_diploid,
            "n_haploid": cohort.n_haploid,
            "sequence_length": cohort.sequence_length,
            "recomb_rate": cohort.recomb_rate,
            "mut_rate": cohort.mut_rate,
            "population_size": cohort.population_size,
        },
        "n_ground_truth": n_ground_truth,
        "sweep_n_values": list(sweep_n_values),
        "n_bootstrap": n_bootstrap,
        "seed": seed,
        "map_mse": map_mse,
        "bootstrap_coverage": bs_cov_global,
        "fig2a_pair_indices": [list(p) for p in pair_indices],
        "posterior_source": "msprime_independent_proxy",
        "pool_layout": {
            "ground_truth": [0, n_ground_truth],
            "posterior": [n_ground_truth, n_ground_truth + max_n],
            "map_single_draw": total_draws - 1,
        },
        "notes": (
            "v1 uses msprime independent draws as proxy SINGER "
            "posterior. Ground-truth pool, posterior pool, and "
            "MAP single draw are DISJOINT (different RNG seeds) so "
            "MSE/coverage panels don't collapse at large N. Real "
            "SINGER integration is a Phase 2 follow-up."),
    }
    metadata_path = output_dir / "fig2_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2))

    return Fig2RunResult(
        cohort=cohort,
        sweep_n_values=list(sweep_n_values),
        output_dir=output_dir,
        n_ground_truth=n_ground_truth,
        truth=truth,
        fig2a_csv=fig2a_csv,
        fig2b_csv=fig2b_csv,
        fig2c_csv=fig2c_csv,
        metadata_path=metadata_path,
    )


def _select_fig2a_pairs(
    truth: np.ndarray, k: int = 12,
) -> List[PairKey]:
    """Select k pairs spanning the kinship-value range for violins.

    Uses upper triangle off-diagonal entries. Picks pairs at evenly-
    spaced quantiles of the truth-value distribution.
    """
    n = truth.shape[0]
    iu = np.triu_indices(n, k=1)  # strictly upper
    if len(iu[0]) == 0:
        return []
    values = truth[iu]
    order = np.argsort(values)
    quantile_idx = np.linspace(0, len(order) - 1, k).round().astype(int)
    pair_list: List[PairKey] = []
    for qi in quantile_idx:
        pos = order[qi]
        pair_list.append((int(iu[0][pos]), int(iu[1][pos])))
    return pair_list


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig2_out")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--n-diploid", type=int, default=50)
    p.add_argument("--sequence-length", type=int, default=50_000)
    p.add_argument("--recomb-rate", type=float, default=1e-5)
    p.add_argument("--mut-rate", type=float, default=1e-4)
    p.add_argument("--population-size", type=int, default=10_000)
    p.add_argument("--n-ground-truth", type=int, default=100,
                   help="Posterior draws averaged to form the truth")
    p.add_argument("--sweep-n", type=int, nargs="+",
                   default=[5, 10, 20, 50, 100],
                   help="N values to sweep for MSE + coverage panels")
    p.add_argument("--n-bootstrap", type=int, default=100)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    cohort = CohortParams(
        n_diploid=args.n_diploid,
        sequence_length=args.sequence_length,
        recomb_rate=args.recomb_rate,
        mut_rate=args.mut_rate,
        population_size=args.population_size,
    )

    def _progress(i, total):
        if not args.quiet:
            print(f"[fig2] draw {i}/{total}", file=sys.stderr)

    result = run_fig2(
        cohort=cohort,
        output_dir=args.output_dir,
        n_ground_truth=args.n_ground_truth,
        sweep_n_values=args.sweep_n,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        progress_callback=None if args.quiet else _progress,
    )

    print(
        f"[fig2] DONE; output → {result.output_dir}",
        file=sys.stderr,
    )
    print(f"  fig2a = {result.fig2a_csv}", file=sys.stderr)
    print(f"  fig2b = {result.fig2b_csv}", file=sys.stderr)
    print(f"  fig2c = {result.fig2c_csv}", file=sys.stderr)
    print(f"  meta  = {result.metadata_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
