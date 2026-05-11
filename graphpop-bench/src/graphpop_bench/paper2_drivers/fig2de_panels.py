"""Paper 2 Fig 2d/2e driver — per-entry posterior SE on a two-
population cohort (P2.5 v1 simulation proxy).

Demonstrates the population-structure signal in Fig 2d/2e
without requiring SINGER or real 1000G data:

- 2-population msprime demographic model (AFR + EUR split at
  ~2000 generations; effective Nes match standard
  Gutenkunst-style estimates).
- N independent msprime draws as proxy SINGER posterior
  (the same Phase-1 Fig 2 strategy, extended to multi-pop).
- Welford aggregation → per-entry mean + variance + SE.
- Pairs stratified by population assignment of the two
  haploids; emits per-pair SE matrix + pop-pair-aggregated
  summary CSVs.

Real-SINGER + 1000G version is P2.5-bis (gated on SINGER
install + data download).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .fig2_panels import welford_aggregate


# ---------------------------------------------------------------------------
# Demographic model + per-draw simulation
# ---------------------------------------------------------------------------

@dataclass
class TwoPopParams:
    """Two-population AFR + EUR demography v1 (illustrative)."""

    afr_size: int = 12_000
    eur_size: int = 4_000
    split_time: int = 2_000              # generations
    n_diploid_per_pop: Dict[str, int] = field(
        default_factory=lambda: {"AFR": 25, "EUR": 25})
    sequence_length: int = 30_000
    recomb_rate: float = 1e-5
    mut_rate: float = 1e-4

    @property
    def total_diploid(self) -> int:
        return sum(self.n_diploid_per_pop.values())

    @property
    def total_haploid(self) -> int:
        return 2 * self.total_diploid


def build_two_pop_demography(params: TwoPopParams):
    """Construct msprime.Demography for a single AFR-EUR split.

    AFR is the ancestral pop; EUR splits off at `split_time`
    generations with a smaller Ne (bottleneck simplification).
    """
    import msprime
    demo = msprime.Demography()
    demo.add_population(name="AFR", initial_size=params.afr_size)
    demo.add_population(name="EUR", initial_size=params.eur_size)
    demo.add_population_split(
        time=params.split_time,
        derived=["EUR"],
        ancestral="AFR",
    )
    return demo


def simulate_draw_two_pop(
    params: TwoPopParams, seed: int,
) -> Tuple[np.ndarray, List[str]]:
    """Run msprime + egrm; return (grm, pop_label_per_haploid).

    `pop_label_per_haploid[k]` ∈ {"AFR", "EUR"} for haploid k.
    The haploid layout follows msprime's sample order (AFR first
    when `samples` is built with AFR before EUR).
    """
    import msprime
    from egrm import varGRM_C

    demo = build_two_pop_demography(params)
    ts = msprime.sim_ancestry(
        samples=params.n_diploid_per_pop,
        demography=demo,
        sequence_length=params.sequence_length,
        recombination_rate=params.recomb_rate,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(
        ts, rate=params.mut_rate,
        model=msprime.BinaryMutationModel(),
        random_seed=seed,
    )
    grm_mat, _, _ = varGRM_C(ts, var=False)
    grm = np.asarray(grm_mat, dtype=np.float64)
    pop_label = _haploid_pop_labels(ts)
    return grm, pop_label


def _haploid_pop_labels(ts) -> List[str]:
    """For each sample node, return the population name."""
    pop_id_to_name = {p.id: p.metadata.get("name", f"pop_{p.id}")
                      for p in ts.populations()}
    # msprime sets the name through Population objects; the
    # metadata key approach is robust to API drift.
    if not pop_id_to_name or all(
        v.startswith("pop_") for v in pop_id_to_name.values()
    ):
        # Fall back: use the table's name attribute directly.
        pop_id_to_name = {p.id: p.name or f"pop_{p.id}"
                          for p in ts.populations()}
    out = []
    for s in ts.samples():
        pop_id = ts.node(s).population
        out.append(pop_id_to_name.get(pop_id, f"pop_{pop_id}"))
    return out


# ---------------------------------------------------------------------------
# Pop-pair stratification
# ---------------------------------------------------------------------------

def stratify_pairs_by_pop(
    pop_labels: Sequence[str], *, include_diagonal: bool = False,
) -> Dict[str, List[Tuple[int, int]]]:
    """Group all (i, j) pairs (i ≤ j) by their pop-pair label.

    Pair label is the alphabetically-sorted join of the two pops
    (e.g. ('EUR', 'AFR') → 'AFR-EUR'). Diagonal pairs are
    optional; by default we emit only off-diagonal pairs (the
    panel statistics use pair distributions).
    """
    n = len(pop_labels)
    buckets: Dict[str, List[Tuple[int, int]]] = {}
    for i in range(n):
        j_start = i if include_diagonal else i + 1
        for j in range(j_start, n):
            a, b = sorted([pop_labels[i], pop_labels[j]])
            key = f"{a}-{b}"
            buckets.setdefault(key, []).append((i, j))
    return buckets


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig2DERunResult:
    output_dir: Path
    fig2d_csv: Path
    fig2e_csv: Path
    metadata_path: Path
    n_posterior: int
    n_ground_truth: int


def run_fig2de(
    *,
    params: TwoPopParams,
    n_ground_truth: int = 100,
    n_posterior: int = 30,
    seed: int = 2026,
    output_dir: Path,
    progress_callback=None,
) -> Fig2DERunResult:
    """End-to-end Fig 2d/2e v1 sweep."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    total_draws = n_ground_truth + n_posterior
    seeds = rng.integers(1, 2**31 - 1, size=total_draws).tolist()

    pop_labels: List[str] | None = None
    draws: List[np.ndarray] = []
    for i, draw_seed in enumerate(seeds):
        grm, pop_label = simulate_draw_two_pop(
            params, int(draw_seed))
        if pop_labels is None:
            pop_labels = pop_label
        else:
            if pop_label != pop_labels:
                raise RuntimeError(
                    "pop_label order changed mid-sweep — "
                    "msprime sample-ordering nondeterminism?")
        draws.append(grm)
        if progress_callback is not None:
            progress_callback(i + 1, total_draws)

    assert pop_labels is not None
    truth = np.mean(draws[:n_ground_truth], axis=0)
    posterior_pool = draws[n_ground_truth:]
    posterior_mean, posterior_var, _ = welford_aggregate(
        posterior_pool)
    posterior_se = np.sqrt(np.maximum(0.0, posterior_var))

    buckets = stratify_pairs_by_pop(pop_labels)

    fig2d_csv = output_dir / "fig2d_panel_data.csv"
    fig2e_csv = output_dir / "fig2e_panel_data.csv"
    _write_fig2d(
        pop_labels, posterior_mean, posterior_se, truth,
        fig2d_csv)
    _write_fig2e(buckets, posterior_se, fig2e_csv)

    metadata = {
        "params": {
            "afr_size": params.afr_size,
            "eur_size": params.eur_size,
            "split_time": params.split_time,
            "n_diploid_per_pop": params.n_diploid_per_pop,
            "sequence_length": params.sequence_length,
            "recomb_rate": params.recomb_rate,
            "mut_rate": params.mut_rate,
        },
        "n_ground_truth": n_ground_truth,
        "n_posterior": n_posterior,
        "seed": seed,
        "pool_layout": {
            "ground_truth": [0, n_ground_truth],
            "posterior": [n_ground_truth, total_draws],
        },
        "notes": (
            "v1: 2-pop AFR+EUR msprime model + independent draws "
            "as proxy posterior. Real-SINGER 1000G version is "
            "P2.5-bis (gated)."),
    }
    metadata_path = output_dir / "fig2de_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2))

    return Fig2DERunResult(
        output_dir=output_dir,
        fig2d_csv=fig2d_csv,
        fig2e_csv=fig2e_csv,
        metadata_path=metadata_path,
        n_posterior=n_posterior,
        n_ground_truth=n_ground_truth,
    )


def _write_fig2d(
    pop_labels: Sequence[str],
    posterior_mean: np.ndarray,
    posterior_se: np.ndarray,
    truth: np.ndarray,
    path: Path,
) -> None:
    """Long-form per-pair CSV — readable by the heatmap renderer."""
    n = len(pop_labels)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "sample_a", "sample_b",
            "pop_a", "pop_b", "pop_pair_label",
            "posterior_mean", "posterior_se",
            "truth", "abs_error",
        ])
        for i in range(n):
            for j in range(i, n):
                a, b = sorted([pop_labels[i], pop_labels[j]])
                label = f"{a}-{b}"
                pm = float(posterior_mean[i, j])
                se = float(posterior_se[i, j])
                tr = float(truth[i, j])
                w.writerow([
                    i, j, pop_labels[i], pop_labels[j], label,
                    f"{pm:.6g}", f"{se:.6g}",
                    f"{tr:.6g}", f"{abs(pm - tr):.6g}",
                ])


def _write_fig2e(
    buckets: Dict[str, List[Tuple[int, int]]],
    posterior_se: np.ndarray,
    path: Path,
) -> None:
    """Per-pop-pair-label summary: count, mean SE, %% above top-quartile."""
    path.parent.mkdir(parents=True, exist_ok=True)
    # Global top-quartile threshold (over all off-diagonal pairs).
    all_pairs = [pq for pairs in buckets.values() for pq in pairs]
    all_se = np.array([
        posterior_se[i, j] for (i, j) in all_pairs])
    top_q_threshold = (float(np.quantile(all_se, 0.75))
                       if all_se.size else 0.0)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pop_pair_label", "n_pairs",
            "mean_se", "median_se", "std_se",
            "frac_above_top_quartile_threshold",
            "top_quartile_threshold",
        ])
        for label in sorted(buckets):
            pairs = buckets[label]
            ses = np.array([
                posterior_se[i, j] for (i, j) in pairs])
            if ses.size == 0:
                continue
            w.writerow([
                label, len(ses),
                f"{ses.mean():.6g}",
                f"{np.median(ses):.6g}",
                f"{ses.std(ddof=1) if ses.size > 1 else 0.0:.6g}",
                f"{(ses > top_q_threshold).mean():.6g}",
                f"{top_q_threshold:.6g}",
            ])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig2_out")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--afr-diploid", type=int, default=25)
    p.add_argument("--eur-diploid", type=int, default=25)
    p.add_argument("--sequence-length", type=int, default=30_000)
    p.add_argument("--n-ground-truth", type=int, default=100)
    p.add_argument("--n-posterior", type=int, default=30)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    params = TwoPopParams(
        n_diploid_per_pop={
            "AFR": args.afr_diploid,
            "EUR": args.eur_diploid,
        },
        sequence_length=args.sequence_length,
    )

    def _progress(i, total):
        if not args.quiet:
            print(f"[fig2de] draw {i}/{total}", file=sys.stderr)

    result = run_fig2de(
        params=params,
        n_ground_truth=args.n_ground_truth,
        n_posterior=args.n_posterior,
        seed=args.seed,
        output_dir=args.output_dir,
        progress_callback=None if args.quiet else _progress,
    )
    print(f"[fig2de] DONE → {result.output_dir}", file=sys.stderr)
    print(f"  fig2d = {result.fig2d_csv}", file=sys.stderr)
    print(f"  fig2e = {result.fig2e_csv}", file=sys.stderr)
    print(f"  meta  = {result.metadata_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
