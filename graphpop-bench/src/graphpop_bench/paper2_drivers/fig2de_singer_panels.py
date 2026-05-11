"""Paper 2 Fig 2d/2e — REAL SINGER posterior (P2.5-bis).

Same panel structure as `fig2de_panels.py` but the posterior pool
comes from **SINGER MCMC posterior samples on a single ARG**
instead of independent msprime draws. The driver:

1. Simulates a 2-pop AFR+EUR msprime cohort (same demography as
   P2.5 v1) — the cohort is the "truth" ARG.
2. Writes the cohort as a VCF.
3. Runs SINGER with N posterior samples + thin.
4. Converts each SINGER output to a `.trees` file via
   `convert_to_tskit`.
5. Computes `egrm.varGRM_C(ts)` on each posterior `.trees`.
6. Aggregates Welford → posterior mean + variance + SE.
7. Stratifies per-pair SE by population pair.

The hypothesis under test: SINGER posterior SE captures
ARG-inference uncertainty, which should be HIGHER for
between-population pairs (deep coalescence + harder to
resolve). v1 found the OPPOSITE polarity under
independent-draw proxy — this driver tests whether real
SINGER flips the polarity.

Real ground-truth: the inferred GRM from the original msprime
ARG itself (computed via `egrm.varGRM_C` on the original
TreeSequence — bypasses SINGER). This serves as the "true"
GRM under the simulated demography for MSE / abs_error.

Compute budget:
- 5+5 diploid = 20 haploid cohort × 30 kb (vary).
- SINGER MCMC: ~5-15 min for N ≈ 20.
- egrm × N + Welford: ~5 min for N = 20.
- Total: ~15-30 min.
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

from ..competitors import SingerRunner
from .fig2_panels import welford_aggregate
from .fig2de_panels import (
    TwoPopParams,
    build_two_pop_demography,
    stratify_pairs_by_pop,
)


# ---------------------------------------------------------------------------
# Simulation + VCF
# ---------------------------------------------------------------------------

def simulate_cohort_for_singer(
    params: TwoPopParams, seed: int,
):
    """Single msprime sim — the "truth" ARG that SINGER will try
    to recover from its mutations. Returns ts + per-haploid pop
    label."""
    import msprime
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
    # tskit's Population row exposes name via metadata in 1.x; the
    # fig2de v1 driver uses the same code path. Fall back to id-based
    # name if both keys are absent.
    pop_id_to_name: dict = {}
    for p in ts.populations():
        meta = getattr(p, "metadata", None) or {}
        name = meta.get("name") if isinstance(meta, dict) else None
        pop_id_to_name[p.id] = name or f"pop_{p.id}"
    pop_labels = [
        pop_id_to_name.get(ts.node(s).population, "pop_?")
        for s in ts.samples()
    ]
    return ts, pop_labels


def write_cohort_vcf(ts, out_path: Path) -> Path:
    """Write the cohort as a VCF SINGER can read.

    Forces contig_id="1" + position ≥ 1.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        kwargs = {
            "contig_id": "1",
            "position_transform": lambda x: np.fmax(1, x),
        }
        if ts.num_individuals == 0:
            kwargs["ploidy"] = 2
        ts.write_vcf(fh, **kwargs)
    return out_path


# ---------------------------------------------------------------------------
# Posterior aggregation
# ---------------------------------------------------------------------------

def aggregate_posterior_from_trees(
    trees_paths: Sequence[Path],
) -> Tuple[np.ndarray, np.ndarray]:
    """For each .trees in the posterior, compute egrm.varGRM_C; aggregate."""
    import tskit
    from egrm import varGRM_C
    draws: List[np.ndarray] = []
    n_samples: int | None = None
    for p in trees_paths:
        ts = tskit.load(str(p))
        if n_samples is None:
            n_samples = ts.num_samples
        elif ts.num_samples != n_samples:
            raise ValueError(
                f"posterior sample size changed mid-aggregation: "
                f"{ts.num_samples} != {n_samples}")
        grm, _, _ = varGRM_C(ts, var=False)
        draws.append(np.asarray(grm, dtype=np.float64))
    if not draws:
        raise ValueError("no posterior .trees files to aggregate")
    posterior_mean, posterior_var, _ = welford_aggregate(draws)
    posterior_se = np.sqrt(np.maximum(0.0, posterior_var))
    return posterior_mean, posterior_se


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig2DESingerResult:
    output_dir: Path
    fig2d_csv: Path
    fig2e_csv: Path
    metadata_path: Path
    n_posterior_converted: int


def run_fig2de_singer(
    *,
    params: TwoPopParams,
    n_posterior: int,
    thin: int,
    Ne: int,
    seed: int,
    output_dir: Path,
    progress_callback=None,
) -> Fig2DESingerResult:
    """End-to-end P2.5-bis pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)

    if progress_callback is not None:
        progress_callback("simulate")
    truth_ts, pop_labels = simulate_cohort_for_singer(params, seed)
    n_haploid = truth_ts.num_samples

    if progress_callback is not None:
        progress_callback("egrm-truth")
    from egrm import varGRM_C
    truth_grm, _, _ = varGRM_C(truth_ts, var=False)
    truth = np.asarray(truth_grm, dtype=np.float64)

    if progress_callback is not None:
        progress_callback("write-vcf")
    vcf_path = output_dir / "cohort.vcf"
    write_cohort_vcf(truth_ts, vcf_path)

    if progress_callback is not None:
        progress_callback("singer-mcmc")
    runner = SingerRunner()
    singer_result = runner.run(
        vcf_path, output_dir / "singer_run",
        Ne=Ne,
        mutation_rate=params.mut_rate,
        start=0, end=params.sequence_length,
        n_samples=n_posterior,
        thin=thin,
        seed=seed,
        graphpop_commit="p2.5-bis",
    )
    if singer_result.n_samples_converted < 2:
        raise RuntimeError(
            f"SINGER produced only {singer_result.n_samples_converted} "
            f"converted posterior trees; need ≥ 2 for Welford")

    if progress_callback is not None:
        progress_callback("egrm-posterior")
    posterior_mean, posterior_se = aggregate_posterior_from_trees(
        singer_result.trees_paths)

    if posterior_se.shape != truth.shape:
        raise ValueError(
            f"shape mismatch: posterior {posterior_se.shape} vs "
            f"truth {truth.shape}; SINGER may have lost samples")

    buckets = stratify_pairs_by_pop(pop_labels)

    if progress_callback is not None:
        progress_callback("write-panels")
    fig2d_csv = output_dir / "fig2d_singer_panel_data.csv"
    fig2e_csv = output_dir / "fig2e_singer_panel_data.csv"
    _write_fig2d(pop_labels, posterior_mean, posterior_se, truth,
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
        "n_posterior_requested": n_posterior,
        "n_posterior_converted": singer_result.n_samples_converted,
        "thin": thin,
        "Ne": Ne,
        "seed": seed,
        "singer_wall_clock_s": singer_result.profiling.wall_clock_s,
        "singer_rss_peak_mb": singer_result.profiling.rss_peak_mb,
        "notes": (
            "P2.5-bis: real SINGER MCMC posterior on a single "
            "msprime ARG. Replaces the independent-draw proxy used "
            "in P2.5 v1. The truth-GRM is computed directly from "
            "the original msprime TreeSequence."),
    }
    metadata_path = output_dir / "fig2de_singer_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2))

    return Fig2DESingerResult(
        output_dir=output_dir,
        fig2d_csv=fig2d_csv,
        fig2e_csv=fig2e_csv,
        metadata_path=metadata_path,
        n_posterior_converted=singer_result.n_samples_converted,
    )


def _write_fig2d(
    pop_labels: Sequence[str],
    posterior_mean: np.ndarray,
    posterior_se: np.ndarray,
    truth: np.ndarray,
    path: Path,
) -> None:
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
    all_pairs = [pq for pairs in buckets.values() for pq in pairs]
    all_se = np.array([
        posterior_se[i, j] for (i, j) in all_pairs])
    top_q_threshold = (float(np.quantile(all_se, 0.75))
                       if all_se.size else 0.0)
    path.parent.mkdir(parents=True, exist_ok=True)
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
    p.add_argument("--output-dir", type=Path,
                   default=DEFAULT_OUTPUT / "singer")
    p.add_argument("--afr-diploid", type=int, default=5)
    p.add_argument("--eur-diploid", type=int, default=5)
    p.add_argument("--sequence-length", type=int, default=30_000)
    p.add_argument("--n-posterior", type=int, default=20)
    p.add_argument("--thin", type=int, default=10)
    p.add_argument("--Ne", type=int, default=10_000)
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

    def _progress(stage: str):
        if not args.quiet:
            print(f"[fig2de-singer] {stage}", file=sys.stderr)

    result = run_fig2de_singer(
        params=params,
        n_posterior=args.n_posterior,
        thin=args.thin,
        Ne=args.Ne,
        seed=args.seed,
        output_dir=args.output_dir,
        progress_callback=None if args.quiet else _progress,
    )
    print(
        f"[fig2de-singer] DONE → {result.output_dir}",
        file=sys.stderr)
    print(
        f"  posterior samples converted: {result.n_posterior_converted}",
        file=sys.stderr)
    print(f"  fig2d = {result.fig2d_csv}", file=sys.stderr)
    print(f"  fig2e = {result.fig2e_csv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
