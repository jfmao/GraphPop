"""Paper 2 Fig 2d/2e real-data driver — R3 Phase 2 real-data.

Per-entry posterior SE on real 1000G EUR + YRI chr22 region via
SINGER posterior N=20. Tests the canonical Deng/Nielsen/Song
claim — between-pop SE > within-pop SE — at full per-population
cohort scale (replacing the P2.5 v1 sim proxy + the P2.5-bis
small-sim-cohort SINGER run).

Reuses:
- Phase-1 / P2.5-bis Welford aggregator + pop stratifier
- P2.5-bis SingerRunner + aggregate_posterior_from_trees
- R0.2 subset_1000g extract_subset_vcf
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
from .fig2de_panels import stratify_pairs_by_pop
from .fig2de_singer_panels import (
    aggregate_posterior_from_trees,
)
from .subset_1000g import (
    PanelRow,
    extract_subset_vcf,
    load_panel,
    samples_for_super_pop,
    samples_for_sub_pop,
    vcf_path_for_chr,
)


# ---------------------------------------------------------------------------
# Pop-label assignment (diploid → 2 haploids per individual)
# ---------------------------------------------------------------------------

def assign_haploid_pop_labels(
    diploid_sample_ids: Sequence[str],
    panel: Sequence[PanelRow],
    use_sub_pop: bool = False,
) -> List[str]:
    """For each diploid sample, emit two haploid pop labels.

    use_sub_pop=False → super_pop (EUR / AFR / etc.)
    use_sub_pop=True  → sub_pop (YRI / GBR / etc.) — finer-grain.
    """
    panel_by_sample = {r.sample: r for r in panel}
    out: List[str] = []
    for sid in diploid_sample_ids:
        row = panel_by_sample.get(sid)
        if row is None:
            raise KeyError(f"sample {sid!r} not in panel")
        label = row.pop if use_sub_pop else row.super_pop
        out.append(label)
        out.append(label)
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig2DERealRunResult:
    output_dir: Path
    fig2d_csv: Path
    fig2e_csv: Path
    metadata_path: Path
    cohort_label: str
    n_haploid: int
    n_posterior_converted: int


def run_fig2de_real(
    *,
    cohort_label: str,
    diploid_sample_ids: Sequence[str],
    pop_labels_per_haploid: Sequence[str],
    chr_label: str,
    start: int,
    end: int,
    n_posterior: int,
    thin: int,
    Ne: int,
    mutation_rate: float,
    output_dir: Path,
    progress_callback=None,
) -> Fig2DERealRunResult:
    """End-to-end R3 pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1 — extract cohort sub-VCF (full chrN; region restriction
    # is passed to SINGER separately).
    if progress_callback is not None:
        progress_callback("extract-cohort")
    cohort_vcf = output_dir / f"{cohort_label}.vcf.gz"
    if not cohort_vcf.exists():
        extract_subset_vcf(
            input_vcf=vcf_path_for_chr(chr_label),
            output_vcf=cohort_vcf,
            sample_ids=diploid_sample_ids,
        )
    # SINGER also needs the uncompressed .vcf with matching prefix.
    cohort_vcf_uncompressed = output_dir / f"{cohort_label}.vcf"
    if not cohort_vcf_uncompressed.exists():
        import subprocess
        subprocess.run(
            ["bash", "-c",
             f"zcat {cohort_vcf} > {cohort_vcf_uncompressed}"],
            check=True)

    # Step 2 — SINGER posterior.
    if progress_callback is not None:
        progress_callback("singer-posterior")
    runner = SingerRunner()
    singer_result = runner.run(
        cohort_vcf_uncompressed,
        output_dir / "singer_run",
        Ne=Ne,
        mutation_rate=mutation_rate,
        start=start, end=end,
        n_samples=n_posterior,
        thin=thin,
        seed=2026,
        graphpop_commit="fig2de-real",
    )
    if singer_result.n_samples_converted < 2:
        raise RuntimeError(
            f"SINGER produced only {singer_result.n_samples_converted} "
            f"converted posterior trees; need ≥ 2 for Welford")

    # Step 3 — egrm per posterior sample + Welford.
    if progress_callback is not None:
        progress_callback("egrm-posterior")
    posterior_mean, posterior_se = aggregate_posterior_from_trees(
        singer_result.trees_paths)
    n_haploid = posterior_mean.shape[0]
    if n_haploid != len(pop_labels_per_haploid):
        raise RuntimeError(
            f"SINGER produced {n_haploid} haploids; pop_labels has "
            f"{len(pop_labels_per_haploid)} — alignment mismatch")

    buckets = stratify_pairs_by_pop(pop_labels_per_haploid)

    fig2d_csv = output_dir / "fig2d_real_panel_data.csv"
    fig2e_csv = output_dir / "fig2e_real_panel_data.csv"
    _write_fig2d(pop_labels_per_haploid, posterior_mean,
                 posterior_se, fig2d_csv)
    _write_fig2e(buckets, posterior_se, fig2e_csv)

    metadata = {
        "cohort_label": cohort_label,
        "n_diploid": len(diploid_sample_ids),
        "n_haploid": n_haploid,
        "chr_label": chr_label,
        "region": [start, end],
        "n_posterior_requested": n_posterior,
        "n_posterior_converted": singer_result.n_samples_converted,
        "thin": thin,
        "Ne": Ne,
        "mutation_rate": mutation_rate,
        "singer_wall_clock_s": singer_result.profiling.wall_clock_s,
        "singer_rss_peak_mb": singer_result.profiling.rss_peak_mb,
        "notes": (
            "R3 Phase 2 real-data: 1000G chr22 EUR+YRI subset, "
            "SINGER posterior, per-entry SE matrix with pop-pair "
            "stratification. Tests canonical Deng/Nielsen/Song "
            "between-pop > within-pop SE claim at real-data scale."),
    }
    metadata_path = output_dir / "fig2de_real_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2))

    return Fig2DERealRunResult(
        output_dir=output_dir,
        fig2d_csv=fig2d_csv, fig2e_csv=fig2e_csv,
        metadata_path=metadata_path,
        cohort_label=cohort_label,
        n_haploid=n_haploid,
        n_posterior_converted=singer_result.n_samples_converted,
    )


def _write_fig2d(
    pop_labels: Sequence[str],
    posterior_mean: np.ndarray,
    posterior_se: np.ndarray,
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
        ])
        for i in range(n):
            for j in range(i, n):
                a, b = sorted([pop_labels[i], pop_labels[j]])
                label = f"{a}-{b}"
                pm = float(posterior_mean[i, j])
                se = float(posterior_se[i, j])
                w.writerow([
                    i, j, pop_labels[i], pop_labels[j], label,
                    f"{pm:.6g}", f"{se:.6g}",
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
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig2_out/real")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--chr", default="22")
    p.add_argument("--start", type=int, default=22_000_000)
    p.add_argument("--end", type=int, default=22_200_000)
    p.add_argument("--n-posterior", type=int, default=20)
    p.add_argument("--thin", type=int, default=5)
    p.add_argument("--Ne", type=int, default=10_000)
    p.add_argument("--mutation-rate", type=float, default=1.25e-8)
    p.add_argument("--cohort-label", default="EUR_YRI")
    p.add_argument("--use-sub-pop-labels", action="store_true",
                   help="Use sub-pop (YRI etc.) instead of super-pop "
                        "labels for stratification")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    panel = load_panel()
    eur = samples_for_super_pop(panel, "EUR")
    yri = samples_for_sub_pop(panel, "YRI")
    diploids = list(eur) + list(yri)
    haploid_labels = assign_haploid_pop_labels(
        diploids, panel, use_sub_pop=args.use_sub_pop_labels)
    print(
        f"[fig2de-real] cohort = {args.cohort_label}: "
        f"EUR {len(eur)} + YRI {len(yri)} = {len(diploids)} diploid",
        file=sys.stderr)

    def _progress(stage):
        print(f"[fig2de-real:{args.cohort_label}] {stage}",
              file=sys.stderr)

    result = run_fig2de_real(
        cohort_label=args.cohort_label,
        diploid_sample_ids=diploids,
        pop_labels_per_haploid=haploid_labels,
        chr_label=args.chr,
        start=args.start, end=args.end,
        n_posterior=args.n_posterior, thin=args.thin,
        Ne=args.Ne, mutation_rate=args.mutation_rate,
        output_dir=args.output_dir,
        progress_callback=_progress,
    )
    print(
        f"[fig2de-real] DONE → {result.output_dir} "
        f"(n_haploid={result.n_haploid}, "
        f"n_posterior_converted={result.n_posterior_converted})",
        file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
