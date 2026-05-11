"""Paper 2 Fig 3f real-data driver — R4 Phase 2 real-data.

Cross-population pathway-h² stability on real 1000G EUR + YRI
chr22 region using ONE posterior sample from R3's existing
SINGER MCMC output (no new MCMC compute).

Reuses Phase-1 v1 helpers verbatim:
- fig3f_panels.simplify_per_pop / remap_child_set
- fig3ab_panels conditional_egrm / lit_child_weight /
  simulate_phenotype / haseman_elston
- fig3e_real_panels position helpers
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

from .fig3ab_panels import (
    _conditional_egrm,
    _lit_child_weight,
    haseman_elston,
    simulate_phenotype,
)
from .fig3f_panels import (
    cross_pop_within_tolerance,
    remap_child_set,
    simplify_per_pop,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_pathway_mutations_in_ts(
    ts, pathway_positions_in_region: set[int],
) -> tuple[list[int], set[int]]:
    """Return (mutation_indices, child_nodes) for sites at
    pathway_positions. Both indexed in the input `ts`.

    pathway_positions_in_region are region-relative (i.e. already
    minus REGION_START so they match the SINGER .trees positions).
    """
    mut_indices: list[int] = []
    children: set[int] = set()
    for site in ts.sites():
        if int(site.position) in pathway_positions_in_region:
            for mutation in site.mutations:
                mut_indices.append(int(mutation.id))
                children.add(int(mutation.node))
    return mut_indices, children


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig3FRealRunResult:
    output_dir: Path
    panel_csv: Path
    summary_path: Path
    rows: list[dict]
    cross_pop_within_tolerance: bool


def run_fig3f_real(
    *,
    ts_path: Path,
    pop_labels_per_haploid: Sequence[str],
    pathway_positions_in_region: set[int],
    n_pheno_replicates: int,
    true_h2: float,
    tolerance: float,
    seed: int,
    output_dir: Path,
    progress_callback=None,
) -> Fig3FRealRunResult:
    """End-to-end R4 v1 pipeline."""
    import tskit

    output_dir.mkdir(parents=True, exist_ok=True)
    if progress_callback is not None:
        progress_callback("load-ts")
    ts = tskit.load(str(ts_path))
    if ts.num_samples != len(pop_labels_per_haploid):
        raise ValueError(
            f"ts.num_samples={ts.num_samples} != "
            f"len(pop_labels)={len(pop_labels_per_haploid)}")

    # Find pathway mutation indices + child nodes in the FULL ts.
    if progress_callback is not None:
        progress_callback("pathway-mut-index")
    _, full_children = find_pathway_mutations_in_ts(
        ts, pathway_positions_in_region)

    # Simplify per pop.
    if progress_callback is not None:
        progress_callback("simplify-per-pop")
    per_pop = simplify_per_pop(ts, pop_labels_per_haploid)

    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    h2_by_pop: Dict[str, List[float]] = {}
    for pop, (ts_pop, _, node_map) in per_pop.items():
        if progress_callback is not None:
            progress_callback(f"compute-grm-{pop}")
        pop_children = remap_child_set(full_children, node_map)
        if not pop_children:
            continue
        weight = _lit_child_weight(pop_children)
        grm_pop, _ = _conditional_egrm(ts_pop, weight)
        grm_pop = np.asarray(grm_pop, dtype=np.float64)

        # Per-pop genotype matrix for phenotype simulation.
        G_pop = ts_pop.genotype_matrix().T.astype(np.float64)
        # Pathway mutations in the simplified ts.
        mut_indices_pop, _ = find_pathway_mutations_in_ts(
            ts_pop, pathway_positions_in_region)
        if not mut_indices_pop:
            continue
        # Map mutation IDs → site/variant column indices.
        site_for_mut: dict[int, int] = {}
        for site in ts_pop.sites():
            for m in site.mutations:
                site_for_mut[int(m.id)] = int(site.id)
        pop_pathway_var_idx = sorted({
            site_for_mut[mi] for mi in mut_indices_pop
            if mi in site_for_mut
        })
        if not pop_pathway_var_idx:
            continue

        h2_by_pop.setdefault(pop, [])
        for rep in range(n_pheno_replicates):
            if progress_callback is not None:
                progress_callback(
                    f"pheno-rep-{pop}-{rep + 1}")
            y_pop, achieved_h2 = simulate_phenotype(
                G_pop, pop_pathway_var_idx, true_h2, rng)
            h2 = haseman_elston(grm_pop, y_pop)
            rows.append({
                "population": pop,
                "pheno_rep": rep,
                "true_h2": achieved_h2,
                "h2_estimate": h2,
                "n_pathway_variants_in_pop": len(
                    pop_pathway_var_idx),
            })
            h2_by_pop[pop].append(h2)

    within, cross_mean, per_pop_mean = (
        cross_pop_within_tolerance(h2_by_pop, tolerance))

    panel_csv = output_dir / "fig3f_real_panel_data.csv"
    summary_path = output_dir / "fig3f_real_summary.json"
    _write_panel(rows, panel_csv)
    summary = {
        "tolerance": tolerance,
        "cross_pop_mean": cross_mean,
        "within_tolerance": within,
        "per_pop_mean": per_pop_mean,
        "ts_path": str(ts_path),
        "notes": (
            "R4 Phase 2 real-data: cross-pop pathway-h² on R3's "
            "SINGER posterior sample 0 (point-estimate-as-MAP). "
            "v1 covers EUR + AFR(via YRI). 5-pop sweep needs a "
            "fresh full-cohort SINGER run (R4-bis)."),
    }
    summary_path.write_text(json.dumps(summary, indent=2))

    return Fig3FRealRunResult(
        output_dir=output_dir, panel_csv=panel_csv,
        summary_path=summary_path, rows=rows,
        cross_pop_within_tolerance=within,
    )


def _write_panel(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "population", "pheno_rep",
            "true_h2", "h2_estimate",
            "n_pathway_variants_in_pop",
        ])
        for r in rows:
            w.writerow([
                r["population"], r["pheno_rep"],
                f"{r['true_h2']:.6g}",
                f"{r['h2_estimate']:.6g}",
                r["n_pathway_variants_in_pop"],
            ])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig3_out/real")
DEFAULT_TS = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/"
    "fig2_out/real/singer_run/ts_0.trees")
DEFAULT_PATHWAY_VARIANTS_TXT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/"
    "fig3_out/real/R-HSA-202733_chr22_variants.txt")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--ts-path", type=Path, default=DEFAULT_TS)
    p.add_argument("--pathway-variants-txt", type=Path,
                   default=DEFAULT_PATHWAY_VARIANTS_TXT)
    p.add_argument("--region-start", type=int, default=22_000_000)
    p.add_argument("--region-end", type=int, default=22_200_000)
    p.add_argument("--n-pheno-replicates", type=int, default=5)
    p.add_argument("--true-h2", type=float, default=0.5)
    p.add_argument("--tolerance", type=float, default=0.05)
    p.add_argument("--seed", type=int, default=2026)
    return p


def main(argv: list[str] | None = None) -> int:
    from .subset_1000g import (
        load_panel, samples_for_super_pop, samples_for_sub_pop,
    )
    from .fig2de_real_panels import assign_haploid_pop_labels

    args = _build_arg_parser().parse_args(argv)

    panel = load_panel()
    eur = samples_for_super_pop(panel, "EUR")
    yri = samples_for_sub_pop(panel, "YRI")
    diploids = list(eur) + list(yri)
    haploid_labels = assign_haploid_pop_labels(
        diploids, panel, use_sub_pop=False)
    print(
        f"[fig3f-real] cohort = EUR {len(eur)} + YRI {len(yri)} "
        f"= {len(diploids)} diploid",
        file=sys.stderr)

    # Build pathway positions (region-relative, matching SINGER ts).
    pathway_positions_abs: set[int] = set()
    with open(args.pathway_variants_txt) as fh:
        for line in fh:
            vid = line.strip()
            parts = vid.split(":")
            if len(parts) < 2:
                continue
            try:
                p = int(parts[1])
            except ValueError:
                continue
            if args.region_start <= p <= args.region_end:
                pathway_positions_abs.add(p)
    pathway_positions_rel = {
        p - args.region_start for p in pathway_positions_abs}
    print(
        f"[fig3f-real] pathway positions in region: "
        f"{len(pathway_positions_rel)}",
        file=sys.stderr)

    def _progress(stage):
        print(f"[fig3f-real] {stage}", file=sys.stderr)

    result = run_fig3f_real(
        ts_path=args.ts_path,
        pop_labels_per_haploid=haploid_labels,
        pathway_positions_in_region=pathway_positions_rel,
        n_pheno_replicates=args.n_pheno_replicates,
        true_h2=args.true_h2,
        tolerance=args.tolerance,
        seed=args.seed,
        output_dir=args.output_dir,
        progress_callback=_progress,
    )
    print(
        f"[fig3f-real] DONE → {result.output_dir} "
        f"(within_tolerance={result.cross_pop_within_tolerance})",
        file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
