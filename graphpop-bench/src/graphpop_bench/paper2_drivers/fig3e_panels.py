"""Paper 2 Fig 3e driver — pathway-conditional GRM: branch vs SNP.

P2.4 v1 (per PLAN_fig3e.md): sim cohort + synthetic
pathway-of-K mutations + GraphPop branch-GRM (conditional_egrm)
vs PLINK pathway-extract GRM, swept over pathway size to expose
the rare-variant regime where SNP-based methods lose signal.

Canonical 1000G + Reactome version is P2.4-bis (gated).
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

from ..competitors import PlinkGrmRunner, parse_grm_bin
# Reuse the inlined conditional_egrm / lit_child_weight from
# fig3ab_panels — they're already verified byte-stable.
from .fig3ab_panels import (
    CohortParams as _BaseCohort,
    _conditional_egrm,
    _lit_child_weight,
    child_nodes_of_mutations,
)
# Reuse the VCF emitter from fig5b_panels (handles diploid layout
# + position-clamp + binary GT).
from .fig5b_panels import write_dosage_vcf


# ---------------------------------------------------------------------------
# Pure helpers — pathway-subset VCF + diploid GRM collapse
# ---------------------------------------------------------------------------

def collapse_haploid_to_diploid_grm(
    haploid_grm: np.ndarray,
) -> np.ndarray:
    """Sum the 4 (haploid × haploid) entries per diploid pair.

    haploid_grm: (n_hap, n_hap), even n_hap, diploid layout (2i,
                 2i+1 share individual i).
    returns:     (n_dipl, n_dipl) where
                 G_dipl[i, j] = sum_{k in {2i,2i+1}, l in {2j,2j+1}} G_hap[k,l].
    """
    n_hap = haploid_grm.shape[0]
    if n_hap % 2 != 0:
        raise ValueError(
            f"haploid_grm must have even n_haploid, got {n_hap}")
    n_dipl = n_hap // 2
    # Reshape to (n_dipl, 2, n_dipl, 2) then sum over the inner axes.
    reshaped = haploid_grm.reshape(n_dipl, 2, n_dipl, 2)
    return reshaped.sum(axis=(1, 3))


def compute_branch_pathway_grm(
    ts, pathway_mutation_indices: Sequence[int],
) -> np.ndarray:
    """Pathway-restricted branch GRM via lit_child_weight + conditional_egrm.

    Returns the n_haploid × n_haploid matrix; caller collapses to
    diploid for the PLINK comparison.
    """
    if len(pathway_mutation_indices) == 0:
        return np.zeros((ts.num_samples, ts.num_samples),
                        dtype=np.float64)
    children = child_nodes_of_mutations(
        ts, pathway_mutation_indices)
    weight_fn = _lit_child_weight(children)
    mat, _ = _conditional_egrm(ts, weight_fn)
    return np.asarray(mat, dtype=np.float64)


def write_pathway_vcf(
    ts, pathway_mutation_indices: Sequence[int],
    n_diploid: int, out_path: Path,
) -> tuple[Path, int]:
    """Emit a binary, polymorphic-only VCF restricted to the pathway.

    Returns (path, n_variants_written). Skips sites that are fixed
    in the cohort (no PLINK signal) or multi-allelic.
    """
    geno_full = ts.genotype_matrix().T   # (n_haploid, n_variant)
    n_var = geno_full.shape[1]
    pathway_set = set(int(m) for m in pathway_mutation_indices)
    if not pathway_set:
        raise ValueError("pathway_mutation_indices is empty")
    keep_mask = np.zeros(n_var, dtype=bool)
    for m in pathway_set:
        if 0 <= m < n_var:
            keep_mask[m] = True
    raw = geno_full[:, keep_mask]
    # Bi-allelic + polymorphic filter.
    polymorphic_mask = (
        (raw.max(axis=0) <= 1)
        & (raw.min(axis=0) >= 0)
        & (raw.sum(axis=0) > 0)
        & (raw.sum(axis=0) < raw.shape[0])
    )
    if not polymorphic_mask.any():
        # Write an empty-variant VCF anyway so PLINK can detect
        # the empty case. We use a 1-site stub with all-zero GT.
        out_path.parent.mkdir(parents=True, exist_ok=True)
        sample_ids = [f"sim_{i:04d}" for i in range(n_diploid)]
        out_path.write_text("\n".join([
            "##fileformat=VCFv4.2",
            "##contig=<ID=1,length=2>",
            ('##FORMAT=<ID=GT,Number=1,Type=String,'
             'Description="Genotype">'),
            ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"
             "FORMAT\t" + "\t".join(sample_ids)),
        ]) + "\n")
        return out_path, 0
    kept_geno = raw[:, polymorphic_mask].astype(np.int8)
    pathway_positions = np.array([
        max(1, int(ts.site(m).position))
        for m in sorted(pathway_set) if 0 <= m < n_var
    ])
    pathway_positions = pathway_positions[
        keep_mask[sorted(pathway_set)] if False else
        np.arange(len(pathway_positions))[
            polymorphic_mask[: len(pathway_positions)]
        ]
    ] if False else None
    # The above branch is convoluted; recompute positions cleanly:
    pathway_positions = []
    sorted_idx = sorted(pathway_set)
    for m in sorted_idx:
        if 0 <= m < n_var:
            pathway_positions.append(max(1, int(ts.site(m).position)))
    pathway_positions = np.array(pathway_positions)
    # Trim positions to the polymorphic subset (positions ordering
    # matches the keep_mask ordering, which was sorted).
    pathway_positions = pathway_positions[polymorphic_mask]
    sample_ids = [f"sim_{i:04d}" for i in range(n_diploid)]
    write_dosage_vcf(kept_geno, pathway_positions, sample_ids, out_path)
    return out_path, int(polymorphic_mask.sum())


# ---------------------------------------------------------------------------
# Density / correlation helpers
# ---------------------------------------------------------------------------

def off_diagonal_pairs(n: int) -> Tuple[np.ndarray, np.ndarray]:
    """Strict-upper-triangle (i < j) indices for an (n, n) matrix."""
    return np.triu_indices(n, k=1)


def density_nnz(
    grm: np.ndarray, threshold: float = 1e-6,
) -> float:
    """Fraction of off-diagonal pairs with `|G_ij| > threshold`."""
    if grm.ndim != 2 or grm.shape[0] != grm.shape[1]:
        raise ValueError(f"grm shape {grm.shape} must be square")
    iu = off_diagonal_pairs(grm.shape[0])
    vals = grm[iu]
    if vals.size == 0:
        return 0.0
    return float((np.abs(vals) > threshold).mean())


def correlate_off_diag(
    grm_a: np.ndarray, grm_b: np.ndarray,
) -> Tuple[float, float]:
    """Pearson + Spearman correlation on off-diagonal pairs.

    Skips pairs where BOTH inputs are exactly zero (uninformative).
    Returns (pearson, spearman) as Python floats. Either may be NaN
    when one of the inputs is constant.
    """
    if grm_a.shape != grm_b.shape:
        raise ValueError(
            f"shape mismatch: {grm_a.shape} vs {grm_b.shape}")
    iu = off_diagonal_pairs(grm_a.shape[0])
    a = grm_a[iu].astype(np.float64)
    b = grm_b[iu].astype(np.float64)
    # Drop double-zeros so a long tail of NaN-pair doesn't dominate.
    nonzero_mask = (np.abs(a) > 0) | (np.abs(b) > 0)
    a = a[nonzero_mask]
    b = b[nonzero_mask]
    if a.size < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan"), float("nan")
    pearson = float(np.corrcoef(a, b)[0, 1])
    # Spearman via rank transform (no scipy dep).
    ra = a.argsort().argsort()
    rb = b.argsort().argsort()
    spearman = float(np.corrcoef(ra, rb)[0, 1])
    return pearson, spearman


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig3ECohortParams:
    afr_diploid: int = 25
    eur_diploid: int = 25
    sequence_length: int = 30_000
    recomb_rate: float = 1e-5
    mut_rate: float = 1e-4

    @property
    def n_diploid(self) -> int:
        return self.afr_diploid + self.eur_diploid

    @property
    def n_haploid(self) -> int:
        return 2 * self.n_diploid


@dataclass
class Fig3ERunResult:
    output_dir: Path
    panel_csv: Path
    summary_path: Path
    rows: List[dict]


def simulate_fig3e_cohort(
    params: Fig3ECohortParams, seed: int,
):
    """msprime 2-pop AFR+EUR + binary mutations; return ts."""
    import msprime
    demo = msprime.Demography()
    demo.add_population(name="AFR", initial_size=12_000)
    demo.add_population(name="EUR", initial_size=4_000)
    demo.add_population_split(
        time=2_000, derived=["EUR"], ancestral="AFR")
    ts = msprime.sim_ancestry(
        samples={"AFR": params.afr_diploid,
                 "EUR": params.eur_diploid},
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
    return ts


def run_fig3e(
    *,
    params: Fig3ECohortParams,
    pathway_sizes: Sequence[int] = (5, 10, 25, 50, 100),
    n_replicates: int = 5,
    threshold: float = 1e-6,
    seed: int = 2026,
    output_dir: Path,
    progress_callback=None,
) -> Fig3ERunResult:
    """End-to-end Fig 3e v1."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    sim_seeds = rng.integers(1, 2**31 - 1, size=n_replicates).tolist()

    rows: List[dict] = []
    total_cells = n_replicates * len(pathway_sizes)
    cell_idx = 0
    for rep, sim_seed in enumerate(sim_seeds):
        ts = simulate_fig3e_cohort(params, int(sim_seed))
        for M in sorted(pathway_sizes):
            cell_idx += 1
            if progress_callback is not None:
                progress_callback(cell_idx, total_cells, rep, M)
            n_var = ts.num_mutations
            if M > n_var:
                rows.append({
                    "pathway_size": M,
                    "replicate": rep,
                    "n_pathway_variants_used": 0,
                    "n_pairs_off_diag": 0,
                    "branch_nnz_frac": float("nan"),
                    "plink_nnz_frac": float("nan"),
                    "branch_minus_plink_nnz": float("nan"),
                    "corr_pearson": float("nan"),
                    "corr_spearman": float("nan"),
                })
                continue
            pathway_idx = list(range(M))
            branch_hap = compute_branch_pathway_grm(ts, pathway_idx)
            branch_dipl = collapse_haploid_to_diploid_grm(branch_hap)

            pathway_dir = output_dir / "plink" / f"rep{rep}_M{M}"
            vcf_path = pathway_dir / "pathway.vcf"
            vcf_path.parent.mkdir(parents=True, exist_ok=True)
            _, n_kept = write_pathway_vcf(
                ts, pathway_idx, params.n_diploid, vcf_path)

            if n_kept == 0:
                # No polymorphic pathway variants — PLINK can't
                # estimate anything; record zero-density.
                plink_dipl = np.zeros(
                    (params.n_diploid, params.n_diploid))
            else:
                runner = PlinkGrmRunner()
                result = runner.run(
                    vcf_path, pathway_dir,
                    input_kind="vcf",
                    extra_args=["--bad-freqs"],
                    seed=rep, graphpop_commit="fig3e",
                )
                plink_dipl, _ = parse_grm_bin(
                    pathway_dir / "plink_grm.grm.bin",
                    pathway_dir / "plink_grm.grm.id",
                )

            branch_nnz = density_nnz(branch_dipl, threshold)
            plink_nnz = density_nnz(plink_dipl, threshold)
            pearson, spearman = correlate_off_diag(
                branch_dipl, plink_dipl)
            iu = off_diagonal_pairs(params.n_diploid)
            n_pairs = int(iu[0].size)
            rows.append({
                "pathway_size": M,
                "replicate": rep,
                "n_pathway_variants_used": int(n_kept),
                "n_pairs_off_diag": n_pairs,
                "branch_nnz_frac": branch_nnz,
                "plink_nnz_frac": plink_nnz,
                "branch_minus_plink_nnz": branch_nnz - plink_nnz,
                "corr_pearson": pearson,
                "corr_spearman": spearman,
            })

    panel_csv = output_dir / "fig3e_panel_data.csv"
    _write_panel(rows, panel_csv)
    summary_path = output_dir / "fig3e_summary.json"
    _write_summary(rows, threshold, summary_path)

    return Fig3ERunResult(
        output_dir=output_dir, panel_csv=panel_csv,
        summary_path=summary_path, rows=rows,
    )


def _write_panel(rows: List[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pathway_size", "replicate", "n_pathway_variants_used",
            "n_pairs_off_diag",
            "branch_nnz_frac", "plink_nnz_frac",
            "branch_minus_plink_nnz",
            "corr_pearson", "corr_spearman",
        ])
        for r in rows:
            w.writerow([
                r["pathway_size"], r["replicate"],
                r["n_pathway_variants_used"],
                r["n_pairs_off_diag"],
                f"{r['branch_nnz_frac']:.6g}",
                f"{r['plink_nnz_frac']:.6g}",
                f"{r['branch_minus_plink_nnz']:.6g}",
                f"{r['corr_pearson']:.6g}",
                f"{r['corr_spearman']:.6g}",
            ])


def _write_summary(
    rows: List[dict], threshold: float, path: Path,
) -> None:
    by_size: Dict[int, List[dict]] = {}
    for r in rows:
        by_size.setdefault(r["pathway_size"], []).append(r)
    summary = {
        "threshold": threshold,
        "by_size": {
            str(M): {
                "n_replicates": len(group),
                "mean_branch_nnz_frac": float(np.mean([
                    r["branch_nnz_frac"] for r in group])),
                "mean_plink_nnz_frac": float(np.mean([
                    r["plink_nnz_frac"] for r in group])),
                "mean_corr_pearson": float(np.nanmean([
                    r["corr_pearson"] for r in group])),
            }
            for M, group in sorted(by_size.items())
        },
        "notes": (
            "v1: msprime 2-pop AFR+EUR + synthetic pathways from "
            "first M mutations. Real-data Fig 3e (1000G + Reactome "
            "+ GTEx) is P2.4-bis."),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig3_out")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--afr-diploid", type=int, default=25)
    p.add_argument("--eur-diploid", type=int, default=25)
    p.add_argument("--sequence-length", type=int, default=30_000)
    p.add_argument("--pathway-sizes", type=int, nargs="+",
                   default=[5, 10, 25, 50, 100])
    p.add_argument("--n-replicates", type=int, default=5)
    p.add_argument("--threshold", type=float, default=1e-6)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    params = Fig3ECohortParams(
        afr_diploid=args.afr_diploid,
        eur_diploid=args.eur_diploid,
        sequence_length=args.sequence_length,
    )

    def _progress(idx, total, rep, M):
        if not args.quiet:
            print(
                f"[fig3e] cell {idx}/{total} rep={rep} M={M}",
                file=sys.stderr)

    result = run_fig3e(
        params=params,
        pathway_sizes=args.pathway_sizes,
        n_replicates=args.n_replicates,
        threshold=args.threshold,
        seed=args.seed,
        output_dir=args.output_dir,
        progress_callback=None if args.quiet else _progress,
    )
    print(f"[fig3e] DONE → {result.output_dir}", file=sys.stderr)
    print(f"  panel  = {result.panel_csv}", file=sys.stderr)
    print(f"  summary = {result.summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
