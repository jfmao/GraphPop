"""Paper 2 Fig 5b real-data driver — R1 Phase 2 real-data.

Cryptic-pair recall on real 1000G chr22 cohorts (replaces the
v1 msprime sim base). Reuses the v1 pure-logic helpers
(`inject_mosaic`, `select_disjoint_pairs`, `measure_recall`,
`write_dosage_vcf`) — only the cohort source changes from
msprime to 1000G VCF.

Strategy per PLAN_fig5b_real.md:
- Subset 1000G chr22 to a cohort (EUR / full 2,504) via
  R0.2's `extract_subset_vcf`.
- Load haploid genotype matrix via cyvcf2 + bi-allelic
  polymorphic filter.
- Inject K disjoint (anchor, target) pairs at IBD fraction f.
- Write injected VCF; run PLINK GRM (H1 wrapper); measure
  recall at threshold = 0.0625.
- Sweep f ∈ {0.0625, 0.125, 0.25, 0.5} × R replicates.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence

import numpy as np

from ..competitors import PlinkGrmRunner, parse_grm_bin
from .fig5b_panels import (
    inject_mosaic,
    measure_recall,
    select_disjoint_pairs,
    write_dosage_vcf,
)
from .subset_1000g import (
    extract_subset_vcf,
    load_panel,
    samples_for_super_pop,
    vcf_path_for_chr,
)


# ---------------------------------------------------------------------------
# Real-cohort genotype loader
# ---------------------------------------------------------------------------

def load_real_cohort_haploids(
    vcf_path: Path,
    sample_ids: Sequence[str] | None = None,
) -> tuple[np.ndarray, np.ndarray, List[str]]:
    """Stream `vcf_path` via cyvcf2; return (haploid_matrix, positions, sids).

    Only bi-allelic polymorphic SNVs are retained (consistent with
    the v1 driver's filter). Haploid layout: rows 2i, 2i+1 belong
    to diploid i. Sample IDs follow VCF header order.

    sample_ids : if provided, RESTRICTS to these samples (caller
                 must have already pre-subsetted; cyvcf2 doesn't
                 filter samples on the fly here — we just verify
                 the VCF header matches).
    """
    import cyvcf2

    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")
    vcf = cyvcf2.VCF(str(vcf_path))
    vcf_samples = list(vcf.samples)
    if sample_ids is not None and list(sample_ids) != vcf_samples:
        raise ValueError(
            f"VCF sample order mismatch: VCF has {len(vcf_samples)} "
            f"samples, requested {len(sample_ids)} — please pre-subset "
            f"the VCF via extract_subset_vcf before loading")

    n_diploid = len(vcf_samples)
    haploids: List[np.ndarray] = []
    positions: List[int] = []
    for variant in vcf:
        # Bi-allelic: REF + 1 ALT.
        if not variant.is_snp:
            continue
        if len(variant.ALT) != 1:
            continue
        # gt_types: 0=hom-ref, 1=het, 2=unknown/missing, 3=hom-alt.
        # We need per-haplotype 0/1 calls; use genotypes (list of
        # [a1, a2, phased]) directly. Skip variants with any missing.
        gts = variant.genotypes
        if any((g[0] == -1 or g[1] == -1) for g in gts):
            continue
        # Skip non-bi-allelic alleles (msprime-style filter).
        if any((g[0] not in (0, 1) or g[1] not in (0, 1))
               for g in gts):
            continue
        row = np.empty(2 * n_diploid, dtype=np.int8)
        for i, g in enumerate(gts):
            row[2 * i] = g[0]
            row[2 * i + 1] = g[1]
        # Polymorphic + non-fixed filter.
        s = int(row.sum())
        if s == 0 or s == row.size:
            continue
        haploids.append(row)
        positions.append(max(1, int(variant.POS)))
    if not haploids:
        raise RuntimeError(
            f"no bi-allelic polymorphic SNPs in {vcf_path}")
    haploid_matrix = np.stack(haploids, axis=1).astype(np.int8)
    # Stack-along-axis=1 gave (n_haploid, n_var); but rows above
    # have n_haploid columns each? Let me redo: each `row` has
    # length 2*n_diploid = n_haploid. We stacked them as columns
    # (axis=1) → (n_haploid, n_var). Good.
    pos_arr = np.array(positions, dtype=np.int64)
    return haploid_matrix, pos_arr, vcf_samples


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig5BRealRunResult:
    output_dir: Path
    panel_csv: Path
    summary_path: Path
    metadata_path: Path
    rows: List[dict]
    cohort_label: str
    n_haploid: int


def run_fig5b_real(
    *,
    cohort_label: str,
    sample_ids: Sequence[str],
    vcf_path: Path,
    output_dir: Path,
    fractions: Sequence[float] = (0.0625, 0.125, 0.25, 0.5),
    n_pairs: int = 50,
    n_replicates: int = 5,
    threshold: float = 0.0625,
    seed: int = 2026,
    bcftools_binary: str = "bcftools",
    progress_callback=None,
) -> Fig5BRealRunResult:
    """End-to-end Fig 5b R1 sweep on a real cohort."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rng_master = np.random.default_rng(seed)

    # Step 1 — extract cohort sub-VCF.
    if progress_callback is not None:
        progress_callback("extract-cohort")
    cohort_vcf = output_dir / f"{cohort_label}.vcf.gz"
    if not cohort_vcf.exists():
        extract_subset_vcf(
            input_vcf=vcf_path,
            output_vcf=cohort_vcf,
            sample_ids=sample_ids,
            bcftools_binary=bcftools_binary,
        )

    # Step 2 — load genotypes via cyvcf2.
    if progress_callback is not None:
        progress_callback("load-genotypes")
    haploids, positions, vcf_samples = load_real_cohort_haploids(
        cohort_vcf, sample_ids=list(sample_ids))
    n_haploid = haploids.shape[0]
    n_dipl = n_haploid // 2

    # Step 3 — sweep.
    rows: List[dict] = []
    cell_idx = 0
    cell_count = n_replicates * len(fractions)
    for rep in range(n_replicates):
        rep_seed = int(rng_master.integers(1, 2**31 - 1))
        for f in fractions:
            cell_idx += 1
            if progress_callback is not None:
                progress_callback(
                    f"cell {cell_idx}/{cell_count} "
                    f"rep={rep} f={f}")
            rng_cell = np.random.default_rng(
                rep_seed + int(f * 1000))
            pairs = select_disjoint_pairs(n_dipl, n_pairs, rng_cell)
            injected = inject_mosaic(
                haploids, pairs, f, rng_cell)

            with tempfile.TemporaryDirectory() as td:
                td_path = Path(td)
                inj_vcf = td_path / f"inj_rep{rep}_f{f}.vcf"
                write_dosage_vcf(
                    injected, positions, vcf_samples, inj_vcf)
                plink_out = (
                    output_dir / "plink"
                    / f"{cohort_label}_rep{rep}_f{f}")
                runner = PlinkGrmRunner()
                runner.run(
                    inj_vcf, plink_out,
                    input_kind="vcf",
                    extra_args=["--bad-freqs"],
                    seed=rep,
                    graphpop_commit=f"fig5b-real-{cohort_label}",
                )

            grm, plink_sample_ids = parse_grm_bin(
                plink_out / "plink_grm.grm.bin",
                plink_out / "plink_grm.grm.id",
            )
            n_recovered, k_inj, recall, phis = measure_recall(
                grm, plink_sample_ids, pairs, threshold)
            rows.append({
                "cohort_label": cohort_label,
                "n_haploid": n_haploid,
                "fraction": f,
                "replicate": rep,
                "k_injected": k_inj,
                "k_recovered": n_recovered,
                "recall": recall,
                "median_phi": (
                    float(np.median(phis)) if phis else 0.0),
            })

    panel_csv = output_dir / "fig5b_real_panel_data.csv"
    summary_path = output_dir / "fig5b_real_summary.json"
    metadata_path = output_dir / "fig5b_real_metadata.json"
    _write_panel(rows, panel_csv)
    _write_summary(rows, threshold, summary_path)
    _write_metadata(
        cohort_label, sample_ids, vcf_path, n_haploid,
        len(positions), fractions, n_pairs, n_replicates,
        threshold, seed, metadata_path)

    return Fig5BRealRunResult(
        output_dir=output_dir, panel_csv=panel_csv,
        summary_path=summary_path,
        metadata_path=metadata_path,
        rows=rows,
        cohort_label=cohort_label,
        n_haploid=n_haploid,
    )


def _write_panel(rows: List[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "cohort_label", "n_haploid",
            "fraction", "replicate", "k_injected",
            "k_recovered", "recall", "median_phi",
        ])
        for r in rows:
            w.writerow([
                r["cohort_label"], r["n_haploid"],
                f"{r['fraction']:.6g}", r["replicate"],
                r["k_injected"], r["k_recovered"],
                f"{r['recall']:.6g}", f"{r['median_phi']:.6g}",
            ])


def _write_summary(
    rows: List[dict], threshold: float, path: Path,
) -> None:
    by_fraction: dict = {}
    for r in rows:
        by_fraction.setdefault(r["fraction"], []).append(r["recall"])
    summary = {
        "threshold": threshold,
        "cohort_label": rows[0]["cohort_label"] if rows else None,
        "n_haploid": rows[0]["n_haploid"] if rows else 0,
        "by_fraction": {
            f"{f:.6g}": {
                "mean_recall": float(np.mean(rs)),
                "std_recall": float(
                    np.std(rs, ddof=1) if len(rs) > 1 else 0.0),
                "n_replicates": len(rs),
            }
            for f, rs in by_fraction.items()
        },
        "notes": (
            "R1 Phase 2 real-data: 1000G chr22 cohort + Bernoulli(f) "
            "mosaic injection + PLINK GRM. Replaces the v1 msprime "
            "sim base; pure-logic helpers reused verbatim."),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2))


def _write_metadata(
    cohort_label: str,
    sample_ids: Sequence[str],
    vcf_path: Path,
    n_haploid: int,
    n_variants: int,
    fractions: Sequence[float],
    n_pairs: int,
    n_replicates: int,
    threshold: float,
    seed: int,
    path: Path,
) -> None:
    metadata = {
        "cohort_label": cohort_label,
        "n_samples_requested": len(sample_ids),
        "n_haploid_after_load": n_haploid,
        "n_variants_after_filter": n_variants,
        "source_vcf": str(vcf_path),
        "fractions": list(fractions),
        "n_pairs_per_replicate": n_pairs,
        "n_replicates": n_replicates,
        "threshold": threshold,
        "seed": seed,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, indent=2))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/"
    "benchmarks/fig5_out/real")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--cohort", choices=["EUR", "AFR", "EAS",
                                        "SAS", "AMR", "all"],
                   default="EUR")
    p.add_argument("--chr", default="22")
    p.add_argument("--fractions", type=float, nargs="+",
                   default=[0.0625, 0.125, 0.25, 0.5])
    p.add_argument("--n-pairs", type=int, default=50)
    p.add_argument("--n-replicates", type=int, default=5)
    p.add_argument("--threshold", type=float, default=0.0625)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    panel = load_panel()
    if args.cohort == "all":
        sample_ids = [r.sample for r in panel]
        label = "all_2504"
    else:
        sample_ids = samples_for_super_pop(panel, args.cohort)
        label = args.cohort
    vcf_path = vcf_path_for_chr(args.chr)

    def _progress(stage):
        if not args.quiet:
            print(f"[fig5b-real:{label}] {stage}", file=sys.stderr)

    result = run_fig5b_real(
        cohort_label=label, sample_ids=sample_ids,
        vcf_path=vcf_path, output_dir=args.output_dir,
        fractions=args.fractions, n_pairs=args.n_pairs,
        n_replicates=args.n_replicates,
        threshold=args.threshold, seed=args.seed,
        progress_callback=None if args.quiet else _progress,
    )
    print(f"[fig5b-real] DONE → {result.output_dir}", file=sys.stderr)
    print(f"  cohort = {result.cohort_label} "
          f"(n_haploid={result.n_haploid})", file=sys.stderr)
    print(f"  panel  = {result.panel_csv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
