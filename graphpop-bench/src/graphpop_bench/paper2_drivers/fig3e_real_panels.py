"""Paper 2 Fig 3e real-data driver — R2 Phase 2 real-data.

Pathway-conditional GRM on real 1000G chr22 + real Reactome
pathway annotations: tests whether branch-GRM retains kinship
signal where PLINK pathway-extract GRM loses it on rare-variant
pathways (the regime P2.4 v1 sim couldn't generate).

Strategy per PLAN_fig3e_real.md:
- Choose one disease-relevant Reactome pathway (default:
  R-HSA-202733 "Cell surface interactions at the vascular wall")
- Resolve to chr22 variant set via has_consequence_edges
- Pre-subset cohort + pathway-only VCFs (via bcftools)
- PLINK GRM on pathway VCF
- SINGER MAP ARG on the cohort's region → conditional_egrm with
  lit_child_weight on pathway mutations
- Compare density + correlation
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence, Set, Tuple

import numpy as np

from ..competitors import PlinkGrmRunner, parse_grm_bin
from .fig3e_panels import (
    collapse_haploid_to_diploid_grm,
    correlate_off_diag,
    density_nnz,
)
from .fig3ab_panels import (
    _conditional_egrm,
    _lit_child_weight,
)
from .pathway_resolver import (
    chr_pathway_variants_by_id,
    load_pathway_names,
    load_pathway_to_genes,
)


# ---------------------------------------------------------------------------
# Variant-ID position parsing
# ---------------------------------------------------------------------------

def parse_variant_position(variant_id: str) -> int | None:
    """`'chr22:1234567:A:G'` → 1234567; else None."""
    parts = variant_id.split(":")
    if len(parts) < 2:
        return None
    try:
        return int(parts[1])
    except ValueError:
        return None


def variant_set_to_positions(
    variant_ids: Set[str],
    in_region: Tuple[int, int] | None = None,
) -> Set[int]:
    """Map a variant-id set → positions. If in_region (start, end)
    is provided, also filter by position range."""
    out: Set[int] = set()
    for vid in variant_ids:
        p = parse_variant_position(vid)
        if p is None:
            continue
        if in_region is not None and not (
            in_region[0] <= p <= in_region[1]
        ):
            continue
        out.add(p)
    return out


# ---------------------------------------------------------------------------
# Branch-GRM with pathway predicate on a real-data SINGER .trees
# ---------------------------------------------------------------------------

def pathway_child_nodes_from_ts(
    ts, pathway_positions: Set[int],
):
    """Find SINGER-converted ts.mutations() at pathway positions →
    their child-node IDs. Used as the lit_child_weight set for
    conditional_egrm.
    """
    children = set()
    for site in ts.sites():
        if int(site.position) in pathway_positions:
            for mutation in site.mutations:
                children.add(int(mutation.node))
    return children


def branch_grm_pathway_from_ts(
    ts, pathway_positions: Set[int],
) -> np.ndarray:
    """conditional_egrm with lit_child_weight on pathway-bearing branches.

    Returns the haploid (n_haploid, n_haploid) matrix.
    """
    children = pathway_child_nodes_from_ts(ts, pathway_positions)
    if not children:
        return np.zeros(
            (ts.num_samples, ts.num_samples), dtype=np.float64)
    weight = _lit_child_weight(children)
    mat, _ = _conditional_egrm(ts, weight)
    return np.asarray(mat, dtype=np.float64)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig3ERealRunResult:
    output_dir: Path
    panel_csv: Path
    metadata_path: Path
    pathway_id: str
    pathway_name: str
    n_pathway_variants_used: int
    branch_nnz_frac: float
    plink_nnz_frac: float
    corr_pearson: float
    corr_spearman: float
    n_diploid: int


def run_fig3e_real(
    *,
    pathway_id: str,
    pathway_vcf: Path,
    region_trees: Path,
    pathway_positions: Set[int],
    cohort_label: str,
    output_dir: Path,
    threshold: float = 0.01,
    plink_extra_args: Sequence[str] = ("--bad-freqs",),
) -> Fig3ERealRunResult:
    """Run PLINK pathway-GRM + branch-pathway-GRM + comparison.

    pathway_vcf : already-extracted pathway-restricted VCF
                  (bcftools + BED).
    region_trees : SINGER .trees on the cohort's region.
    pathway_positions : set of pathway-variant positions to
                        restrict the branch GRM to.
    """
    import tskit
    output_dir.mkdir(parents=True, exist_ok=True)

    # PLINK GRM on pathway VCF.
    plink_dir = output_dir / "plink"
    plink_dir.mkdir(parents=True, exist_ok=True)
    runner = PlinkGrmRunner()
    plink_result = runner.run(
        pathway_vcf, plink_dir,
        input_kind="vcf",
        extra_args=list(plink_extra_args),
        seed=0,
        graphpop_commit=f"fig3e-real-{pathway_id}",
    )
    plink_grm, plink_sample_ids = parse_grm_bin(
        plink_dir / "plink_grm.grm.bin",
        plink_dir / "plink_grm.grm.id",
    )

    # SINGER → tskit .trees → branch-pathway-GRM.
    ts = tskit.load(str(region_trees))
    branch_grm_hap = branch_grm_pathway_from_ts(
        ts, pathway_positions)
    branch_grm = collapse_haploid_to_diploid_grm(branch_grm_hap)

    # If sample counts differ, surface clearly (e.g. SINGER drops
    # samples on multi-allelic / fixed sites).
    if branch_grm.shape != plink_grm.shape:
        raise RuntimeError(
            f"GRM shape mismatch: branch {branch_grm.shape} vs "
            f"PLINK {plink_grm.shape}; sample-set alignment "
            f"between SINGER VCF and PLINK VCF needs reconciliation"
        )

    # Density + correlation metrics.
    branch_nnz = density_nnz(branch_grm, threshold)
    plink_nnz = density_nnz(plink_grm, threshold)
    pearson, spearman = correlate_off_diag(branch_grm, plink_grm)
    n_dipl = plink_grm.shape[0]

    pathway_names = load_pathway_names()
    pathway_name = pathway_names.get(pathway_id, "(unknown)")

    panel_csv = output_dir / "fig3e_real_panel_data.csv"
    metadata_path = output_dir / "fig3e_real_metadata.json"
    _write_panel(
        pathway_id, pathway_name, len(pathway_positions),
        branch_nnz, plink_nnz, pearson, spearman, n_dipl,
        panel_csv,
    )
    metadata = {
        "pathway_id": pathway_id,
        "pathway_name": pathway_name,
        "cohort_label": cohort_label,
        "n_diploid": n_dipl,
        "n_pathway_positions_in_region": len(pathway_positions),
        "threshold": threshold,
        "pathway_vcf": str(pathway_vcf),
        "region_trees": str(region_trees),
        "plink_wall_clock_s": plink_result.profiling.wall_clock_s,
        "plink_rss_peak_mb": plink_result.profiling.rss_peak_mb,
        "notes": (
            "R2 Phase 2 real-data: 1000G chr22 + Reactome "
            "pathway annotations. SINGER MAP ARG on a chr22 "
            "region + pathway-restricted PLINK GRM + branch-GRM "
            "via conditional_egrm. Region-restricted to keep "
            "single-CPU compute tractable."),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2))

    return Fig3ERealRunResult(
        output_dir=output_dir,
        panel_csv=panel_csv,
        metadata_path=metadata_path,
        pathway_id=pathway_id,
        pathway_name=pathway_name,
        n_pathway_variants_used=len(pathway_positions),
        branch_nnz_frac=branch_nnz,
        plink_nnz_frac=plink_nnz,
        corr_pearson=pearson,
        corr_spearman=spearman,
        n_diploid=n_dipl,
    )


def _write_panel(
    pathway_id: str,
    pathway_name: str,
    n_pathway_variants: int,
    branch_nnz: float,
    plink_nnz: float,
    pearson: float,
    spearman: float,
    n_diploid: int,
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pathway_id", "pathway_name", "n_pathway_variants",
            "n_diploid",
            "branch_nnz_frac", "plink_nnz_frac",
            "corr_pearson", "corr_spearman",
        ])
        w.writerow([
            pathway_id, pathway_name, n_pathway_variants, n_diploid,
            f"{branch_nnz:.6g}",
            f"{plink_nnz:.6g}",
            f"{pearson:.6g}",
            f"{spearman:.6g}",
        ])
