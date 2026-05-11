"""Paper 2 Fig 5b driver — cryptic-pair recall (simulation proxy).

v1 strategy (per PLAN_fig5b.md):

- msprime base cohort,
- inject K (anchor, target) pairs by per-variant Bernoulli(f)
  mosaicking of the target's haploids onto the anchor's,
- run PLINK GRM (H1 wrapper) on the mutated cohort,
- recall = fraction of injected pairs whose GRM[anchor, target]
  exceeds the 2nd-degree threshold (0.0625).

Sweep f ∈ {0.0625, 0.125, 0.25, 0.5} × replicates → recall curve.

Real-data version (1000G chr22 + Neo4j ingest) is P2.3-bis,
gated on user OK for Neo4j + data download.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence, Tuple

import numpy as np

from ..competitors import PlinkGrmRunner, parse_grm_bin


# ---------------------------------------------------------------------------
# Pure-logic helpers — mosaic injection + recall
# ---------------------------------------------------------------------------

def inject_mosaic(
    haploid_matrix: np.ndarray,           # (n_haploid, n_variants)
    pairs: Sequence[Tuple[int, int]],     # (anchor_dipl, target_dipl) pairs
    fraction: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Inject K synthetic cryptic-relative pairs at controlled IBD.

    For each (anchor_diploid_idx, target_diploid_idx) pair, both
    haploids of the target are independently overwritten on a
    per-variant Bernoulli(`fraction`) draw from the corresponding
    haploid of the anchor; the remaining (1-fraction) keeps the
    original target haploid. The returned matrix has the same
    shape as the input.

    The diploid layout convention: haploid `2*i` and `2*i + 1`
    belong to diploid sample `i`.
    """
    if not 0.0 <= fraction <= 1.0:
        raise ValueError(
            f"fraction must be in [0, 1], got {fraction!r}")
    if haploid_matrix.ndim != 2:
        raise ValueError(
            f"haploid_matrix must be 2-D, got {haploid_matrix.shape}")
    n_hap, n_var = haploid_matrix.shape
    if n_hap % 2 != 0:
        raise ValueError(
            f"n_haploid must be even (diploid layout); got {n_hap}")

    out = haploid_matrix.copy()
    n_dipl = n_hap // 2
    for anchor, target in pairs:
        if not (0 <= anchor < n_dipl and 0 <= target < n_dipl):
            raise ValueError(
                f"pair ({anchor}, {target}) out of range "
                f"for {n_dipl} diploids")
        if anchor == target:
            raise ValueError(
                f"pair ({anchor}, {target}) cannot self-anchor")
        for hap_offset in (0, 1):
            anchor_hap = haploid_matrix[2 * anchor + hap_offset]
            target_hap = out[2 * target + hap_offset]
            mask = rng.random(n_var) < fraction
            target_hap[mask] = anchor_hap[mask]
    return out


def measure_recall(
    grm: np.ndarray, sample_ids: Sequence[str],
    pairs: Sequence[Tuple[int, int]],
    threshold: float,
) -> Tuple[int, int, float, List[float]]:
    """Per-pair recall against the 0.0625 (or other) threshold.

    Returns ``(n_recovered, k_injected, recall, per_pair_phis)``.

    Sample IDs are matched to the diploid index via
    ``sample_ids[i]`` corresponding to diploid `i` (PLINK writes
    them in FID/IID order matching the source VCF).
    """
    n_recovered = 0
    phis: List[float] = []
    n_total = len(pairs)
    if n_total == 0:
        return 0, 0, 0.0, []
    for anchor, target in pairs:
        if anchor >= len(sample_ids) or target >= len(sample_ids):
            raise IndexError(
                f"pair ({anchor}, {target}) out of range for "
                f"{len(sample_ids)} samples")
        phi = float(grm[anchor, target])
        phis.append(phi)
        if phi > threshold:
            n_recovered += 1
    return n_recovered, n_total, n_recovered / n_total, phis


def select_disjoint_pairs(
    n_diploid: int, k: int, rng: np.random.Generator,
) -> List[Tuple[int, int]]:
    """Return K disjoint (anchor, target) diploid-index pairs.

    Each diploid index appears in at most one pair to ensure the
    injection doesn't double-count signal.
    """
    if 2 * k > n_diploid:
        raise ValueError(
            f"cannot select {k} disjoint pairs from {n_diploid} "
            f"diploids (need ≥ 2k)")
    permutation = rng.permutation(n_diploid)
    selected = permutation[: 2 * k]
    return [
        (int(selected[2 * i]), int(selected[2 * i + 1]))
        for i in range(k)
    ]


# ---------------------------------------------------------------------------
# msprime + VCF emission
# ---------------------------------------------------------------------------

@dataclass
class CohortParams:
    n_diploid: int = 150
    sequence_length: int = 30_000
    recomb_rate: float = 1e-5
    mut_rate: float = 1e-4
    population_size: int = 10_000


def simulate_cohort_haploids(
    cohort: CohortParams, seed: int,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Run msprime; return (haploid_matrix, positions, sample_ids).

    haploid_matrix is (n_haploid, n_variants) int8 with diploid
    layout (haploid 2i, 2i+1 belong to diploid i). The matrix is
    restricted to bi-allelic sites (max allele index ≤ 1) so the
    downstream VCF emitter can use plain `0|1` GT calls.
    """
    import msprime
    ts = msprime.sim_ancestry(
        samples=cohort.n_diploid,
        sequence_length=cohort.sequence_length,
        recombination_rate=cohort.recomb_rate,
        population_size=cohort.population_size,
        random_seed=seed,
    )
    # Use BinaryMutationModel so every site has alleles ∈ {0, 1};
    # the default JC69 model produces 4-allele sites that confuse
    # PLINK's bi-allelic VCF expectation.
    ts = msprime.sim_mutations(
        ts, rate=cohort.mut_rate,
        model=msprime.BinaryMutationModel(),
        random_seed=seed)
    raw = ts.genotype_matrix().T  # (haploid, variant)
    # Keep only bi-allelic sites — under a binary model some sites
    # can still receive a second mutation that maps back to 0 (no
    # variation), so filter those out too.
    biallelic_mask = (raw.max(axis=0) <= 1) & (raw.min(axis=0) >= 0) & (raw.sum(axis=0) > 0) & (raw.sum(axis=0) < raw.shape[0])
    if not biallelic_mask.any():
        raise RuntimeError(
            "no bi-allelic polymorphic sites in this msprime sim; "
            "try a different seed / mut_rate")
    geno = raw[:, biallelic_mask].astype(np.int8)
    all_positions = np.fromiter(
        (max(1, int(s.position)) for s in ts.sites()),
        dtype=np.int64, count=ts.num_sites)
    positions = all_positions[biallelic_mask]
    sample_ids = [f"sim_{i:04d}" for i in range(cohort.n_diploid)]
    return geno, positions, sample_ids


def write_dosage_vcf(
    haploid_matrix: np.ndarray, positions: np.ndarray,
    sample_ids: Sequence[str], out_path: Path,
) -> Path:
    """Emit a minimal-but-valid VCF that PLINK 2.0 can read.

    The genotype block uses phased haplotype calls (`0|1`); positions
    are forced ≥ 1 to satisfy the VCF spec.
    """
    n_hap, n_var = haploid_matrix.shape
    n_dipl = n_hap // 2
    if len(sample_ids) != n_dipl:
        raise ValueError(
            f"sample_ids size {len(sample_ids)} != n_dipl {n_dipl}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID=1,length={int(positions.max()) + 1}>\n")
        fh.write(
            "##FORMAT=<ID=GT,Number=1,Type=String,"
            "Description=\"Genotype\">\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for sid in sample_ids:
            fh.write(f"\t{sid}")
        fh.write("\n")
        for vi in range(n_var):
            pos = int(positions[vi])
            row_parts = [
                "1", str(pos), f"rs{vi}", "A", "G",
                ".", "PASS", ".", "GT",
            ]
            for di in range(n_dipl):
                h0 = int(haploid_matrix[2 * di, vi])
                h1 = int(haploid_matrix[2 * di + 1, vi])
                row_parts.append(f"{h0}|{h1}")
            fh.write("\t".join(row_parts) + "\n")
    return out_path


def grm_diploid_matrix(
    grm_bin: Path, grm_id: Path,
) -> Tuple[np.ndarray, List[str]]:
    """Thin wrapper around `parse_grm_bin` — returns (grm, sample_ids).

    The H1 parser already returns the full symmetric matrix; this
    helper exists to keep the driver readable + give us a single
    place to swap parsers if needed.
    """
    grm, sample_ids = parse_grm_bin(grm_bin, grm_id)
    return grm, sample_ids


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig5BRunResult:
    output_dir: Path
    panel_csv: Path
    summary_path: Path
    rows: List[dict]


def run_fig5b(
    *,
    cohort: CohortParams,
    fractions: Sequence[float],
    n_pairs: int = 20,
    n_replicates: int = 5,
    threshold: float = 0.0625,
    seed: int = 2026,
    output_dir: Path,
    progress_callback=None,
) -> Fig5BRunResult:
    """End-to-end Fig 5b v1 sweep."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rng_master = np.random.default_rng(seed)

    rows: List[dict] = []
    cell_count = n_replicates * len(fractions)
    cell_idx = 0
    for rep in range(n_replicates):
        rep_seed = int(rng_master.integers(1, 2**31 - 1))
        base_geno, positions, sample_ids = simulate_cohort_haploids(
            cohort, rep_seed)
        for f in fractions:
            cell_idx += 1
            if progress_callback is not None:
                progress_callback(cell_idx, cell_count, rep, f)
            rng_cell = np.random.default_rng(
                int(rng_master.integers(1, 2**31 - 1)))
            pairs = select_disjoint_pairs(
                cohort.n_diploid, n_pairs, rng_cell)
            injected = inject_mosaic(
                base_geno, pairs, f, rng_cell)

            with tempfile.TemporaryDirectory() as td:
                td_path = Path(td)
                vcf_path = td_path / f"cohort_rep{rep}_f{f}.vcf"
                write_dosage_vcf(
                    injected, positions, sample_ids, vcf_path)
                plink_out = output_dir / "plink" / f"rep{rep}_f{f}"
                runner = PlinkGrmRunner()
                result = runner.run(
                    vcf_path, plink_out,
                    input_kind="vcf",
                    extra_args=["--bad-freqs"],
                    seed=rep,
                    graphpop_commit="fig5b",
                )

            grm, plink_sample_ids = grm_diploid_matrix(
                plink_out / "plink_grm.grm.bin",
                plink_out / "plink_grm.grm.id",
            )
            # PLINK preserves the VCF's sample order.
            n_recovered, k_injected, recall, phis = measure_recall(
                grm, plink_sample_ids, pairs, threshold)
            rows.append({
                "fraction": f,
                "replicate": rep,
                "k_injected": k_injected,
                "k_recovered": n_recovered,
                "recall": recall,
                "median_phi": float(np.median(phis)) if phis else 0.0,
            })

    panel_csv = output_dir / "fig5b_panel_data.csv"
    summary_path = output_dir / "fig5b_summary.json"
    _write_panel(rows, panel_csv)
    _write_summary(rows, threshold, summary_path)

    return Fig5BRunResult(
        output_dir=output_dir,
        panel_csv=panel_csv,
        summary_path=summary_path,
        rows=rows,
    )


def _write_panel(rows: List[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "fraction", "replicate", "k_injected",
            "k_recovered", "recall", "median_phi",
        ])
        for r in rows:
            w.writerow([
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
            "v1 simulation proxy: msprime cohort + Bernoulli(f) "
            "mosaic injection + PLINK GRM. Real-data 1000G chr22 + "
            "Neo4j-ingest version is P2.3-bis (gated)."),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig5_out")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--n-diploid", type=int, default=150)
    p.add_argument("--sequence-length", type=int, default=30_000)
    p.add_argument("--fractions", type=float, nargs="+",
                   default=[0.0625, 0.125, 0.25, 0.5])
    p.add_argument("--n-pairs", type=int, default=20)
    p.add_argument("--n-replicates", type=int, default=5)
    p.add_argument("--threshold", type=float, default=0.0625)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    cohort = CohortParams(
        n_diploid=args.n_diploid,
        sequence_length=args.sequence_length,
    )

    def _progress(idx, total, rep, f):
        if not args.quiet:
            print(
                f"[fig5b] cell {idx}/{total} "
                f"rep={rep} fraction={f}",
                file=sys.stderr)

    result = run_fig5b(
        cohort=cohort,
        fractions=args.fractions,
        n_pairs=args.n_pairs,
        n_replicates=args.n_replicates,
        threshold=args.threshold,
        seed=args.seed,
        output_dir=args.output_dir,
        progress_callback=None if args.quiet else _progress,
    )
    print(f"[fig5b] DONE → {result.output_dir}", file=sys.stderr)
    print(f"  panel  = {result.panel_csv}", file=sys.stderr)
    print(f"  summary = {result.summary_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
