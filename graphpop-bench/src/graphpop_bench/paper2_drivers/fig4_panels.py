"""Paper 2 Fig 4 driver — sim-only biobank-scaling sweep (G4).

Sweeps `n_diploid` × replicates × competitor wrapper; collects
the wall-clock + RSS-peak receipts produced by the H3/H4
wrappers' `profile_command` subprocess isolation; emits a tidy
per-(N, tool, replicate) panel CSV.

v1 scope (Phase 1, sim-only, on a dev box):
- Competitors: **tskit** (H3) + **egrm** (H4). Both implement
  the branch GRM and ship with subprocess-isolated profiling.
- PLINK / KING / GraphPop curves are out-of-scope until those
  binaries (or a parameterised Java fixture) are wired up. The
  driver still detects them via `is_available()` so the figure
  caption can honestly explain absent curves.
- Sweep: n_diploid ∈ {25, 50, 100, 250}; haploid ∈ {50..500};
  3 replicates per cell. ~15 min wall-clock on a dev box.

The driver is pure orchestration: it never measures wall-clock
itself — measurement is delegated to the wrapper layer where
`profile_command` runs the inner Python subprocess under
`/usr/bin/time -v` (or rusage fallback) for sandbox-isolated
RSS. Receipts are the audit trail.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

from ..competitors import (
    EgrmRunner,
    KingRunner,
    PlinkGrmRunner,
    TskitBranchGrmRunner,
)


# ---------------------------------------------------------------------------
# Tool registry
# ---------------------------------------------------------------------------

# Each entry maps a `tool` label → (RunnerClass, extra_run_kwargs).
# Some wrappers accept .trees (tskit, egrm); the SNP-level
# wrappers want VCF/BED instead. Driver routes accordingly.
TOOL_TREESEQ: Dict[str, Tuple[Any, Dict[str, Any]]] = {
    "tskit": (TskitBranchGrmRunner, {}),
    "egrm":  (EgrmRunner, {"compute_var": False}),
}

TOOL_VCF: Dict[str, Tuple[Any, Dict[str, Any]]] = {
    "plink": (PlinkGrmRunner, {}),
    "king":  (KingRunner, {"mode": "kinship"}),
}

ALL_TOOLS = list(TOOL_TREESEQ.keys()) + list(TOOL_VCF.keys())


def is_tool_available(tool: str) -> bool:
    """Whether the chosen tool's wrapper can run end-to-end."""
    if tool in TOOL_TREESEQ:
        Runner, _ = TOOL_TREESEQ[tool]
        return bool(Runner.is_available())
    if tool in TOOL_VCF:
        Runner, _ = TOOL_VCF[tool]
        return bool(Runner.is_available())
    raise ValueError(f"unknown tool: {tool!r}")


# ---------------------------------------------------------------------------
# Pure-logic helpers (no msprime / wrapper imports required)
# ---------------------------------------------------------------------------

@dataclass
class CellResult:
    """One (N, tool, replicate) row of the scaling sweep."""

    n_diploid: int
    n_haploid: int
    sequence_length: int
    tool: str
    replicate: int
    wall_clock_s: float
    rss_peak_mb: float
    user_cpu_s: float
    system_cpu_s: float
    exit_code: int
    backend: str
    receipt_path: str


def aggregate_by_tool(
    rows: Sequence[CellResult],
) -> Dict[Tuple[int, str], Dict[str, float]]:
    """Compute mean + std-of-mean per (n_haploid, tool) cell."""
    buckets: Dict[Tuple[int, str], List[CellResult]] = {}
    for r in rows:
        buckets.setdefault((r.n_haploid, r.tool), []).append(r)

    out: Dict[Tuple[int, str], Dict[str, float]] = {}
    for key, group in buckets.items():
        ws = [g.wall_clock_s for g in group]
        rs = [g.rss_peak_mb for g in group]
        n = len(group)
        out[key] = {
            "n": n,
            "wall_mean": sum(ws) / n,
            "wall_std": _sample_std(ws),
            "rss_mean": sum(rs) / n,
            "rss_std": _sample_std(rs),
        }
    return out


def _sample_std(xs: List[float]) -> float:
    n = len(xs)
    if n < 2:
        return 0.0
    m = sum(xs) / n
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))


def log_log_fit(
    xs: Sequence[float], ys: Sequence[float],
) -> tuple[float, float]:
    """Least-squares power-law fit: `log y = slope · log x + intercept`.

    Returns ``(slope, intercept)``. Requires at least 2 strictly
    positive points; raises otherwise.
    """
    if len(xs) != len(ys):
        raise ValueError(f"xs/ys length mismatch: {len(xs)} != {len(ys)}")
    if len(xs) < 2:
        raise ValueError("need ≥ 2 points for log-log fit")
    xs_log = [math.log(x) for x in xs]
    ys_log = [math.log(y) for y in ys]
    n = len(xs_log)
    mean_x = sum(xs_log) / n
    mean_y = sum(ys_log) / n
    num = sum((xs_log[i] - mean_x) * (ys_log[i] - mean_y) for i in range(n))
    den = sum((xs_log[i] - mean_x) ** 2 for i in range(n))
    if den == 0:
        raise ValueError("xs collinear after log; cannot fit slope")
    slope = num / den
    intercept = mean_y - slope * mean_x
    return slope, intercept


def write_panel_csv(rows: Sequence[CellResult], path: Path) -> None:
    """Emit the tidy panel CSV for fig4_figures consumption."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "n_diploid", "n_haploid", "sequence_length",
            "tool", "replicate",
            "wall_clock_s", "rss_peak_mb",
            "user_cpu_s", "system_cpu_s",
            "exit_code", "backend", "receipt_path",
        ])
        for r in rows:
            w.writerow([
                r.n_diploid, r.n_haploid, r.sequence_length,
                r.tool, r.replicate,
                f"{r.wall_clock_s:.6g}", f"{r.rss_peak_mb:.6g}",
                f"{r.user_cpu_s:.6g}", f"{r.system_cpu_s:.6g}",
                r.exit_code, r.backend, r.receipt_path,
            ])


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

@dataclass
class CohortParams:
    n_diploid: int
    sequence_length: int = 30_000
    recomb_rate: float = 1e-5
    mut_rate: float = 1e-4
    population_size: int = 10_000

    @property
    def n_haploid(self) -> int:
        return self.n_diploid * 2


def simulate_cohort_trees(
    cohort: CohortParams, seed: int, out_path: Path,
) -> Path:
    """Write a .trees file for one cohort + seed. Returns the path."""
    import msprime
    ts = msprime.sim_ancestry(
        samples=cohort.n_diploid,
        sequence_length=cohort.sequence_length,
        recombination_rate=cohort.recomb_rate,
        population_size=cohort.population_size,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=cohort.mut_rate,
                                random_seed=seed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ts.dump(str(out_path))
    return out_path


def write_vcf_for_cohort(trees_path: Path, vcf_path: Path) -> Path:
    """Convert a `.trees` file to a VCF for SNP-level wrappers.

    Uses `tskit.TreeSequence.write_vcf`. msprime sims via
    `sim_ancestry(samples=...)` carry diploid individual rows
    in the TreeSequence; in that case tskit infers ploidy from
    the individuals and forbids passing `ploidy=` explicitly.
    """
    import tskit
    ts = tskit.load(str(trees_path))
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    with open(vcf_path, "w") as fh:
        # contig_id="1", increment any zero positions to 1 to
        # satisfy VCF spec (msprime may place a variant at pos 0).
        kwargs = {
            "contig_id": "1",
            "position_transform": lambda x: __import__(
                "numpy").fmax(1, x),
        }
        if ts.num_individuals == 0:
            kwargs["ploidy"] = 2
        ts.write_vcf(fh, **kwargs)
    return vcf_path


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig4SweepConfig:
    n_diploid_values: List[int] = field(
        default_factory=lambda: [25, 50, 100, 250])
    n_replicates: int = 3
    sequence_length: int = 30_000
    recomb_rate: float = 1e-5
    mut_rate: float = 1e-4
    population_size: int = 10_000
    seed: int = 2026


@dataclass
class Fig4RunResult:
    output_dir: Path
    panel_csv: Path
    skip_path: Path
    metadata_path: Path
    rows: List[CellResult]


def run_fig4(
    *, config: Fig4SweepConfig,
    output_dir: Path,
    tools: Sequence[str] = ("tskit", "egrm"),
    progress_callback=None,
) -> Fig4RunResult:
    """End-to-end Fig 4 scaling sweep."""
    output_dir.mkdir(parents=True, exist_ok=True)
    receipts_dir = output_dir / "receipts"
    receipts_dir.mkdir(parents=True, exist_ok=True)
    trees_dir = output_dir / "trees_cache"
    trees_dir.mkdir(parents=True, exist_ok=True)

    skipped: Dict[str, str] = {}
    available_tools: List[str] = []
    for t in tools:
        if not is_tool_available(t):
            skipped[t] = "wrapper.is_available() returned False"
        else:
            available_tools.append(t)

    rows: List[CellResult] = []
    cell_count = (
        len(config.n_diploid_values)
        * config.n_replicates * len(available_tools)
    )
    cell_idx = 0
    for n_diploid in config.n_diploid_values:
        cohort = CohortParams(
            n_diploid=n_diploid,
            sequence_length=config.sequence_length,
            recomb_rate=config.recomb_rate,
            mut_rate=config.mut_rate,
            population_size=config.population_size,
        )
        for rep in range(config.n_replicates):
            cohort_seed = config.seed + 10_000 * n_diploid + rep
            trees_path = (
                trees_dir
                / f"cohort_n{n_diploid}_r{rep}.trees"
            )
            if not trees_path.exists():
                simulate_cohort_trees(cohort, cohort_seed,
                                       trees_path)

            # Per-replicate VCF cache (produced lazily on first
            # SNP-level wrapper invocation; skipped if all tools
            # are tree-seq-native).
            vcf_path = (
                trees_dir / f"cohort_n{n_diploid}_r{rep}.vcf"
            )
            for tool in available_tools:
                cell_idx += 1
                if progress_callback is not None:
                    progress_callback(
                        cell_idx, cell_count,
                        n_diploid, tool, rep)
                cell_outdir = (
                    receipts_dir
                    / f"{tool}_n{n_diploid}_r{rep}"
                )
                if tool in TOOL_VCF and not vcf_path.exists():
                    write_vcf_for_cohort(trees_path, vcf_path)
                row = _run_one_cell(
                    tool=tool,
                    cohort=cohort,
                    trees_path=trees_path,
                    vcf_path=vcf_path,
                    output_dir=cell_outdir,
                    replicate=rep,
                )
                rows.append(row)

    panel_csv = output_dir / "fig4_panel_data.csv"
    write_panel_csv(rows, panel_csv)

    skip_path = output_dir / "fig4_skipped.json"
    skip_path.write_text(json.dumps(
        {"skipped": skipped,
         "ran": available_tools,
         "n_cells": cell_count},
        indent=2))

    metadata = {
        "config": {
            "n_diploid_values": list(config.n_diploid_values),
            "n_replicates": config.n_replicates,
            "sequence_length": config.sequence_length,
            "recomb_rate": config.recomb_rate,
            "mut_rate": config.mut_rate,
            "population_size": config.population_size,
            "seed": config.seed,
        },
        "tools_requested": list(tools),
        "tools_ran": available_tools,
        "tools_skipped": skipped,
        "panel_csv": str(panel_csv),
        "n_rows": len(rows),
        "notes": (
            "v1: tskit + egrm scaling only; PLINK/KING/GraphPop "
            "curves deferred (binaries / parameterised Java "
            "fixture not wired)."),
    }
    metadata_path = output_dir / "fig4_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2))

    return Fig4RunResult(
        output_dir=output_dir,
        panel_csv=panel_csv,
        skip_path=skip_path,
        metadata_path=metadata_path,
        rows=rows,
    )


def _run_one_cell(
    *,
    tool: str,
    cohort: CohortParams,
    trees_path: Path,
    vcf_path: Path,
    output_dir: Path,
    replicate: int,
) -> CellResult:
    """Run one wrapper on one cohort; return its profiling row."""
    output_dir.mkdir(parents=True, exist_ok=True)
    if tool in TOOL_TREESEQ:
        Runner, extra_kwargs = TOOL_TREESEQ[tool]
        runner = Runner()
        if tool == "tskit":
            result = runner.run(
                trees_path, output_dir,
                seed=replicate, graphpop_commit="fig4")
        elif tool == "egrm":
            result = runner.run(
                trees_path, output_dir,
                seed=replicate, graphpop_commit="fig4",
                **extra_kwargs)
        else:
            raise ValueError(f"no dispatch for tool {tool!r}")
        prof = result.profiling
    elif tool in TOOL_VCF:
        Runner, extra_kwargs = TOOL_VCF[tool]
        runner = Runner()
        if tool == "plink":
            result = runner.run(
                vcf_path, output_dir,
                input_kind="vcf",
                extra_args=["--bad-freqs"],  # tiny cohorts → PLINK gate-off
                seed=replicate, graphpop_commit="fig4",
                **extra_kwargs)
        elif tool == "king":
            result = runner.run(
                vcf_path, output_dir,
                input_kind="vcf",
                seed=replicate, graphpop_commit="fig4",
                **extra_kwargs)
        else:
            raise ValueError(f"no dispatch for tool {tool!r}")
        prof = result.profiling
    else:
        raise ValueError(f"unknown tool: {tool!r}")
    return CellResult(
        n_diploid=cohort.n_diploid,
        n_haploid=cohort.n_haploid,
        sequence_length=cohort.sequence_length,
        tool=tool, replicate=replicate,
        wall_clock_s=prof.wall_clock_s,
        rss_peak_mb=prof.rss_peak_mb,
        user_cpu_s=prof.user_cpu_s,
        system_cpu_s=prof.system_cpu_s,
        exit_code=prof.exit_code,
        backend=prof.backend,
        receipt_path=str(output_dir / "receipt.json"),
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig4_out")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--n-diploid", type=int, nargs="+",
                   default=[25, 50, 100, 250])
    p.add_argument("--n-replicates", type=int, default=3)
    p.add_argument("--sequence-length", type=int, default=30_000)
    p.add_argument("--tools", nargs="+",
                   default=["tskit", "egrm"],
                   choices=ALL_TOOLS)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    config = Fig4SweepConfig(
        n_diploid_values=args.n_diploid,
        n_replicates=args.n_replicates,
        sequence_length=args.sequence_length,
        seed=args.seed,
    )

    def _progress(idx, total, n_diploid, tool, rep):
        if not args.quiet:
            print(
                f"[fig4] {idx}/{total} "
                f"n_diploid={n_diploid} tool={tool} rep={rep}",
                file=sys.stderr)

    result = run_fig4(
        config=config, output_dir=args.output_dir,
        tools=args.tools,
        progress_callback=None if args.quiet else _progress,
    )
    print(f"[fig4] DONE → {result.output_dir}", file=sys.stderr)
    print(f"  panel = {result.panel_csv}", file=sys.stderr)
    print(f"  skip  = {result.skip_path}", file=sys.stderr)
    print(f"  meta  = {result.metadata_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
