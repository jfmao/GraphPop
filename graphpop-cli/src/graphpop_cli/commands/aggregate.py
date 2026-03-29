"""graphpop aggregate — aggregate results and generate summary tables."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import click

from ..cli import pass_ctx


@click.command()
@click.option("--results-dir", "-d", type=click.Path(exists=True), required=True,
              help="Directory with per-procedure TSV results (from run-all)")
@click.option("--json-results", "-j", type=click.Path(exists=True),
              help="JSON results file (from run-all)")
@click.option("--output-dir", "-o", type=click.Path(), default="graphpop_tables",
              help="Output directory for summary tables")
@pass_ctx
def aggregate(ctx, results_dir, json_results, output_dir):
    """Aggregate per-population results into summary tables.

    Reads TSV results from a run-all output directory and produces
    publication-ready summary tables:

    \b
      population_summary.tsv   — per-pop diversity, theta, Tajima's D, Fis
      fst_matrix.tsv           — pairwise Fst matrix
      pinpis.tsv               — piN/piS ratios (if conditioned results exist)
      selection_peaks.tsv      — top iHS/XP-EHH/nSL peaks per population
      roh_summary.tsv          — per-pop FROH statistics
    """
    results_path = Path(results_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load JSON results if provided
    all_results = {}
    if json_results:
        with open(json_results) as f:
            all_results = json.load(f)
        click.echo(f"Loaded {len(all_results)} results from JSON")

    # --- Table 1: Population Summary ---
    diversity_dir = results_path / "diversity"
    if diversity_dir.exists():
        click.echo("Generating population_summary.tsv...")
        pop_stats = _aggregate_single_row_tsv(diversity_dir)
        _write_summary(out_dir / "population_summary.tsv", pop_stats,
                       ["population", "chr", "pi", "theta_w", "tajima_d",
                        "het_exp", "het_obs", "fis", "n_variants", "n_segregating"])

    # --- Table 2: Fst Matrix ---
    divergence_dir = results_path / "divergence"
    if divergence_dir.exists():
        click.echo("Generating fst_matrix.tsv...")
        div_stats = _aggregate_single_row_tsv(divergence_dir)
        _write_summary(out_dir / "fst_matrix.tsv", div_stats,
                       ["pop1", "pop2", "chr", "fst_hudson", "fst_wc", "dxy", "da"])

    # --- Table 3: ROH Summary ---
    roh_dir = results_path / "roh"
    if roh_dir.exists():
        click.echo("Generating roh_summary.tsv...")
        roh_data = _aggregate_multi_row_tsv(roh_dir)
        # Compute per-population means
        pop_roh = {}
        for rec in roh_data:
            pop = rec.get("population", rec.get("file_pop", "unknown"))
            if pop not in pop_roh:
                pop_roh[pop] = {"n_samples": 0, "total_froh": 0.0,
                                "total_n_roh": 0, "max_froh": 0.0}
            pop_roh[pop]["n_samples"] += 1
            froh = float(rec.get("froh", 0))
            pop_roh[pop]["total_froh"] += froh
            pop_roh[pop]["total_n_roh"] += int(rec.get("n_roh", 0))
            pop_roh[pop]["max_froh"] = max(pop_roh[pop]["max_froh"], froh)

        rows = []
        for pop, s in sorted(pop_roh.items()):
            rows.append({
                "population": pop,
                "n_samples": s["n_samples"],
                "mean_froh": f"{s['total_froh'] / s['n_samples']:.6f}",
                "mean_n_roh": f"{s['total_n_roh'] / s['n_samples']:.1f}",
                "max_froh": f"{s['max_froh']:.6f}",
            })
        _write_dict_tsv(out_dir / "roh_summary.tsv", rows)

    # --- Table 4: Selection Peaks ---
    for proc in ("ihs", "nsl", "xpehh"):
        proc_dir = results_path / proc
        if proc_dir.exists():
            click.echo(f"Generating {proc}_peaks.tsv...")
            peaks = _extract_peaks(proc_dir, proc, top_n=100)
            _write_dict_tsv(out_dir / f"{proc}_peaks.tsv", peaks)

    # --- Table 5: Garud's H Sweep Windows ---
    garud_dir = results_path / "garud_h"
    if garud_dir.exists():
        click.echo("Generating sweep_windows.tsv...")
        sweeps = _extract_sweep_windows(garud_dir, h12_threshold=0.1)
        _write_dict_tsv(out_dir / "sweep_windows.tsv", sweeps)

    click.echo(f"\nSummary tables written to {out_dir}/")
    for f in sorted(out_dir.glob("*.tsv")):
        n_lines = sum(1 for _ in open(f)) - 1
        click.echo(f"  {f.name}: {n_lines} rows")


def _aggregate_single_row_tsv(directory: Path) -> list[dict]:
    """Read TSV files with single data row, extract pop/chr from filename."""
    rows = []
    for tsv in sorted(directory.glob("*.tsv")):
        parts = tsv.stem.split("_")
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for rec in reader:
                # Infer pop and chr from filename: POP_CHR.tsv
                if len(parts) >= 2:
                    rec["population"] = "_".join(parts[:-1])
                    rec["chr"] = parts[-1]
                elif "vs" in tsv.stem:
                    # Pairwise: POP1_vs_POP2_CHR.tsv
                    vs_idx = parts.index("vs")
                    rec["pop1"] = "_".join(parts[:vs_idx])
                    rec["pop2"] = "_".join(parts[vs_idx + 1:-1])
                    rec["chr"] = parts[-1]
                rows.append(rec)
    return rows


def _aggregate_multi_row_tsv(directory: Path) -> list[dict]:
    """Read TSV files with multiple data rows."""
    rows = []
    for tsv in sorted(directory.glob("*.tsv")):
        parts = tsv.stem.split("_")
        pop = "_".join(parts[:-1]) if len(parts) >= 2 else parts[0]
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for rec in reader:
                rec["file_pop"] = pop
                rows.append(rec)
    return rows


def _extract_peaks(directory: Path, stat_name: str,
                   top_n: int = 100) -> list[dict]:
    """Extract top peaks from per-variant result files."""
    all_variants = []
    for tsv in sorted(directory.glob("*.tsv")):
        parts = tsv.stem.split("_")
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for rec in reader:
                score = rec.get(stat_name, rec.get(f"{stat_name}_unstd", "0"))
                try:
                    rec["abs_score"] = abs(float(score))
                except (ValueError, TypeError):
                    rec["abs_score"] = 0
                rec["source_file"] = tsv.stem
                all_variants.append(rec)

    all_variants.sort(key=lambda r: r["abs_score"], reverse=True)
    return all_variants[:top_n]


def _extract_sweep_windows(directory: Path,
                           h12_threshold: float = 0.1) -> list[dict]:
    """Extract windows exceeding H12 threshold."""
    sweeps = []
    for tsv in sorted(directory.glob("*.tsv")):
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for rec in reader:
                try:
                    if float(rec.get("h12", 0)) >= h12_threshold:
                        sweeps.append(rec)
                except (ValueError, TypeError):
                    pass
    sweeps.sort(key=lambda r: float(r.get("h12", 0)), reverse=True)
    return sweeps


def _write_summary(path: Path, rows: list[dict], columns: list[str]):
    """Write summary table with specified columns."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_dict_tsv(path: Path, rows: list[dict]):
    """Write list of dicts as TSV."""
    if not rows:
        with open(path, "w") as f:
            f.write("# No results\n")
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
