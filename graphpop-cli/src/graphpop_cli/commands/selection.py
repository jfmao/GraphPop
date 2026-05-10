"""graphpop selection — ARG-aware selection scans (M8)."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def selection():
    """ARG-aware selection scans (Phase 5).

    Subcommands:

      allele-age-scan      Per-variant z-score of allele age | freq
                           (sweep signal: too young for frequency).
      branch-outlier-scan  Per-window z-score of total branch length
                           (Speidel et al. 2019 sweep detector).
    """


@selection.command("allele-age-scan")
@click.argument("run_id")
@click.option("--n-freq-bins", type=int, default=20, show_default=True,
              help="Number of logit-spaced frequency bins")
@click.option("--min-freq", type=float, default=0.05, show_default=True,
              help="Lowest frequency to bin")
@click.option("--max-freq", type=float, default=0.95, show_default=True,
              help="Highest frequency to bin")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def allele_age_scan(ctx, run_id, n_freq_bins, min_freq, max_freq,
                    output_path, fmt):
    """Per-variant allele-age z-score conditional on frequency (M8).

    Stratifies all variants by derived-allele frequency, computes
    per-bin (mean, sd) of log-age, and emits a z-score per variant.
    Negative z = candidate sweep.
    """
    opts: dict[str, object] = {
        "n_freq_bins": n_freq_bins,
        "min_freq": min_freq,
        "max_freq": max_freq,
    }
    cypher = build_cypher(
        "graphpop.selection.allele_age_scan",
        [f"'{run_id}'"],
        options=opts,
        yield_cols=["variant_id", "freq", "n_carriers", "n_samples",
                    "age_midpoint", "log_age", "bin_index", "bin_n",
                    "bin_mean_log_age", "bin_sd_log_age", "z_score",
                    "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "selection-allele-age-scan",
                  {"run_id": run_id, "n_freq_bins": n_freq_bins})


@selection.command("branch-outlier-scan")
@click.argument("run_id")
@click.argument("sample_ids")
@click.option("--window-size", type=int, default=10_000, show_default=True,
              help="Window size in bp")
@click.option("--step", type=int, default=None,
              help="Step size in bp (default: window-size)")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def branch_outlier_scan(ctx, run_id, sample_ids, window_size, step,
                        output_path, fmt):
    """Per-window branch-length outlier scan (M8).

    SAMPLE_IDS is comma-separated. Strong negative z-score = sweep
    candidate (lineages pulled tight in the focal set).
    """
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    opts: dict[str, object] = {"window_size": window_size}
    if step is not None:
        opts["step"] = step

    cypher = build_cypher(
        "graphpop.selection.branch_outlier_scan",
        [f"'{run_id}'", _cypher_str_list(sids)],
        options=opts,
        yield_cols=["start", "end", "total_branch_length", "mean_total",
                    "sd_total", "z_score", "n_samples", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "selection-branch-outlier-scan",
                  {"run_id": run_id, "n_samples": len(sids),
                   "window_size": window_size})


def _cypher_str_list(items: list[str]) -> str:
    inner = ", ".join(f"'{s}'" for s in items)
    return f"[{inner}]"
