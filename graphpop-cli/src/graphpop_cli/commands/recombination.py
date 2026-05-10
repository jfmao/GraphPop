"""graphpop recombination — recombination-map inference (M11)."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def recombination():
    """Recombination-rate inference (Phase 5).

    Subcommands:

      arg-breakpoints  Per-window ρ from ARG breakpoint density.
      ld-decay         Per-window ρ from Hudson-Kaplan moment estimator
                       on pairwise r² over :CARRIES.
    """


@recombination.command("arg-breakpoints")
@click.argument("run_id")
@click.option("--window-size", type=int, default=10_000, show_default=True,
              help="Window size in bp")
@click.option("--step", type=int, default=None,
              help="Step size in bp (default: window-size)")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def arg_breakpoints(ctx, run_id, window_size, step, output_path, fmt):
    """ARG-derived per-window recombination rate (M11)."""
    opts: dict[str, object] = {"window_size": window_size}
    if step is not None:
        opts["step"] = step
    cypher = build_cypher(
        "graphpop.recombination.arg_breakpoints",
        [f"'{run_id}'"],
        options=opts,
        yield_cols=["start", "end", "n_breakpoints", "n_marginal_trees",
                    "total_branch_length", "rho_per_bp", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "recombination-arg-breakpoints",
                  {"run_id": run_id, "window_size": window_size})


@recombination.command("ld-decay")
@click.argument("sample_ids")
@click.option("--window-size", type=int, default=10_000, show_default=True,
              help="Window size in bp")
@click.option("--step", type=int, default=None,
              help="Step size in bp (default: window-size)")
@click.option("--min-maf", type=float, default=0.05, show_default=True,
              help="Minimum derived-allele frequency to include a variant")
@click.option("--max-pair-distance", type=int, default=5_000, show_default=True,
              help="Maximum pair distance (bp) for r² computation")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def ld_decay(ctx, sample_ids, window_size, step, min_maf,
             max_pair_distance, output_path, fmt):
    """LD-decay-based per-window recombination rate (M11).

    SAMPLE_IDS is a comma-separated list; v1 requires ≤ 63 samples
    (uses a packed bitmask carrier representation).
    """
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    opts: dict[str, object] = {
        "window_size": window_size,
        "min_maf": min_maf,
        "max_pair_distance": max_pair_distance,
    }
    if step is not None:
        opts["step"] = step
    cypher = build_cypher(
        "graphpop.recombination.ld_decay",
        [_cypher_str_list(sids)],
        options=opts,
        yield_cols=["start", "end", "n_variant_pairs", "mean_r2",
                    "mean_pair_distance", "rho_per_bp", "n_samples",
                    "method", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "recombination-ld-decay",
                  {"n_samples": len(sids), "window_size": window_size})


def _cypher_str_list(items: list[str]) -> str:
    inner = ", ".join(f"'{s}'" for s in items)
    return f"[{inner}]"
