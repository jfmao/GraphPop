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


@recombination.command("ldhat-mcmc")
@click.argument("sample_ids")
@click.option("--window-size", type=int, default=10_000, show_default=True)
@click.option("--step", type=int, default=None)
@click.option("--max-pair-distance", type=int, default=5_000, show_default=True)
@click.option("--min-maf", type=float, default=0.05, show_default=True)
@click.option("--n-iter", type=int, default=5_000, show_default=True)
@click.option("--burn-in", type=int, default=1_000, show_default=True)
@click.option("--prop-sd", type=float, default=0.5, show_default=True,
              help="M-H proposal std on log10(ρ)")
@click.option("--sigma", type=float, default=0.1, show_default=True,
              help="Likelihood Gaussian std on r² residuals")
@click.option("--seed", type=int, default=42, show_default=True)
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def ldhat_mcmc(ctx, sample_ids, window_size, step, max_pair_distance,
               min_maf, n_iter, burn_in, prop_sd, sigma, seed,
               output_path, fmt):
    """Per-window LDhat-style MCMC posterior on ρ (M13.A)."""
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    opts: dict[str, object] = {
        "window_size": window_size,
        "max_pair_distance": max_pair_distance,
        "min_maf": min_maf,
        "n_iter": n_iter,
        "burn_in": burn_in,
        "prop_sd": prop_sd,
        "sigma": sigma,
        "seed": seed,
    }
    if step is not None:
        opts["step"] = step
    cypher = build_cypher(
        "graphpop.recombination.ldhat_mcmc",
        [_cypher_str_list(sids)],
        options=opts,
        yield_cols=["start", "end", "n_variant_pairs",
                    "rho_posterior_mean", "rho_lower_2_5",
                    "rho_upper_97_5", "n_iter", "n_accepted",
                    "n_samples", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "recombination-ldhat-mcmc",
                  {"n_samples": len(sids), "n_iter": n_iter})


@recombination.command("hmm-smooth")
@click.argument("sample_ids")
@click.option("--window-size", type=int, default=10_000, show_default=True)
@click.option("--step", type=int, default=None)
@click.option("--max-pair-distance", type=int, default=5_000, show_default=True)
@click.option("--n-states", type=int, default=20, show_default=True)
@click.option("--state-log-lo", type=float, default=-10.0, show_default=True)
@click.option("--state-log-hi", type=float, default=-4.0, show_default=True)
@click.option("--emission-sd", type=float, default=0.5, show_default=True)
@click.option("--switch-rate", type=float, default=0.1, show_default=True)
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def hmm_smooth(ctx, sample_ids, window_size, step, max_pair_distance,
               n_states, state_log_lo, state_log_hi, emission_sd,
               switch_rate, output_path, fmt):
    """Pyrho-style HMM smoothing of per-window ρ (M13.B)."""
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    opts: dict[str, object] = {
        "window_size": window_size,
        "max_pair_distance": max_pair_distance,
        "n_states": n_states,
        "state_log_lo": state_log_lo,
        "state_log_hi": state_log_hi,
        "emission_sd": emission_sd,
        "switch_rate": switch_rate,
    }
    if step is not None:
        opts["step"] = step
    cypher = build_cypher(
        "graphpop.recombination.hmm_smooth",
        [_cypher_str_list(sids)],
        options=opts,
        yield_cols=["start", "end", "rho_per_bp", "rho_smoothed",
                    "hmm_state", "n_variant_pairs", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "recombination-hmm-smooth",
                  {"n_samples": len(sids), "n_states": n_states})


@recombination.command("stratified-ld-decay")
@click.argument("sample_ids")
@click.option("--stratify-by", default="population", show_default=True,
              type=click.Choice(["population", "sex"]))
@click.option("--window-size", type=int, default=10_000, show_default=True)
@click.option("--step", type=int, default=None)
@click.option("--max-pair-distance", type=int, default=5_000, show_default=True)
@click.option("--min-maf", type=float, default=0.05, show_default=True)
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def stratified_ld_decay(ctx, sample_ids, stratify_by, window_size, step,
                         max_pair_distance, min_maf, output_path, fmt):
    """Stratified Hudson-Kaplan ρ-map (M13.C)."""
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    opts: dict[str, object] = {
        "stratify_by": stratify_by,
        "window_size": window_size,
        "max_pair_distance": max_pair_distance,
        "min_maf": min_maf,
    }
    if step is not None:
        opts["step"] = step
    cypher = build_cypher(
        "graphpop.recombination.stratified_ld_decay",
        [_cypher_str_list(sids)],
        options=opts,
        yield_cols=["start", "end", "stratum", "n_variant_pairs",
                    "mean_r2", "mean_pair_distance", "rho_per_bp",
                    "n_samples", "stratify_by", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt,
                  "recombination-stratified-ld-decay",
                  {"n_samples": len(sids), "stratify_by": stratify_by})


@recombination.command("hotspots")
@click.argument("sample_ids")
@click.option("--method", default="ld_decay", show_default=True,
              type=click.Choice(["ld_decay"]))
@click.option("--window-size", type=int, default=10_000, show_default=True)
@click.option("--step", type=int, default=None)
@click.option("--max-pair-distance", type=int, default=5_000, show_default=True)
@click.option("--min-maf", type=float, default=0.05, show_default=True)
@click.option("--fdr-q", type=float, default=0.05, show_default=True)
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def hotspots(ctx, sample_ids, method, window_size, step,
              max_pair_distance, min_maf, fdr_q, output_path, fmt):
    """Recombination-hotspot detection with BH FDR (M13.D)."""
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    opts: dict[str, object] = {
        "method": method,
        "window_size": window_size,
        "max_pair_distance": max_pair_distance,
        "min_maf": min_maf,
        "fdr_q": fdr_q,
    }
    if step is not None:
        opts["step"] = step
    cypher = build_cypher(
        "graphpop.recombination.hotspots",
        [_cypher_str_list(sids)],
        options=opts,
        yield_cols=["start", "end", "rho_per_bp", "z_score", "p_value",
                    "adj_p_value", "is_hotspot", "n_variant_pairs",
                    "method", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "recombination-hotspots",
                  {"n_samples": len(sids), "fdr_q": fdr_q})
