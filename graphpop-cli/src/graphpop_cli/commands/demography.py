"""graphpop demography — population-size trajectory inference (M7)."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def demography():
    """Demographic inference (Phase 5).

    Subcommands:

      ne-trajectory  Closed-form Ne(t) by inverting graphpop.arg.coalescence_rate.
    """


@demography.command("ne-trajectory")
@click.argument("run_id")
@click.argument("sample_ids", required=False, default="")
@click.option("--time-bins", required=True,
              help="Comma-separated monotonically increasing bin edges, "
                   "e.g. 0,0.25,0.5,1,2,1e9")
@click.option("--ploidy", type=int, default=2, show_default=True,
              help="Sample ploidy (Ne = 1 / (ploidy * rate))")
@click.option("--population", default=None,
              help="Resolve sample list from :Sample.population instead of "
                   "an explicit SAMPLE_IDS argument")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def ne_trajectory(ctx, run_id, sample_ids, time_bins, ploidy, population,
                  output_path, fmt):
    """Closed-form per-bin Ne(t) trajectory (M7).

    SAMPLE_IDS is a comma-separated list (e.g. hap_0,hap_1,...,hap_19).
    Pass an empty string and use --population to resolve via
    :Sample.population.
    """
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    bins = [float(x) for x in time_bins.split(",")]
    opts: dict[str, object] = {"time_bins": bins, "ploidy": ploidy}
    if population:
        opts["population"] = population
    if not sids and not population:
        raise click.UsageError(
            "Either SAMPLE_IDS or --population must be provided")

    cypher = build_cypher(
        "graphpop.demography.ne_trajectory",
        [f"'{run_id}'", _cypher_str_list(sids)],
        options=opts,
        yield_cols=["time_lo", "time_hi", "n_coalescent_events",
                    "lineage_pair_time", "rate", "ne", "ne_se",
                    "flag", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "demography-ne-trajectory",
                  {"run_id": run_id, "n_samples": len(sids), "ploidy": ploidy})


def _cypher_str_list(items: list[str]) -> str:
    inner = ", ".join(f"'{s}'" for s in items)
    return f"[{inner}]"
