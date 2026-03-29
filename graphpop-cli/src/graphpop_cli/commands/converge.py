"""graphpop converge — find regions where multiple selection statistics converge."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..formatters import format_output


# Statistics stored on Variant nodes vs GenomicWindow nodes
VARIANT_STATS = {"ihs", "xpehh", "nsl"}
WINDOW_STATS = {"h12", "fst", "pi", "tajima_d"}


@click.command("converge")
@click.option("--stats", required=True,
              help="Comma-separated statistic names (ihs, xpehh, nsl, h12, fst, pi, tajima_d)")
@click.option("--thresholds", required=True,
              help="Comma-separated threshold values (matched positionally to stats)")
@click.option("--chr", "chromosome", help="Chromosome (optional, all if not specified)")
@click.option("--pop", "population", required=True, help="Population name")
@click.option("--pop2", help="Second population (for xpehh, fst)")
@click.option("--window", type=int, default=0,
              help="Aggregate into windows of this size (default: per-variant)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--limit", type=int, default=10000, help="Maximum rows (default: 10000)")
@pass_ctx
def converge(ctx, stats, thresholds, chromosome, population, pop2,
             window, output_path, fmt, limit):
    """Find genomic regions where multiple statistics exceed thresholds.

    Identifies convergent selection signals by requiring multiple persisted
    statistics to simultaneously exceed user-defined thresholds.

    \b
    For variant-based stats (ihs, xpehh, nsl): queries Variant nodes.
    For window-based stats (h12, fst, pi, tajima_d): queries GenomicWindow nodes.
    If mixed, runs two queries and merges results by position.

    \b
    Examples:
      graphpop converge --stats ihs,xpehh --thresholds 2.0,2.0 --pop EUR --pop2 AFR
      graphpop converge --stats ihs,nsl --thresholds 2.0,2.0 --chr chr22 --pop EUR -o conv.tsv
      graphpop converge --stats h12,fst --thresholds 0.3,0.5 --pop GJ-tmp --window 100000
      graphpop converge --stats ihs,xpehh,h12,fst --thresholds 2.0,2.0,0.3,0.5 --pop EUR --pop2 AFR
    """
    stat_list = [s.strip() for s in stats.split(",")]
    thresh_list = [float(t.strip()) for t in thresholds.split(",")]

    if len(stat_list) != len(thresh_list):
        click.echo("Error: --stats and --thresholds must have the same number of items.", err=True)
        raise SystemExit(1)

    requested_variant = [s for s in stat_list if s in VARIANT_STATS]
    requested_window = [s for s in stat_list if s in WINDOW_STATS]
    unknown = [s for s in stat_list if s not in VARIANT_STATS and s not in WINDOW_STATS]
    if unknown:
        click.echo(f"Warning: unknown statistics ignored: {unknown}", err=True)

    stat_thresh = dict(zip(stat_list, thresh_list))

    results = []

    # --- Query variant-based stats ---
    if requested_variant:
        where_parts = []
        if chromosome:
            where_parts.append(f"v.chr = '{chromosome}'")

        return_cols = [
            "v.variantId AS variant_id",
            "v.chr AS chr",
            "v.pos AS pos",
        ]

        for stat in requested_variant:
            prop = _prop_name(stat, population, pop2)
            if prop is None:
                click.echo(f"Warning: skipping {stat} (need --pop2 for xpehh)", err=True)
                continue
            thresh = stat_thresh[stat]
            where_parts.append(f"abs(v.{prop}) >= {thresh}")
            return_cols.append(f"v.{prop} AS {stat}")

        # Join to gene annotation
        cypher = (
            f"MATCH (v:Variant) "
            f"WHERE {' AND '.join(where_parts)} "
            f"OPTIONAL MATCH (v)-[:HAS_CONSEQUENCE]->(g:Gene) "
            f"RETURN DISTINCT {', '.join(return_cols)}, "
            f"g.symbol AS gene "
            f"ORDER BY v.pos LIMIT {limit}"
        )
        variant_records = ctx.run(cypher)
        results.extend(variant_records)

    # --- Query window-based stats ---
    if requested_window:
        where_parts = [f"w.population = '{population}'"]
        if chromosome:
            where_parts.append(f"w.chr = '{chromosome}'")

        return_cols = [
            "w.windowId AS window_id",
            "w.chr AS chr",
            "w.start AS start",
            "w.end AS end",
        ]

        for stat in requested_window:
            prop = stat
            thresh = stat_thresh[stat]
            if stat == "fst" and pop2:
                prop = f"fst_{population}_{pop2}"
            if stat in ("h12",):
                where_parts.append(f"w.{prop} >= {thresh}")
            elif stat in ("tajima_d",):
                # Tajima's D: extreme negative = selection
                where_parts.append(f"w.{prop} <= -{thresh}")
            else:
                where_parts.append(f"w.{prop} >= {thresh}")
            return_cols.append(f"w.{prop} AS {stat}")

        cypher = (
            f"MATCH (w:GenomicWindow) "
            f"WHERE {' AND '.join(where_parts)} "
            f"RETURN {', '.join(return_cols)} "
            f"ORDER BY w.start LIMIT {limit}"
        )
        window_records = ctx.run(cypher)
        results.extend(window_records)

    if not results:
        click.echo("No convergent signals found with the given thresholds.", err=True)
        return

    click.echo(f"Found {len(results)} convergent records.", err=True)
    format_output(results, output_path, fmt, "converge",
                  {"stats": stats, "thresholds": thresholds,
                   "chr": chromosome, "pop": population, "pop2": pop2})


def _prop_name(stat: str, population: str, pop2: str | None) -> str | None:
    """Build the Neo4j property name for a given statistic."""
    if stat in ("ihs", "nsl"):
        return f"{stat}_{population}"
    elif stat == "xpehh":
        if not pop2:
            return None
        return f"xpehh_{population}_{pop2}"
    return stat
