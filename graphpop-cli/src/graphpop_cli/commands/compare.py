"""graphpop compare — compare statistics between two populations."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.command("compare")
@click.argument("pop1")
@click.argument("pop2")
@click.argument("chr", metavar="CHR")
@click.option("--stat", required=True,
              type=click.Choice(["pi", "theta_w", "tajima_d", "fst", "ihs"]),
              help="Statistic to compare")
@click.option("--window-size", type=int, default=100000,
              help="Sliding window size for comparison (default: 100000)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--limit", type=int, default=10000, help="Maximum rows (default: 10000)")
@pass_ctx
def compare(ctx, pop1, pop2, chr, stat, window_size, output_path, fmt, limit):
    """Compare statistics between two populations across a chromosome.

    For window-based stats (pi, theta_w, tajima_d, fst), queries GenomicWindow
    nodes for both populations and joins by window position to compute delta.

    For ihs, queries Variant nodes with ihs_{POP1} and ihs_{POP2} properties
    and computes the per-variant difference.

    Output columns: window_start, window_end, stat_pop1, stat_pop2, delta, abs_delta.

    \b
    Examples:
      graphpop compare EUR AFR chr22 --stat pi -o delta_pi.tsv
      graphpop compare GJ-tmp GJ-trp Chr1 --stat fst -o delta.tsv
      graphpop compare EUR EAS chr22 --stat ihs -o ihs_diff.tsv
    """
    if stat == "ihs":
        records = _compare_variant_stat(ctx, pop1, pop2, chr, stat, limit)
    else:
        records = _compare_window_stat(ctx, pop1, pop2, chr, stat, window_size, limit)

    if not records:
        click.echo(f"No comparison data found for {stat} on {chr} "
                    f"({pop1} vs {pop2}).", err=True)
        return

    click.echo(f"Found {len(records)} comparison rows.", err=True)
    format_output(records, output_path, fmt, "compare",
                  {"pop1": pop1, "pop2": pop2, "chr": chr, "stat": stat})


def _compare_window_stat(ctx, pop1, pop2, chr, stat, window_size, limit):
    """Compare window-based statistics between two populations."""
    prop = stat

    cypher = (
        f"MATCH (w1:GenomicWindow) "
        f"WHERE w1.chr = '{chr}' AND w1.population = '{pop1}' "
        f"AND w1.{prop} IS NOT NULL "
        f"WITH w1 "
        f"MATCH (w2:GenomicWindow) "
        f"WHERE w2.chr = '{chr}' AND w2.population = '{pop2}' "
        f"AND w2.start = w1.start AND w2.end = w1.end "
        f"AND w2.{prop} IS NOT NULL "
        f"RETURN w1.start AS window_start, "
        f"w1.end AS window_end, "
        f"w1.{prop} AS {stat}_{pop1}, "
        f"w2.{prop} AS {stat}_{pop2}, "
        f"(w1.{prop} - w2.{prop}) AS delta, "
        f"abs(w1.{prop} - w2.{prop}) AS abs_delta "
        f"ORDER BY w1.start LIMIT {limit}"
    )
    return ctx.run(cypher)


def _compare_variant_stat(ctx, pop1, pop2, chr, stat, limit):
    """Compare variant-based statistics (ihs) between two populations."""
    prop1 = f"{stat}_{pop1}"
    prop2 = f"{stat}_{pop2}"

    cypher = (
        f"MATCH (v:Variant) "
        f"WHERE v.chr = '{chr}' "
        f"AND v.{prop1} IS NOT NULL AND v.{prop2} IS NOT NULL "
        f"RETURN v.pos AS pos, "
        f"v.variantId AS variantId, "
        f"v.{prop1} AS {stat}_{pop1}, "
        f"v.{prop2} AS {stat}_{pop2}, "
        f"(v.{prop1} - v.{prop2}) AS delta, "
        f"abs(v.{prop1} - v.{prop2}) AS abs_delta "
        f"ORDER BY v.pos LIMIT {limit}"
    )
    return ctx.run(cypher)
