"""graphpop filter — query persisted results with annotation filters."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.command("filter")
@click.argument("statistic", type=click.Choice([
    "ihs", "xpehh", "nsl", "fst", "pi", "tajima_d", "h12",
]))
@click.argument("chr")
@click.argument("population")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--consequence", help="Filter by VEP consequence type (e.g., missense_variant)")
@click.option("--pathway", help="Filter by pathway name")
@click.option("--gene", help="Filter by gene name or ID")
@click.option("--min-score", type=float, help="Minimum absolute score")
@click.option("--max-score", type=float, help="Maximum absolute score")
@click.option("--pop2", help="Second population (for xpehh)")
@click.option("--limit", type=int, default=10000, help="Maximum rows (default: 10000)")
@pass_ctx
def filter_results(ctx, statistic, chr, population, output_path, fmt,
                   consequence, pathway, gene, min_score, max_score, pop2, limit):
    """Query persisted statistics with annotation-based filters.

    This command retrieves already-computed statistics (iHS, XP-EHH, nSL, etc.)
    from graph nodes and filters them by functional annotation. It is the
    recommended way to perform conditioned analysis for haplotype-based
    statistics, which must be computed genome-wide first and then filtered.

    \b
    Workflow:
      1. Compute statistics:  graphpop ihs chr1 EUR --persist
      2. Filter by annotation: graphpop filter ihs chr1 EUR --consequence missense_variant

    \b
    Examples:
      graphpop filter ihs chr1 EUR --consequence missense_variant -o ihs_missense.tsv
      graphpop filter xpehh chr1 EUR --pop2 AFR --pathway "Cardiac repolarization"
      graphpop filter nsl chr1 GJ-tmp --gene GW5 --min-score 2.0
      graphpop filter h12 chr1 GJ-tmp --consequence missense_variant
    """
    # Build the property name for this statistic
    if statistic == "xpehh" and pop2:
        prop = f"xpehh_{population}_{pop2}"
        prop_unstd = f"xpehh_unstd_{population}_{pop2}"
    elif statistic == "xpehh":
        # Try to find any xpehh property
        prop = f"xpehh_{population}_*"
        click.echo("Warning: --pop2 not specified; will search for any XP-EHH involving this population.", err=True)
        prop = None
    elif statistic in ("ihs", "nsl"):
        prop = f"{statistic}_{population}"
        prop_unstd = f"{statistic}_unstd_{population}"
    elif statistic in ("fst", "pi", "tajima_d", "h12"):
        prop = statistic
        prop_unstd = None
    else:
        prop = statistic
        prop_unstd = None

    # Build Cypher query
    if statistic in ("ihs", "xpehh", "nsl"):
        # Per-variant statistics stored on Variant nodes
        match_clause = "MATCH (v:Variant)"
        where_parts = [f"v.chr = '{chr}'"]
        if prop:
            where_parts.append(f"v.{prop} IS NOT NULL")

        # Annotation join
        if consequence:
            match_clause += "-[:HAS_CONSEQUENCE]->(hc)"
            where_parts.append(f"hc.consequence = '{consequence}'")
        if pathway:
            match_clause += "-[:HAS_CONSEQUENCE]->(:Gene)-[:IN_PATHWAY]->(pw:Pathway)"
            where_parts.append(f"pw.name CONTAINS '{pathway}'")
        if gene:
            match_clause += "-[:HAS_CONSEQUENCE]->(g:Gene)"
            where_parts.append(f"(g.geneId = '{gene}' OR g.symbol = '{gene}')")

        if min_score is not None and prop:
            where_parts.append(f"abs(v.{prop}) >= {min_score}")
        if max_score is not None and prop:
            where_parts.append(f"abs(v.{prop}) <= {max_score}")

        return_cols = [
            "v.variantId AS variant_id",
            "v.pos AS pos",
        ]
        if prop:
            return_cols.append(f"v.{prop} AS {statistic}")
        if prop_unstd:
            return_cols.append(f"v.{prop_unstd} AS {statistic}_unstd")
        if consequence:
            return_cols.append("hc.consequence AS consequence")
            return_cols.append("hc.impact AS impact")
        if gene:
            return_cols.append("g.symbol AS gene")

        cypher = (
            f"{match_clause} "
            f"WHERE {' AND '.join(where_parts)} "
            f"RETURN DISTINCT {', '.join(return_cols)} "
            f"ORDER BY v.pos LIMIT {limit}"
        )

    elif statistic == "h12":
        # Garud's H stored on GenomicWindow nodes
        match_clause = "MATCH (w:GenomicWindow)"
        where_parts = [
            f"w.chr = '{chr}'",
            f"w.population = '{population}'",
        ]
        if min_score is not None:
            where_parts.append(f"w.h12 >= {min_score}")

        cypher = (
            f"{match_clause} "
            f"WHERE {' AND '.join(where_parts)} "
            f"RETURN w.windowId AS window_id, w.chr AS chr, "
            f"w.start AS start, w.end AS end, "
            f"w.h12 AS h12, w.h2_h1 AS h2_h1, w.hap_diversity AS hap_div "
            f"ORDER BY w.h12 DESC LIMIT {limit}"
        )

    else:
        # Window-level statistics (fst, pi, tajima_d)
        match_clause = "MATCH (w:GenomicWindow)"
        where_parts = [
            f"w.chr = '{chr}'",
            f"w.population = '{population}'",
        ]
        if min_score is not None:
            where_parts.append(f"w.{prop} >= {min_score}")
        if max_score is not None:
            where_parts.append(f"w.{prop} <= {max_score}")

        cypher = (
            f"{match_clause} "
            f"WHERE {' AND '.join(where_parts)} "
            f"RETURN w.windowId AS window_id, w.start AS start, w.end AS end, "
            f"w.{prop} AS {statistic}, w.n_variants AS n_variants "
            f"ORDER BY w.start LIMIT {limit}"
        )

    records = ctx.run(cypher)

    if not records:
        click.echo(f"No results found for {statistic} on {chr}/{population} "
                    f"with given filters.", err=True)
        return

    click.echo(f"Found {len(records)} records.", err=True)
    format_output(records, output_path, fmt, "filter",
                  {"statistic": statistic, "chr": chr, "population": population,
                   "consequence": consequence, "pathway": pathway, "gene": gene})
