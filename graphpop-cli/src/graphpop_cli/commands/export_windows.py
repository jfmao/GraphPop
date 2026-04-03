"""graphpop export-windows — batch export GenomicWindow nodes to TSV."""
from __future__ import annotations

from pathlib import Path

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.command("export-windows")
@click.argument("chr", required=False)
@click.argument("population", required=False)
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--min-pi", type=float, help="Minimum pi filter")
@click.option("--max-pi", type=float, help="Maximum pi filter")
@click.option("--min-fst", type=float, help="Minimum Fst filter")
@click.option("--min-tajima-d", type=float, help="Minimum Tajima's D filter")
@click.option("--max-tajima-d", type=float, help="Maximum Tajima's D filter")
@click.option("--run-id", help="Filter by specific run ID")
@click.option("--limit", type=int, help="Maximum number of windows to return")
@pass_ctx
def export_windows(ctx, chr, population, output_path, fmt,
                   min_pi, max_pi, min_fst, min_tajima_d, max_tajima_d,
                   run_id, limit):
    """Export GenomicWindow nodes from the graph as TSV.

    Query persisted genome scan results (GenomicWindow nodes) with optional
    filters. Without arguments, exports all windows. With CHR and POPULATION,
    exports windows for that combination.

    \b
    Examples:
      graphpop export-windows                          # all windows
      graphpop export-windows chr22 EUR -o windows.tsv # specific region
      graphpop export-windows --min-fst 0.5            # high-Fst windows
      graphpop export-windows chr1 AFR --max-tajima-d -2  # negative Tajima's D
    """
    # Build Cypher query with parameterized values to prevent injection
    where_clauses = []
    params: dict = {}

    if chr:
        where_clauses.append("w.chr = $chr")
        params["chr"] = chr
    if population:
        where_clauses.append("w.population = $population")
        params["population"] = population
    if min_pi is not None:
        where_clauses.append(f"w.pi >= {min_pi}")
    if max_pi is not None:
        where_clauses.append(f"w.pi <= {max_pi}")
    if min_fst is not None:
        where_clauses.append(f"w.fst >= {min_fst}")
    if min_tajima_d is not None:
        where_clauses.append(f"w.tajima_d >= {min_tajima_d}")
    if max_tajima_d is not None:
        where_clauses.append(f"w.tajima_d <= {max_tajima_d}")
    if run_id:
        where_clauses.append("w.run_id = $run_id")
        params["run_id"] = run_id

    where = " AND ".join(where_clauses) if where_clauses else "TRUE"
    limit_clause = " LIMIT $limit" if limit else ""
    if limit:
        params["limit"] = limit

    cypher = (
        f"MATCH (w:GenomicWindow) WHERE {where} "
        f"RETURN w.windowId AS window_id, w.chr AS chr, "
        f"w.start AS start, w.end AS end, "
        f"w.population AS population, w.run_id AS run_id, "
        f"w.n_variants AS n_variants, w.n_segregating AS n_segregating, "
        f"w.pi AS pi, w.theta_w AS theta_w, w.tajima_d AS tajima_d, "
        f"w.fst AS fst, w.fst_wc AS fst_wc, w.dxy AS dxy, "
        f"w.pbs AS pbs, w.fay_wu_h AS fay_wu_h "
        f"ORDER BY w.chr, w.start"
        f"{limit_clause}"
    )

    records = ctx.run(cypher, params)

    if not records:
        click.echo("No windows found matching criteria.", err=True)
        return

    click.echo(f"Exporting {len(records)} windows...", err=True)
    format_output(records, output_path, fmt, "export-windows",
                  {"chr": chr, "population": population})
