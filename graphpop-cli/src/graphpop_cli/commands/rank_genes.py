"""graphpop rank-genes — rank genes by composite selection evidence."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.command("rank-genes")
@click.option("--pop", "population", required=True, help="Population name")
@click.option("--pop2", help="Second population (for xpehh)")
@click.option("--chr", "chromosome", help="Restrict to chromosome")
@click.option("--top", type=int, default=100, help="Number of top genes (default: 100)")
@click.option("--sort-by", "sort_by", default="composite",
              type=click.Choice(["composite", "max_abs_ihs", "max_abs_xpehh",
                                  "max_h12", "mean_fst", "n_high_impact"]),
              help="Sort criterion (default: composite)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def rank_genes(ctx, population, pop2, chromosome, top, sort_by,
               output_path, fmt):
    """Rank genes by composite selection evidence.

    For each gene computes:
      - max_abs_ihs:   max |iHS| across variants in the gene
      - max_abs_xpehh: max |XP-EHH| (if --pop2 provided)
      - max_h12:       max H12 from overlapping GenomicWindow
      - mean_fst:      mean Fst from GenomicWindow overlapping the gene
      - n_high_impact: count of HIGH-impact variants

    Composite score = sum of per-stat percentile ranks (higher = stronger signal).

    \b
    Examples:
      graphpop rank-genes --pop EUR --top 50 -o top_genes.tsv
      graphpop rank-genes --pop GJ-tmp --pop2 GJ-trop --chr Chr01 --sort-by max_abs_ihs
      graphpop rank-genes --pop EUR --pop2 AFR --sort-by mean_fst --format json
    """
    # Dynamic property names cannot be parameterized — kept as f-strings.
    ihs_prop = f"ihs_{population}"
    xpehh_prop = f"xpehh_{population}_{pop2}" if pop2 else None

    # Build parameterized chromosome filter strings and params dict.
    chr_filter_v = "AND v.chr = $chromosome" if chromosome else ""
    chr_filter_g = "AND g.chr = $chromosome" if chromosome else ""

    params: dict = {}
    if chromosome:
        params["chromosome"] = chromosome

    # --- Query: per-gene iHS, high-impact count ---
    cypher_variant = f"""
    MATCH (v:Variant)-[hc:HAS_CONSEQUENCE]->(g:Gene)
    WHERE v.{ihs_prop} IS NOT NULL {chr_filter_v}
    WITH g,
         MAX(abs(v.{ihs_prop})) AS max_abs_ihs,
         {'MAX(abs(v.' + xpehh_prop + ')) AS max_abs_xpehh,' if xpehh_prop else ''}
         SUM(CASE WHEN hc.impact = 'HIGH' THEN 1 ELSE 0 END) AS n_high_impact,
         COUNT(DISTINCT v) AS n_variants
    RETURN g.symbol AS gene,
           g.geneId AS gene_id,
           g.chr AS chr,
           g.start AS gene_start,
           g.end AS gene_end,
           max_abs_ihs,
           {'max_abs_xpehh,' if xpehh_prop else ''}
           n_high_impact,
           n_variants
    """
    variant_records = ctx.run(cypher_variant, params)

    if not variant_records:
        # Fallback: try without iHS requirement
        click.echo("No iHS data found; querying genes by annotation only.", err=True)
        cypher_fallback = f"""
        MATCH (v:Variant)-[hc:HAS_CONSEQUENCE]->(g:Gene)
        WHERE TRUE {chr_filter_v}
        WITH g,
             SUM(CASE WHEN hc.impact = 'HIGH' THEN 1 ELSE 0 END) AS n_high_impact,
             COUNT(DISTINCT v) AS n_variants
        RETURN g.symbol AS gene,
               g.geneId AS gene_id,
               g.chr AS chr,
               g.start AS gene_start,
               g.end AS gene_end,
               0.0 AS max_abs_ihs,
               n_high_impact,
               n_variants
        """
        variant_records = ctx.run(cypher_fallback, params)

    if not variant_records:
        click.echo("No genes found.", err=True)
        return

    # Build gene lookup
    gene_data = {}
    for rec in variant_records:
        gene = rec.get("gene") or rec.get("gene_id")
        if not gene:
            continue
        gene_data[gene] = {
            "gene": gene,
            "gene_id": rec.get("gene_id", ""),
            "chr": rec.get("chr", ""),
            "gene_start": rec.get("gene_start", 0),
            "gene_end": rec.get("gene_end", 0),
            "max_abs_ihs": rec.get("max_abs_ihs", 0) or 0,
            "max_abs_xpehh": rec.get("max_abs_xpehh", 0) or 0,
            "n_high_impact": rec.get("n_high_impact", 0) or 0,
            "n_variants": rec.get("n_variants", 0) or 0,
            "max_h12": 0.0,
            "mean_fst": 0.0,
        }

    # --- Query: window-based stats (H12, Fst) overlapping genes ---
    click.echo("Querying window-based statistics...", err=True)
    window_params: dict = {
        "gene_names": list(gene_data.keys()),
        "population": population,
    }
    if chromosome:
        window_params["chromosome"] = chromosome

    cypher_windows = f"""
    MATCH (g:Gene)
    WHERE g.symbol IN $gene_names {chr_filter_g}
    MATCH (w:GenomicWindow)
    WHERE w.chr = g.chr AND w.population = $population
      AND w.start <= g.end AND w.end >= g.start
    WITH g, MAX(w.h12) AS max_h12, AVG(w.fst) AS mean_fst
    RETURN g.symbol AS gene, max_h12, mean_fst
    """
    try:
        window_records = ctx.run(cypher_windows, window_params)
        for rec in window_records:
            gene = rec.get("gene")
            if gene in gene_data:
                gene_data[gene]["max_h12"] = rec.get("max_h12", 0) or 0
                gene_data[gene]["mean_fst"] = rec.get("mean_fst", 0) or 0
    except SystemExit:
        click.echo("Warning: could not query GenomicWindow stats.", err=True)

    # --- Compute composite score (sum of ranks) ---
    genes = list(gene_data.values())
    stat_cols = ["max_abs_ihs", "max_abs_xpehh", "max_h12", "mean_fst", "n_high_impact"]

    for col in stat_cols:
        vals = sorted(set(g[col] for g in genes))
        rank_map = {v: i for i, v in enumerate(vals)}
        n = max(len(vals) - 1, 1)
        for g in genes:
            g[f"_rank_{col}"] = rank_map[g[col]] / n if n > 0 else 0

    for g in genes:
        g["composite"] = sum(g[f"_rank_{col}"] for col in stat_cols)

    # Clean up internal rank columns
    for g in genes:
        for col in stat_cols:
            del g[f"_rank_{col}"]

    # Sort
    if sort_by == "composite":
        genes.sort(key=lambda g: g["composite"], reverse=True)
    else:
        genes.sort(key=lambda g: g.get(sort_by, 0), reverse=True)

    # Top N
    genes = genes[:top]

    click.echo(f"Ranked {len(genes)} genes by {sort_by}.", err=True)
    format_output(genes, output_path, fmt, "rank-genes",
                  {"pop": population, "pop2": pop2, "chr": chromosome,
                   "sort_by": sort_by, "top": top})
