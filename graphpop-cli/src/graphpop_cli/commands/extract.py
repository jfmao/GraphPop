"""graphpop extract — extract variants, samples, and genotypes from the graph."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.group()
def extract():
    """Extract data from the graph: variants, samples, or genotypes.

    \b
    Subcommands:
      variants   Query Variant nodes with flexible filters
      samples    Query Sample nodes for a population
      genotypes  Extract sample x variant dosage matrix for a region

    \b
    Examples:
      graphpop extract variants --chr chr22 --pop EUR --consequence missense_variant -o variants.tsv
      graphpop extract samples --pop EUR -o samples.tsv
      graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 --pop EUR -o geno.tsv
    """
    pass


@extract.command("variants")
@click.option("--chr", "chromosome", help="Chromosome filter")
@click.option("--start", type=int, help="Region start position")
@click.option("--end", type=int, help="Region end position")
@click.option("--pop", "population", help="Population name (for AF lookup)")
@click.option("--min-af", type=float, help="Minimum allele frequency")
@click.option("--max-af", type=float, help="Maximum allele frequency")
@click.option("--consequence", help="VEP consequence type (e.g., missense_variant)")
@click.option("--pathway", help="Pathway name (substring match)")
@click.option("--gene", help="Gene symbol or ID")
@click.option("--fields", default="variantId,pos,ref,alt,af",
              help="Comma-separated fields to return (default: variantId,pos,ref,alt,af)")
@click.option("--limit", type=int, default=10000, help="Maximum rows (default: 10000)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def extract_variants(ctx, chromosome, start, end, population, min_af, max_af,
                     consequence, pathway, gene, fields, limit, output_path, fmt):
    """Query Variant nodes with optional filters.

    Builds Cypher dynamically based on provided filters. Use --fields to
    select which variant properties to return.

    \b
    Examples:
      graphpop extract variants --chr chr22 --pop EUR --consequence missense_variant -o variants.tsv
      graphpop extract variants --chr chr22 --start 16000000 --end 17000000 --fields pos,ref,alt,af,fst,ihs -o region.tsv
      graphpop extract variants --gene KCNE1 --pop EUR -o kcne1_variants.tsv
    """
    field_list = [f.strip() for f in fields.split(",")]

    # Build MATCH clause with optional annotation joins
    match_clause = "MATCH (v:Variant)"
    where_parts = []

    if consequence:
        match_clause += "-[:HAS_CONSEQUENCE]->(hc)"
        where_parts.append(f"hc.consequence = '{consequence}'")
    if gene:
        if "-[:HAS_CONSEQUENCE]->" not in match_clause:
            match_clause += "-[:HAS_CONSEQUENCE]->(g:Gene)"
        else:
            match_clause = match_clause.replace("->(hc)", "->(hc)-[:ON_GENE]->(g:Gene)")
            # Simpler: use separate MATCH
            match_clause = "MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)"
            if consequence:
                match_clause = "MATCH (v:Variant)-[hc_rel:HAS_CONSEQUENCE]->(g:Gene)"
                where_parts = [p for p in where_parts if "hc.consequence" not in p]
                where_parts.append(f"hc_rel.consequence = '{consequence}'")
        where_parts.append(f"(g.symbol = '{gene}' OR g.geneId = '{gene}')")
    if pathway:
        if "Gene" not in match_clause:
            match_clause += "-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(pw:Pathway)"
        else:
            match_clause += "-[:IN_PATHWAY]->(pw:Pathway)"
        where_parts.append(f"pw.name CONTAINS '{pathway}'")

    if chromosome:
        where_parts.append(f"v.chr = '{chromosome}'")
    if start is not None:
        where_parts.append(f"v.pos >= {start}")
    if end is not None:
        where_parts.append(f"v.pos <= {end}")

    # AF filtering: if population is given, look up index in pop_ids array
    if population and (min_af is not None or max_af is not None):
        # Use dynamic property access for population-specific AF
        if min_af is not None:
            where_parts.append(
                f"ANY(i IN range(0, size(v.pop_ids)-1) "
                f"WHERE v.pop_ids[i] = '{population}' AND v.af[i] >= {min_af})"
            )
        if max_af is not None:
            where_parts.append(
                f"ANY(i IN range(0, size(v.pop_ids)-1) "
                f"WHERE v.pop_ids[i] = '{population}' AND v.af[i] <= {max_af})"
            )

    # Build RETURN clause from requested fields
    return_cols = []
    for f in field_list:
        if f == "af" and population:
            return_cols.append(
                f"[i IN range(0, size(v.pop_ids)-1) "
                f"WHERE v.pop_ids[i] = '{population}' | v.af[i]][0] AS af_{population}"
            )
        elif f in ("gene", "gene_symbol") and "Gene" in match_clause:
            return_cols.append("g.symbol AS gene")
        elif f == "consequence" and consequence:
            col = "hc.consequence" if "hc)" in match_clause else "hc_rel.consequence"
            return_cols.append(f"{col} AS consequence")
        else:
            return_cols.append(f"v.{f} AS {f}")

    where_str = " AND ".join(where_parts) if where_parts else "true"
    cypher = (
        f"{match_clause} "
        f"WHERE {where_str} "
        f"RETURN DISTINCT {', '.join(return_cols)} "
        f"ORDER BY v.pos LIMIT {limit}"
    )

    records = ctx.run(cypher)
    if not records:
        click.echo("No variants found with given filters.", err=True)
        return

    click.echo(f"Found {len(records)} variants.", err=True)
    format_output(records, output_path, fmt, "extract variants",
                  {"chr": chromosome, "start": start, "end": end,
                   "pop": population, "consequence": consequence,
                   "pathway": pathway, "gene": gene})


@extract.command("samples")
@click.option("--pop", "population", required=True, help="Population name")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def extract_samples(ctx, population, output_path, fmt):
    """Query Sample nodes for a population.

    Returns sampleId, population, and packed_index for each sample. If
    population-level summary stats (e.g., FROH) are available, they are
    included.

    \b
    Examples:
      graphpop extract samples --pop EUR -o samples.tsv
      graphpop extract samples --pop GJ-tmp --format json
    """
    cypher = f"""
    MATCH (s:Sample)
    WHERE s.population = '{population}'
    OPTIONAL MATCH (p:Population {{name: '{population}'}})
    RETURN s.sampleId AS sampleId,
           s.population AS population,
           s.packed_index AS packed_index,
           s.froh AS froh,
           p.n_samples AS pop_n_samples,
           p.mean_froh AS pop_mean_froh
    ORDER BY s.packed_index
    """
    records = ctx.run(cypher)
    if not records:
        click.echo(f"No samples found for population '{population}'.", err=True)
        return

    click.echo(f"Found {len(records)} samples for {population}.", err=True)
    format_output(records, output_path, fmt, "extract samples",
                  {"population": population})


@extract.command("genotypes")
@click.option("--chr", "chromosome", required=True, help="Chromosome")
@click.option("--start", type=int, required=True, help="Region start position")
@click.option("--end", type=int, required=True, help="Region end position")
@click.option("--pop", "population", required=True, help="Population name")
@click.option("--format-gt", "gt_format", default="dosage",
              type=click.Choice(["dosage", "gt", "raw"]),
              help="Output format: dosage (0/1/2), gt (0/0, 0/1, 1/1), or raw (hex gt_packed)")
@click.option("--limit", type=int, default=1000, help="Maximum variants (default: 1000)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def extract_genotypes(ctx, chromosome, start, end, population, gt_format,
                      limit, output_path, fmt):
    """Extract genotype data for a region and population.

    Queries CARRIES edges between Sample and Variant nodes to build a
    sample x variant matrix. For large regions, consider using --limit
    to cap the number of variants.

    Note: gt_packed decoding requires bit operations. With --format-gt raw,
    the raw gt_packed byte array is returned as a hex string per variant.
    With dosage or gt mode, individual CARRIES edges are queried instead.

    \b
    Examples:
      graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 --pop EUR -o geno.tsv
      graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 --pop EUR --format-gt raw -o geno_raw.tsv
    """
    if gt_format == "raw":
        # Return per-variant summary with raw gt_packed as hex
        cypher = (
            f"MATCH (v:Variant) "
            f"WHERE v.chr = '{chromosome}' AND v.pos >= {start} AND v.pos <= {end} "
            f"RETURN v.variantId AS variantId, v.pos AS pos, v.ref AS ref, v.alt AS alt, "
            f"v.gt_packed AS gt_packed_hex, "
            f"[i IN range(0, size(v.pop_ids)-1) "
            f"WHERE v.pop_ids[i] = '{population}' | v.af[i]][0] AS af "
            f"ORDER BY v.pos LIMIT {limit}"
        )
        records = ctx.run(cypher)
        if not records:
            click.echo("No variants found in region.", err=True)
            return
        click.echo(f"Found {len(records)} variants (raw gt_packed mode).", err=True)
        format_output(records, output_path, fmt, "extract genotypes",
                      {"chr": chromosome, "start": start, "end": end,
                       "pop": population, "format": gt_format})
    else:
        # Query CARRIES edges for individual genotypes
        gt_label = "c.gt" if gt_format == "dosage" else (
            "CASE c.gt WHEN 1 THEN '0/1' WHEN 2 THEN '1/1' ELSE '0/0' END"
        )
        cypher = (
            f"MATCH (s:Sample)-[c:CARRIES]->(v:Variant) "
            f"WHERE s.population = '{population}' "
            f"AND v.chr = '{chromosome}' AND v.pos >= {start} AND v.pos <= {end} "
            f"RETURN s.sampleId AS sampleId, v.variantId AS variantId, "
            f"v.pos AS pos, {gt_label} AS genotype "
            f"ORDER BY v.pos, s.sampleId "
            f"LIMIT {limit * 100}"
        )
        records = ctx.run(cypher)
        if not records:
            click.echo("No genotype data found. CARRIES edges may not exist "
                       "for this region/population.", err=True)
            return
        click.echo(f"Found {len(records)} genotype entries.", err=True)
        format_output(records, output_path, fmt, "extract genotypes",
                      {"chr": chromosome, "start": start, "end": end,
                       "pop": population, "format": gt_format})
