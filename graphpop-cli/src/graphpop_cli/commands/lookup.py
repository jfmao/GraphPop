"""graphpop lookup — query genes, pathways, variants, and regions in the graph."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.group()
def lookup():
    """Look up genes, pathways, variants, or genomic regions.

    \b
    Subcommands:
      gene      Look up a gene by symbol or ID
      pathway   Look up a pathway by name
      variant   Look up a variant by ID
      region    Look up genes and stats in a genomic region

    \b
    Examples:
      graphpop lookup gene KCNE1
      graphpop lookup pathway "Cardiac repolarization"
      graphpop lookup variant chr22:16050075:A:G
      graphpop lookup region chr6 9000000 9600000
    """
    pass


@lookup.command("gene")
@click.argument("gene_name")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def lookup_gene(ctx, gene_name, output_path, fmt):
    """Look up a gene: variant count, consequences, pathways, and selection stats.

    GENE_NAME can be a gene symbol (e.g., KCNE1) or gene ID (e.g., ENSG00000180509).

    \b
    Examples:
      graphpop lookup gene KCNE1
      graphpop lookup gene GW5 -o gw5_info.tsv
      graphpop lookup gene ENSG00000180509 --format json
    """
    cypher = """
    MATCH (g:Gene)
    WHERE g.symbol = $gene_name OR g.geneId = $gene_name
    OPTIONAL MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g)
    OPTIONAL MATCH (g)-[:IN_PATHWAY]->(pw:Pathway)
    WITH g, v, COLLECT(DISTINCT pw.name) AS pathways
    RETURN g.symbol AS gene,
           g.geneId AS gene_id,
           g.chr AS chr,
           g.start AS start,
           g.end AS end,
           v.variantId AS variant_id,
           v.pos AS pos,
           v.ref AS ref,
           v.alt AS alt,
           pathways,
           CASE WHEN v IS NOT NULL THEN [k IN keys(v) WHERE k STARTS WITH 'ihs_' | k + '=' + toString(v[k])] ELSE [] END AS ihs_scores,
           CASE WHEN v IS NOT NULL THEN [k IN keys(v) WHERE k STARTS WITH 'xpehh_' | k + '=' + toString(v[k])] ELSE [] END AS xpehh_scores
    ORDER BY v.pos
    """
    records = ctx.run(cypher, {"gene_name": gene_name})

    if not records:
        click.echo(f"Gene '{gene_name}' not found in the graph.", err=True)
        return

    click.echo(f"Found {len(records)} variants for gene {gene_name}.", err=True)
    format_output(records, output_path, fmt, "lookup gene",
                  {"gene": gene_name})


@lookup.command("pathway")
@click.argument("pw_name")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def lookup_pathway(ctx, pw_name, output_path, fmt):
    """Look up a pathway: member genes and variant counts.

    PW_NAME is matched as a substring (CONTAINS) against pathway names.

    \b
    Examples:
      graphpop lookup pathway "Cardiac repolarization"
      graphpop lookup pathway "starch" -o starch_pathway.tsv
    """
    cypher = """
    MATCH (pw:Pathway)
    WHERE pw.name CONTAINS $pw_name
    OPTIONAL MATCH (g:Gene)-[:IN_PATHWAY]->(pw)
    OPTIONAL MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g)
    WITH pw, g, COUNT(DISTINCT v) AS variant_count
    RETURN pw.name AS pathway,
           pw.pathwayId AS pathway_id,
           g.symbol AS gene,
           g.geneId AS gene_id,
           g.chr AS chr,
           g.start AS gene_start,
           g.end AS gene_end,
           variant_count
    ORDER BY pw.name, g.symbol
    """
    records = ctx.run(cypher, {"pw_name": pw_name})

    if not records:
        click.echo(f"No pathways matching '{pw_name}' found.", err=True)
        return

    click.echo(f"Found {len(records)} gene entries across matching pathways.", err=True)
    format_output(records, output_path, fmt, "lookup pathway",
                  {"pathway": pw_name})


@lookup.command("variant")
@click.argument("var_id")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def lookup_variant(ctx, var_id, output_path, fmt):
    """Full annotation for a single variant.

    VAR_ID format: chr:pos:ref:alt (e.g., chr22:16050075:A:G).

    \b
    Examples:
      graphpop lookup variant chr22:16050075:A:G
      graphpop lookup variant Chr01:12345:A:T --format json
    """
    cypher = """
    MATCH (v:Variant {variantId: $var_id})
    OPTIONAL MATCH (v)-[:HAS_CONSEQUENCE]->(g:Gene)
    OPTIONAL MATCH (g)-[:IN_PATHWAY]->(pw:Pathway)
    RETURN v AS variant_props,
           g.symbol AS gene,
           g.geneId AS gene_id,
           COLLECT(DISTINCT pw.name) AS pathways
    """
    records = ctx.run(cypher, {"var_id": var_id})

    if not records:
        click.echo(f"Variant '{var_id}' not found.", err=True)
        return

    # Flatten variant properties into columns
    flat_records = []
    for rec in records:
        row = {}
        vprops = rec.get("variant_props", {})
        if vprops:
            for k, v in vprops.items():
                row[k] = v
        row["gene"] = rec.get("gene")
        row["gene_id"] = rec.get("gene_id")
        row["pathways"] = rec.get("pathways", [])
        flat_records.append(row)

    format_output(flat_records, output_path, fmt, "lookup variant",
                  {"variant_id": var_id})


@lookup.command("region")
@click.argument("chr")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def lookup_region(ctx, chr, start, end, output_path, fmt):
    """Genes and summary stats in a genomic region.

    Returns per-gene variant counts and mean allele frequencies in the region.

    \b
    Examples:
      graphpop lookup region chr6 9000000 9600000
      graphpop lookup region chr22 16000000 17000000 -o region.tsv
    """
    cypher = """
    MATCH (v:Variant)
    WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end
    OPTIONAL MATCH (v)-[:HAS_CONSEQUENCE]->(g:Gene)
    WITH g, COUNT(DISTINCT v) AS variant_count,
         MIN(v.pos) AS min_pos, MAX(v.pos) AS max_pos
    RETURN COALESCE(g.symbol, 'intergenic') AS gene,
           g.geneId AS gene_id,
           g.start AS gene_start,
           g.end AS gene_end,
           variant_count,
           min_pos,
           max_pos
    ORDER BY min_pos
    """
    records = ctx.run(cypher, {"chr": chr, "start": start, "end": end})

    if not records:
        click.echo(f"No variants found in {chr}:{start}-{end}.", err=True)
        return

    click.echo(f"Found {len(records)} genes/regions in {chr}:{start}-{end}.", err=True)
    format_output(records, output_path, fmt, "lookup region",
                  {"chr": chr, "start": start, "end": end})
