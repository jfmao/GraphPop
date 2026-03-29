"""graphpop neighbors -- explore the graph neighborhood around a gene."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.command()
@click.argument("gene")
@click.option("--hops", default=1, type=click.IntRange(1, 3),
              help="Number of hops to traverse (default: 1, max: 3)")
@click.option("--via", default="IN_PATHWAY",
              type=click.Choice(["IN_PATHWAY", "LD", "HAS_GO_TERM"],
                                case_sensitive=False),
              help="Relationship type to traverse (default: IN_PATHWAY)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def neighbors(ctx, gene, hops, via, output_path, fmt):
    """Explore the graph neighborhood around a gene.

    Traverses shared pathways, LD edges, or GO terms to find related genes.

    \b
    Examples:
      graphpop neighbors KCNE1 -o neighbors.tsv
      graphpop neighbors KCNE1 --hops 2 -o neighbors_2hop.tsv
      graphpop neighbors GW5 --via HAS_GO_TERM --format json
    """
    via = via.upper()

    if via == "IN_PATHWAY":
        cypher = _pathway_query(gene, hops)
    elif via == "LD":
        cypher = _ld_query(gene, hops)
    elif via == "HAS_GO_TERM":
        cypher = _go_query(gene, hops)
    else:
        click.echo(f"Unsupported --via type: {via}", err=True)
        raise SystemExit(1)

    records = ctx.run(cypher)

    if not records:
        click.echo(f"No neighbors found for gene '{gene}' via {via} "
                    f"({hops} hop(s)).", err=True)
        return

    click.echo(f"Found {len(records)} neighbor(s) for {gene} via {via} "
                f"({hops} hop(s)).", err=True)
    format_output(records, output_path, fmt, "neighbors",
                  {"gene": gene, "hops": hops, "via": via})


def _pathway_query(gene: str, hops: int) -> str:
    """Build pathway-based neighbor query."""
    if hops == 1:
        return (
            f"MATCH (g1:Gene)-[:IN_PATHWAY]->(p:Pathway)<-[:IN_PATHWAY]-(g2:Gene) "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') AND g1 <> g2 "
            f"RETURN DISTINCT g2.symbol AS gene, p.name AS shared_pathway, "
            f"g2.chr AS chr, g2.start AS start, g2.end AS end "
            f"ORDER BY gene"
        )
    elif hops == 2:
        return (
            f"MATCH (g1:Gene)-[:IN_PATHWAY]->(p1:Pathway)<-[:IN_PATHWAY]-(g2:Gene)"
            f"-[:IN_PATHWAY]->(p2:Pathway)<-[:IN_PATHWAY]-(g3:Gene) "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') "
            f"AND g1 <> g2 AND g1 <> g3 AND g2 <> g3 "
            f"RETURN DISTINCT g3.symbol AS gene, "
            f"g2.symbol AS via_gene, p1.name AS pathway_1, p2.name AS pathway_2, "
            f"g3.chr AS chr, g3.start AS start, g3.end AS end "
            f"ORDER BY gene"
        )
    else:  # hops == 3
        return (
            f"MATCH path = (g1:Gene)"
            f"(-[:IN_PATHWAY]->(:Pathway)<-[:IN_PATHWAY]-(:Gene)){{3}} "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') "
            f"WITH g1, last(nodes(path)) AS g_end, "
            f"[n IN nodes(path) WHERE 'Pathway' IN labels(n) | n.name] AS pws, "
            f"[n IN nodes(path) WHERE 'Gene' IN labels(n) | n.symbol] AS genes "
            f"WHERE g1 <> g_end "
            f"RETURN DISTINCT g_end.symbol AS gene, "
            f"g_end.chr AS chr, g_end.start AS start, g_end.end AS end, "
            f"pws AS pathways, genes AS via_genes "
            f"ORDER BY gene "
            f"LIMIT 500"
        )


def _ld_query(gene: str, hops: int) -> str:
    """Build LD-based neighbor query."""
    if hops == 1:
        return (
            f"MATCH (g1:Gene)<-[:HAS_CONSEQUENCE]-(v1:Variant)"
            f"-[ld:LD]-(v2:Variant)-[:HAS_CONSEQUENCE]->(g2:Gene) "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') AND g1 <> g2 "
            f"RETURN DISTINCT g2.symbol AS gene, "
            f"max(ld.r2) AS max_r2, g2.chr AS chr "
            f"ORDER BY max_r2 DESC"
        )
    elif hops == 2:
        return (
            f"MATCH (g1:Gene)<-[:HAS_CONSEQUENCE]-(v1:Variant)"
            f"-[:LD]-(v2:Variant)-[:LD]-(v3:Variant)"
            f"-[:HAS_CONSEQUENCE]->(g2:Gene) "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') AND g1 <> g2 "
            f"RETURN DISTINCT g2.symbol AS gene, g2.chr AS chr "
            f"ORDER BY gene "
            f"LIMIT 500"
        )
    else:  # hops == 3
        return (
            f"MATCH (g1:Gene)<-[:HAS_CONSEQUENCE]-(v1:Variant)"
            f"-[:LD]-(v2:Variant)-[:LD]-(v3:Variant)-[:LD]-(v4:Variant)"
            f"-[:HAS_CONSEQUENCE]->(g2:Gene) "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') AND g1 <> g2 "
            f"RETURN DISTINCT g2.symbol AS gene, g2.chr AS chr "
            f"ORDER BY gene "
            f"LIMIT 500"
        )


def _go_query(gene: str, hops: int) -> str:
    """Build GO term-based neighbor query."""
    if hops == 1:
        return (
            f"MATCH (g1:Gene)-[:HAS_GO_TERM]->(go:GOTerm)<-[:HAS_GO_TERM]-(g2:Gene) "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') AND g1 <> g2 "
            f"RETURN DISTINCT g2.symbol AS gene, go.name AS shared_go_term, "
            f"go.goId AS go_id, g2.chr AS chr "
            f"ORDER BY gene"
        )
    elif hops == 2:
        return (
            f"MATCH (g1:Gene)-[:HAS_GO_TERM]->(:GOTerm)<-[:HAS_GO_TERM]-(g2:Gene)"
            f"-[:HAS_GO_TERM]->(go2:GOTerm)<-[:HAS_GO_TERM]-(g3:Gene) "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') "
            f"AND g1 <> g2 AND g1 <> g3 AND g2 <> g3 "
            f"RETURN DISTINCT g3.symbol AS gene, "
            f"g2.symbol AS via_gene, go2.name AS go_term, "
            f"g3.chr AS chr "
            f"ORDER BY gene "
            f"LIMIT 500"
        )
    else:  # hops == 3
        return (
            f"MATCH path = (g1:Gene)"
            f"(-[:HAS_GO_TERM]->(:GOTerm)<-[:HAS_GO_TERM]-(:Gene)){{3}} "
            f"WHERE (g1.symbol = '{gene}' OR g1.geneId = '{gene}') "
            f"WITH g1, last(nodes(path)) AS g_end, "
            f"[n IN nodes(path) WHERE 'GOTerm' IN labels(n) | n.name] AS terms, "
            f"[n IN nodes(path) WHERE 'Gene' IN labels(n) | n.symbol] AS genes "
            f"WHERE g1 <> g_end "
            f"RETURN DISTINCT g_end.symbol AS gene, "
            f"g_end.chr AS chr, "
            f"terms AS go_terms, genes AS via_genes "
            f"ORDER BY gene "
            f"LIMIT 500"
        )
