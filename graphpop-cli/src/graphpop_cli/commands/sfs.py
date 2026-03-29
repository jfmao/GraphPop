"""graphpop sfs — site frequency spectrum."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command()
@click.argument("chr")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.argument("population")
@click.option("--unfolded", is_flag=True, default=False,
              help="Compute unfolded SFS (requires ancestral allele)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--consequence", help="Filter by VEP consequence type")
@click.option("--pathway", help="Filter by pathway name")
@click.option("--gene", help="Filter by gene name")
@click.option("--min-af", type=float, help="Minimum allele frequency")
@click.option("--max-af", type=float, help="Maximum allele frequency")
@pass_ctx
def sfs(ctx, chr, start, end, population, unfolded, output_path, fmt,
         consequence, pathway, gene, min_af, max_af):
    """Compute the site frequency spectrum."""
    opts = build_options_map(consequence=consequence, pathway=pathway, gene=gene,
                             min_af=min_af, max_af=max_af)
    cypher = build_cypher(
        "graphpop.sfs",
        [f"'{chr}'", str(start), str(end), f"'{population}'",
         "true" if unfolded else "false"],
        options=opts if opts else None,
        yield_cols=["sfs", "n_variants", "max_ac", "n_polarized"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "sfs",
                  {"chr": chr, "start": start, "end": end, "pop": population,
                   "unfolded": unfolded})
