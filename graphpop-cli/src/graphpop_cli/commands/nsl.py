"""graphpop nsl — number of segregating sites by length."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command()
@click.argument("chr")
@click.argument("population")
@click.option("--min-af", type=float, help="Minimum allele frequency filter")
@click.option("--persist", is_flag=True, default=False,
              help="Write nSL scores to Variant nodes")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def nsl(ctx, chr, population, min_af, persist, output_path, fmt):
    """Compute number of segregating sites by length (nSL)."""
    opts = build_options_map(min_af=min_af, persist=persist)
    cypher = build_cypher(
        "graphpop.nsl",
        [f"'{chr}'", f"'{population}'"],
        options=opts if opts else None,
        yield_cols=["variantId", "pos", "af", "nsl_unstd", "nsl"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "nsl",
                  {"chr": chr, "pop": population, "min_af": min_af,
                   "persist": persist})
