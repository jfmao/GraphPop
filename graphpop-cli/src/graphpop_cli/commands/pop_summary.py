"""graphpop pop-summary — whole-chromosome population summary statistics."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command("pop-summary")
@click.argument("chr")
@click.argument("population")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--consequence", help="Filter by VEP consequence type")
@click.option("--pathway", help="Filter by pathway name")
@click.option("--gene", help="Filter by gene name")
@pass_ctx
def pop_summary(ctx, chr, population, output_path, fmt,
                consequence, pathway, gene):
    """Compute population summary statistics for a chromosome."""
    opts = build_options_map(consequence=consequence, pathway=pathway, gene=gene)
    cypher = build_cypher(
        "graphpop.pop_summary",
        [f"'{chr}'", f"'{population}'"],
        options=opts if opts else None,
        yield_cols=["pi", "theta_w", "tajima_d", "fay_wu_h", "mean_he", "mean_ho",
                     "mean_fis", "n_variants", "n_segregating", "n_polarized"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "pop-summary",
                  {"chr": chr, "pop": population})
