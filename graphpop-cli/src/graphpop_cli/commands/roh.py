"""graphpop roh — runs of homozygosity."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command()
@click.argument("chr")
@click.argument("population")
@click.option("--method", type=click.Choice(["hmm", "window"]), default="hmm",
              help="ROH detection method (default: hmm)")
@click.option("--min-length", type=int, help="Minimum ROH length in bp")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def roh(ctx, chr, population, method, min_length, output_path, fmt):
    """Detect runs of homozygosity (ROH) per sample."""
    opts = build_options_map(method=method, min_length=min_length)
    cypher = build_cypher(
        "graphpop.roh",
        [f"'{chr}'", f"'{population}'"],
        options=opts if opts else None,
        yield_cols=["sampleId", "n_roh", "total_length", "froh",
                     "mean_length", "max_length"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "roh",
                  {"chr": chr, "pop": population, "method": method,
                   "min_length": min_length})
