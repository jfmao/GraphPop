"""graphpop ld — linkage disequilibrium (r2, D')."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command()
@click.argument("chr")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.argument("population")
@click.argument("max_dist", type=int)
@click.argument("threshold", type=float)
@click.option("--persist", is_flag=True, default=False,
              help="Write LD edges to the graph")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--min-af", type=float, help="Minimum allele frequency")
@pass_ctx
def ld(ctx, chr, start, end, population, max_dist, threshold, persist,
       output_path, fmt, min_af):
    """Compute pairwise linkage disequilibrium (r2 and D')."""
    opts = build_options_map(min_af=min_af, write_edges=persist)
    cypher = build_cypher(
        "graphpop.ld",
        [f"'{chr}'", str(start), str(end), f"'{population}'",
         str(max_dist), str(threshold)],
        options=opts if opts else None,
        yield_cols=["variant1", "variant2", "r2", "dprime", "distance"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "ld",
                  {"chr": chr, "start": start, "end": end, "pop": population,
                   "max_dist": max_dist, "threshold": threshold, "persist": persist})
