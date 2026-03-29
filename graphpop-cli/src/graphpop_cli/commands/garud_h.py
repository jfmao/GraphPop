"""graphpop garud-h — Garud's H statistics for haplotype homozygosity."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command("garud-h")
@click.argument("chr")
@click.argument("population")
@click.argument("window_size", type=int)
@click.argument("step_size", type=int)
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--min-af", type=float, help="Minimum allele frequency")
@pass_ctx
def garud_h(ctx, chr, population, window_size, step_size, output_path, fmt, min_af):
    """Compute Garud's H1, H12, H2/H1 in sliding windows."""
    opts = build_options_map(min_af=min_af)
    cypher = build_cypher(
        "graphpop.garud_h",
        [f"'{chr}'", f"'{population}'", str(window_size), str(step_size)],
        options=opts if opts else None,
        yield_cols=["chr", "start", "end", "population", "h1", "h12", "h2_h1",
                     "hap_diversity", "n_haplotypes", "n_variants"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "garud-h",
                  {"chr": chr, "pop": population, "window": window_size,
                   "step": step_size})
