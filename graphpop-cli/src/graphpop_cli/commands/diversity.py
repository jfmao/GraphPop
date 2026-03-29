"""graphpop diversity — nucleotide diversity, theta, Tajima's D, Fay & Wu's H."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command()
@click.argument("chr")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.argument("population")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--consequence", help="Filter by VEP consequence type")
@click.option("--pathway", help="Filter by pathway name")
@click.option("--gene", help="Filter by gene name")
@click.option("--min-af", type=float, help="Minimum allele frequency")
@click.option("--max-af", type=float, help="Maximum allele frequency")
@pass_ctx
def diversity(ctx, chr, start, end, population, output_path, fmt,
              consequence, pathway, gene, min_af, max_af):
    """Compute nucleotide diversity, theta_W, Tajima's D, Fay & Wu's H."""
    opts = build_options_map(consequence=consequence, pathway=pathway, gene=gene,
                             min_af=min_af, max_af=max_af)
    cypher = build_cypher(
        "graphpop.diversity",
        [f"'{chr}'", str(start), str(end), f"'{population}'"],
        options=opts if opts else None,
        yield_cols=["pi", "theta_w", "tajima_d", "fay_wu_h", "fay_wu_h_norm",
                     "het_exp", "het_obs", "fis", "n_variants", "n_segregating",
                     "n_polarized"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "diversity",
                  {"chr": chr, "start": start, "end": end, "pop": population})
