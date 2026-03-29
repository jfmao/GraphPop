"""graphpop xpehh — cross-population extended haplotype homozygosity."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command()
@click.argument("chr")
@click.argument("pop1")
@click.argument("pop2")
@click.option("--min-af", type=float, help="Minimum allele frequency filter")
@click.option("--persist", is_flag=True, default=False,
              help="Write XP-EHH scores to Variant nodes")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def xpehh(ctx, chr, pop1, pop2, min_af, persist, output_path, fmt):
    """Compute cross-population extended haplotype homozygosity (XP-EHH)."""
    opts = build_options_map(min_af=min_af, persist=persist)
    cypher = build_cypher(
        "graphpop.xpehh",
        [f"'{chr}'", f"'{pop1}'", f"'{pop2}'"],
        options=opts if opts else None,
        yield_cols=["variantId", "pos", "af_pop1", "af_pop2",
                     "xpehh_unstd", "xpehh"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "xpehh",
                  {"chr": chr, "pop1": pop1, "pop2": pop2,
                   "min_af": min_af, "persist": persist})
