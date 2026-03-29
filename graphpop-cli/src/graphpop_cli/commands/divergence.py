"""graphpop divergence — Fst, Dxy, Da, PBS."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command()
@click.argument("chr")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.argument("pop1")
@click.argument("pop2")
@click.option("--pop3", help="Third population for PBS computation")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--consequence", help="Filter by VEP consequence type")
@click.option("--pathway", help="Filter by pathway name")
@click.option("--gene", help="Filter by gene name")
@click.option("--min-af", type=float, help="Minimum allele frequency")
@pass_ctx
def divergence(ctx, chr, start, end, pop1, pop2, pop3, output_path, fmt,
               consequence, pathway, gene, min_af):
    """Compute Hudson Fst, W&C Fst, Dxy, Da, and optionally PBS."""
    opts = build_options_map(consequence=consequence, pathway=pathway, gene=gene,
                             min_af=min_af)
    positional = [f"'{chr}'", str(start), str(end), f"'{pop1}'", f"'{pop2}'"]
    if pop3:
        positional.append(f"'{pop3}'")
    cypher = build_cypher(
        "graphpop.divergence", positional,
        options=opts if opts else None,
        yield_cols=["fst_hudson", "fst_wc", "dxy", "da", "pbs", "n_variants"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "divergence",
                  {"chr": chr, "pop1": pop1, "pop2": pop2, "pop3": pop3})
