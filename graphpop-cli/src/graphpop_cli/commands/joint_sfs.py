"""graphpop joint-sfs — joint site frequency spectrum between two populations."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command("joint-sfs")
@click.argument("chr")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.argument("pop1")
@click.argument("pop2")
@click.option("--unfolded", is_flag=True, default=False,
              help="Compute unfolded joint SFS (requires ancestral allele)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--consequence", help="Filter by VEP consequence type")
@click.option("--pathway", help="Filter by pathway name")
@click.option("--gene", help="Filter by gene name")
@click.option("--min-af", type=float, help="Minimum allele frequency")
@pass_ctx
def joint_sfs(ctx, chr, start, end, pop1, pop2, unfolded, output_path, fmt,
              consequence, pathway, gene, min_af):
    """Compute the joint site frequency spectrum between two populations."""
    opts = build_options_map(consequence=consequence, pathway=pathway, gene=gene,
                             min_af=min_af)
    cypher = build_cypher(
        "graphpop.joint_sfs",
        [f"'{chr}'", str(start), str(end), f"'{pop1}'", f"'{pop2}'",
         "true" if unfolded else "false"],
        options=opts if opts else None,
        yield_cols=["joint_sfs", "n_variants", "max_ac1", "max_ac2", "dim1", "dim2"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "joint-sfs",
                  {"chr": chr, "start": start, "end": end,
                   "pop1": pop1, "pop2": pop2, "unfolded": unfolded})
