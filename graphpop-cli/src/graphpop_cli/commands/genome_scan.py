"""graphpop genome-scan — sliding-window genome scan."""
import click
from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


@click.command("genome-scan")
@click.argument("chr")
@click.argument("population")
@click.argument("window_size", type=int)
@click.argument("step_size", type=int)
@click.option("--pop2", help="Second population for Fst/Dxy/PBS")
@click.option("--persist", is_flag=True, default=False,
              help="Persist window results to graph (default behavior)")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@click.option("--consequence", help="Filter by VEP consequence type")
@click.option("--pathway", help="Filter by pathway name")
@click.option("--gene", help="Filter by gene name")
@click.option("--min-af", type=float, help="Minimum allele frequency")
@pass_ctx
def genome_scan(ctx, chr, population, window_size, step_size, pop2, persist,
                output_path, fmt, consequence, pathway, gene, min_af):
    """Run a sliding-window genome scan (pi, theta, Tajima's D, Fst, etc.)."""
    opts = build_options_map(consequence=consequence, pathway=pathway, gene=gene,
                             min_af=min_af)
    positional = [f"'{chr}'", f"'{population}'", str(window_size), str(step_size)]
    if pop2:
        positional.append(f"'{pop2}'")
    cypher = build_cypher(
        "graphpop.genome_scan", positional,
        options=opts if opts else None,
        yield_cols=["window_id", "chr", "start", "end", "population",
                     "n_variants", "n_segregating", "pi", "theta_w", "tajima_d",
                     "fst", "fst_wc", "dxy", "pbs", "fay_wu_h"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "genome-scan",
                  {"chr": chr, "pop": population, "window": window_size,
                   "step": step_size, "pop2": pop2})
