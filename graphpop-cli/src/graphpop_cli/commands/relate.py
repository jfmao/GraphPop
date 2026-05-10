"""graphpop relate — relationship classification + family clustering (M5)."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def relate():
    """Per-pair relationship classification and family clustering.

    Subcommands:

      classify  Manichaikul 2010 categorical relationships from kinship.
      families  Connected-components family clustering on a relatedness graph.
    """


@relate.command("classify")
@click.argument("source")
@click.option("--min-phi", type=float, default=None,
              help="Minimum phi to emit a row (default 0.0442 = 3rd-degree)")
@click.option("--persist/--no-persist", default=True,
              help="Write :RELATIVE edges (default: persist)")
@click.option("--identical", type=float, default=None,
              help="Override identical / MZ-twin cutoff (default 0.354)")
@click.option("--first-degree", type=float, default=None,
              help="Override first-degree cutoff (default 0.177)")
@click.option("--second-degree", type=float, default=None,
              help="Override second-degree cutoff (default 0.0884)")
@click.option("--third-degree", type=float, default=None,
              help="Override third-degree cutoff (default 0.0442)")
@click.option("--ibs0-threshold", type=float, default=None,
              help="IBS0-fraction cutoff that splits parent-child from "
                   "full-sibling (default 0.0050)")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def classify(ctx, source, min_phi, persist, identical, first_degree,
             second_degree, third_degree, ibs0_threshold, output_path, fmt):
    """Per-pair relationship classification (Manichaikul 2010).

    SOURCE selects the pairwise signal:

      \b
      king          KING-robust kinship + IBS0 from :KINSHIP edges
      <ibd_source>  any IBD source name from :IBD_SEGMENT (e.g. hap_ibd,
                    arg_derived) — phi is total IBD length / (2 × genome)

    Writes idempotent (:Sample)-[:RELATIVE {relationship, degree, phi,
    ibs0_frac, source}]->(:Sample) edges per source.
    """
    opts: dict[str, object] = {"persist": persist}
    if min_phi is not None:
        opts["min_phi"] = min_phi
    if identical is not None:
        opts["identical"] = identical
    if first_degree is not None:
        opts["first_degree"] = first_degree
    if second_degree is not None:
        opts["second_degree"] = second_degree
    if third_degree is not None:
        opts["third_degree"] = third_degree
    if ibs0_threshold is not None:
        opts["ibs0_threshold"] = ibs0_threshold

    cypher = build_cypher(
        "graphpop.relate.classify",
        [f"'{source}'"],
        options=opts,
        yield_cols=["sample_a", "sample_b", "relationship", "degree",
                    "phi", "ibs0_frac", "source"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "relate-classify",
                  {"source": source})


@relate.command("families")
@click.argument("source")
@click.option("--max-degree", type=int, default=None,
              help="Edge predicate: link pairs with :RELATIVE.degree <= N "
                   "(default 2 = extended family)")
@click.option("--min-total-ibd-bp", type=int, default=None,
              help="Edge predicate (mutually exclusive): link pairs with "
                   "summed :IBD_SEGMENT length >= N bp")
@click.option("--min-phi", type=float, default=None,
              help="Edge predicate (mutually exclusive): link pairs with "
                   ":RELATIVE.phi >= F")
@click.option("--persist/--no-persist", default=True,
              help="Write :Family + :IN_FAMILY edges (default: persist)")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def families(ctx, source, max_degree, min_total_ibd_bp, min_phi,
             persist, output_path, fmt):
    """Family clustering via connected components.

    Pick exactly one edge predicate; default is --max-degree=2. Each
    connected component is one extended family; isolated samples form
    singleton families.
    """
    n_predicates = sum(x is not None
                       for x in (max_degree, min_total_ibd_bp, min_phi))
    if n_predicates > 1:
        raise click.UsageError(
            "--max-degree, --min-total-ibd-bp, and --min-phi are mutually "
            "exclusive; pick one.")

    opts: dict[str, object] = {"persist": persist}
    if max_degree is not None:
        opts["max_degree"] = max_degree
    if min_total_ibd_bp is not None:
        opts["min_total_ibd_bp"] = min_total_ibd_bp
    if min_phi is not None:
        opts["min_phi"] = min_phi

    cypher = build_cypher(
        "graphpop.relate.families",
        [f"'{source}'"],
        options=opts,
        yield_cols=["sample_id", "family_id", "family_size", "method"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "relate-families",
                  {"source": source})
