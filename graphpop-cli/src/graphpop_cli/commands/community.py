"""graphpop community — Louvain modularity detection (M9)."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def community():
    """Community detection on the relatedness graph (M9).

    Subcommands:

      louvain  Pure-Java Louvain modularity on :RELATIVE edges.
    """


@community.command("louvain")
@click.argument("source")
@click.option("--edge-weight", default="unit", show_default=True,
              type=click.Choice(["unit", "phi", "degree"]),
              help="Edge weighting: 'unit' (1.0), 'phi' (RELATIVE.phi), "
                   "or 'degree' (1 / (1 + degree); closer relatives weigh more)")
@click.option("--seed", type=int, default=42, show_default=True,
              help="RNG seed for deterministic node-iteration order")
@click.option("--persist/--no-persist", default=True,
              help="Write :IN_COMMUNITY edges + :Community nodes")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def louvain(ctx, source, edge_weight, seed, persist, output_path, fmt):
    """Louvain modularity community detection (M9).

    SOURCE is the :RELATIVE.source from M5 (e.g. 'king', 'hap_ibd').
    """
    opts: dict[str, object] = {
        "edge_weight": edge_weight,
        "seed": seed,
        "persist": persist,
    }
    cypher = build_cypher(
        "graphpop.community.louvain",
        [f"'{source}'"],
        options=opts,
        yield_cols=["sample_id", "community_id", "modularity",
                    "n_communities", "source"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "community-louvain",
                  {"source": source, "edge_weight": edge_weight})
