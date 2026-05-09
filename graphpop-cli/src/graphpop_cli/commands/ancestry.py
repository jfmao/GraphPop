"""graphpop ancestry — local-ancestry painting ingest and management (M4.B)."""
from __future__ import annotations

from pathlib import Path

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.group()
def ancestry():
    """Local-ancestry painting (:HAS_ANCESTRY) management.

    Subcommands:

      ingest                Ingest per-TreeNode painting from a TSV.
      ingest-from-samples   Propagate sample-level ancestry to all
                            TreeNodes via majority-vote DFS.
      list                  List every painting currently stored.
      delete                Remove a painting by (run_id, painter).
    """


@ancestry.command("ingest")
@click.argument("painting_tsv", type=click.Path(exists=True, dir_okay=False))
@click.option("--run-id", required=True,
              help=":ARGRun.runId to attach the painting to")
@click.option("--painter", default="manual", show_default=True,
              help="Identifier of the painting source (e.g. 'manual', 'rfmix')")
@click.option("--replace", is_flag=True,
              help="Delete any prior painting with the same (run_id, painter) "
                   "before inserting")
@pass_ctx
def ingest(ctx, painting_tsv, run_id, painter, replace):
    """Ingest a per-TreeNode painting TSV.

    TSV columns (no header required):

      tskit_node_id<TAB>population_id[<TAB>posterior_prob]

    posterior_prob defaults to 1.0 when omitted. Lines starting with
    '#' and blank lines are skipped.
    """
    from graphpop_import.ancestry_ingester import (
        AncestryIngester, PaintingRow,
    )

    rows: list[PaintingRow] = []
    with open(painting_tsv) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise click.ClickException(
                    f"line {lineno}: expected at least 2 TSV columns, got: {raw!r}"
                )
            try:
                node = int(parts[0])
            except ValueError as exc:
                raise click.ClickException(
                    f"line {lineno}: first column must be an int; got {parts[0]!r}"
                ) from exc
            prob = float(parts[2]) if len(parts) >= 3 else 1.0
            rows.append(PaintingRow(node, parts[1], prob))

    ing = AncestryIngester(ctx.driver, database=ctx.database)
    summary = ing.ingest(run_id, rows, painter=painter, replace=replace)
    click.echo(
        f"Ingested {summary.n_edges} :HAS_ANCESTRY edges over "
        f"{summary.n_painted_nodes} TreeNodes ({summary.n_populations} "
        f"populations) for run_id={run_id}, painter={painter}.",
        err=True,
    )


@ancestry.command("ingest-from-samples")
@click.argument("sample_tsv", type=click.Path(exists=True, dir_okay=False))
@click.option("--run-id", required=True,
              help=":ARGRun.runId whose TreeNodes will be painted")
@click.option("--painter", default="majority_vote", show_default=True,
              help="Identifier of the painting source")
@click.option("--replace", is_flag=True,
              help="Delete any prior painting with the same (run_id, painter)")
@pass_ctx
def ingest_from_samples(ctx, sample_tsv, run_id, painter, replace):
    """Propagate sample-level ancestry to every TreeNode in the run.

    Sample TSV columns:

      sample_id<TAB>population_id

    Internal TreeNodes get the argmax-of-descendants label, with
    posterior_prob = majority fraction. Ties tip to the alphabetically
    smallest label.
    """
    from graphpop_import.ancestry_ingester import AncestryIngester

    sample_anc: dict[str, str] = {}
    with open(sample_tsv) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise click.ClickException(
                    f"line {lineno}: expected at least 2 TSV columns; got {raw!r}"
                )
            sample_anc[parts[0]] = parts[1]

    ing = AncestryIngester(ctx.driver, database=ctx.database)
    summary = ing.ingest_from_samples(
        run_id, sample_anc, painter=painter, replace=replace)
    click.echo(
        f"Propagated {len(sample_anc)} sample labels -> "
        f"{summary.n_edges} :HAS_ANCESTRY edges over "
        f"{summary.n_painted_nodes} TreeNodes ({summary.n_populations} "
        f"populations) for run_id={run_id}, painter={painter}.",
        err=True,
    )


@ancestry.command("list")
@click.option("--run-id", default=None,
              help="Filter to a single :ARGRun (default: all runs)")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def list_(ctx, run_id, output_path, fmt):
    """List every (run_id, painter) painting currently in the database."""
    from graphpop_import.ancestry_ingester import AncestryIngester

    ing = AncestryIngester(ctx.driver, database=ctx.database)
    paintings = ing.list_paintings(run_id=run_id)
    records = [
        {
            "run_id": p.run_id,
            "painter": p.painter,
            "n_painted_nodes": p.n_painted_nodes,
            "n_edges": p.n_edges,
            "n_populations": p.n_populations,
        }
        for p in paintings
    ]
    format_output(records, output_path, fmt, "ancestry-list", {})


@ancestry.command("delete")
@click.option("--run-id", required=True)
@click.option("--painter", required=True)
@click.option("--yes", is_flag=True,
              help="Skip the confirmation prompt")
@pass_ctx
def delete(ctx, run_id, painter, yes):
    """Delete every :HAS_ANCESTRY edge for a (run_id, painter) pair."""
    if not yes:
        click.confirm(
            f"Delete painting run_id='{run_id}', painter='{painter}'?",
            abort=True,
        )
    from graphpop_import.ancestry_ingester import AncestryIngester

    ing = AncestryIngester(ctx.driver, database=ctx.database)
    n = ing.delete_painting(run_id, painter)
    click.echo(f"Deleted {n} :HAS_ANCESTRY edges.", err=True)
