"""graphpop arg — ingest, list, and delete ARG runs (M4.A)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.group()
def arg():
    """Ancestral recombination graph ingest and management.

    Subcommands:

      ingest  Ingest a tskit TreeSequence into the GraphPop ARG layer.
      list    List every :ARGRun currently stored.
      delete  Remove a single :ARGRun and all its derived data.

    The ARG layer (M4.A) is additive: it expects the target database to
    already contain :Variant, :Sample, :Population nodes from a prior
    bulk import.
    """


@arg.command("ingest")
@click.argument("tree_sequence", type=click.Path(exists=True, dir_okay=False))
@click.option("--run-id", required=True,
              help="Globally unique identifier for this run")
@click.option("--source",
              type=click.Choice(
                  ["tsinfer", "tsdate", "singer", "relate", "msprime"],
                  case_sensitive=False),
              default="tsinfer",
              show_default=True,
              help="Inferer that produced the ARG")
@click.option("--sample-id-map", type=click.Path(exists=True, dir_okay=False),
              default=None,
              help="TSV mapping tskit sample index -> :Sample.sampleId; "
                   "first column is the integer index, second is the sample ID. "
                   "Defaults to identity mapping (sample_{i}).")
@click.option("--params-json", default=None,
              help="JSON dict of free-form parameters to record on :ARGRun")
@click.option("--batch-size", type=int, default=10_000, show_default=True,
              help="Rows per UNWIND batch")
@pass_ctx
def ingest(ctx, tree_sequence, run_id, source, sample_id_map, params_json,
           batch_size):
    """Ingest a tskit TreeSequence file into Neo4j as :ARGRun + :TreeNode + ARG edges.

    The ingest is one-shot and additive: it never touches existing
    :Variant or :Sample nodes. To replace a prior run with the same
    run-id, run ``graphpop arg delete <run-id>`` first.
    """
    try:
        import tskit  # type: ignore
    except ImportError as exc:
        raise click.ClickException(
            "tskit is not installed. Install with "
            "`pip install graphpop-cli[arg]` or `pip install tskit`."
        ) from exc

    from graphpop_import.arg_importer import ARGIngester

    params: dict | None = None
    if params_json:
        try:
            params = json.loads(params_json)
        except json.JSONDecodeError as exc:
            raise click.ClickException(
                f"--params-json is not valid JSON: {exc}"
            ) from exc

    sample_id_callable = None
    if sample_id_map:
        mapping: dict[int, str] = {}
        with open(sample_id_map) as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    raise click.ClickException(
                        f"--sample-id-map line {lineno} must have two TSV "
                        f"columns; got: {raw!r}"
                    )
                try:
                    idx = int(parts[0])
                except ValueError as exc:
                    raise click.ClickException(
                        f"--sample-id-map line {lineno} first column must "
                        f"be an integer; got: {parts[0]!r}"
                    ) from exc
                mapping[idx] = parts[1]

        def sample_id_callable(i: int) -> str:  # type: ignore[no-redef]
            try:
                return mapping[i]
            except KeyError as exc:
                raise click.ClickException(
                    f"--sample-id-map has no entry for sample index {i}"
                ) from exc

    treeseq = tskit.load(tree_sequence)
    click.echo(
        f"Loaded {tree_sequence}: {treeseq.num_samples} samples, "
        f"{treeseq.num_trees} trees, {treeseq.num_edges} edges, "
        f"{treeseq.num_mutations} mutations.",
        err=True,
    )

    ingester = ARGIngester(
        ctx.driver, database=ctx.database, batch_size=batch_size
    )
    summary = ingester.ingest(
        treeseq,
        run_id=run_id,
        source=source.lower(),
        params=params,
        sample_id_map=sample_id_callable,
    )
    click.echo(
        f"Ingested ARGRun {summary.run_id}: {summary.n_nodes} TreeNodes, "
        f"{summary.n_edges} PARENT_OF edges, {summary.n_mutations} mutations.",
        err=True,
    )


@arg.command("list")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def list_runs(ctx, output_path, fmt):
    """List every :ARGRun in the database."""
    from graphpop_import.arg_importer import ARGIngester

    ingester = ARGIngester(ctx.driver, database=ctx.database)
    runs = ingester.list_runs()
    records = [
        {
            "run_id": r.run_id,
            "source": r.source,
            "n_samples": r.n_samples,
            "sequence_length": r.sequence_length,
            "n_trees": r.n_trees,
            "n_nodes": r.n_nodes,
            "n_edges": r.n_edges,
            "n_mutations": r.n_mutations,
            "created_at": r.created_at.isoformat() if r.created_at else None,
        }
        for r in runs
    ]
    format_output(records, output_path, fmt, "arg-list", {})


@arg.command("delete")
@click.argument("run_id")
@click.option("--yes", is_flag=True,
              help="Skip the confirmation prompt")
@pass_ctx
def delete(ctx, run_id, yes):
    """Delete an :ARGRun and every node/relationship that references it."""
    if not yes:
        click.confirm(
            f"Delete ARGRun '{run_id}' and all derived TreeNodes/edges?",
            abort=True,
        )
    from graphpop_import.arg_importer import ARGIngester

    ingester = ARGIngester(ctx.driver, database=ctx.database)
    n = ingester.delete_run(run_id)
    click.echo(f"Deleted run {run_id}: {n} TreeNodes removed.", err=True)
