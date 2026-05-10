"""graphpop arg — ingest, list, delete ARG runs (M4.A) +
ARG-derived statistics (M6: tmrca, branch-diversity, coalescence-rate,
allele-age)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def arg():
    """Ancestral recombination graph ingest and management.

    Subcommands:

      ingest             Ingest a tskit TreeSequence.
      list               List every :ARGRun currently stored.
      delete             Remove a single :ARGRun and all its derived data.
      tmrca              Pairwise TMRCA at a position (or window-mean).
      branch-diversity   Branch-mode pi (tskit-equivalent).
      coalescence-rate   Per-time-bin coalescence rate.
      allele-age         Time bracket for a variant's :MUTATED_ON edge.

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


@arg.command("tmrca")
@click.argument("run_id")
@click.argument("sample_a")
@click.argument("sample_b")
@click.option("--position", type=int, default=None,
              help="Single-position TMRCA at this bp (default: span-weighted "
                   "mean over the full sequence)")
@click.option("--window-start", type=int, default=None,
              help="Window start (bp) for window-mean TMRCA")
@click.option("--window-end", type=int, default=None,
              help="Window end (bp) for window-mean TMRCA")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def tmrca(ctx, run_id, sample_a, sample_b, position, window_start, window_end,
          output_path, fmt):
    """TMRCA between two samples in an ARG (M6).

    Three modes: --position for a single bp; --window-start/--window-end
    for a window-mean TMRCA; neither argument set for genome-wide
    span-weighted mean.
    """
    opts: dict[str, object] = {}
    if position is not None:
        opts["position"] = position
    if window_start is not None:
        opts["window_start"] = window_start
    if window_end is not None:
        opts["window_end"] = window_end

    cypher = build_cypher(
        "graphpop.arg.tmrca",
        [f"'{run_id}'", f"'{sample_a}'", f"'{sample_b}'"],
        options=opts if opts else None,
        yield_cols=["sample_a", "sample_b", "position", "tmrca",
                    "mean_tmrca", "mrca_node_id", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "arg-tmrca", {"run_id": run_id})


@arg.command("branch-diversity")
@click.argument("run_id")
@click.argument("sample_ids")
@click.option("--mode", default="pi", show_default=True,
              type=click.Choice(["pi"]),
              help="Diversity mode (only 'pi' supported in v1)")
@click.option("--windows", default=None,
              help="Comma-separated bp boundaries, e.g. 0,25000,50000 → "
                   "two windows. Default: single window over full sequence.")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def branch_diversity(ctx, run_id, sample_ids, mode, windows, output_path, fmt):
    """Branch-mode diversity (M6). SAMPLE_IDS is a comma-separated list."""
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    opts: dict[str, object] = {"mode": mode}
    if windows:
        opts["windows"] = [int(x) for x in windows.split(",")]

    cypher = build_cypher(
        "graphpop.arg.branch_diversity",
        [f"'{run_id}'", _cypher_str_list(sids)],
        options=opts,
        yield_cols=["start", "end", "branch_pi", "n_samples", "mode", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "arg-branch-diversity",
                  {"run_id": run_id, "n_samples": len(sids)})


@arg.command("coalescence-rate")
@click.argument("run_id")
@click.argument("sample_ids")
@click.option("--time-bins", required=True,
              help="Comma-separated monotonically increasing bin edges, "
                   "e.g. 0,100,1000,10000,100000")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def coalescence_rate(ctx, run_id, sample_ids, time_bins, output_path, fmt):
    """Per-time-bin coalescence rate (M6). SAMPLE_IDS is comma-separated."""
    sids = [s.strip() for s in sample_ids.split(",") if s.strip()]
    bins = [float(x) for x in time_bins.split(",")]
    opts: dict[str, object] = {"time_bins": bins}

    cypher = build_cypher(
        "graphpop.arg.coalescence_rate",
        [f"'{run_id}'", _cypher_str_list(sids)],
        options=opts,
        yield_cols=["time_lo", "time_hi", "n_coalescent_events",
                    "lineage_pair_time", "rate", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "arg-coalescence-rate",
                  {"run_id": run_id, "n_samples": len(sids)})


@arg.command("allele-age")
@click.argument("run_id")
@click.argument("variant_id")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def allele_age(ctx, run_id, variant_id, output_path, fmt):
    """Time bracket for a variant's :MUTATED_ON edge (M6)."""
    cypher = build_cypher(
        "graphpop.arg.allele_age",
        [f"'{run_id}'", f"'{variant_id}'"],
        yield_cols=["variant_id", "child_node_id", "parent_node_id",
                    "child_time", "parent_time", "midpoint_time",
                    "n_carriers", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "arg-allele-age",
                  {"run_id": run_id, "variant_id": variant_id})


def _cypher_str_list(items: list[str]) -> str:
    """Render a Python str list as a Cypher list literal."""
    inner = ", ".join(f"'{s}'" for s in items)
    return f"[{inner}]"
