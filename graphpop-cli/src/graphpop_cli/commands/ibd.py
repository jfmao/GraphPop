"""graphpop ibd — IBD-segment ingest, ARG-derived caller, kinship aggregator (M4.3)."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def ibd():
    """Identity-by-descent segments and IBD-derived kinship.

    Subcommands:

      ingest    Ingest external IBD-caller output (hap-IBD-style TSV).
      from-arg  Derive segments from a :ARGRun via :PARENT_OF traversal.
      kinship   Browning-style kinship from total IBD length per pair.
      list      List ingested IBD sources.
      delete    Remove all segments for a given source.
    """


@ibd.command("ingest")
@click.argument("tsv_path", type=click.Path(exists=True, dir_okay=False))
@click.option("--source", default="hap_ibd", show_default=True,
              help="Identifier for the IBD source")
@click.option("--has-header", is_flag=True,
              help="Skip the first line of the TSV")
@click.option("--replace", is_flag=True,
              help="Delete existing segments with the same source first")
@pass_ctx
def ingest(ctx, tsv_path, source, has_header, replace):
    """Ingest a hap-IBD-style 6-column TSV.

    TSV columns (no header by default):

      sample_a<TAB>sample_b<TAB>chr<TAB>start_bp<TAB>end_bp<TAB>length_cM

    length_cM is optional. Lines starting with '#' are skipped.
    """
    from graphpop_import.ibd_ingester import IBDIngester

    ing = IBDIngester(ctx.driver, database=ctx.database)
    summary = ing.ingest_tsv(tsv_path, source=source,
                               has_header=has_header, replace=replace)
    click.echo(
        f"Ingested {summary.n_segments} :IBD_SEGMENT edges over "
        f"{summary.n_pairs} pairs for source={source}.",
        err=True)


@ibd.command("from-arg")
@click.argument("run_id")
@click.option("--max-tmrca", type=float, default=None,
              help="Drop segments whose MRCA time exceeds this (generations)")
@click.option("--min-length-bp", type=int, default=0,
              help="Drop segments shorter than this")
@click.option("--chr", default="chr1", show_default=True,
              help="Chromosome label to record on emitted edges")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def from_arg(ctx, run_id, max_tmrca, min_length_bp, chr, output_path, fmt):
    """Derive IBD segments from an :ARGRun (graph-native; equivalent to
    tskit.TreeSequence.ibd_segments).

    Writes :IBD_SEGMENT edges with source='arg_derived' and streams
    one row per segment back. Idempotent per run_id.
    """
    opts: dict[str, object] = {"chr": chr}
    if max_tmrca is not None:
        opts["max_tmrca"] = max_tmrca
    if min_length_bp:
        opts["min_length_bp"] = min_length_bp

    cypher = build_cypher(
        "graphpop.ibd.from_arg",
        [f"'{run_id}'"],
        options=opts if opts else None,
        yield_cols=["sample_a", "sample_b", "chr", "start", "end",
                    "length_bp", "mrca_node_id", "tmrca", "source", "runId"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "ibd-from-arg",
                  {"run_id": run_id})


@ibd.command("kinship")
@click.argument("source")
@click.option("--min-length-bp", type=int, default=0,
              help="Skip pairs with total IBD length below this")
@click.option("--total-genome-cm", type=float, default=None,
              help="Override the total genome length (cM) for the cM-mode "
                   "denominator")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def kinship(ctx, source, min_length_bp, total_genome_cm, output_path, fmt):
    """Browning-style kinship from total IBD length per pair.

    Uses length_cM when present on every segment; otherwise falls back
    to length_bp with method='ibd_bp'.
    """
    opts: dict[str, object] = {}
    if min_length_bp:
        opts["min_length_bp"] = min_length_bp
    if total_genome_cm is not None:
        opts["total_genome_cM"] = total_genome_cm

    cypher = build_cypher(
        "graphpop.ibd.kinship",
        [f"'{source}'"],
        options=opts if opts else None,
        yield_cols=["sample_a", "sample_b", "phi", "ibs0", "het_het",
                    "n_snp", "n_aa_min", "method"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "ibd-kinship",
                  {"source": source})


@ibd.command("list")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def list_(ctx, output_path, fmt):
    """List every IBD source currently in the database."""
    from graphpop_import.ibd_ingester import IBDIngester

    ing = IBDIngester(ctx.driver, database=ctx.database)
    sources = ing.list_sources()
    records = [
        {"source": s.source, "n_edges": s.n_edges, "n_pairs": s.n_pairs}
        for s in sources
    ]
    format_output(records, output_path, fmt, "ibd-list", {})


@ibd.command("delete")
@click.option("--source", required=True)
@click.option("--yes", is_flag=True, help="Skip confirmation prompt")
@pass_ctx
def delete(ctx, source, yes):
    """Delete every :IBD_SEGMENT edge for the given source."""
    if not yes:
        click.confirm(f"Delete all :IBD_SEGMENT edges for source='{source}'?",
                       abort=True)
    from graphpop_import.ibd_ingester import IBDIngester

    ing = IBDIngester(ctx.driver, database=ctx.database)
    n = ing.delete_source(source)
    click.echo(f"Deleted {n} :IBD_SEGMENT edges.", err=True)
