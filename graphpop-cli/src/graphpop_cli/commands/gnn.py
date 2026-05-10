"""graphpop gnn — GNN embeddings (M10).

Thin facade over graphpop-gnn (training pipeline) + the Java
:Sample.embedding query procedures (knn, cluster).
"""
from __future__ import annotations

import subprocess
import sys

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def gnn():
    """GNN embeddings on the :Sample × :Variant graph (M10).

    Subcommands:

      export   Dump :CARRIES bipartite to .npz (for training).
      train    GraphSAGE + InfoNCE training (requires graphpop-gnn[gpu]).
      embed    Persist trained model embeddings as :Sample.embedding.
      knn      Top-k nearest neighbours by cosine similarity.
      cluster  k-means clustering over :Sample.embedding.
    """


@gnn.command("export")
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True))
@click.option("--min-samples", type=int, default=1, show_default=True)
@click.option("--limit-variants", type=int, default=None)
@pass_ctx
def export(ctx, output, min_samples, limit_variants):
    """Delegate to graphpop-gnn export (torch-free)."""
    args = [
        sys.executable, "-m", "graphpop_gnn.cli", "export",
        "--output", output,
        "--min-samples", str(min_samples),
        "--database", ctx.database,
    ]
    if limit_variants is not None:
        args += ["--limit-variants", str(limit_variants)]
    raise SystemExit(subprocess.call(args))


@gnn.command("train")
@click.option("--input", "input_path", required=True,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True))
@click.option("--epochs", type=int, default=50, show_default=True)
@click.option("--batch-size", type=int, default=256, show_default=True)
@click.option("--lr", type=float, default=1e-3, show_default=True)
@click.option("--out-dim", type=int, default=32, show_default=True)
@click.option("--seed", type=int, default=42, show_default=True)
def train(input_path, output, epochs, batch_size, lr, out_dim, seed):
    """Delegate to graphpop-gnn train (requires [gpu] extras)."""
    args = [
        sys.executable, "-m", "graphpop_gnn.cli", "train",
        "--input", input_path,
        "--output", output,
        "--epochs", str(epochs),
        "--batch-size", str(batch_size),
        "--lr", str(lr),
        "--out-dim", str(out_dim),
        "--seed", str(seed),
    ]
    raise SystemExit(subprocess.call(args))


@gnn.command("embed")
@click.option("--model", "model_path", required=True,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--persist/--no-persist", default=True)
@pass_ctx
def embed(ctx, model_path, persist):
    """Persist trained embeddings as :Sample.embedding."""
    args = [
        sys.executable, "-m", "graphpop_gnn.cli", "embed",
        "--model", model_path,
        "--database", ctx.database,
        "--persist" if persist else "--no-persist",
    ]
    raise SystemExit(subprocess.call(args))


@gnn.command("knn")
@click.argument("sample_id")
@click.option("--k", type=int, default=10, show_default=True)
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def knn(ctx, sample_id, k, output_path, fmt):
    """Top-k nearest neighbours by cosine over :Sample.embedding."""
    cypher = build_cypher(
        "graphpop.embedding.knn",
        [f"'{sample_id}'", str(k)],
        yield_cols=["query_sample_id", "neighbor_sample_id",
                    "cosine_similarity", "rank"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "embedding-knn",
                  {"sample_id": sample_id, "k": k})


@gnn.command("cluster")
@click.argument("method", type=click.Choice(["kmeans"]))
@click.option("--k", type=int, required=True,
              help="Number of clusters")
@click.option("--seed", type=int, default=42, show_default=True)
@click.option("--max-iter", type=int, default=100, show_default=True)
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def cluster(ctx, method, k, seed, max_iter, output_path, fmt):
    """k-means clustering over :Sample.embedding."""
    opts = {"seed": seed, "max_iter": max_iter}
    cypher = build_cypher(
        "graphpop.embedding.cluster",
        [f"'{method}'", str(k)],
        options=opts,
        yield_cols=["sample_id", "cluster_id", "distance_to_centroid",
                    "n_clusters", "method"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "embedding-cluster",
                  {"method": method, "k": k})
