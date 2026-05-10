"""graphpop-gnn CLI: export, train, embed.

Lazy imports torch / torch_geometric so non-training subcommands
work without the [gpu] extras installed.
"""
from __future__ import annotations

import os
import sys

import click


@click.group()
def main():
    """GNN embeddings for GraphPop."""


def _build_driver():
    from neo4j import GraphDatabase

    uri = os.environ.get("GRAPHPOP_URI", "bolt://localhost:7687")
    user = os.environ.get("GRAPHPOP_USER", "neo4j")
    password = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
    return GraphDatabase.driver(uri, auth=(user, password))


@main.command()
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True),
              help="Output .npz path")
@click.option("--min-samples", type=int, default=1, show_default=True,
              help="Drop variants with fewer than this many CARRIES edges")
@click.option("--limit-variants", type=int, default=None,
              help="Cap the number of variants (for smoke runs)")
@click.option("--database", default="neo4j", show_default=True)
def export(output: str, min_samples: int,
            limit_variants: int | None, database: str) -> None:
    """Pull the :Sample × :Variant bipartite graph from Neo4j → .npz."""
    from .exporter import GraphExporter

    driver = _build_driver()
    try:
        ex = GraphExporter(driver, database=database)
        graph = ex.export(min_samples=min_samples,
                            limit_variants=limit_variants)
        graph.save(output)
        click.echo(
            f"Exported {graph.n_samples} samples × {graph.n_variants} "
            f"variants × {graph.n_edges} edges → {output}",
            err=True,
        )
    finally:
        driver.close()


@main.command()
@click.option("--input", "input_path", required=True,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True),
              help="Output .pt model path")
@click.option("--epochs", type=int, default=50, show_default=True)
@click.option("--batch-size", type=int, default=256, show_default=True)
@click.option("--lr", type=float, default=1e-3, show_default=True)
@click.option("--out-dim", type=int, default=32, show_default=True)
@click.option("--seed", type=int, default=42, show_default=True)
def train(input_path: str, output: str, epochs: int, batch_size: int,
          lr: float, out_dim: int, seed: int) -> None:
    """Train GraphSAGE + InfoNCE on the exported graph."""
    try:
        import torch  # noqa: F401  (presence check)
    except ImportError as exc:
        raise click.ClickException(
            "Training requires torch + torch-geometric; install with "
            "`pip install graphpop-gnn[gpu]`."
        ) from exc

    from .exporter import ExportedGraph
    from .train import TrainConfig, train as train_fn
    import numpy as np
    import torch

    graph = ExportedGraph.load(input_path)
    cfg = TrainConfig(
        epochs=epochs, batch_size=batch_size, lr=lr,
        out_dim=out_dim, seed=seed,
    )
    model, embeddings = train_fn(graph, cfg)
    torch.save(
        {
            "state_dict": model.state_dict(),
            "config": {
                "n_samples": graph.n_samples,
                "n_variants": graph.n_variants,
                "init_dim": cfg.init_dim,
                "hidden_dim": cfg.hidden_dim,
                "out_dim": cfg.out_dim,
            },
            "sample_ids": graph.sample_ids,
            "embeddings": embeddings,
        },
        output,
    )
    click.echo(
        f"Trained {epochs} epochs over {graph.n_samples} samples; "
        f"saved model + embeddings to {output}",
        err=True,
    )


@main.command()
@click.option("--model", "model_path", required=True,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--persist/--no-persist", default=True,
              help="Write embeddings as :Sample.embedding properties")
@click.option("--database", default="neo4j", show_default=True)
def embed(model_path: str, persist: bool, database: str) -> None:
    """Persist a trained model's embeddings to Neo4j."""
    try:
        import torch  # noqa: F401
    except ImportError as exc:
        raise click.ClickException(
            "Loading model checkpoints requires torch; install with "
            "`pip install graphpop-gnn[gpu]`."
        ) from exc

    import torch
    from .persist import EmbeddingPersister

    ckpt = torch.load(model_path, map_location="cpu", weights_only=False)
    sample_ids = ckpt["sample_ids"]
    embeddings = ckpt["embeddings"]

    if not persist:
        click.echo(
            f"Loaded {len(sample_ids)} embeddings × dim "
            f"{embeddings.shape[1]} (not persisted; --no-persist)",
            err=True,
        )
        return
    driver = _build_driver()
    try:
        per = EmbeddingPersister(driver, database=database)
        n = per.persist(sample_ids, embeddings)
        click.echo(
            f"Persisted {n} :Sample.embedding properties.", err=True)
    finally:
        driver.close()


if __name__ == "__main__":
    main()
