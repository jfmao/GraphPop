"""Self-supervised GraphSAGE training loop.

Imports torch on module load; not safe to import from a torch-less
environment. Use lazy imports in CLI / orchestrator code.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from torch import optim

from .exporter import ExportedGraph
from .model import SampleEmbeddingModel, info_nce_loss


@dataclass
class TrainConfig:
    epochs: int = 50
    batch_size: int = 256
    lr: float = 1e-3
    init_dim: int = 32
    hidden_dim: int = 64
    out_dim: int = 32
    temperature: float = 0.1
    n_negatives: int = 16
    seed: int = 42
    device: str = "cuda" if torch.cuda.is_available() else "cpu"


def build_edge_index(graph: ExportedGraph) -> torch.Tensor:
    """Build the undirected bipartite edge_index for PyG.

    Variants are indexed in `[n_samples, n_samples + n_variants)`.
    Returns shape `(2, 2 * n_edges)` (each undirected edge stored
    twice).
    """
    n = graph.n_samples
    src = np.concatenate([graph.sample_idx, n + graph.variant_idx])
    dst = np.concatenate([n + graph.variant_idx, graph.sample_idx])
    return torch.tensor(np.stack([src, dst]), dtype=torch.long)


def sample_co_carrier_pairs(
    graph: ExportedGraph,
    n_pairs: int,
    rng: np.random.Generator,
) -> torch.Tensor:
    """Build positive pairs: samples that co-carry the same variant.

    Picks random variants weighted by carrier-count, then picks
    a random pair of carriers. Returns `(n_pairs, 2)` tensor of
    sample indices.
    """
    if graph.n_edges == 0 or n_pairs <= 0:
        return torch.zeros((0, 2), dtype=torch.long)
    # Group sample_idx by variant_idx.
    order = np.argsort(graph.variant_idx, kind="stable")
    s_sorted = graph.sample_idx[order]
    v_sorted = graph.variant_idx[order]
    # Find run lengths per variant.
    boundaries = np.concatenate([
        [0], np.where(np.diff(v_sorted) != 0)[0] + 1, [len(v_sorted)]
    ])
    runs = []
    for i in range(len(boundaries) - 1):
        lo, hi = boundaries[i], boundaries[i + 1]
        if hi - lo >= 2:
            runs.append((lo, hi))
    if not runs:
        return torch.zeros((0, 2), dtype=torch.long)
    counts = np.array([hi - lo for lo, hi in runs], dtype=np.int64)
    weights = counts.astype(np.float64)
    weights = weights / weights.sum()
    chosen_runs = rng.choice(len(runs), size=n_pairs, p=weights)
    pairs = np.empty((n_pairs, 2), dtype=np.int64)
    for k, ri in enumerate(chosen_runs):
        lo, hi = runs[ri]
        a, b = rng.choice(np.arange(lo, hi), size=2, replace=False)
        pairs[k, 0] = s_sorted[a]
        pairs[k, 1] = s_sorted[b]
    return torch.tensor(pairs, dtype=torch.long)


def train(graph: ExportedGraph, cfg: TrainConfig | None = None
          ) -> tuple[SampleEmbeddingModel, np.ndarray]:
    """Train and return (model, sample_embeddings_numpy)."""
    cfg = cfg or TrainConfig()
    torch.manual_seed(cfg.seed)
    rng = np.random.default_rng(cfg.seed)

    device = torch.device(cfg.device)
    model = SampleEmbeddingModel(
        n_samples=graph.n_samples,
        n_variants=graph.n_variants,
        init_dim=cfg.init_dim,
        hidden_dim=cfg.hidden_dim,
        out_dim=cfg.out_dim,
    ).to(device)
    edge_index = build_edge_index(graph).to(device)
    optim_ = optim.Adam(model.parameters(), lr=cfg.lr)

    model.train()
    for epoch in range(cfg.epochs):
        positive_pairs = sample_co_carrier_pairs(
            graph, cfg.batch_size, rng).to(device)
        emb = model.sample_embeddings(edge_index)
        loss = info_nce_loss(
            emb, positive_pairs,
            n_negatives=cfg.n_negatives,
            temperature=cfg.temperature,
        )
        optim_.zero_grad()
        loss.backward()
        optim_.step()

    model.eval()
    with torch.no_grad():
        sample_emb = model.sample_embeddings(edge_index).cpu().numpy()
    return model, sample_emb
