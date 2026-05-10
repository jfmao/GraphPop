"""GraphSAGE / GAT model definitions.

Imports torch and torch_geometric on module load; tests should skip
this module when those deps are missing. The CLI imports lazily so
non-training commands stay usable without [gpu] installed.
"""
from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn
from torch_geometric.nn import GraphSAGE, SAGEConv


class SampleEmbeddingModel(nn.Module):
    """2-layer GraphSAGE on the bipartite Sample × Variant graph.

    Input: per-node feature tensor (one-hot or learned embedding
    initialisation; we use a learnable tensor of shape
    `(n_samples + n_variants, init_dim)`).

    Output: per-sample embedding of dimension `out_dim`.
    """

    def __init__(
        self,
        n_samples: int,
        n_variants: int,
        init_dim: int = 32,
        hidden_dim: int = 64,
        out_dim: int = 32,
    ):
        super().__init__()
        self.n_samples = n_samples
        self.n_variants = n_variants
        self.init_embedding = nn.Embedding(n_samples + n_variants, init_dim)
        self.conv1 = SAGEConv(init_dim, hidden_dim)
        self.conv2 = SAGEConv(hidden_dim, out_dim)

    def forward(self, edge_index: torch.Tensor) -> torch.Tensor:
        x = self.init_embedding.weight
        h = F.relu(self.conv1(x, edge_index))
        h = self.conv2(h, edge_index)
        return h

    def sample_embeddings(self, edge_index: torch.Tensor) -> torch.Tensor:
        """Return per-sample embeddings (first `n_samples` rows)."""
        all_emb = self(edge_index)
        return all_emb[: self.n_samples]


def info_nce_loss(
    embeddings: torch.Tensor,
    positive_pairs: torch.Tensor,
    n_negatives: int = 16,
    temperature: float = 0.1,
) -> torch.Tensor:
    """Contrastive InfoNCE loss over positive pairs.

    embeddings : `(n_samples, dim)` per-sample embedding tensor
    positive_pairs : `(P, 2)` integer tensor; row i = (anchor_i,
                      positive_i) of sample indices that should be
                      pulled together (e.g., samples that co-carry
                      a rare variant).
    n_negatives : number of in-batch negatives per anchor.
    """
    if positive_pairs.numel() == 0:
        return torch.zeros((), device=embeddings.device)
    anchors = embeddings[positive_pairs[:, 0]]
    positives = embeddings[positive_pairs[:, 1]]
    n = embeddings.shape[0]
    neg_idx = torch.randint(
        0, n, (positive_pairs.shape[0], n_negatives),
        device=embeddings.device,
    )
    negatives = embeddings[neg_idx]  # (P, n_negatives, dim)

    a = F.normalize(anchors, dim=-1)
    p = F.normalize(positives, dim=-1)
    n_n = F.normalize(negatives, dim=-1)

    pos_logit = (a * p).sum(dim=-1, keepdim=True) / temperature
    neg_logits = torch.einsum("pd,pnd->pn", a, n_n) / temperature
    logits = torch.cat([pos_logit, neg_logits], dim=-1)  # (P, 1+N)
    labels = torch.zeros(logits.shape[0], dtype=torch.long,
                          device=embeddings.device)
    return F.cross_entropy(logits, labels)
