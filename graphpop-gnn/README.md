# graphpop-gnn

GNN-based sample embeddings for GraphPop (M10). Trains
self-supervised graph embeddings on the `:Sample` × `:Variant`
bipartite graph using PyTorch Geometric, persists them as
`:Sample.embedding` (vector property), and exposes Java procedures
(`graphpop.embedding.knn`, `graphpop.embedding.cluster`) for query.

## Install

The package itself has light deps (neo4j-driver, numpy, click).
Training and inference require the optional `[gpu]` extra:

```sh
pip install graphpop-gnn                     # query / persist only
pip install "graphpop-gnn[gpu]"             # training pipeline (torch + PyG)
```

## CLI

```sh
graphpop-gnn export --output cohort.npz                 # dump :CARRIES bipartite
graphpop-gnn train --input cohort.npz --output model.pt # GraphSAGE + InfoNCE
graphpop-gnn embed --model model.pt --persist           # write :Sample.embedding
```

Or via the main `graphpop` CLI: `graphpop gnn …`.

## Defaults

- 2-layer GraphSAGE
- Hidden dim 64, output dim 32
- InfoNCE self-supervised contrastive loss
- Adam, lr 1e-3
- 50 epochs, batch size 256

Smoke tests run on CPU on a 100-variant fixture; full training
requires the RTX 4090.
