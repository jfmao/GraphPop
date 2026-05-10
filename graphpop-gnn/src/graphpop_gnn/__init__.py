"""GNN embeddings for GraphPop (M10)."""

from .exporter import GraphExporter, ExportedGraph
from .persist import EmbeddingPersister

__all__ = ["GraphExporter", "ExportedGraph", "EmbeddingPersister"]
__version__ = "0.2.0.dev0"
