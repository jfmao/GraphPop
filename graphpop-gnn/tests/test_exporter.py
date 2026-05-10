"""Tests for the torch-free exporter / persist modules.

These tests use a tiny in-memory mock Neo4j driver — they verify
the data-shape contract without requiring a live Neo4j or torch.
"""
from __future__ import annotations

import numpy as np
import pytest

from graphpop_gnn import EmbeddingPersister, ExportedGraph, GraphExporter


class _MockResult:
    def __init__(self, rows):
        self._rows = rows
        self._idx = 0

    def __iter__(self):
        for row in self._rows:
            yield row

    def single(self):
        return self._rows[0] if self._rows else None


class _MockSession:
    def __init__(self, query_handlers):
        self._handlers = query_handlers

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def run(self, cypher, params=None):
        for matcher, handler in self._handlers:
            if matcher in cypher:
                return _MockResult(handler(params or {}))
        return _MockResult([])


class _MockDriver:
    def __init__(self, query_handlers):
        self._handlers = query_handlers

    def session(self, database="neo4j"):
        return _MockSession(self._handlers)


def test_exporter_collapses_samples_variants_carries_into_packed_arrays():
    handlers = [
        ("MATCH (s:Sample) ", lambda p: [
            {"sid": "S1", "pop": "EUR"},
            {"sid": "S2", "pop": "EUR"},
            {"sid": "S3", "pop": "AFR"},
        ]),
        ("MATCH (v:Variant)", lambda p: [
            {"vid": "v1"},
            {"vid": "v2"},
        ]),
        ("MATCH (s:Sample)-[r:CARRIES]->(v:Variant)", lambda p: [
            {"sid": "S1", "vid": "v1", "gt": 1},
            {"sid": "S2", "vid": "v1", "gt": 2},
            {"sid": "S3", "vid": "v2", "gt": 1},
        ]),
    ]
    ex = GraphExporter(_MockDriver(handlers))
    g = ex.export()
    assert g.n_samples == 3
    assert g.n_variants == 2
    assert g.n_edges == 3
    assert list(g.sample_ids) == ["S1", "S2", "S3"]
    assert list(g.variant_ids) == ["v1", "v2"]
    assert list(g.sample_idx) == [0, 1, 2]
    assert list(g.variant_idx) == [0, 0, 1]
    assert list(g.gt) == [1, 2, 1]
    assert g.populations == ["EUR", "EUR", "AFR"]


def test_exporter_round_trip_save_load(tmp_path):
    g = ExportedGraph(
        sample_ids=["A", "B"],
        variant_ids=["x", "y", "z"],
        sample_idx=np.array([0, 0, 1], dtype=np.int64),
        variant_idx=np.array([0, 1, 2], dtype=np.int64),
        gt=np.array([1, 2, 1], dtype=np.int8),
        populations=["EUR", "AFR"],
    )
    path = tmp_path / "g.npz"
    g.save(path)
    g2 = ExportedGraph.load(path)
    assert g2.sample_ids == g.sample_ids
    assert g2.variant_ids == g.variant_ids
    np.testing.assert_array_equal(g2.sample_idx, g.sample_idx)
    np.testing.assert_array_equal(g2.variant_idx, g.variant_idx)
    np.testing.assert_array_equal(g2.gt, g.gt)
    assert g2.populations == g.populations


def test_persister_pushes_embeddings_in_batches():
    persisted = []
    cleared = []

    def carries_handler(p):
        return [{"c": len(p["batch"])}]

    def clear_handler(p):
        return [{"c": 5}]

    handlers = [
        ("MATCH (s:Sample {sampleId:", carries_handler),
        ("UNWIND $batch", carries_handler),
        ("REMOVE s.embedding", clear_handler),
    ]
    drv = _MockDriver(handlers)
    per = EmbeddingPersister(drv)
    emb = np.array([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]], dtype=np.float32)
    n = per.persist(["A", "B", "C"], emb, batch_size=2)
    # First batch (2 rows) + second (1 row) = 3 total.
    assert n == 3


def test_persister_validates_shape():
    per = EmbeddingPersister(_MockDriver([]))
    with pytest.raises(ValueError):
        per.persist(["A"], np.array([1.0, 2.0]))  # not 2-D
    with pytest.raises(ValueError):
        per.persist(["A", "B"], np.array([[1.0]]))  # length mismatch


def test_persister_clear_runs_remove_query():
    handlers = [
        ("REMOVE s.embedding", lambda p: [{"c": 7}]),
    ]
    per = EmbeddingPersister(_MockDriver(handlers))
    assert per.clear() == 7


def test_exported_graph_empty_sane_shapes():
    g = ExportedGraph()
    assert g.n_samples == 0
    assert g.n_variants == 0
    assert g.n_edges == 0
