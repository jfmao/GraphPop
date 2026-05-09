"""Tests for graphpop_import.ancestry_ingester (mocked driver)."""
from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from graphpop_import.ancestry_ingester import (
    AncestryIngester,
    PaintingRow,
    _Topology,
    _propagate_majority,
    _chunked,
)


# ---- Recording mock helpers ------------------------------------------------


class _RecordingTx:
    def __init__(self, read_handlers=None):
        self.calls = []
        self._read_handlers = read_handlers or []

    def run(self, cypher, **kwargs):
        self.calls.append((cypher, kwargs))
        for prefix, rows in self._read_handlers:
            if cypher.startswith(prefix):
                return _Result(rows)
        return _Result([])


class _Result:
    def __init__(self, rows):
        self._rows = rows

    def __iter__(self):
        for r in self._rows:
            yield _Record(r)

    def single(self):
        return _Record(self._rows[0]) if self._rows else None


class _Record(dict):
    pass


class _Session:
    def __init__(self, read_handlers=None):
        self.tx = _RecordingTx(read_handlers)

    def __enter__(self): return self
    def __exit__(self, *args): return False

    def execute_write(self, callback, *args, **kwargs):
        return callback(self.tx, *args, **kwargs)

    def execute_read(self, callback, *args, **kwargs):
        return callback(self.tx, *args, **kwargs)


def _mock_driver(read_handlers=None):
    driver = MagicMock()
    session = _Session(read_handlers)
    driver.session.return_value = session
    return driver, session


# ---- Direct ingest --------------------------------------------------------


def test_ingest_creates_population_nodes_and_edges():
    driver, session = _mock_driver()
    ing = AncestryIngester(driver)
    rows = [
        PaintingRow(0, "EUR", 1.0),
        PaintingRow(1, "EUR", 1.0),
        PaintingRow(2, "AFR", 1.0),
    ]
    summary = ing.ingest("r0", rows, painter="manual")

    pop_calls = [c for c in session.tx.calls if c[0].startswith("MERGE (p:Population")]
    assert {c[1]["pop"] for c in pop_calls} == {"EUR", "AFR"}

    edge_calls = [c for c in session.tx.calls
                  if c[0].startswith("UNWIND $rows AS r")]
    assert len(edge_calls) == 1
    assert len(edge_calls[0][1]["rows"]) == 3
    assert edge_calls[0][1]["runId"] == "r0"
    assert edge_calls[0][1]["painter"] == "manual"

    assert summary.n_edges == 3
    assert summary.n_painted_nodes == 3
    assert summary.n_populations == 2


def test_ingest_default_replace_false_raises_on_conflict():
    driver, session = _mock_driver(
        read_handlers=[
            ("MATCH ()-[r:HAS_ANCESTRY", [{"c": 5}]),
        ]
    )
    ing = AncestryIngester(driver)
    with pytest.raises(RuntimeError, match="already exists"):
        ing.ingest("r0", [PaintingRow(0, "EUR")], painter="manual")


def test_ingest_replace_true_deletes_existing():
    driver, session = _mock_driver(
        read_handlers=[
            ("MATCH ()-[r:HAS_ANCESTRY", [{"n_deleted": 5}]),
        ]
    )
    ing = AncestryIngester(driver)
    ing.ingest("r0", [PaintingRow(0, "EUR")],
               painter="manual", replace=True)

    delete_calls = [c for c in session.tx.calls
                    if "DELETE r RETURN n_deleted" in c[0]]
    assert len(delete_calls) == 1


def test_ingest_probabilistic_painting():
    driver, session = _mock_driver()
    ing = AncestryIngester(driver)
    ing.ingest("r0", [
        PaintingRow(0, "EUR", 0.6),
        PaintingRow(0, "AFR", 0.4),
    ], painter="rfmix")

    edge_calls = [c for c in session.tx.calls
                  if c[0].startswith("UNWIND $rows AS r")]
    rows = edge_calls[0][1]["rows"]
    probs = sorted(r["posterior_prob"] for r in rows)
    assert probs == [0.4, 0.6]


# ---- Sample-level propagation ---------------------------------------------


def test_propagate_majority_balanced_two_ancestries():
    """4-leaf balanced tree (fractional painting preserves the partition):
            6
           / \\
          4   5
         /\\   /\\
        0 1  2 3
    Samples 0,1 = EUR; 2,3 = AFR. Internal nodes 4 = EUR×1.0,
    5 = AFR×1.0, 6 = EUR×0.5 + AFR×0.5 (two edges)."""
    topology = _Topology(
        nodes=[
            {"nodeId": 0, "is_sample": True},
            {"nodeId": 1, "is_sample": True},
            {"nodeId": 2, "is_sample": True},
            {"nodeId": 3, "is_sample": True},
            {"nodeId": 4, "is_sample": False},
            {"nodeId": 5, "is_sample": False},
            {"nodeId": 6, "is_sample": False},
        ],
        edges=[(4, 0), (4, 1), (5, 2), (5, 3), (6, 4), (6, 5)],
        sample_to_node={"hap_0": 0, "hap_1": 1, "hap_2": 2, "hap_3": 3},
    )
    sample_anc = {"hap_0": "EUR", "hap_1": "EUR", "hap_2": "AFR", "hap_3": "AFR"}
    rows = _propagate_majority(topology, sample_anc)

    by_node: dict[int, dict[str, float]] = {}
    for r in rows:
        by_node.setdefault(r.tskit_node_id, {})[r.population_id] = r.posterior_prob

    assert by_node[0] == {"EUR": 1.0}
    assert by_node[1] == {"EUR": 1.0}
    assert by_node[2] == {"AFR": 1.0}
    assert by_node[3] == {"AFR": 1.0}
    assert by_node[4] == {"EUR": 1.0}
    assert by_node[5] == {"AFR": 1.0}
    # Root has 2 EUR + 2 AFR -> two edges, each at 0.5
    assert by_node[6] == {"EUR": 0.5, "AFR": 0.5}


def test_propagate_majority_skewed():
    """Three samples, all EUR -> internal node also EUR with prob 1.0."""
    topology = _Topology(
        nodes=[
            {"nodeId": 0, "is_sample": True},
            {"nodeId": 1, "is_sample": True},
            {"nodeId": 2, "is_sample": True},
            {"nodeId": 3, "is_sample": False},
        ],
        edges=[(3, 0), (3, 1), (3, 2)],
        sample_to_node={"a": 0, "b": 1, "c": 2},
    )
    sample_anc = {"a": "EUR", "b": "EUR", "c": "EUR"}
    rows = _propagate_majority(topology, sample_anc)
    by_id = {r.tskit_node_id: r for r in rows}
    assert by_id[3].population_id == "EUR"
    assert by_id[3].posterior_prob == 1.0


def test_propagate_majority_skips_unmapped_samples():
    """A sample not in the ancestry dict must not contribute votes."""
    topology = _Topology(
        nodes=[
            {"nodeId": 0, "is_sample": True},
            {"nodeId": 1, "is_sample": True},
            {"nodeId": 2, "is_sample": False},
        ],
        edges=[(2, 0), (2, 1)],
        sample_to_node={"hap_0": 0},  # hap_1 not exposed
    )
    rows = _propagate_majority(topology, {"hap_0": "EUR"})
    by_id = {r.tskit_node_id: r for r in rows}
    # Internal node sees only one descendant labelled EUR.
    assert 0 in by_id and by_id[0].population_id == "EUR"
    assert 2 in by_id and by_id[2].population_id == "EUR"
    # Sample hap_1 has no painting -> no row.
    assert 1 not in by_id


def test_propagate_majority_empty_input_yields_no_rows():
    topology = _Topology(
        nodes=[{"nodeId": 0, "is_sample": True}],
        edges=[],
        sample_to_node={},
    )
    assert _propagate_majority(topology, {}) == []


# ---- list / delete --------------------------------------------------------


def test_list_paintings_all():
    canned = [
        {"runId": "r0", "painter": "manual",
         "n_nodes": 7, "n_edges": 7, "n_pops": 2},
        {"runId": "r1", "painter": "rfmix",
         "n_nodes": 5, "n_edges": 6, "n_pops": 3},
    ]
    driver, _ = _mock_driver(
        read_handlers=[("MATCH (n:TreeNode)-[r:HAS_ANCESTRY]->", canned)])
    ing = AncestryIngester(driver)
    out = ing.list_paintings()
    assert len(out) == 2
    assert out[0].run_id == "r0" and out[0].n_edges == 7
    assert out[1].run_id == "r1"


def test_delete_painting_returns_count():
    driver, _ = _mock_driver(
        read_handlers=[("MATCH ()-[r:HAS_ANCESTRY", [{"n_deleted": 12}])])
    ing = AncestryIngester(driver)
    assert ing.delete_painting("r0", "manual") == 12


# ---- chunked --------------------------------------------------------------


def test_chunked_basic():
    assert list(_chunked(range(7), 3)) == [[0, 1, 2], [3, 4, 5], [6]]
