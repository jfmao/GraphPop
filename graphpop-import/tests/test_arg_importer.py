"""Tests for graphpop_import.arg_importer.

Uses real msprime simulations as fixtures (msprime is in the [arg]
extras) and a hand-rolled session/tx mock to verify the Cypher emitted
without requiring a live Neo4j.
"""
from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime
from typing import Any
from unittest.mock import MagicMock

import pytest

msprime = pytest.importorskip("msprime")
tskit = pytest.importorskip("tskit")

from graphpop_import.arg_importer import ARGIngester, ARGRunSummary, _chunked


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tiny_ts():
    """A 4-sample (2-diploid), 100bp, no-recombination tree sequence."""
    ts = msprime.sim_ancestry(
        samples=2,
        sequence_length=100,
        recombination_rate=0,
        random_seed=42,
    )
    ts = msprime.sim_mutations(ts, rate=1e-2, random_seed=42)
    return ts


@pytest.fixture
def small_ts():
    """A 20-sample (10-diploid), 50kb, recombining tree sequence with mutations.

    Bumped from 5 kb / mu=1e-5 in step 1 of M4.1 step 3: at 50 kb / mu=1e-4
    we get >100 mutations and multiple marginal trees, sufficient for eGRM
    validation downstream.
    """
    ts = msprime.sim_ancestry(
        samples=10,
        sequence_length=50_000,
        recombination_rate=1e-5,
        random_seed=43,
    )
    ts = msprime.sim_mutations(ts, rate=1e-4, random_seed=43)
    return ts


# ---------------------------------------------------------------------------
# Recording session/tx mocks
# ---------------------------------------------------------------------------


class _RecordingTx:
    """Mock tx that records every (cypher, kwargs) pair that runs."""

    def __init__(self, read_handlers: dict[str, list[dict]] | None = None):
        self.calls: list[tuple[str, dict[str, Any]]] = []
        self._read_handlers = read_handlers or {}

    def run(self, cypher: str, **kwargs: Any):
        # accept the alternative driver-style positional rows= via **kwargs
        self.calls.append((cypher, kwargs))
        # If a test supplied a canned reader for this cypher, return rows.
        for prefix, rows in self._read_handlers.items():
            if cypher.startswith(prefix):
                return _RecordingResult(rows)
        return _RecordingResult([])


class _RecordingResult:
    def __init__(self, rows: list[dict]):
        self._rows = rows

    def __iter__(self):
        for r in self._rows:
            yield _RecordingRecord(r)

    def single(self):
        return _RecordingRecord(self._rows[0]) if self._rows else None


class _RecordingRecord(dict):
    """Subclass of dict so callers can use record["key"] like Neo4j Records."""
    pass


class _RecordingSession:
    def __init__(self, read_handlers: dict[str, list[dict]] | None = None):
        self.tx = _RecordingTx(read_handlers)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def execute_write(self, callback, *args, **kwargs):
        return callback(self.tx, *args, **kwargs)

    def execute_read(self, callback, *args, **kwargs):
        return callback(self.tx, *args, **kwargs)


def _mock_driver(read_handlers: dict[str, list[dict]] | None = None):
    driver = MagicMock()
    session = _RecordingSession(read_handlers)
    driver.session.return_value = session
    return driver, session


# ---------------------------------------------------------------------------
# helpers: chunking
# ---------------------------------------------------------------------------


def test_chunked_splits_evenly():
    assert list(_chunked(range(5), 2)) == [[0, 1], [2, 3], [4]]


def test_chunked_handles_empty():
    assert list(_chunked([], 100)) == []


# ---------------------------------------------------------------------------
# ingest() orchestration
# ---------------------------------------------------------------------------


def test_ingest_creates_arg_run_node_with_summary_props(tiny_ts):
    driver, session = _mock_driver()
    ingester = ARGIngester(driver)

    summary = ingester.ingest(
        tiny_ts, run_id="t01", source="msprime", params={"seed": 42}
    )

    assert isinstance(summary, ARGRunSummary)
    assert summary.run_id == "t01"
    assert summary.source == "msprime"
    assert summary.n_samples == int(tiny_ts.num_samples)
    assert summary.n_trees == int(tiny_ts.num_trees)
    assert summary.n_nodes == int(tiny_ts.num_nodes)

    # First N statements are arg_schema constraints+indices, then the
    # CREATE (r:ARGRun ...). Find the ARGRun create call and check kwargs.
    arg_run_calls = [
        c for c in session.tx.calls if c[0].startswith("CREATE (r:ARGRun")
    ]
    assert len(arg_run_calls) == 1
    cypher, kwargs = arg_run_calls[0]
    assert kwargs["runId"] == "t01"
    assert kwargs["source"] == "msprime"
    assert kwargs["n_samples"] == summary.n_samples
    assert kwargs["n_trees"] == summary.n_trees
    assert json.loads(kwargs["params"]) == {"seed": 42}
    # created_at must be ISO 8601
    datetime.fromisoformat(kwargs["created_at"])


def test_ingest_runs_schema_setup_first(tiny_ts):
    driver, session = _mock_driver()
    ingester = ARGIngester(driver)
    ingester.ingest(tiny_ts, run_id="t01", source="msprime")

    # The first CREATE must be a CONSTRAINT (schema setup), not an
    # ARGRun creation.
    first_create = next(
        c for c in session.tx.calls if c[0].startswith("CREATE")
    )
    assert "CONSTRAINT" in first_create[0]


def test_ingest_emits_one_tree_node_per_tskit_node(tiny_ts):
    driver, session = _mock_driver()
    ARGIngester(driver).ingest(tiny_ts, run_id="t01", source="msprime")

    tree_node_calls = [
        c for c in session.tx.calls
        if c[0].startswith("UNWIND $rows AS r CREATE (:TreeNode")
    ]
    total_rows = sum(len(c[1]["rows"]) for c in tree_node_calls)
    assert total_rows == int(tiny_ts.num_nodes)


def test_ingest_emits_one_parent_of_per_edge(tiny_ts):
    driver, session = _mock_driver()
    ARGIngester(driver).ingest(tiny_ts, run_id="t01", source="msprime")

    parent_of_calls = [
        c for c in session.tx.calls
        if "PARENT_OF" in c[0] and c[0].startswith("UNWIND $rows AS r")
    ]
    total = sum(len(c[1]["rows"]) for c in parent_of_calls)
    assert total == int(tiny_ts.num_edges)

    # Spot-check that runId/start/end are passed correctly
    sample_row = parent_of_calls[0][1]["rows"][0]
    assert sample_row["runId"] == "t01"
    assert isinstance(sample_row["start"], int)
    assert isinstance(sample_row["end"], int)
    assert sample_row["parent_id"].startswith("t01:")
    assert sample_row["child_id"].startswith("t01:")


def test_ingest_emits_one_represents_per_sample(tiny_ts):
    driver, session = _mock_driver()
    ARGIngester(driver).ingest(tiny_ts, run_id="t01", source="msprime")

    rep_calls = [
        c for c in session.tx.calls
        if "REPRESENTS" in c[0] and c[0].startswith("UNWIND $rows AS r")
    ]
    total = sum(len(c[1]["rows"]) for c in rep_calls)
    assert total == int(tiny_ts.num_samples)

    haplotypes_seen = set()
    for c in rep_calls:
        for row in c[1]["rows"]:
            haplotypes_seen.add(row["haplotype"])
            assert row["sampleId"].startswith("sample_")
    assert haplotypes_seen == {0, 1}, "expected both hap 0 and hap 1 present"


def test_ingest_uses_custom_sample_id_map(tiny_ts):
    driver, session = _mock_driver()
    ingester = ARGIngester(driver)
    ingester.ingest(
        tiny_ts,
        run_id="t01",
        source="msprime",
        sample_id_map=lambda i: f"NA{i:04d}",
    )

    rep_rows = [
        row
        for c in session.tx.calls
        if "REPRESENTS" in c[0] and c[0].startswith("UNWIND $rows AS r")
        for row in c[1]["rows"]
    ]
    sample_ids = sorted({row["sampleId"] for row in rep_rows})
    assert sample_ids == ["NA0000", "NA0001"]


def test_ingest_emits_mutated_on_per_mutation(tiny_ts):
    driver, session = _mock_driver()
    ARGIngester(driver).ingest(tiny_ts, run_id="t01", source="msprime")

    mut_calls = [
        c for c in session.tx.calls
        if "MUTATED_ON" in c[0] and c[0].startswith("UNWIND $rows AS r")
    ]
    total = sum(len(c[1]["rows"]) for c in mut_calls)
    assert total == int(tiny_ts.num_mutations)


def test_ingest_runs_with_recombination(small_ts):
    """A recombining tree sequence (multiple trees) ingests cleanly."""
    driver, session = _mock_driver()
    summary = ARGIngester(driver).ingest(
        small_ts, run_id="r01", source="msprime"
    )
    assert summary.n_trees > 1, "fixture must actually recombine"
    # All UNWIND batches stayed under the default batch_size cap.
    for cypher, kwargs in session.tx.calls:
        if cypher.startswith("UNWIND $rows AS r"):
            assert len(kwargs["rows"]) <= 10_000


# ---------------------------------------------------------------------------
# list_runs / delete_run
# ---------------------------------------------------------------------------


def test_list_runs_returns_summaries():
    canned_row = {
        "runId": "r01",
        "source": "msprime",
        "n_samples": 4,
        "sequence_length": 100,
        "n_trees": 1,
        "n_nodes": 7,
        "n_edges": 6,
        "n_mutations": 3,
        "created_at": datetime(2026, 5, 9, 12, 0, 0),
    }
    driver, session = _mock_driver(
        read_handlers={"MATCH (r:ARGRun)": [canned_row]}
    )
    runs = ARGIngester(driver).list_runs()
    assert len(runs) == 1
    assert runs[0].run_id == "r01"
    assert runs[0].n_nodes == 7


def test_delete_run_calls_detach_delete():
    canned = [{"n_nodes": 7}]
    driver, session = _mock_driver(
        read_handlers={"MATCH (n:TreeNode {runId": canned}
    )
    n = ARGIngester(driver).delete_run("r01")
    assert n == 7

    # We should see at least one MATCH ... DETACH DELETE for the ARGRun.
    delete_calls = [c for c in session.tx.calls if "DETACH DELETE" in c[0]]
    assert any("(r:ARGRun" in c[0] for c in delete_calls)


# ---------------------------------------------------------------------------
# Run-id isolation: two ingests don't mix
# ---------------------------------------------------------------------------


def test_two_runs_use_distinct_tree_node_id_prefixes(tiny_ts):
    driver, session = _mock_driver()
    ingester = ARGIngester(driver)
    ingester.ingest(tiny_ts, run_id="A", source="msprime")
    ingester.ingest(tiny_ts, run_id="B", source="msprime")

    tree_node_rows = [
        row
        for c in session.tx.calls
        if c[0].startswith("UNWIND $rows AS r CREATE (:TreeNode")
        for row in c[1]["rows"]
    ]
    prefixes = defaultdict(int)
    for row in tree_node_rows:
        prefixes[row["treeNodeId"].split(":")[0]] += 1
    assert prefixes["A"] == int(tiny_ts.num_nodes)
    assert prefixes["B"] == int(tiny_ts.num_nodes)
