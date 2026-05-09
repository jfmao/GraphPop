"""Round-trip integration tests for graphpop_import.arg_importer.

Skipped by default. Set ``GRAPHPOP_NEO4J_URI`` (and optionally
``GRAPHPOP_NEO4J_USER`` / ``GRAPHPOP_NEO4J_PASSWORD``) to a *test*
Neo4j instance to enable. Tests use a unique run_id prefix and clean
up after themselves, but do not isolate themselves from other writes
to the same database.

These checks are the headline correctness gates from the M4.A plan:

* round-trip on msprime fixture (n=20, 5 kb)
* :PARENT_OF count matches tskit edges
* :TreeNode count matches tskit nodes
* :ARGRun summary matches tskit metadata
* TMRCA spot-check against tskit at random positions
* multi-posterior coexistence (two runs do not interfere)
* delete_run cleans up cleanly

Run manually:

::

    GRAPHPOP_NEO4J_URI=bolt://localhost:7687 \\
    GRAPHPOP_NEO4J_USER=neo4j \\
    GRAPHPOP_NEO4J_PASSWORD=graphpop \\
    pytest graphpop-import/tests/test_arg_importer_integration.py -v
"""
from __future__ import annotations

import os
import uuid

import pytest

msprime = pytest.importorskip("msprime")
tskit = pytest.importorskip("tskit")
neo4j = pytest.importorskip("neo4j")

from graphpop_import.arg_importer import ARGIngester


pytestmark = pytest.mark.skipif(
    "GRAPHPOP_NEO4J_URI" not in os.environ,
    reason="set GRAPHPOP_NEO4J_URI to enable Neo4j integration tests",
)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def driver():
    from neo4j import GraphDatabase

    uri = os.environ["GRAPHPOP_NEO4J_URI"]
    user = os.environ.get("GRAPHPOP_NEO4J_USER", "neo4j")
    password = os.environ.get("GRAPHPOP_NEO4J_PASSWORD", "graphpop")
    drv = GraphDatabase.driver(uri, auth=(user, password))
    drv.verify_connectivity()
    yield drv
    drv.close()


@pytest.fixture
def run_id():
    return f"test_{uuid.uuid4().hex[:8]}"


@pytest.fixture
def cleanup_runs(driver):
    """Track runs created during the test and tear them down."""
    created: list[str] = []
    yield created
    if created:
        ingester = ARGIngester(driver)
        for rid in created:
            try:
                ingester.delete_run(rid)
            except Exception:  # pragma: no cover
                pass


@pytest.fixture
def small_ts():
    ts = msprime.sim_ancestry(
        samples=10,
        sequence_length=5_000,
        recombination_rate=1e-4,
        random_seed=43,
    )
    ts = msprime.sim_mutations(ts, rate=1e-5, random_seed=43)
    return ts


# ---------------------------------------------------------------------------
# round-trip
# ---------------------------------------------------------------------------


def test_round_trip_counts_match(driver, run_id, cleanup_runs, small_ts):
    cleanup_runs.append(run_id)
    ingester = ARGIngester(driver)
    summary = ingester.ingest(small_ts, run_id=run_id, source="msprime")

    with driver.session() as s:
        n_tree_nodes = s.run(
            "MATCH (n:TreeNode {runId: $rid}) RETURN count(n) AS c",
            rid=run_id,
        ).single()["c"]
        n_parent_of = s.run(
            "MATCH ()-[r:PARENT_OF {runId: $rid}]->() RETURN count(r) AS c",
            rid=run_id,
        ).single()["c"]
        n_represents = s.run(
            "MATCH (:TreeNode {runId: $rid})-[r:REPRESENTS]->(:Sample) "
            "RETURN count(r) AS c",
            rid=run_id,
        ).single()["c"]

    assert n_tree_nodes == int(small_ts.num_nodes)
    assert n_parent_of == int(small_ts.num_edges)
    # REPRESENTS may be < num_samples if the test DB doesn't contain
    # matching :Sample nodes; the count is bounded by num_samples.
    assert 0 <= n_represents <= int(small_ts.num_samples)
    assert summary.n_nodes == int(small_ts.num_nodes)


def test_arg_run_node_summary(driver, run_id, cleanup_runs, small_ts):
    cleanup_runs.append(run_id)
    ARGIngester(driver).ingest(small_ts, run_id=run_id, source="msprime")
    with driver.session() as s:
        rec = s.run(
            "MATCH (r:ARGRun {runId: $rid}) RETURN r", rid=run_id
        ).single()
    assert rec is not None
    props = dict(rec["r"])
    assert props["source"] == "msprime"
    assert props["n_nodes"] == int(small_ts.num_nodes)
    assert props["n_edges"] == int(small_ts.num_edges)
    assert props["n_mutations"] == int(small_ts.num_mutations)


def test_tmrca_spot_check_matches_tskit(driver, run_id, cleanup_runs,
                                        small_ts):
    """For 5 random pairs at the midpoint, Cypher MRCA matches tskit."""
    import random

    cleanup_runs.append(run_id)
    ARGIngester(driver).ingest(small_ts, run_id=run_id, source="msprime")

    rng = random.Random(0)
    samples = list(small_ts.samples())
    midpoint = small_ts.sequence_length // 2
    tree = small_ts.at(midpoint)
    pairs = [
        tuple(rng.sample(samples, 2))
        for _ in range(5)
    ]

    with driver.session() as s:
        for a, b in pairs:
            mrca_tskit = tree.mrca(a, b)
            tmrca_tskit = small_ts.node(mrca_tskit).time
            # Cypher MRCA: walk up :PARENT_OF from each leaf, collecting
            # ancestors whose interval contains the midpoint, find the
            # lowest-time common ancestor.
            cy = """
            MATCH path_a = (a:TreeNode {treeNodeId: $aid})
                          <-[:PARENT_OF*]-(anc_a:TreeNode {runId: $rid})
            MATCH path_b = (b:TreeNode {treeNodeId: $bid})
                          <-[:PARENT_OF*]-(anc_b:TreeNode {runId: $rid})
            WHERE anc_a = anc_b
              AND ALL(r IN relationships(path_a) WHERE r.start <= $pos AND r.end > $pos)
              AND ALL(r IN relationships(path_b) WHERE r.start <= $pos AND r.end > $pos)
            RETURN min(anc_a.time) AS tmrca
            """
            rec = s.run(
                cy,
                aid=f"{run_id}:{a}",
                bid=f"{run_id}:{b}",
                rid=run_id,
                pos=midpoint,
            ).single()
            tmrca_neo4j = rec["tmrca"]
            assert tmrca_neo4j == pytest.approx(tmrca_tskit, rel=1e-9), (
                f"pair ({a},{b}) at pos {midpoint}: "
                f"tskit={tmrca_tskit}, neo4j={tmrca_neo4j}"
            )


def test_multi_posterior_runs_do_not_interfere(driver, run_id, cleanup_runs,
                                               small_ts):
    """Two ingests with different runIds yield independent counts."""
    rid_a = run_id + "_A"
    rid_b = run_id + "_B"
    cleanup_runs.extend([rid_a, rid_b])
    ingester = ARGIngester(driver)

    ingester.ingest(small_ts, run_id=rid_a, source="msprime")
    ingester.ingest(small_ts, run_id=rid_b, source="msprime")

    with driver.session() as s:
        for rid in (rid_a, rid_b):
            n_nodes = s.run(
                "MATCH (n:TreeNode {runId: $rid}) RETURN count(n) AS c",
                rid=rid,
            ).single()["c"]
            assert n_nodes == int(small_ts.num_nodes)


def test_delete_run_removes_only_target(driver, run_id, cleanup_runs,
                                        small_ts):
    rid_a = run_id + "_A"
    rid_b = run_id + "_B"
    cleanup_runs.extend([rid_a, rid_b])
    ingester = ARGIngester(driver)

    ingester.ingest(small_ts, run_id=rid_a, source="msprime")
    ingester.ingest(small_ts, run_id=rid_b, source="msprime")
    n = ingester.delete_run(rid_a)
    assert n == int(small_ts.num_nodes)

    with driver.session() as s:
        n_a = s.run(
            "MATCH (n:TreeNode {runId: $rid}) RETURN count(n) AS c",
            rid=rid_a,
        ).single()["c"]
        n_b = s.run(
            "MATCH (n:TreeNode {runId: $rid}) RETURN count(n) AS c",
            rid=rid_b,
        ).single()["c"]
    assert n_a == 0, f"{rid_a} should be cleaned up"
    assert n_b == int(small_ts.num_nodes), f"{rid_b} must be untouched"


def test_indexes_present_after_ingest(driver, run_id, cleanup_runs,
                                      small_ts):
    cleanup_runs.append(run_id)
    ARGIngester(driver).ingest(small_ts, run_id=run_id, source="msprime")
    with driver.session() as s:
        names = {r["name"] for r in s.run("SHOW INDEXES YIELD name")}
    for required in ("tree_node_run", "tree_node_run_sample",
                     "parent_of_run_start", "parent_of_run_end",
                     "represents_run", "mutated_on_run", "has_ancestry_run"):
        assert required in names, f"missing index: {required}"
