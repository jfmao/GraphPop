"""Local-ancestry painting ingest into the GraphPop ARG layer.

Populates the ``:HAS_ANCESTRY`` edges reserved by M4.A:

    (:TreeNode)-[:HAS_ANCESTRY {runId, posterior_prob, painter}]->(:Population)

Two entry-points:

* :meth:`AncestryIngester.ingest` — direct, one row per
  ``(tskit_node_id, population_id, posterior_prob)``.
* :meth:`AncestryIngester.ingest_from_samples` — sample-level
  ancestry plus majority-vote DFS propagation to every TreeNode in
  the run.

Closes the schema-only-in-M4.A gap so that
``graphpop.kinship.branch_grm_by_ancestry`` (M4.1 step 7) returns
non-empty matrices.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, Iterable, Iterator, Mapping, Optional

if TYPE_CHECKING:  # pragma: no cover
    from neo4j import Driver, ManagedTransaction


@dataclass(frozen=True)
class PaintingRow:
    """Single per-TreeNode ancestry assignment."""

    tskit_node_id: int
    population_id: str
    posterior_prob: float = 1.0


@dataclass
class PaintingSummary:
    run_id: str
    painter: str
    n_painted_nodes: int
    n_edges: int
    n_populations: int
    created_at: Optional[datetime] = None


@dataclass
class AncestryIngestSummary:
    run_id: str
    painter: str
    n_painted_nodes: int
    n_edges: int
    n_populations: int
    created_at: datetime = field(
        default_factory=lambda: datetime.now(timezone.utc))


class AncestryIngester:
    """Additive Cypher ingest of ``:HAS_ANCESTRY`` edges."""

    def __init__(
        self,
        driver: "Driver",
        *,
        database: str = "neo4j",
        batch_size: int = 10_000,
    ) -> None:
        self._driver = driver
        self._database = database
        self._batch_size = batch_size

    # ---- public API --------------------------------------------------

    def ingest(
        self,
        run_id: str,
        painting: Iterable[PaintingRow],
        *,
        painter: str = "manual",
        replace: bool = False,
    ) -> AncestryIngestSummary:
        """Direct ingest from per-TreeNode rows."""
        rows = list(painting)
        populations = sorted({row.population_id for row in rows})

        with self._driver.session(database=self._database) as session:
            if replace:
                session.execute_write(
                    self._delete_painter, run_id, painter)
            else:
                session.execute_write(
                    self._fail_if_painter_exists, run_id, painter)
            session.execute_write(
                self._ensure_population_nodes, populations)
            for batch in _chunked(rows, self._batch_size):
                session.execute_write(
                    self._create_has_ancestry_batch,
                    run_id, painter, batch)

        return AncestryIngestSummary(
            run_id=run_id,
            painter=painter,
            n_painted_nodes=len({row.tskit_node_id for row in rows}),
            n_edges=len(rows),
            n_populations=len(populations),
        )

    def ingest_from_samples(
        self,
        run_id: str,
        sample_ancestry: Mapping[str, str],
        *,
        painter: str = "majority_vote",
        replace: bool = False,
    ) -> AncestryIngestSummary:
        """Propagate sample-level ancestry to every TreeNode in the run.

        Algorithm:
        1. Load ARG topology + sample-id mapping via Cypher.
        2. Annotate sample-flagged TreeNodes with their per-sample
           ancestry label.
        3. DFS post-order from leaves: each internal node's painting
           is the argmax of descendant-sample counts; tie → smallest
           label alphabetically. ``posterior_prob`` is the majority
           fraction (count / total descendants).

        Notes:
        * Sample TreeNodes that have no entry in ``sample_ancestry``
          are skipped (and so are their up-tree contributions).
        * For ARGs whose painting differs per genomic interval, this
          single-painting-per-node approximation may be too coarse;
          use :meth:`ingest` directly for per-interval painting.
        """
        with self._driver.session(database=self._database) as session:
            topology = session.execute_read(self._load_topology, run_id)

        rows = _propagate_majority(topology, sample_ancestry)
        return self.ingest(run_id, rows, painter=painter, replace=replace)

    def list_paintings(
        self, run_id: Optional[str] = None
    ) -> list[PaintingSummary]:
        with self._driver.session(database=self._database) as session:
            return session.execute_read(self._list_paintings, run_id)

    def delete_painting(self, run_id: str, painter: str) -> int:
        with self._driver.session(database=self._database) as session:
            return session.execute_write(
                self._delete_painter, run_id, painter)

    # ---- internal callbacks (each takes a tx) ------------------------

    @staticmethod
    def _ensure_population_nodes(
        tx: "ManagedTransaction", populations: list[str]
    ) -> None:
        for pop in populations:
            tx.run(
                "MERGE (p:Population {populationId: $pop})",
                pop=pop)

    @staticmethod
    def _create_has_ancestry_batch(
        tx: "ManagedTransaction", run_id: str, painter: str,
        batch: list[PaintingRow]
    ) -> None:
        rows = [
            {
                "tskit_node_id": row.tskit_node_id,
                "population_id": row.population_id,
                "posterior_prob": float(row.posterior_prob),
            }
            for row in batch
        ]
        tx.run(
            "UNWIND $rows AS r "
            "MATCH (n:TreeNode {treeNodeId: $runId + ':' + toString(r.tskit_node_id)}) "
            "MATCH (p:Population {populationId: r.population_id}) "
            "CREATE (n)-[:HAS_ANCESTRY {"
            "runId: $runId, posterior_prob: r.posterior_prob, painter: $painter"
            "}]->(p)",
            rows=rows, runId=run_id, painter=painter)

    @staticmethod
    def _fail_if_painter_exists(
        tx: "ManagedTransaction", run_id: str, painter: str
    ) -> None:
        result = tx.run(
            "MATCH ()-[r:HAS_ANCESTRY {runId: $runId, painter: $painter}]->() "
            "RETURN count(r) AS c LIMIT 1",
            runId=run_id, painter=painter)
        record = result.single()
        if record is not None and record["c"] > 0:
            raise RuntimeError(
                f"Painting already exists for run_id='{run_id}', "
                f"painter='{painter}'. Pass replace=True to overwrite."
            )

    @staticmethod
    def _delete_painter(
        tx: "ManagedTransaction", run_id: str, painter: str
    ) -> int:
        result = tx.run(
            "MATCH ()-[r:HAS_ANCESTRY {runId: $runId, painter: $painter}]->() "
            "WITH count(r) AS n_deleted "
            "MATCH ()-[r:HAS_ANCESTRY {runId: $runId, painter: $painter}]->() "
            "DELETE r RETURN n_deleted",
            runId=run_id, painter=painter)
        record = result.single()
        return int(record["n_deleted"]) if record is not None else 0

    @staticmethod
    def _list_paintings(
        tx: "ManagedTransaction", run_id: Optional[str]
    ) -> list[PaintingSummary]:
        if run_id is None:
            cypher = (
                "MATCH (n:TreeNode)-[r:HAS_ANCESTRY]->(p:Population) "
                "WITH r.runId AS runId, r.painter AS painter, "
                "count(DISTINCT n) AS n_nodes, count(r) AS n_edges, "
                "count(DISTINCT p) AS n_pops "
                "RETURN runId, painter, n_nodes, n_edges, n_pops "
                "ORDER BY runId, painter"
            )
            params: dict[str, Any] = {}
        else:
            cypher = (
                "MATCH (n:TreeNode {runId: $runId})"
                "-[r:HAS_ANCESTRY {runId: $runId}]->(p:Population) "
                "WITH r.painter AS painter, count(DISTINCT n) AS n_nodes, "
                "count(r) AS n_edges, count(DISTINCT p) AS n_pops "
                "RETURN $runId AS runId, painter, n_nodes, n_edges, n_pops "
                "ORDER BY painter"
            )
            params = {"runId": run_id}

        out: list[PaintingSummary] = []
        for record in tx.run(cypher, **params):
            out.append(PaintingSummary(
                run_id=record["runId"],
                painter=record["painter"],
                n_painted_nodes=int(record["n_nodes"]),
                n_edges=int(record["n_edges"]),
                n_populations=int(record["n_pops"]),
            ))
        return out

    @staticmethod
    def _load_topology(
        tx: "ManagedTransaction", run_id: str
    ) -> "_Topology":
        nodes: list[dict[str, Any]] = []
        for r in tx.run(
            "MATCH (n:TreeNode {runId: $runId}) "
            "RETURN n.nodeId AS nodeId, n.is_sample AS is_sample "
            "ORDER BY n.nodeId",
            runId=run_id,
        ):
            nodes.append({
                "nodeId": int(r["nodeId"]),
                "is_sample": bool(r["is_sample"]),
            })

        edges: list[tuple[int, int]] = []
        for r in tx.run(
            "MATCH (p:TreeNode {runId: $runId})-[rel:PARENT_OF "
            "{runId: $runId}]->(c:TreeNode) "
            "RETURN p.nodeId AS parent, c.nodeId AS child",
            runId=run_id,
        ):
            edges.append((int(r["parent"]), int(r["child"])))

        sample_to_node: dict[str, int] = {}
        for r in tx.run(
            "MATCH (n:TreeNode {runId: $runId})"
            "-[rel:REPRESENTS {runId: $runId}]->(s:Sample) "
            "RETURN s.sampleId AS sampleId, n.nodeId AS nodeId",
            runId=run_id,
        ):
            sample_to_node[str(r["sampleId"])] = int(r["nodeId"])

        return _Topology(
            nodes=nodes, edges=edges, sample_to_node=sample_to_node)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


@dataclass
class _Topology:
    nodes: list[dict[str, Any]]
    edges: list[tuple[int, int]]   # (parent, child) tskit ids
    sample_to_node: dict[str, int]  # sampleId -> tskit nodeId


def _propagate_majority(
    topology: _Topology, sample_ancestry: Mapping[str, str]
) -> list[PaintingRow]:
    """Majority-vote propagation from sample leaves up to every node.

    Multi-parent ARG nodes (different parents across different
    intervals) are handled by computing descendants across all edges:
    a node's descendant-sample set is the union of all reachable
    sample-flagged nodes. For deterministic output we DFS by tskit
    nodeId order.
    """
    # children[parent_id] = list of child_ids
    children: dict[int, list[int]] = {}
    for parent, child in topology.edges:
        children.setdefault(parent, []).append(child)

    is_sample = {n["nodeId"]: n["is_sample"] for n in topology.nodes}

    # Sample tskit nodeId -> ancestry label (None if missing).
    node_to_anc: dict[int, str] = {}
    for sample_id, anc in sample_ancestry.items():
        nid = topology.sample_to_node.get(sample_id)
        if nid is None:
            continue
        node_to_anc[nid] = anc

    # Descendant-ancestry counters per node, computed by post-order DFS.
    desc_counts: dict[int, Counter] = {}

    def dfs(n: int) -> Counter:
        if n in desc_counts:
            return desc_counts[n]
        if is_sample.get(n, False):
            anc = node_to_anc.get(n)
            counter: Counter = Counter()
            if anc is not None:
                counter[anc] = 1
            desc_counts[n] = counter
            return counter
        counter = Counter()
        for c in children.get(n, []):
            counter.update(dfs(c))
        desc_counts[n] = counter
        return counter

    for node in topology.nodes:
        dfs(node["nodeId"])

    # Emit ALL labels with their fractional probabilities — preserves the
    # ancestry partition (Σ_a prob(a) = 1) for the downstream
    # branch_grm_by_ancestry decomposition. Internal nodes with mixed
    # descendants yield multiple :HAS_ANCESTRY edges per node.
    rows: list[PaintingRow] = []
    for node in topology.nodes:
        nid = node["nodeId"]
        counts = desc_counts[nid]
        if not counts:
            continue
        total = sum(counts.values())
        # Deterministic ordering by label name.
        for label in sorted(counts):
            rows.append(PaintingRow(
                tskit_node_id=nid,
                population_id=label,
                posterior_prob=counts[label] / total,
            ))
    return rows


def _chunked(iterable: Iterable[Any], size: int) -> Iterator[list[Any]]:
    batch: list[Any] = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch
