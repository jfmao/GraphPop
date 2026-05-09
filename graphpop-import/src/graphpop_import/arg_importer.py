"""Online Cypher ingest of tskit TreeSequences into the GraphPop ARG layer.

Public API: :class:`ARGIngester`. See ``M4.A`` plan in
``/home/jfmao/.claude/plans/`` and section 13 of
``docs/GraphPop_Compiled_Holistic_Design.md`` for the schema this writer
emits into.

This is the first online-Cypher importer in GraphPop. It is *additive*:
the target database already contains :Variant, :Sample, :Population
nodes from a prior bulk import; we attach an ARG layer on top without
disturbing them.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, Callable, Iterable, Iterator, Literal

from . import arg_schema

if TYPE_CHECKING:  # only for type hints; tskit is an optional dep
    import tskit
    from neo4j import Driver, ManagedTransaction


Source = Literal["tsinfer", "tsdate", "singer", "relate", "msprime"]


@dataclass
class ARGRunSummary:
    """Summary of an ingested ARG run (one tskit TreeSequence)."""

    run_id: str
    source: str
    n_samples: int
    sequence_length: int
    n_trees: int
    n_nodes: int
    n_edges: int
    n_mutations: int
    created_at: datetime


class ARGIngester:
    """Ingest tskit TreeSequences as :ARGRun + :TreeNode + relationships.

    Parameters
    ----------
    driver:
        An open neo4j ``Driver``. The caller owns the driver's lifecycle.
    database:
        Target database name. Defaults to ``"neo4j"``.
    batch_size:
        Maximum rows per ``UNWIND`` batch. Default 10,000 rows -- a
        sweet spot for online Cypher with the project's 4 GB heap.
    """

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

    # ---- public API ----------------------------------------------------

    def ingest(
        self,
        treeseq: "tskit.TreeSequence",
        *,
        run_id: str,
        source: Source,
        params: dict | None = None,
        sample_id_map: Callable[[int], str] | None = None,
    ) -> ARGRunSummary:
        """Ingest one tskit TreeSequence as :ARGRun + ARG layer relationships.

        The full sequence -- :ARGRun, :TreeNode batch, :PARENT_OF batch,
        :REPRESENTS batch, :MUTATED_ON batch -- runs as a sequence of
        ``session.execute_write(...)`` calls, one per logical batch.
        Driver-managed retries cover transient failures.

        Parameters
        ----------
        treeseq:
            The tree sequence to ingest. A *single* sample of an ARG
            posterior; for SINGER's M=100 posterior, call ``ingest()``
            M times with M distinct ``run_id`` strings.
        run_id:
            Globally unique identifier for this run. Used as the prefix
            of every ``treeNodeId`` and as the ``runId`` payload on
            every relationship.
        source:
            Which inferer produced the ARG. Surfaced on the :ARGRun
            node so downstream procedures can filter.
        params:
            Free-form parameter dict (e.g. recombination map, mutation
            rate) -- serialised to JSON and stored on :ARGRun.
        sample_id_map:
            Maps tskit's local sample IDs (``0``-indexed integers) to
            existing :Sample.sampleId values. Defaults to the identity
            mapping ``i -> "sample_{i}"``.
        """
        if sample_id_map is None:
            sample_id_map = lambda i: f"sample_{i}"

        summary = self._build_summary(
            treeseq, run_id=run_id, source=source, params=params
        )

        with self._driver.session(database=self._database) as session:
            session.execute_write(arg_schema.create_indices)
            session.execute_write(self._create_run_node, summary, params or {})
            session.execute_write(
                self._ingest_tree_nodes, summary.run_id, treeseq
            )
            session.execute_write(
                self._ingest_parent_of_edges, summary.run_id, treeseq
            )
            session.execute_write(
                self._ingest_represents_edges,
                summary.run_id,
                treeseq,
                sample_id_map,
            )
            session.execute_write(
                self._ingest_mutated_on_edges, summary.run_id, treeseq
            )

        return summary

    def list_runs(self) -> list[ARGRunSummary]:
        """Return a summary of every :ARGRun currently in the database."""
        with self._driver.session(database=self._database) as session:
            records = session.execute_read(self._read_runs)
        return [
            ARGRunSummary(
                run_id=r["runId"],
                source=r["source"],
                n_samples=r["n_samples"],
                sequence_length=r["sequence_length"],
                n_trees=r["n_trees"],
                n_nodes=r["n_nodes"],
                n_edges=r["n_edges"],
                n_mutations=r["n_mutations"],
                created_at=r["created_at"].to_native()
                if hasattr(r["created_at"], "to_native")
                else r["created_at"],
            )
            for r in records
        ]

    def delete_run(self, run_id: str) -> int:
        """Delete every node and relationship belonging to ``run_id``.

        Returns the number of :TreeNode nodes deleted (a proxy for
        "did anything happen"). Other runs in the database are left
        untouched.
        """
        with self._driver.session(database=self._database) as session:
            return session.execute_write(self._delete_run, run_id)

    # ---- internal: summary ---------------------------------------------

    @staticmethod
    def _build_summary(
        treeseq: "tskit.TreeSequence",
        *,
        run_id: str,
        source: str,
        params: dict | None,
    ) -> ARGRunSummary:
        return ARGRunSummary(
            run_id=run_id,
            source=source,
            n_samples=int(treeseq.num_samples),
            sequence_length=int(treeseq.sequence_length),
            n_trees=int(treeseq.num_trees),
            n_nodes=int(treeseq.num_nodes),
            n_edges=int(treeseq.num_edges),
            n_mutations=int(treeseq.num_mutations),
            created_at=datetime.now(timezone.utc),
        )

    # ---- internal: write callbacks (each takes a tx) -------------------

    @staticmethod
    def _create_run_node(
        tx: "ManagedTransaction",
        summary: ARGRunSummary,
        params: dict,
    ) -> None:
        tx.run(
            "CREATE (r:ARGRun {"
            "runId: $runId, source: $source, params: $params, "
            "n_samples: $n_samples, sequence_length: $sequence_length, "
            "n_trees: $n_trees, n_nodes: $n_nodes, n_edges: $n_edges, "
            "n_mutations: $n_mutations, created_at: datetime($created_at)"
            "})",
            runId=summary.run_id,
            source=summary.source,
            params=json.dumps(params, sort_keys=True),
            n_samples=summary.n_samples,
            sequence_length=summary.sequence_length,
            n_trees=summary.n_trees,
            n_nodes=summary.n_nodes,
            n_edges=summary.n_edges,
            n_mutations=summary.n_mutations,
            created_at=summary.created_at.isoformat(),
        )

    @staticmethod
    def _ingest_tree_nodes(
        tx: "ManagedTransaction",
        run_id: str,
        treeseq: "tskit.TreeSequence",
    ) -> None:
        # Avoid importing tskit's NODE_IS_SAMPLE flag at module import time.
        try:
            import tskit  # type: ignore

            sample_flag = tskit.NODE_IS_SAMPLE
        except ImportError:  # pragma: no cover - tskit is an optional dep
            sample_flag = 1

        rows = [
            {
                "treeNodeId": f"{run_id}:{i}",
                "runId": run_id,
                "nodeId": int(i),
                "time": float(node.time),
                "is_sample": bool(node.flags & sample_flag),
                "flags": int(node.flags),
            }
            for i, node in enumerate(treeseq.nodes())
        ]
        for batch in _chunked(rows, 10_000):
            tx.run(
                "UNWIND $rows AS r "
                "CREATE (:TreeNode {"
                "treeNodeId: r.treeNodeId, runId: r.runId, "
                "nodeId: r.nodeId, time: r.time, "
                "is_sample: r.is_sample, flags: r.flags"
                "})",
                rows=batch,
            )

    @staticmethod
    def _ingest_parent_of_edges(
        tx: "ManagedTransaction",
        run_id: str,
        treeseq: "tskit.TreeSequence",
    ) -> None:
        rows = [
            {
                "parent_id": f"{run_id}:{int(edge.parent)}",
                "child_id": f"{run_id}:{int(edge.child)}",
                "runId": run_id,
                "start": int(edge.left),
                "end": int(edge.right),
            }
            for edge in treeseq.edges()
        ]
        for batch in _chunked(rows, 10_000):
            tx.run(
                "UNWIND $rows AS r "
                "MATCH (parent:TreeNode {treeNodeId: r.parent_id}) "
                "MATCH (child:TreeNode {treeNodeId: r.child_id}) "
                "CREATE (parent)-[:PARENT_OF {"
                "runId: r.runId, start: r.start, end: r.end"
                "}]->(child)",
                rows=batch,
            )

    @staticmethod
    def _ingest_represents_edges(
        tx: "ManagedTransaction",
        run_id: str,
        treeseq: "tskit.TreeSequence",
        sample_id_map: Callable[[int], str],
    ) -> None:
        # tskit ``treeseq.samples()`` returns the local node IDs of every
        # sample-flagged node; haplotype ordering is implicit (typically
        # 2*i and 2*i+1 for diploid sample i).
        samples: list[int] = list(treeseq.samples())
        rows: list[dict[str, Any]] = []
        for sample_index, tskit_node_id in enumerate(samples):
            # Default identity haplotype ordering: even -> hap 0, odd -> hap 1.
            haplotype = sample_index % 2
            sample_diploid_index = sample_index // 2
            rows.append(
                {
                    "treeNodeId": f"{run_id}:{int(tskit_node_id)}",
                    "sampleId": sample_id_map(sample_diploid_index),
                    "runId": run_id,
                    "haplotype": haplotype,
                }
            )
        for batch in _chunked(rows, 10_000):
            tx.run(
                "UNWIND $rows AS r "
                "MATCH (n:TreeNode {treeNodeId: r.treeNodeId}) "
                "MATCH (s:Sample {sampleId: r.sampleId}) "
                "CREATE (n)-[:REPRESENTS {"
                "runId: r.runId, haplotype: r.haplotype"
                "}]->(s)",
                rows=batch,
            )

    @staticmethod
    def _ingest_mutated_on_edges(
        tx: "ManagedTransaction",
        run_id: str,
        treeseq: "tskit.TreeSequence",
    ) -> None:
        # Build (variantId, child_node_id, parent_node_id, derived_state)
        # tuples. We assume a Variant is uniquely identified by its
        # genomic position (matching the existing GraphPop convention);
        # if multiple alts exist at the same site, the test fixture must
        # encode them as distinct sites.
        rows: list[dict[str, Any]] = []
        edges_by_child: dict[int, list[Any]] = {}
        for edge in treeseq.edges():
            edges_by_child.setdefault(int(edge.child), []).append(edge)

        for site in treeseq.sites():
            position = int(site.position)
            for mutation in site.mutations:
                child = int(mutation.node)
                # Find the edge whose interval covers this site for this
                # child. tskit guarantees uniqueness of the parent at
                # each (child, position) pair.
                parent_node = -1
                for edge in edges_by_child.get(child, ()):
                    if int(edge.left) <= position < int(edge.right):
                        parent_node = int(edge.parent)
                        break
                rows.append(
                    {
                        "position": position,
                        "treeNodeId": f"{run_id}:{child}",
                        "runId": run_id,
                        "parent_node_id": parent_node,
                        "derived_state": str(mutation.derived_state),
                    }
                )
        for batch in _chunked(rows, 10_000):
            tx.run(
                "UNWIND $rows AS r "
                "MATCH (v:Variant) WHERE v.pos = r.position "
                "MATCH (n:TreeNode {treeNodeId: r.treeNodeId}) "
                "CREATE (v)-[:MUTATED_ON {"
                "runId: r.runId, parent_node_id: r.parent_node_id, "
                "derived_state: r.derived_state"
                "}]->(n)",
                rows=batch,
            )

    @staticmethod
    def _read_runs(tx: "ManagedTransaction") -> list[dict[str, Any]]:
        result = tx.run(
            "MATCH (r:ARGRun) RETURN r.runId AS runId, r.source AS source, "
            "r.n_samples AS n_samples, r.sequence_length AS sequence_length, "
            "r.n_trees AS n_trees, r.n_nodes AS n_nodes, "
            "r.n_edges AS n_edges, r.n_mutations AS n_mutations, "
            "r.created_at AS created_at "
            "ORDER BY r.created_at"
        )
        return [dict(record) for record in result]

    @staticmethod
    def _delete_run(tx: "ManagedTransaction", run_id: str) -> int:
        # Three explicit-tx steps. CALL { ... } IN TRANSACTIONS requires an
        # implicit (auto-commit) transaction, so we keep delete in a single
        # explicit tx -- fine for 1000G-class scale (a few million edges).
        result = tx.run(
            "MATCH (n:TreeNode {runId: $runId}) RETURN count(n) AS n_nodes",
            runId=run_id,
        )
        n_nodes = int(result.single()["n_nodes"])
        tx.run(
            "MATCH (n:TreeNode {runId: $runId}) DETACH DELETE n",
            runId=run_id,
        )
        tx.run("MATCH (r:ARGRun {runId: $runId}) DETACH DELETE r", runId=run_id)
        return n_nodes


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _chunked(iterable: Iterable[Any], size: int) -> Iterator[list[Any]]:
    """Yield consecutive ``size``-element chunks from ``iterable``."""
    batch: list[Any] = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch
