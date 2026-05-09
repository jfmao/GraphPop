"""Schema declarations for the GraphPop ARG layer (M4.A).

Constraints and indices that must exist in Neo4j before any ARG ingest
runs. All statements use ``IF NOT EXISTS`` and are safe to re-run.

Mirrors the existing pattern in
``data/processed/rice_csv/post_import_indexes.cypher`` and
``scripts/post-deploy-setup.sh``.

Schema (matches docs/GraphPop_Compiled_Holistic_Design.md s 13):

* ``(:ARGRun {runId})`` — one node per inferred ARG (e.g. SINGER posterior
  samples each get their own run).
* ``(:TreeNode {treeNodeId, runId, nodeId, time, is_sample, flags})``.
* ``(:TreeNode)-[:PARENT_OF {runId, start, end}]->(:TreeNode)``.
* ``(:TreeNode)-[:REPRESENTS {runId, haplotype}]->(:Sample)``.
* ``(:Variant)-[:MUTATED_ON {runId, parent_node_id, derived_state}]
  ->(:TreeNode)``.
* ``(:TreeNode)-[:HAS_ANCESTRY {runId, posterior_prob, painter}]
  ->(:Population)`` -- schema reserved; populated by a separate
  post-ingest tool.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from neo4j import ManagedTransaction


# ---------------------------------------------------------------------------
# Schema statements
# ---------------------------------------------------------------------------

CONSTRAINTS: tuple[str, ...] = (
    "CREATE CONSTRAINT arg_run_id IF NOT EXISTS "
    "FOR (r:ARGRun) REQUIRE r.runId IS UNIQUE",
    "CREATE CONSTRAINT tree_node_id IF NOT EXISTS "
    "FOR (n:TreeNode) REQUIRE n.treeNodeId IS UNIQUE",
)

NODE_INDICES: tuple[str, ...] = (
    "CREATE INDEX tree_node_run IF NOT EXISTS "
    "FOR (n:TreeNode) ON (n.runId)",
    "CREATE INDEX tree_node_run_sample IF NOT EXISTS "
    "FOR (n:TreeNode) ON (n.runId, n.is_sample)",
)

# Range indices on the :PARENT_OF interval. Mirrors the gw_range index on
# :GenomicWindow (chr, start, end). Two single-purpose indices on
# (runId, start) and (runId, end) cover the typical traversal: filter by
# runId first, then range-scan by genomic position.
RELATIONSHIP_INDICES: tuple[str, ...] = (
    "CREATE INDEX parent_of_run_start IF NOT EXISTS "
    "FOR ()-[r:PARENT_OF]-() ON (r.runId, r.start)",
    "CREATE INDEX parent_of_run_end IF NOT EXISTS "
    "FOR ()-[r:PARENT_OF]-() ON (r.runId, r.end)",
    "CREATE INDEX represents_run IF NOT EXISTS "
    "FOR ()-[r:REPRESENTS]-() ON (r.runId)",
    "CREATE INDEX mutated_on_run IF NOT EXISTS "
    "FOR ()-[r:MUTATED_ON]-() ON (r.runId)",
    "CREATE INDEX has_ancestry_run IF NOT EXISTS "
    "FOR ()-[r:HAS_ANCESTRY]-() ON (r.runId)",
)


def all_statements() -> tuple[str, ...]:
    """Return every constraint + index statement in deterministic order.

    The order is constraints first (so duplicate-id ingests fail fast),
    then node indices, then relationship indices. Useful for tests that
    want to assert exactly what gets emitted, and for callers that want
    to bulk-execute the schema setup.
    """
    return CONSTRAINTS + NODE_INDICES + RELATIONSHIP_INDICES


def create_indices(tx: "ManagedTransaction") -> None:
    """Create every constraint and index for the ARG layer.

    Idempotent: every statement uses ``IF NOT EXISTS``. Intended to be
    called inside ``session.execute_write(...)``.
    """
    for stmt in all_statements():
        tx.run(stmt)
