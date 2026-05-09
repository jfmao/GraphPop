"""Tests for graphpop_import.arg_schema."""
from __future__ import annotations

from unittest.mock import MagicMock

from graphpop_import import arg_schema


def test_all_statements_count_matches_named_groups():
    """Sanity: every named group is included in all_statements()."""
    expected = (
        len(arg_schema.CONSTRAINTS)
        + len(arg_schema.NODE_INDICES)
        + len(arg_schema.RELATIONSHIP_INDICES)
    )
    assert len(arg_schema.all_statements()) == expected


def test_constraints_use_if_not_exists_and_unique():
    for stmt in arg_schema.CONSTRAINTS:
        assert "IF NOT EXISTS" in stmt
        assert "REQUIRE" in stmt and "IS UNIQUE" in stmt


def test_indices_use_if_not_exists():
    for stmt in arg_schema.NODE_INDICES + arg_schema.RELATIONSHIP_INDICES:
        assert "IF NOT EXISTS" in stmt
        assert "CREATE INDEX" in stmt


def test_constraint_targets():
    """The two unique constraints must target ARGRun.runId and TreeNode.treeNodeId."""
    joined = " ".join(arg_schema.CONSTRAINTS)
    assert "(r:ARGRun) REQUIRE r.runId IS UNIQUE" in joined
    assert "(n:TreeNode) REQUIRE n.treeNodeId IS UNIQUE" in joined


def test_node_indices_target_tree_node_run_id():
    joined = " ".join(arg_schema.NODE_INDICES)
    assert "(n:TreeNode) ON (n.runId)" in joined
    assert "(n:TreeNode) ON (n.runId, n.is_sample)" in joined


def test_relationship_indices_cover_all_arg_edges():
    """Every ARG-layer relationship type must have a runId index."""
    joined = " ".join(arg_schema.RELATIONSHIP_INDICES)
    for rel in ("PARENT_OF", "REPRESENTS", "MUTATED_ON", "HAS_ANCESTRY"):
        assert f"[r:{rel}]" in joined, f"no relationship index for :{rel}"


def test_parent_of_uses_composite_run_start_run_end_range():
    """Range traversal needs (runId, start) and (runId, end)."""
    joined = " ".join(arg_schema.RELATIONSHIP_INDICES)
    assert "[r:PARENT_OF]-() ON (r.runId, r.start)" in joined
    assert "[r:PARENT_OF]-() ON (r.runId, r.end)" in joined


def test_create_indices_invokes_tx_run_for_every_statement():
    """create_indices() must call tx.run() once per statement, in order."""
    tx = MagicMock()
    arg_schema.create_indices(tx)

    assert tx.run.call_count == len(arg_schema.all_statements())
    actual = [call.args[0] for call in tx.run.call_args_list]
    assert tuple(actual) == arg_schema.all_statements()


def test_create_indices_is_idempotent_when_replayed():
    """Re-running the function emits the same statements (driver/server
    handles the IF NOT EXISTS at runtime; this guards us against
    accidentally adding state that would mutate between calls)."""
    tx1 = MagicMock()
    tx2 = MagicMock()
    arg_schema.create_indices(tx1)
    arg_schema.create_indices(tx2)
    first = [call.args[0] for call in tx1.run.call_args_list]
    second = [call.args[0] for call in tx2.run.call_args_list]
    assert first == second
