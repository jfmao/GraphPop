"""Unit tests for Paper 2 Fig 5 driver helpers.

Covers Cypher template emission + ecosystem introspection.
No matplotlib / msprime / Java dep — reads only the codebase.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from graphpop_bench.paper2_drivers import fig5_panels as f5


# ---------------------------------------------------------------------------
# Cypher template
# ---------------------------------------------------------------------------

def test_build_cypher_template_emits_files(tmp_path):
    path = f5.build_cypher_template(tmp_path)
    assert path.exists()
    text = path.read_text()
    # Key tokens.
    for token in [
        "graphpop.kinship.branch_grm_by_ancestry",
        "restrict_to_pathway",
        "ancestry",
        "time_window",
        "graphpop.relate.classify",
    ]:
        assert token in text, f"expected {token!r} in template"
    # Annotation JSON is alongside.
    json_path = tmp_path / "fig5a_cypher_annotations.json"
    assert json_path.exists()
    annot = json.loads(json_path.read_text())
    assert len(annot) >= 5
    for row in annot:
        assert "token" in row and "explanation" in row


# ---------------------------------------------------------------------------
# Java procedure introspection
# ---------------------------------------------------------------------------

def test_find_java_procedures_finds_known_names():
    """Live introspection: known procedure names must appear."""
    procs = f5.find_java_procedures(f5.DEFAULT_JAVA_SRC_ROOT)
    names = {p["name"] for p in procs}
    expected_subset = {
        "graphpop.kinship.branch_grm",
        "graphpop.kinship.branch_grm_posterior",
        "graphpop.kinship.branch_grm_by_ancestry",
        "graphpop.relate.classify",
        "graphpop.kinship.king",
        "graphpop.diversity",
    }
    missing = expected_subset - names
    assert not missing, (
        f"missing procedures in introspection result: {missing}")


def test_find_java_procedures_returns_empty_on_missing_root(tmp_path):
    """Robustness: nonexistent root returns []."""
    assert f5.find_java_procedures(tmp_path / "nope") == []


def test_find_java_procedures_attaches_module():
    procs = f5.find_java_procedures(f5.DEFAULT_JAVA_SRC_ROOT)
    by_name = {p["name"]: p for p in procs}
    # branch_grm lives under pairwise/.
    bg = by_name.get("graphpop.kinship.branch_grm")
    assert bg is not None
    assert bg["module"] == "pairwise"
    # diversity lives at the procedures/ root.
    div = by_name.get("graphpop.diversity")
    assert div is not None


# ---------------------------------------------------------------------------
# Wrapper / driver / CLI listings
# ---------------------------------------------------------------------------

def test_list_competitor_wrappers_has_all_phase1():
    wrappers = f5.list_competitor_wrappers()
    labels = {w["label"] for w in wrappers}
    assert labels == {
        "plink_grm", "king", "tskit_branch_grm", "egrm", "s_ldsc"}


def test_list_paper2_drivers_contains_figs():
    drivers = f5.list_paper2_drivers()
    fig_labels = {d["figure"] for d in drivers}
    # Expect at least fig1de, fig2, fig3ab, fig3c, fig4.
    expected = {"fig1de", "fig2", "fig3ab", "fig3c", "fig4"}
    missing = expected - fig_labels
    assert not missing, f"missing driver figures: {missing}"


def test_list_cli_commands_includes_run_subcommands():
    cmds = f5.list_cli_commands()
    assert "graphpop-bench profile" in cmds
    # All five Phase-1 wrappers appear as `run` subcommands.
    for label in ("plink_grm", "king", "tskit_branch_grm",
                  "egrm", "s_ldsc"):
        assert any(f"run {label}" in c for c in cmds), (
            f"no CLI entry for run {label}")


# ---------------------------------------------------------------------------
# introspect_ecosystem composite
# ---------------------------------------------------------------------------

def test_introspect_ecosystem_returns_all_keys():
    eco = f5.introspect_ecosystem()
    assert set(eco.keys()) == {
        "java_procedures", "competitor_wrappers",
        "paper2_drivers", "cli_commands"}
    assert len(eco["java_procedures"]) > 10
    assert len(eco["competitor_wrappers"]) == 5
    assert len(eco["paper2_drivers"]) >= 5
    assert len(eco["cli_commands"]) >= 10


def test_write_ecosystem_json_roundtrip(tmp_path):
    eco = {"java_procedures": [{"name": "x", "mode": "READ",
                                "file": "x.java", "module": "test"}],
           "competitor_wrappers": [], "paper2_drivers": [],
           "cli_commands": []}
    p = f5.write_ecosystem_json(eco, tmp_path)
    assert p.exists()
    back = json.loads(p.read_text())
    assert back == eco
