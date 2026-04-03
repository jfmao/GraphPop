"""Tests for graphpop_cli.config."""
from __future__ import annotations

from graphpop_cli.config import build_cypher, build_options_map


# ---------------------------------------------------------------------------
# build_options_map
# ---------------------------------------------------------------------------

def test_build_options_map_empty():
    assert build_options_map() == {}


def test_build_options_map_filters():
    opts = build_options_map(consequence="missense", min_af=0.01)
    assert opts == {"consequence": "missense", "min_af": 0.01}


def test_build_options_map_none_excluded():
    opts = build_options_map(pathway=None)
    assert "pathway" not in opts


# ---------------------------------------------------------------------------
# build_cypher
# ---------------------------------------------------------------------------

def test_build_cypher_simple():
    cypher = build_cypher("graphpop.diversity",
                          ["'chr1'", "1", "100", "'EUR'"])
    assert cypher == "CALL graphpop.diversity('chr1', 1, 100, 'EUR')"


def test_build_cypher_with_options():
    cypher = build_cypher("graphpop.diversity",
                          ["'chr1'"],
                          options={"gene": "KCNE1"})
    assert "gene: 'KCNE1'" in cypher
    assert cypher.startswith("CALL graphpop.diversity('chr1', {")


def test_build_cypher_with_yield():
    cypher = build_cypher("graphpop.diversity",
                          ["'chr1'"],
                          yield_cols=["pi", "theta_w"])
    assert cypher.endswith("YIELD pi, theta_w")


def test_build_cypher_bool_option():
    cypher = build_cypher("graphpop.diversity",
                          ["'chr1'"],
                          options={"persist": True})
    assert "persist: true" in cypher


def test_build_cypher_list_option():
    cypher = build_cypher("graphpop.diversity",
                          ["'chr1'"],
                          options={"samples": ["S1", "S2"]})
    assert "samples: ['S1', 'S2']" in cypher
