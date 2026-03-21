"""Unit tests for MCP server parameter translation (no Neo4j required)."""

import json
import pytest

from graphpop_mcp.server import _build_options


class TestBuildOptions:
    """Test the _build_options helper that constructs Neo4j options maps."""

    def test_empty_when_all_none(self):
        opts = _build_options()
        assert opts == {}

    def test_single_filter(self):
        opts = _build_options(min_af=0.05)
        assert opts == {"min_af": 0.05}

    def test_multiple_filters(self):
        opts = _build_options(min_af=0.05, max_af=0.95, variant_type="SNP")
        assert opts == {"min_af": 0.05, "max_af": 0.95, "variant_type": "SNP"}

    def test_consequence_filter(self):
        opts = _build_options(consequence="missense_variant")
        assert opts == {"consequence": "missense_variant"}

    def test_pathway_filter(self):
        opts = _build_options(pathway="R-HSA-73857")
        assert opts == {"pathway": "R-HSA-73857"}

    def test_samples_list(self):
        opts = _build_options(samples=["HG00096", "HG00097"])
        assert opts == {"samples": ["HG00096", "HG00097"]}

    def test_two_pop_samples(self):
        opts = _build_options(
            samples1=["HG00096"], samples2=["NA12878"],
            pop2="EUR",
        )
        assert opts == {
            "samples1": ["HG00096"],
            "samples2": ["NA12878"],
            "pop2": "EUR",
        }

    def test_pop3_for_pbs(self):
        opts = _build_options(pop3="SAS")
        assert opts == {"pop3": "SAS"}

    def test_extra_dict_merged(self):
        opts = _build_options(min_af=0.01, extra={"run_id": "test123"})
        assert opts == {"min_af": 0.01, "run_id": "test123"}

    def test_none_values_excluded(self):
        opts = _build_options(
            min_af=None, max_af=0.5, consequence=None, variant_type="SNP",
        )
        assert opts == {"max_af": 0.5, "variant_type": "SNP"}


class TestToolUniverseRegistry:
    """Test ToolUniverse wrapper dispatch."""

    def test_registry_has_12_tools(self):
        from tooluniverse import graphpop_tool
        assert len(graphpop_tool.TOOL_REGISTRY) == 12

    def test_unknown_tool_returns_error(self):
        from tooluniverse import graphpop_tool
        result = graphpop_tool.run("GraphPop_Unknown", {})
        assert result["success"] is False
        assert "Unknown tool" in result["result"]

    def test_all_tool_names_match_json(self):
        """Verify all JSON-defined tool names have implementations."""
        import pathlib
        from tooluniverse import graphpop_tool

        json_path = pathlib.Path(__file__).parent.parent / "tooluniverse" / "graphpop_tools.json"
        with open(json_path) as f:
            tools_json = json.load(f)

        json_names = {t["name"] for t in tools_json}
        registry_names = set(graphpop_tool.TOOL_REGISTRY.keys())
        assert json_names == registry_names, (
            f"Mismatch: JSON has {json_names - registry_names}, "
            f"registry has {registry_names - json_names}"
        )


class TestToolUniverseBuildOptions:
    """Test ToolUniverse _build_options helper."""

    def test_extracts_present_keys(self):
        from tooluniverse.graphpop_tool import _build_options
        args = {"chr": "22", "pop": "AFR", "min_af": 0.05, "consequence": "missense_variant"}
        opts = _build_options(args, ["min_af", "consequence", "pathway"])
        assert opts == {"min_af": 0.05, "consequence": "missense_variant"}

    def test_skips_none_keys(self):
        from tooluniverse.graphpop_tool import _build_options
        args = {"chr": "22", "pop": "AFR", "min_af": None}
        opts = _build_options(args, ["min_af", "consequence"])
        assert opts == {}
