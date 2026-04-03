"""Unit tests for MCP-native tool handler functions (no Neo4j required).

Tests the 9 tool functions that are NOT in the ToolUniverse TOOL_REGISTRY:
  graphpop_status, graphpop_query, graphpop_converge, graphpop_rank_genes,
  graphpop_lookup_gene, graphpop_lookup_pathway, graphpop_lookup_region,
  graphpop_filter, graphpop_inventory.

All Neo4j interactions are mocked via _run_procedure.
"""

import json
from unittest.mock import patch

import pytest

from graphpop_mcp.server import (
    graphpop_converge,
    graphpop_filter,
    graphpop_inventory,
    graphpop_lookup_gene,
    graphpop_lookup_pathway,
    graphpop_lookup_region,
    graphpop_query,
    graphpop_rank_genes,
    graphpop_status,
)

# All tools call _run_procedure, so we patch it consistently.
PATCH_TARGET = "graphpop_mcp.server._run_procedure"


# ---------------------------------------------------------------------------
# graphpop_status
# ---------------------------------------------------------------------------

class TestGraphpopStatus:

    @patch(PATCH_TARGET)
    def test_returns_counts_and_populations(self, mock_run):
        mock_run.side_effect = [
            # First call: counts_cypher
            [{"variants": 1100, "samples": 2504, "populations": 5,
              "genes": 300, "pathways": 50, "chromosomes": 1}],
            # Second call: pops_cypher
            [{"name": "AFR"}, {"name": "EUR"}, {"name": "SAS"}],
        ]
        result = json.loads(graphpop_status())
        assert result["counts"]["variants"] == 1100
        assert result["populations"] == ["AFR", "EUR", "SAS"]
        assert mock_run.call_count == 2

    @patch(PATCH_TARGET)
    def test_empty_database(self, mock_run):
        mock_run.side_effect = [
            [],   # no counts
            [],   # no populations
        ]
        result = json.loads(graphpop_status())
        assert result["counts"] == {}
        assert result["populations"] == []

    @patch(PATCH_TARGET)
    def test_cypher_contains_expected_labels(self, mock_run):
        mock_run.side_effect = [
            [{"variants": 0, "samples": 0, "populations": 0,
              "genes": 0, "pathways": 0, "chromosomes": 0}],
            [],
        ]
        graphpop_status()
        counts_cypher = mock_run.call_args_list[0][0][0]
        for label in ("Variant", "Sample", "Population", "Gene", "Pathway", "Chromosome"):
            assert label in counts_cypher


# ---------------------------------------------------------------------------
# graphpop_query
# ---------------------------------------------------------------------------

class TestGraphpopQuery:

    @patch(PATCH_TARGET)
    def test_basic_query(self, mock_run):
        mock_run.return_value = [{"n": 42}]
        result = json.loads(graphpop_query("RETURN 42 AS n"))
        assert result == [{"n": 42}]
        mock_run.assert_called_once_with("RETURN 42 AS n", None)

    @patch(PATCH_TARGET)
    def test_query_with_params(self, mock_run):
        mock_run.return_value = [{"name": "BRCA1"}]
        params_json = json.dumps({"gene": "BRCA1"})
        result = json.loads(graphpop_query(
            "MATCH (g:Gene {symbol: $gene}) RETURN g.symbol AS name",
            params=params_json,
        ))
        assert result == [{"name": "BRCA1"}]
        call_args = mock_run.call_args
        assert call_args[1] is None or call_args[0][1] == {"gene": "BRCA1"}

    @patch(PATCH_TARGET)
    def test_empty_result(self, mock_run):
        mock_run.return_value = []
        result = json.loads(graphpop_query("MATCH (n:Nothing) RETURN n"))
        assert result == []

    @patch(PATCH_TARGET)
    def test_params_none_by_default(self, mock_run):
        mock_run.return_value = []
        graphpop_query("RETURN 1")
        mock_run.assert_called_once_with("RETURN 1", None)


# ---------------------------------------------------------------------------
# graphpop_converge
# ---------------------------------------------------------------------------

class TestGraphpopConverge:

    @patch(PATCH_TARGET)
    def test_basic_convergence(self, mock_run):
        mock_run.return_value = [
            {"variant": "22:17000000:A:G", "pos": 17000000, "chr": "22", "gene": "KCNE1"},
        ]
        result = json.loads(graphpop_converge(
            pop="AFR", stats="ihs,nsl", thresholds="2.0,2.0", chr="22",
        ))
        assert len(result) == 1
        assert result[0]["gene"] == "KCNE1"

    @patch(PATCH_TARGET)
    def test_cypher_contains_abs_filter(self, mock_run):
        mock_run.return_value = []
        graphpop_converge(pop="EUR", stats="ihs", thresholds="3.0")
        cypher = mock_run.call_args[0][0]
        assert "abs(v.ihs_EUR) >= 3.0" in cypher

    @patch(PATCH_TARGET)
    def test_xpehh_with_pop2(self, mock_run):
        mock_run.return_value = []
        graphpop_converge(pop="AFR", stats="xpehh", thresholds="2.5", pop2="EUR")
        cypher = mock_run.call_args[0][0]
        assert "xpehh_AFR_EUR" in cypher

    @patch(PATCH_TARGET)
    def test_chr_filter_in_cypher(self, mock_run):
        mock_run.return_value = []
        graphpop_converge(pop="AFR", stats="ihs", thresholds="2.0", chr="22")
        cypher = mock_run.call_args[0][0]
        assert "v.chr = '22'" in cypher

    @patch(PATCH_TARGET)
    def test_limit_in_cypher(self, mock_run):
        mock_run.return_value = []
        graphpop_converge(pop="AFR", stats="ihs", thresholds="2.0", limit=50)
        cypher = mock_run.call_args[0][0]
        assert "LIMIT 50" in cypher

    def test_no_variant_stats_returns_error(self):
        # Window-only stats (h12, fst) don't generate variant-level WHERE clauses
        result = json.loads(graphpop_converge(
            pop="AFR", stats="h12,fst", thresholds="0.3,0.1",
        ))
        assert "error" in result

    @patch(PATCH_TARGET)
    def test_empty_result(self, mock_run):
        mock_run.return_value = []
        result = json.loads(graphpop_converge(
            pop="AFR", stats="ihs", thresholds="99.0",
        ))
        assert result == []


# ---------------------------------------------------------------------------
# graphpop_rank_genes
# ---------------------------------------------------------------------------

class TestGraphpopRankGenes:

    @patch(PATCH_TARGET)
    def test_basic_ranking(self, mock_run):
        mock_run.return_value = [
            {"gene": "KCNE1", "chr": "22", "max_ihs": 4.2, "n_high": 2, "n_variants": 15},
            {"gene": "LDLR", "chr": "22", "max_ihs": 3.8, "n_high": 1, "n_variants": 10},
        ]
        result = json.loads(graphpop_rank_genes(pop="AFR"))
        assert len(result) == 2
        assert result[0]["gene"] == "KCNE1"

    @patch(PATCH_TARGET)
    def test_cypher_uses_pop_property(self, mock_run):
        mock_run.return_value = []
        graphpop_rank_genes(pop="EUR")
        cypher = mock_run.call_args[0][0]
        assert "ihs_EUR" in cypher

    @patch(PATCH_TARGET)
    def test_chr_filter(self, mock_run):
        mock_run.return_value = []
        graphpop_rank_genes(pop="AFR", chr="1")
        cypher = mock_run.call_args[0][0]
        assert "v.chr = '1'" in cypher

    @patch(PATCH_TARGET)
    def test_top_limit(self, mock_run):
        mock_run.return_value = []
        graphpop_rank_genes(pop="AFR", top=10)
        cypher = mock_run.call_args[0][0]
        assert "LIMIT 10" in cypher

    @patch(PATCH_TARGET)
    def test_empty_result(self, mock_run):
        mock_run.return_value = []
        result = json.loads(graphpop_rank_genes(pop="AFR"))
        assert result == []


# ---------------------------------------------------------------------------
# graphpop_lookup_gene
# ---------------------------------------------------------------------------

class TestGraphpopLookupGene:

    @patch(PATCH_TARGET)
    def test_found_gene(self, mock_run):
        mock_run.return_value = [{
            "symbol": "KCNE1", "geneId": "ENSG00000180509",
            "chr": "22", "start": 37000000, "end": 37050000,
            "n_variants": 45, "impacts": ["HIGH", "MODERATE"],
            "consequences": ["missense_variant", "synonymous_variant"],
            "pathways": ["Cardiac conduction"],
        }]
        result = json.loads(graphpop_lookup_gene("KCNE1"))
        assert result[0]["symbol"] == "KCNE1"
        assert result[0]["n_variants"] == 45

    @patch(PATCH_TARGET)
    def test_gene_not_found(self, mock_run):
        mock_run.return_value = []
        result = json.loads(graphpop_lookup_gene("FAKEGENE"))
        assert result == []

    @patch(PATCH_TARGET)
    def test_cypher_checks_symbol_and_id(self, mock_run):
        mock_run.return_value = []
        graphpop_lookup_gene("BRCA1")
        cypher = mock_run.call_args[0][0]
        assert "g.symbol = 'BRCA1'" in cypher
        assert "g.geneId = 'BRCA1'" in cypher


# ---------------------------------------------------------------------------
# graphpop_lookup_pathway
# ---------------------------------------------------------------------------

class TestGraphpopLookupPathway:

    @patch(PATCH_TARGET)
    def test_found_pathway(self, mock_run):
        mock_run.return_value = [{
            "pathway": "Signal Transduction",
            "pathwayId": "R-HSA-162582",
            "genes": [
                {"gene": "KCNE1", "n_variants": 10},
                {"gene": "KCNQ1", "n_variants": 20},
            ],
        }]
        result = json.loads(graphpop_lookup_pathway("Signal"))
        assert len(result) == 1
        assert result[0]["pathway"] == "Signal Transduction"

    @patch(PATCH_TARGET)
    def test_pathway_not_found(self, mock_run):
        mock_run.return_value = []
        result = json.loads(graphpop_lookup_pathway("NonexistentPathway"))
        assert result == []

    @patch(PATCH_TARGET)
    def test_cypher_uses_case_insensitive_contains(self, mock_run):
        mock_run.return_value = []
        graphpop_lookup_pathway("cardiac")
        cypher = mock_run.call_args[0][0]
        assert "toLower" in cypher
        assert "CONTAINS" in cypher


# ---------------------------------------------------------------------------
# graphpop_lookup_region
# ---------------------------------------------------------------------------

class TestGraphpopLookupRegion:

    @patch(PATCH_TARGET)
    def test_basic_lookup(self, mock_run):
        mock_run.return_value = [
            {"gene": "KCNE1", "chr": "22", "start": 37000000, "end": 37050000, "n_variants": 15},
        ]
        result = json.loads(graphpop_lookup_region("22", 37000000, 37050000))
        assert len(result) == 1
        assert result[0]["gene"] == "KCNE1"

    @patch(PATCH_TARGET)
    def test_empty_region(self, mock_run):
        mock_run.return_value = []
        result = json.loads(graphpop_lookup_region("22", 1, 100))
        assert result == []

    @patch(PATCH_TARGET)
    def test_cypher_contains_region_bounds(self, mock_run):
        mock_run.return_value = []
        graphpop_lookup_region("22", 100000, 200000)
        cypher = mock_run.call_args[0][0]
        assert "v.chr = '22'" in cypher
        assert "v.pos >= 100000" in cypher
        assert "v.pos <= 200000" in cypher


# ---------------------------------------------------------------------------
# graphpop_filter
# ---------------------------------------------------------------------------

class TestGraphpopFilter:

    @patch(PATCH_TARGET)
    def test_basic_ihs_filter(self, mock_run):
        mock_run.return_value = [
            {"variant": "22:17000000:A:G", "pos": 17000000, "ihs": 3.5},
        ]
        result = json.loads(graphpop_filter(
            statistic="ihs", chr="22", pop="AFR", min_score=2.0,
        ))
        assert len(result) == 1
        assert result[0]["ihs"] == 3.5

    @patch(PATCH_TARGET)
    def test_ihs_cypher_property_name(self, mock_run):
        mock_run.return_value = []
        graphpop_filter(statistic="ihs", chr="22", pop="EUR")
        cypher = mock_run.call_args[0][0]
        assert "ihs_EUR" in cypher

    @patch(PATCH_TARGET)
    def test_xpehh_with_pop2(self, mock_run):
        mock_run.return_value = []
        graphpop_filter(statistic="xpehh", chr="22", pop="AFR", pop2="EUR")
        cypher = mock_run.call_args[0][0]
        assert "xpehh_AFR_EUR" in cypher

    @patch(PATCH_TARGET)
    def test_consequence_filter(self, mock_run):
        mock_run.return_value = []
        graphpop_filter(
            statistic="ihs", chr="22", pop="AFR", consequence="missense_variant",
        )
        cypher = mock_run.call_args[0][0]
        assert "HAS_CONSEQUENCE" in cypher
        assert "missense_variant" in cypher

    @patch(PATCH_TARGET)
    def test_pathway_filter(self, mock_run):
        mock_run.return_value = []
        graphpop_filter(
            statistic="nsl", chr="22", pop="AFR", pathway="Signal",
        )
        cypher = mock_run.call_args[0][0]
        assert "IN_PATHWAY" in cypher
        assert "Signal" in cypher

    @patch(PATCH_TARGET)
    def test_gene_filter(self, mock_run):
        mock_run.return_value = []
        graphpop_filter(
            statistic="ihs", chr="22", pop="AFR", gene="KCNE1",
        )
        cypher = mock_run.call_args[0][0]
        assert "g.symbol = 'KCNE1'" in cypher

    @patch(PATCH_TARGET)
    def test_min_score_filter(self, mock_run):
        mock_run.return_value = []
        graphpop_filter(
            statistic="ihs", chr="22", pop="AFR", min_score=3.0,
        )
        cypher = mock_run.call_args[0][0]
        assert "abs(v.ihs_AFR) >= 3.0" in cypher

    @patch(PATCH_TARGET)
    def test_limit(self, mock_run):
        mock_run.return_value = []
        graphpop_filter(
            statistic="ihs", chr="22", pop="AFR", limit=500,
        )
        cypher = mock_run.call_args[0][0]
        assert "LIMIT 500" in cypher

    @patch(PATCH_TARGET)
    def test_empty_result(self, mock_run):
        mock_run.return_value = []
        result = json.loads(graphpop_filter(
            statistic="ihs", chr="22", pop="AFR",
        ))
        assert result == []


# ---------------------------------------------------------------------------
# graphpop_inventory
# ---------------------------------------------------------------------------

class TestGraphpopInventory:

    @patch(PATCH_TARGET)
    def test_full_inventory(self, mock_run):
        # The function makes 8 calls:
        # 1. populations, 2. chromosomes, 3-8. node counts for 6 labels, 9. props
        mock_run.side_effect = [
            # populations
            [{"pop": "AFR", "n": 661}, {"pop": "EUR", "n": 503}],
            # chromosomes
            [{"chr": "22", "length": 50818468}],
            # 6 node count calls: Variant, Sample, Gene, Pathway, GOTerm, GenomicWindow
            [{"cnt": 1100000}],
            [{"cnt": 2504}],
            [{"cnt": 500}],
            [{"cnt": 100}],
            [{"cnt": 200}],
            [{"cnt": 5000}],
            # props from sample variant
            [{"props": [
                "variantId", "chr", "pos", "ref", "alt",
                "ihs_AFR", "ihs_EUR", "xpehh_AFR_EUR", "nsl_AFR",
                "ancestral_allele",
            ]}],
        ]
        result = json.loads(graphpop_inventory())
        assert len(result["populations"]) == 2
        assert result["n_variant"] == 1100000
        assert result["n_sample"] == 2504
        assert "ihs_AFR" in result["persisted_ihs"]
        assert "xpehh_AFR_EUR" in result["persisted_xpehh"]
        assert result["has_ancestral_allele"] is True

    @patch(PATCH_TARGET)
    def test_empty_database_inventory(self, mock_run):
        mock_run.side_effect = [
            [],   # populations
            [],   # chromosomes
            # 6 node counts all empty
            [], [], [], [], [], [],
            # no props (empty db)
            [],
        ]
        result = json.loads(graphpop_inventory())
        assert result["populations"] == []
        assert result["chromosomes"] == []
        assert result["n_variant"] == 0

    @patch(PATCH_TARGET)
    def test_inventory_checks_all_labels(self, mock_run):
        """Verify that inventory queries all 6 expected node labels."""
        mock_run.side_effect = [
            [],  # populations
            [],  # chromosomes
            [{"cnt": 0}], [{"cnt": 0}], [{"cnt": 0}],
            [{"cnt": 0}], [{"cnt": 0}], [{"cnt": 0}],
            [],  # props
        ]
        graphpop_inventory()
        # Calls: populations, chromosomes, 6 labels, props = 9 total
        assert mock_run.call_count == 9
        all_cypher = " ".join(call[0][0] for call in mock_run.call_args_list)
        for label in ["Variant", "Sample", "Gene", "Pathway", "GOTerm", "GenomicWindow"]:
            assert label in all_cypher

    @patch(PATCH_TARGET)
    def test_no_persisted_stats(self, mock_run):
        """When variant has no selection stat properties."""
        mock_run.side_effect = [
            [{"pop": "AFR", "n": 100}],  # populations
            [],  # chromosomes
            [{"cnt": 100}], [{"cnt": 50}], [{"cnt": 10}],
            [{"cnt": 5}], [{"cnt": 0}], [{"cnt": 0}],
            [{"props": ["variantId", "chr", "pos", "ref", "alt"]}],
        ]
        result = json.loads(graphpop_inventory())
        assert result["persisted_ihs"] == []
        assert result["persisted_xpehh"] == []
        assert result["persisted_nsl"] == []
        assert result["has_ancestral_allele"] is False
