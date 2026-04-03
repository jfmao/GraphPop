"""Tests for GraphPop CLI commands.

Each test invokes a CLI command via Click's CliRunner with a mocked
GraphPopContext.run() method, so no live Neo4j connection is needed.
We verify:
  1. exit_code == 0
  2. mock_run was called the expected number of times
  3. The generated Cypher contains expected patterns
  4. Parameterized queries pass the correct parameter keys/values
"""
from __future__ import annotations

import json

import pytest
from unittest.mock import patch, MagicMock, PropertyMock

from click.testing import CliRunner

from graphpop_cli.cli import main


@pytest.fixture
def cli_runner():
    return CliRunner()


def _invoke(runner, args, mock_records):
    """Invoke a CLI command with mocked Neo4j results.

    Patches GraphPopContext.run to return *mock_records* and prevents
    any real Neo4j driver from being created.
    """
    with patch('graphpop_cli.cli.GraphPopContext.run',
               return_value=mock_records) as mock_run, \
         patch('graphpop_cli.cli.GraphPopContext.driver',
               new_callable=PropertyMock,
               return_value=MagicMock()):
        result = runner.invoke(main, args)
    return result, mock_run


# ---------------------------------------------------------------------------
# 1. diversity command
# ---------------------------------------------------------------------------

class TestDiversity:

    MOCK_DIVERSITY = [{
        "pi": 0.00123, "theta_w": 0.00145, "tajima_d": -0.5,
        "fay_wu_h": 0.0, "fay_wu_h_norm": 0.0,
        "het_exp": 0.32, "het_obs": 0.30, "fis": 0.0625,
        "n_variants": 1000, "n_segregating": 800, "n_polarized": 750,
    }]

    def test_diversity_basic(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR"],
            self.MOCK_DIVERSITY,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        assert "CALL graphpop.diversity" in cypher
        assert "'chr22'" in cypher
        assert "'EUR'" in cypher
        assert "YIELD" in cypher

    def test_diversity_with_options(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR",
             "--consequence", "missense_variant", "--gene", "KCNE1"],
            self.MOCK_DIVERSITY,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        assert "CALL graphpop.diversity" in cypher
        assert "consequence" in cypher
        assert "'missense_variant'" in cypher
        assert "gene" in cypher
        assert "'KCNE1'" in cypher

    def test_diversity_empty_result(self, cli_runner):
        """Empty result should still exit 0 (formatter handles empty list)."""
        result, mock_run = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR"],
            [],
        )
        assert result.exit_code == 0

    def test_diversity_json_output(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR", "--format", "json"],
            self.MOCK_DIVERSITY,
        )
        assert result.exit_code == 0
        parsed = json.loads(result.output)
        assert isinstance(parsed, list)
        assert parsed[0]["pi"] == 0.00123

    def test_diversity_csv_output(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR", "--format", "csv"],
            self.MOCK_DIVERSITY,
        )
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        # CSV header + 1 data row
        assert len(lines) == 2
        assert "," in lines[0]


# ---------------------------------------------------------------------------
# 2. sfs command
# ---------------------------------------------------------------------------

class TestSFS:

    MOCK_SFS = [{
        "sfs": [100, 50, 30, 20, 10],
        "n_variants": 1000,
        "max_ac": 10,
        "n_polarized": 0,
    }]

    def test_sfs_folded(self, cli_runner):
        """Default invocation (no --unfolded flag) passes 'false' for unfolded."""
        result, mock_run = _invoke(
            cli_runner,
            ["sfs", "chr22", "1", "50000000", "EUR"],
            self.MOCK_SFS,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        assert "graphpop.sfs" in cypher
        assert "false" in cypher  # unfolded=false by default

    def test_sfs_unfolded(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["sfs", "chr22", "1", "50000000", "EUR", "--unfolded"],
            self.MOCK_SFS,
        )
        assert result.exit_code == 0
        cypher = mock_run.call_args[0][0]
        assert "true" in cypher  # unfolded=true


# ---------------------------------------------------------------------------
# 3. lookup gene subcommand
# ---------------------------------------------------------------------------

class TestLookupGene:

    MOCK_GENE = [{
        "gene": "KCNE1", "gene_id": "ENSG00000180509",
        "chr": "chr21", "start": 35000, "end": 36000,
        "variant_id": "chr21:35500:A:G", "pos": 35500,
        "ref": "A", "alt": "G",
        "pathways": ["Cardiac"],
        "ihs_scores": [], "xpehh_scores": [],
    }]

    def test_lookup_gene(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["lookup", "gene", "KCNE1"],
            self.MOCK_GENE,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        params = mock_run.call_args[0][1]
        assert "$gene_name" in cypher
        assert params["gene_name"] == "KCNE1"

    def test_lookup_gene_not_found(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["lookup", "gene", "FAKE_GENE"],
            [],
        )
        assert result.exit_code == 0
        assert "not found" in result.stderr


# ---------------------------------------------------------------------------
# 4. lookup region subcommand
# ---------------------------------------------------------------------------

class TestLookupRegion:

    MOCK_REGION = [{
        "gene": "KCNE1", "gene_id": "ENSG00000180509",
        "gene_start": 9050000, "gene_end": 9100000,
        "variant_count": 42, "min_pos": 9050000, "max_pos": 9099999,
    }]

    def test_lookup_region(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["lookup", "region", "chr6", "9000000", "9600000"],
            self.MOCK_REGION,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        params = mock_run.call_args[0][1]
        assert "$chr" in cypher
        assert "$start" in cypher
        assert "$end" in cypher
        assert params["chr"] == "chr6"
        assert params["start"] == 9000000
        assert params["end"] == 9600000

    def test_lookup_region_empty(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["lookup", "region", "chrX", "1", "100"],
            [],
        )
        assert result.exit_code == 0
        assert "No variants found" in result.stderr


# ---------------------------------------------------------------------------
# 5. extract variants subcommand
# ---------------------------------------------------------------------------

class TestExtractVariants:

    MOCK_VARIANTS = [{
        "variantId": "chr22:100:A:G", "pos": 100,
        "ref": "A", "alt": "G", "af_EUR": 0.15,
    }]

    def test_extract_variants(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["extract", "variants", "--chr", "chr22",
             "--pop", "EUR", "--gene", "KCNE1"],
            self.MOCK_VARIANTS,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        params = mock_run.call_args[0][1]
        assert "$chromosome" in cypher
        assert params["chromosome"] == "chr22"
        assert "$gene" in cypher
        assert params["gene"] == "KCNE1"

    def test_extract_variants_empty(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["extract", "variants", "--chr", "chr99", "--pop", "ZZZ"],
            [],
        )
        assert result.exit_code == 0
        assert "No variants found" in result.stderr


# ---------------------------------------------------------------------------
# 6. extract samples subcommand
# ---------------------------------------------------------------------------

class TestExtractSamples:

    MOCK_SAMPLES = [{
        "sampleId": "HG00096", "population": "EUR",
        "packed_index": 0, "froh": 0.05,
        "pop_n_samples": 503, "pop_mean_froh": 0.04,
    }]

    def test_extract_samples(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["extract", "samples", "--pop", "EUR"],
            self.MOCK_SAMPLES,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        params = mock_run.call_args[0][1]
        assert "$population" in cypher
        assert params["population"] == "EUR"

    def test_extract_samples_empty(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["extract", "samples", "--pop", "NOPOP"],
            [],
        )
        assert result.exit_code == 0
        assert "No samples found" in result.stderr


# ---------------------------------------------------------------------------
# 7. compare command
# ---------------------------------------------------------------------------

class TestCompare:

    MOCK_COMPARE = [{
        "window_start": 0, "window_end": 100000,
        "pi_EUR": 0.001, "pi_AFR": 0.002,
        "delta": -0.001, "abs_delta": 0.001,
    }]

    def test_compare_pi(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["compare", "EUR", "AFR", "chr22", "--stat", "pi"],
            self.MOCK_COMPARE,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        params = mock_run.call_args[0][1]
        assert "$chr" in cypher
        assert "$pop1" in cypher
        assert "$pop2" in cypher
        assert params["chr"] == "chr22"
        assert params["pop1"] == "EUR"
        assert params["pop2"] == "AFR"
        # Should reference the pi property
        assert ".pi" in cypher

    def test_compare_empty(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["compare", "EUR", "AFR", "chr22", "--stat", "pi"],
            [],
        )
        assert result.exit_code == 0
        assert "No comparison data" in result.stderr


# ---------------------------------------------------------------------------
# 8. neighbors command
# ---------------------------------------------------------------------------

class TestNeighbors:

    MOCK_NEIGHBORS = [{
        "gene": "KCNE2", "shared_pathway": "Cardiac",
        "chr": "chr21", "start": 40000, "end": 41000,
    }]

    def test_neighbors_pathway(self, cli_runner):
        result, mock_run = _invoke(
            cli_runner,
            ["neighbors", "KCNE1"],
            self.MOCK_NEIGHBORS,
        )
        assert result.exit_code == 0, result.output + (result.stderr or "")
        assert mock_run.call_count == 1

        cypher = mock_run.call_args[0][0]
        params = mock_run.call_args[0][1]
        assert "$gene" in cypher
        assert params["gene"] == "KCNE1"
        # Default --via is IN_PATHWAY
        assert "IN_PATHWAY" in cypher

    def test_neighbors_empty(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["neighbors", "FAKE_GENE"],
            [],
        )
        assert result.exit_code == 0
        assert "No neighbors found" in result.stderr

    def test_neighbors_via_ld(self, cli_runner):
        mock_ld = [{"gene": "KCNE2", "max_r2": 0.85, "chr": "chr21"}]
        result, mock_run = _invoke(
            cli_runner,
            ["neighbors", "KCNE1", "--via", "LD"],
            mock_ld,
        )
        assert result.exit_code == 0
        cypher = mock_run.call_args[0][0]
        assert "LD" in cypher
        assert mock_run.call_args[0][1]["gene"] == "KCNE1"


# ---------------------------------------------------------------------------
# 9. Error handling
# ---------------------------------------------------------------------------

class TestErrorHandling:

    def test_diversity_empty_exits_zero(self, cli_runner):
        """Diversity with empty results still exits 0."""
        result, _ = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR"],
            [],
        )
        assert result.exit_code == 0

    def test_neighbors_empty_stderr_message(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["neighbors", "NONEXISTENT"],
            [],
        )
        assert result.exit_code == 0
        assert "No neighbors found" in result.stderr


# ---------------------------------------------------------------------------
# 10. Output formats
# ---------------------------------------------------------------------------

class TestOutputFormats:

    MOCK_DATA = [{
        "pi": 0.00123, "theta_w": 0.00145, "tajima_d": -0.5,
        "fay_wu_h": 0.0, "fay_wu_h_norm": 0.0,
        "het_exp": 0.32, "het_obs": 0.30, "fis": 0.0625,
        "n_variants": 1000, "n_segregating": 800, "n_polarized": 750,
    }]

    def test_diversity_json_output(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR", "--format", "json"],
            self.MOCK_DATA,
        )
        assert result.exit_code == 0
        parsed = json.loads(result.output)
        assert isinstance(parsed, list)
        assert len(parsed) == 1
        assert "pi" in parsed[0]

    def test_diversity_csv_output(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR", "--format", "csv"],
            self.MOCK_DATA,
        )
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert len(lines) >= 2  # header + data
        assert "," in lines[0]
        # All expected column names present in header
        for col in ("pi", "theta_w", "tajima_d"):
            assert col in lines[0]

    def test_diversity_tsv_output(self, cli_runner):
        result, _ = _invoke(
            cli_runner,
            ["diversity", "chr22", "1", "50000000", "EUR", "--format", "tsv"],
            self.MOCK_DATA,
        )
        assert result.exit_code == 0
        lines = result.output.strip().split("\n")
        assert len(lines) >= 2
        assert "\t" in lines[0]
