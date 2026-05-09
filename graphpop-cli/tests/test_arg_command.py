"""Tests for ``graphpop arg ingest|list|delete``.

Driver and tskit are mocked so no live Neo4j or real tree-sequence
file is needed.
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, PropertyMock, patch

import pytest
from click.testing import CliRunner

from graphpop_cli.cli import main


@pytest.fixture
def cli_runner():
    return CliRunner()


def _patch_driver():
    """Patch GraphPopContext.driver so no real connection is attempted."""
    return patch(
        "graphpop_cli.cli.GraphPopContext.driver",
        new_callable=PropertyMock,
        return_value=MagicMock(),
    )


# ---------------------------------------------------------------------------
# graphpop arg ingest
# ---------------------------------------------------------------------------


def test_arg_ingest_invokes_ingester(cli_runner, tmp_path: Path):
    fixture = tmp_path / "fake.trees"
    fixture.write_bytes(b"\x00")  # any non-empty file; tskit.load is mocked

    summary = MagicMock(
        run_id="r01",
        n_nodes=7,
        n_edges=6,
        n_mutations=3,
    )

    with _patch_driver(), \
         patch("graphpop_import.arg_importer.ARGIngester") as mock_cls, \
         patch("tskit.load") as mock_load:
        mock_load.return_value = MagicMock(
            num_samples=4, num_trees=1, num_edges=6, num_mutations=3
        )
        mock_ingester = MagicMock()
        mock_ingester.ingest.return_value = summary
        mock_cls.return_value = mock_ingester

        result = cli_runner.invoke(
            main,
            [
                "arg", "ingest", str(fixture),
                "--run-id", "r01",
                "--source", "msprime",
                "--params-json", '{"seed": 42}',
            ],
        )

    assert result.exit_code == 0, result.output + (result.stderr or "")
    mock_load.assert_called_once_with(str(fixture))
    mock_ingester.ingest.assert_called_once()
    kwargs = mock_ingester.ingest.call_args.kwargs
    assert kwargs["run_id"] == "r01"
    assert kwargs["source"] == "msprime"
    assert kwargs["params"] == {"seed": 42}


def test_arg_ingest_rejects_invalid_params_json(cli_runner, tmp_path: Path):
    fixture = tmp_path / "fake.trees"
    fixture.write_bytes(b"\x00")

    with _patch_driver(), \
         patch("graphpop_import.arg_importer.ARGIngester"), \
         patch("tskit.load"):
        result = cli_runner.invoke(
            main,
            [
                "arg", "ingest", str(fixture),
                "--run-id", "r01",
                "--params-json", "{not valid json",
            ],
        )

    assert result.exit_code != 0
    assert "valid JSON" in (result.output + (result.stderr or ""))


def test_arg_ingest_uses_sample_id_map_file(cli_runner, tmp_path: Path):
    fixture = tmp_path / "fake.trees"
    fixture.write_bytes(b"\x00")
    map_file = tmp_path / "map.tsv"
    map_file.write_text("# tskit_index\tsample_id\n0\tNA0001\n1\tNA0002\n")

    with _patch_driver(), \
         patch("graphpop_import.arg_importer.ARGIngester") as mock_cls, \
         patch("tskit.load") as mock_load:
        mock_load.return_value = MagicMock(
            num_samples=4, num_trees=1, num_edges=6, num_mutations=0
        )
        mock_ingester = MagicMock()
        mock_ingester.ingest.return_value = MagicMock(
            run_id="r01", n_nodes=7, n_edges=6, n_mutations=0
        )
        mock_cls.return_value = mock_ingester

        result = cli_runner.invoke(
            main,
            [
                "arg", "ingest", str(fixture),
                "--run-id", "r01",
                "--sample-id-map", str(map_file),
            ],
        )

    assert result.exit_code == 0, result.output + (result.stderr or "")
    sample_id_callable = mock_ingester.ingest.call_args.kwargs["sample_id_map"]
    assert sample_id_callable(0) == "NA0001"
    assert sample_id_callable(1) == "NA0002"


# ---------------------------------------------------------------------------
# graphpop arg list
# ---------------------------------------------------------------------------


def test_arg_list_emits_tsv(cli_runner):
    runs = [
        MagicMock(
            run_id="r01",
            source="msprime",
            n_samples=4,
            sequence_length=100,
            n_trees=1,
            n_nodes=7,
            n_edges=6,
            n_mutations=3,
            created_at=datetime(2026, 5, 9, 12, 0, 0),
        ),
    ]

    with _patch_driver(), \
         patch("graphpop_import.arg_importer.ARGIngester") as mock_cls:
        mock_ingester = MagicMock()
        mock_ingester.list_runs.return_value = runs
        mock_cls.return_value = mock_ingester

        result = cli_runner.invoke(main, ["arg", "list"])

    assert result.exit_code == 0, result.output + (result.stderr or "")
    out = result.output
    assert "r01" in out
    assert "msprime" in out


# ---------------------------------------------------------------------------
# graphpop arg delete
# ---------------------------------------------------------------------------


def test_arg_delete_with_yes_flag_runs_delete_run(cli_runner):
    with _patch_driver(), \
         patch("graphpop_import.arg_importer.ARGIngester") as mock_cls:
        mock_ingester = MagicMock()
        mock_ingester.delete_run.return_value = 7
        mock_cls.return_value = mock_ingester

        result = cli_runner.invoke(main, ["arg", "delete", "r01", "--yes"])

    assert result.exit_code == 0, result.output + (result.stderr or "")
    mock_ingester.delete_run.assert_called_once_with("r01")


def test_arg_delete_aborts_without_yes_flag(cli_runner):
    with _patch_driver(), \
         patch("graphpop_import.arg_importer.ARGIngester") as mock_cls:
        mock_ingester = MagicMock()
        mock_cls.return_value = mock_ingester

        # No --yes and no input → click.confirm aborts (exit 1)
        result = cli_runner.invoke(main, ["arg", "delete", "r01"], input="n\n")

    assert result.exit_code != 0
    mock_ingester.delete_run.assert_not_called()
