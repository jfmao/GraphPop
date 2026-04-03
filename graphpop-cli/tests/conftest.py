"""Shared fixtures for GraphPop CLI tests."""
import pytest
from unittest.mock import MagicMock
from click.testing import CliRunner


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_ctx():
    """A mocked GraphPopContext whose run() returns canned records."""
    ctx = MagicMock()
    ctx.run.return_value = []
    ctx.cfg = {"uri": "bolt://localhost:7687", "user": "neo4j",
               "password": "neo4j", "database": "neo4j"}
    ctx.database = "neo4j"
    return ctx
