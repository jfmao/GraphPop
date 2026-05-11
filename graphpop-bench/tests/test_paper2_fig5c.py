"""Unit tests for Paper 2 Fig 5c driver helpers."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from graphpop_bench.paper2_drivers import fig5c_panels as f5c


# ---------------------------------------------------------------------------
# count_lines
# ---------------------------------------------------------------------------

def test_count_lines_total_and_code(tmp_path):
    """Mix of blanks + comments + code → expected counts."""
    p = tmp_path / "x.sh"
    p.write_text(
        "#!/usr/bin/env bash\n"
        "# this is a comment\n"
        "\n"
        "set -euo pipefail\n"
        "echo hello\n"
        "   # leading-whitespace comment\n"
        "\n"
        "exit 0\n"
    )
    total, code = f5c.count_lines(p)
    assert total == 8
    # set, echo, exit = 3 code lines
    assert code == 3


def test_count_lines_cypher_comments(tmp_path):
    p = tmp_path / "q.cypher"
    p.write_text(
        "// header comment\n"
        "MATCH (v:Variant)\n"
        "\n"
        "// inline comment\n"
        "RETURN v\n"
    )
    total, code = f5c.count_lines(p)
    assert total == 5
    assert code == 2


def test_count_lines_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="file not found"):
        f5c.count_lines(tmp_path / "nope.txt")


# ---------------------------------------------------------------------------
# Pipeline introspection
# ---------------------------------------------------------------------------

def test_introspect_pipeline_reference_returns_seven_stages():
    rows = f5c.introspect_pipeline_reference()
    assert len(rows) == 7
    stages = [r.stage for r in rows]
    assert stages == [s.stage for s in f5c.PIPELINE_STAGES]
    # Every row has positive code lines.
    for r in rows:
        assert r.code_lines > 0


def test_introspect_pipeline_reference_includes_known_files():
    rows = f5c.introspect_pipeline_reference()
    files = {r.file for r in rows}
    expected = {
        "01_qc_and_filter.sh", "02_compute_grm.sh",
        "03_compute_kinship.sh", "04_admixture.sh",
        "05_local_ancestry.sh", "06_annotate_pathway.py",
        "07_join_and_filter.py",
    }
    assert files == expected


def test_introspect_pipeline_reference_pipeline_loc_ge_expected():
    """Sanity: total pipeline code-line count is non-trivial."""
    rows = f5c.introspect_pipeline_reference()
    total_code = sum(r.code_lines for r in rows)
    assert total_code >= 60, (
        f"pipeline LOC unexpectedly low: {total_code}")


# ---------------------------------------------------------------------------
# Cypher template row
# ---------------------------------------------------------------------------

def test_cypher_template_row_against_committed(tmp_path):
    """If fig5a_cypher_template.txt is missing, raise clearly."""
    with pytest.raises(FileNotFoundError, match="file not found"):
        f5c.cypher_template_row(tmp_path / "nope.txt")


def test_cypher_template_row_with_synthetic_file(tmp_path):
    p = tmp_path / "q.cypher"
    p.write_text(
        "// header\n"
        "CALL graphpop.kinship.branch_grm()\n"
        "YIELD x\n"
        "RETURN x\n"
    )
    row = f5c.cypher_template_row(p)
    assert row.source == "graphpop"
    assert row.code_lines == 3
    assert row.runtime_s > 0


# ---------------------------------------------------------------------------
# CSV + summary emit
# ---------------------------------------------------------------------------

def test_write_panel_csv_schema(tmp_path):
    rows = f5c.introspect_pipeline_reference()
    csv_path = tmp_path / "panel.csv"
    f5c.write_panel_csv(rows, csv_path)
    with open(csv_path) as fh:
        reader = csv.DictReader(fh)
        out = list(reader)
    assert len(out) == 7
    assert set(out[0].keys()) >= {
        "source", "stage", "tool", "file",
        "total_lines", "code_lines",
        "runtime_s", "runtime_source",
    }


def test_write_summary_json_aggregates(tmp_path):
    """Summary JSON sums LOC + runtime + file count per source."""
    pipeline_rows = f5c.introspect_pipeline_reference()
    fake_cypher = f5c.PanelRow(
        source="graphpop", stage="composed_query",
        tool="GraphPop Cypher", file="x.cypher",
        total_lines=20, code_lines=15,
        runtime_s=30, runtime_source="estimate")
    rows = pipeline_rows + [fake_cypher]
    summary_path = tmp_path / "s.json"
    f5c.write_summary_json(rows, summary_path)
    summary = json.loads(summary_path.read_text())
    assert summary["pipeline"]["n_files"] == 7
    assert summary["graphpop"]["n_files"] == 1
    assert summary["graphpop"]["code_lines"] == 15
    assert summary["graphpop"]["runtime_s"] == 30
    # Pipeline runtime sums the literature estimates.
    expected_runtime = sum(
        s.runtime_s for s in f5c.PIPELINE_STAGES)
    assert summary["pipeline"]["runtime_s"] == expected_runtime
