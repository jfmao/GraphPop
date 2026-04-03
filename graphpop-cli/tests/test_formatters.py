"""Tests for graphpop_cli.formatters."""
from __future__ import annotations

import io
import json

from graphpop_cli.formatters import (
    format_output,
    write_csv,
    write_json,
    write_tsv,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _two_records():
    """Return a pair of records used across multiple tests."""
    return [
        {"variant": "chr1:100:A:T", "count": 42, "af": 0.123456789},
        {"variant": "chr1:200:G:C", "count": 7, "af": 0.00034},
    ]


# ---------------------------------------------------------------------------
# write_tsv
# ---------------------------------------------------------------------------

def test_write_tsv_basic():
    buf = io.StringIO()
    write_tsv(_two_records(), buf)
    lines = buf.getvalue().splitlines()

    # Header line
    assert lines[0] == "variant\tcount\taf"
    # Data lines
    assert lines[1].split("\t") == ["chr1:100:A:T", "42", "0.123457"]
    assert lines[2].split("\t") == ["chr1:200:G:C", "7", "0.00034"]


def test_write_tsv_none_becomes_na():
    buf = io.StringIO()
    write_tsv([{"name": "x", "value": None}], buf)
    data_line = buf.getvalue().splitlines()[1]
    assert "NA" in data_line


def test_write_tsv_list_joined():
    buf = io.StringIO()
    write_tsv([{"pops": ["EUR", "AFR", "EAS"], "n": 3}], buf)
    data_line = buf.getvalue().splitlines()[1]
    assert "EUR,AFR,EAS" in data_line


# ---------------------------------------------------------------------------
# write_csv
# ---------------------------------------------------------------------------

def test_write_csv_basic():
    buf = io.StringIO()
    write_csv(_two_records(), buf)
    lines = buf.getvalue().splitlines()

    assert lines[0] == "variant,count,af"
    # csv.DictWriter will quote fields only when needed
    assert "chr1:100:A:T" in lines[1]
    assert "chr1:200:G:C" in lines[2]


# ---------------------------------------------------------------------------
# write_json
# ---------------------------------------------------------------------------

def test_write_json_basic():
    buf = io.StringIO()
    write_json(_two_records(), buf)
    parsed = json.loads(buf.getvalue())

    assert isinstance(parsed, list)
    assert len(parsed) == 2
    assert parsed[0]["variant"] == "chr1:100:A:T"


def test_write_json_roundtrip():
    records = _two_records()
    buf = io.StringIO()
    write_json(records, buf)
    parsed = json.loads(buf.getvalue())

    assert parsed == records


# ---------------------------------------------------------------------------
# format_output
# ---------------------------------------------------------------------------

def test_format_output_tsv_to_file(tmp_path):
    out_file = tmp_path / "out.tsv"
    records = _two_records()
    format_output(records, str(out_file), fmt="tsv",
                  command="diversity", args={"chrom": "chr1", "pop": "EUR"})

    content = out_file.read_text()
    # Header comment should be present when writing to a file
    assert content.startswith("# graphpop diversity")
    assert "chrom: chr1" in content
    # Data should follow
    assert "variant\tcount\taf" in content


def test_format_output_stdout(capsys):
    records = _two_records()
    format_output(records, None, fmt="tsv", command="diversity",
                  args={"chrom": "chr1"})

    captured = capsys.readouterr().out
    # No header comment when writing to stdout
    assert not captured.startswith("#")
    assert "variant\tcount\taf" in captured
