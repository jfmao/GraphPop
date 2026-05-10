"""Tests for graphpop_import.ibd_ingester (mocked driver)."""
from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from graphpop_import.ibd_ingester import (
    IBDIngester,
    IBDSegmentRow,
    _chunked,
)


# ---- Recording mock helpers ----------------------------------------------


class _RecordingTx:
    def __init__(self, read_handlers=None):
        self.calls = []
        self._read_handlers = read_handlers or []

    def run(self, cypher, **kwargs):
        self.calls.append((cypher, kwargs))
        for prefix, rows in self._read_handlers:
            if cypher.startswith(prefix):
                return _Result(rows)
        return _Result([])


class _Result:
    def __init__(self, rows): self._rows = rows
    def __iter__(self):
        for r in self._rows:
            yield _Record(r)
    def single(self):
        return _Record(self._rows[0]) if self._rows else None


class _Record(dict): pass


class _Session:
    def __init__(self, read_handlers=None):
        self.tx = _RecordingTx(read_handlers)
    def __enter__(self): return self
    def __exit__(self, *args): return False
    def execute_write(self, cb, *a, **k): return cb(self.tx, *a, **k)
    def execute_read(self, cb, *a, **k): return cb(self.tx, *a, **k)


def _mock_driver(read_handlers=None):
    driver = MagicMock()
    driver.session.return_value = _Session(read_handlers)
    return driver, driver.session.return_value


# ---- IBDSegmentRow basics ------------------------------------------------


def test_segment_row_canonical_swaps():
    s = IBDSegmentRow("Z", "A", "chr1", 0, 100)
    c = s.canonical()
    assert c.sample_a == "A" and c.sample_b == "Z"


def test_segment_row_length_bp():
    s = IBDSegmentRow("A", "B", "chr1", 100, 250)
    assert s.length_bp() == 150


# ---- ingest direct path ---------------------------------------------------


def test_ingest_creates_edges_and_returns_summary():
    driver, session = _mock_driver()
    ing = IBDIngester(driver)
    rows = [
        IBDSegmentRow("S0", "S1", "chr1", 0, 30_000_000, length_cM=30.0),
        IBDSegmentRow("S0", "S2", "chr1", 50_000_000, 70_000_000),
    ]
    summary = ing.ingest(rows, source="hap_ibd")

    edge_calls = [c for c in session.tx.calls
                  if c[0].startswith("UNWIND $rows AS r")]
    assert len(edge_calls) == 1
    payload = edge_calls[0][1]["rows"]
    assert len(payload) == 2
    assert payload[0]["length_bp"] == 30_000_000
    assert payload[0]["length_cM"] == 30.0
    assert payload[1]["length_cM"] is None
    assert summary.n_segments == 2
    assert summary.n_pairs == 2


def test_ingest_default_replace_false_raises_on_conflict():
    driver, _ = _mock_driver(read_handlers=[
        ("MATCH ()-[r:IBD_SEGMENT", [{"c": 5}])])
    ing = IBDIngester(driver)
    with pytest.raises(RuntimeError, match="already exist"):
        ing.ingest([IBDSegmentRow("S0", "S1", "chr1", 0, 100)],
                    source="hap_ibd")


def test_ingest_replace_true_calls_delete_first():
    driver, session = _mock_driver(read_handlers=[
        ("MATCH ()-[r:IBD_SEGMENT", [{"n_deleted": 5}])])
    ing = IBDIngester(driver)
    ing.ingest([IBDSegmentRow("S0", "S1", "chr1", 0, 100)],
                source="hap_ibd", replace=True)
    assert any("DELETE r RETURN n_deleted" in c[0]
                for c in session.tx.calls)


def test_ingest_canonicalises_sample_order():
    driver, session = _mock_driver()
    ing = IBDIngester(driver)
    ing.ingest([IBDSegmentRow("Z", "A", "chr1", 0, 100)], source="t")
    edge_calls = [c for c in session.tx.calls
                  if c[0].startswith("UNWIND $rows AS r")]
    payload = edge_calls[0][1]["rows"][0]
    assert payload["sample_a"] == "A"
    assert payload["sample_b"] == "Z"


# ---- TSV parser ----------------------------------------------------------


def test_ingest_tsv_5col(tmp_path):
    p = tmp_path / "ibd.tsv"
    p.write_text("S0\tS1\tchr1\t0\t30000000\nS0\tS2\tchr1\t50000000\t70000000\n")
    driver, session = _mock_driver()
    ing = IBDIngester(driver)
    summary = ing.ingest_tsv(str(p), source="hap_ibd")
    assert summary.n_segments == 2
    edge_rows = [c for c in session.tx.calls
                 if c[0].startswith("UNWIND $rows AS r")][0][1]["rows"]
    # 5-col -> length_cM is None
    for r in edge_rows:
        assert r["length_cM"] is None


def test_ingest_tsv_6col_with_cM(tmp_path):
    p = tmp_path / "ibd.tsv"
    p.write_text("S0\tS1\tchr1\t0\t30000000\t30.5\n")
    driver, session = _mock_driver()
    ing = IBDIngester(driver)
    ing.ingest_tsv(str(p), source="hap_ibd")
    edge_rows = [c for c in session.tx.calls
                 if c[0].startswith("UNWIND $rows AS r")][0][1]["rows"]
    assert edge_rows[0]["length_cM"] == 30.5


def test_ingest_tsv_with_header(tmp_path):
    p = tmp_path / "ibd.tsv"
    p.write_text("sample_a\tsample_b\tchr\tstart\tend\tcM\n"
                  "S0\tS1\tchr1\t0\t100\t1.0\n")
    driver, _ = _mock_driver()
    ing = IBDIngester(driver)
    summary = ing.ingest_tsv(str(p), source="hap_ibd", has_header=True)
    assert summary.n_segments == 1


def test_ingest_tsv_skips_blank_and_comment_lines(tmp_path):
    p = tmp_path / "ibd.tsv"
    p.write_text("# this is a comment\n"
                  "\n"
                  "S0\tS1\tchr1\t0\t100\n")
    driver, _ = _mock_driver()
    ing = IBDIngester(driver)
    summary = ing.ingest_tsv(str(p))
    assert summary.n_segments == 1


def test_ingest_tsv_rejects_short_lines(tmp_path):
    p = tmp_path / "ibd.tsv"
    p.write_text("S0\tS1\tchr1\t100\n")
    driver, _ = _mock_driver()
    ing = IBDIngester(driver)
    with pytest.raises(ValueError, match=">=5 TSV columns|≥5 TSV columns"):
        ing.ingest_tsv(str(p))


# ---- list / delete -------------------------------------------------------


def test_list_sources_returns_per_source_summary():
    canned = [
        {"source": "hap_ibd", "n_edges": 100, "n_pairs": 30},
        {"source": "arg_derived", "n_edges": 642, "n_pairs": 190},
    ]
    driver, _ = _mock_driver(read_handlers=[
        ("MATCH (a:Sample)-[r:IBD_SEGMENT]->", canned)])
    ing = IBDIngester(driver)
    out = ing.list_sources()
    assert len(out) == 2
    assert out[0].source == "hap_ibd"
    assert out[1].n_edges == 642


def test_delete_source_returns_count():
    driver, _ = _mock_driver(read_handlers=[
        ("MATCH ()-[r:IBD_SEGMENT", [{"n_deleted": 7}])])
    ing = IBDIngester(driver)
    assert ing.delete_source("hap_ibd") == 7


# ---- helpers -------------------------------------------------------------


def test_chunked():
    assert list(_chunked(range(7), 3)) == [[0, 1, 2], [3, 4, 5], [6]]
