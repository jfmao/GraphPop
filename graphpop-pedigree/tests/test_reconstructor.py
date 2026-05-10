"""Tests for PRIMUS-lite pedigree reconstruction."""
from __future__ import annotations

from graphpop_pedigree import (
    PedExporter,
    PedigreeReconstructor,
    PedigreeRow,
)


def _by_id(rows):
    return {r.sample_id: r for r in rows}


def test_trio_with_sex_metadata():
    # Father (M), mother (F), child (unknown).
    samples = [
        ("F1", "dad", 1),
        ("F1", "mom", 2),
        ("F1", "kid", 0),
    ]
    edges = [
        ("dad", "kid", "parent_child", 1),
        ("mom", "kid", "parent_child", 1),
    ]
    rec = PedigreeReconstructor(samples, edges)
    rows = _by_id(rec.reconstruct_family("F1"))

    assert rows["dad"].father_id == "0"
    assert rows["dad"].mother_id == "0"
    assert rows["mom"].father_id == "0"
    assert rows["kid"].father_id == "dad"
    assert rows["kid"].mother_id == "mom"


def test_full_sibling_propagation():
    # Two parents + two children identified as full siblings.
    samples = [
        ("F1", "dad", 1),
        ("F1", "mom", 2),
        ("F1", "kid_a", 1),
        ("F1", "kid_b", 2),
    ]
    edges = [
        ("dad", "kid_a", "parent_child", 1),
        ("mom", "kid_a", "parent_child", 1),
        ("kid_a", "kid_b", "full_sibling", 1),
    ]
    rec = PedigreeReconstructor(samples, edges)
    rows = _by_id(rec.reconstruct_family("F1"))

    # kid_b inherits dad+mom from kid_a via sib propagation.
    assert rows["kid_b"].father_id == "dad"
    assert rows["kid_b"].mother_id == "mom"


def test_three_generation_lineage():
    # Grandparents → parent → grandchild.
    samples = [
        ("F1", "gp_dad", 1),
        ("F1", "gp_mom", 2),
        ("F1", "parent", 1),
        ("F1", "kid", 0),
        ("F1", "spouse", 2),
    ]
    edges = [
        ("gp_dad", "parent", "parent_child", 1),
        ("gp_mom", "parent", "parent_child", 1),
        ("parent", "kid", "parent_child", 1),
        ("spouse", "kid", "parent_child", 1),
    ]
    rec = PedigreeReconstructor(samples, edges)
    rows = _by_id(rec.reconstruct_family("F1"))

    assert rows["parent"].father_id == "gp_dad"
    assert rows["parent"].mother_id == "gp_mom"
    assert rows["kid"].father_id == "parent"
    assert rows["kid"].mother_id == "spouse"


def test_singletons_with_no_edges():
    samples = [
        ("F1", "alone1", 0),
        ("F1", "alone2", 0),
    ]
    rec = PedigreeReconstructor(samples, [])
    rows = _by_id(rec.reconstruct_family("F1"))
    assert rows["alone1"].father_id == "0"
    assert rows["alone1"].mother_id == "0"
    assert rows["alone2"].father_id == "0"


def test_unknown_family_returns_empty():
    rec = PedigreeReconstructor([], [])
    assert rec.reconstruct_family("missing") == []


def test_reconstruct_all_iterates_every_family():
    samples = [
        ("FA", "a1", 1),
        ("FA", "a2", 0),
        ("FB", "b1", 2),
    ]
    edges = [("a1", "a2", "parent_child", 1)]
    rec = PedigreeReconstructor(samples, edges)
    all_rows = rec.reconstruct_all()
    assert {r.family_id for r in all_rows} == {"FA", "FB"}


def test_unknown_sex_falls_back_to_arbitrary_slot():
    samples = [
        ("F1", "p1", 0),
        ("F1", "p2", 0),
        ("F1", "kid", 0),
    ]
    edges = [
        ("p1", "kid", "parent_child", 1),
        ("p2", "kid", "parent_child", 1),
    ]
    rec = PedigreeReconstructor(samples, edges)
    rows = _by_id(rec.reconstruct_family("F1"))
    # Both parents recorded; lex-first gets father slot.
    assert {rows["kid"].father_id, rows["kid"].mother_id} == {"p1", "p2"}
    assert rows["kid"].father_id == "p1"


def test_ped_exporter_writes_six_column_format(tmp_path):
    rows = [
        PedigreeRow("F1", "dad", "0", "0", 1),
        PedigreeRow("F1", "mom", "0", "0", 2),
        PedigreeRow("F1", "kid", "dad", "mom", 0),
    ]
    out = tmp_path / "out.ped"
    n = PedExporter.write(rows, out)
    assert n == 3
    text = out.read_text()
    lines = text.strip().split("\n")
    assert len(lines) == 3
    fields = lines[2].split("\t")
    assert fields == ["F1", "kid", "dad", "mom", "0", "-9"]


def test_ped_exporter_to_string_round_trip():
    rows = [PedigreeRow("F", "S", "0", "0", 1)]
    text = PedExporter.to_string(rows)
    assert text.strip().split("\t") == ["F", "S", "0", "0", "1", "-9"]
