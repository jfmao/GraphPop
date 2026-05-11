"""Unit tests for the Paper 2 Fig 1d/1e driver pure-logic helpers.

Covers JSON-matrix loading, TSV pair-map loading, matrix→pair-map
flattening, join + rel-err computation, and CSV emission. Has no
Java / Maven / egrm dependency — the integration steps are
manual one-shots.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from graphpop_bench.paper2_drivers import fig1de_panels as driver


# ---------------------------------------------------------------------------
# load_egrm_json_matrix
# ---------------------------------------------------------------------------

def test_load_egrm_json_matrix_roundtrip(tmp_path):
    blob = {
        "n_samples": 3,
        "sample_ids": [0, 1, 2],
        "matrix": [
            [0.5, 0.1, -0.1],
            [0.1, 0.5, 0.2],
            [-0.1, 0.2, 0.5],
        ],
    }
    p = tmp_path / "ref.json"
    p.write_text(json.dumps(blob))
    sample_ids, matrix = driver.load_egrm_json_matrix(p)
    assert sample_ids == [0, 1, 2]
    assert matrix[1][2] == pytest.approx(0.2)


def test_load_egrm_json_matrix_shape_mismatch(tmp_path):
    blob = {
        "n_samples": 3,
        "sample_ids": [0, 1, 2],
        "matrix": [[0.5, 0.1], [0.1, 0.5]],   # 2x2, not 3x3
    }
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(blob))
    with pytest.raises(ValueError, match="shape mismatch"):
        driver.load_egrm_json_matrix(p)


def test_load_egrm_json_matrix_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="egrm JSON reference"):
        driver.load_egrm_json_matrix(tmp_path / "no.json")


# ---------------------------------------------------------------------------
# load_pair_tsv
# ---------------------------------------------------------------------------

def test_load_pair_tsv_canonicalises_pairs(tmp_path):
    tsv = (
        "sample_a\tsample_b\tkinship\n"
        "0\t1\t0.25\n"
        "2\t0\t0.10\n"     # out-of-order pair; key should canonicalise
        "1\t1\t0.50\n"     # diagonal
    )
    p = tmp_path / "x.tsv"
    p.write_text(tsv)
    pairs = driver.load_pair_tsv(p)
    assert pairs[(0, 1)] == pytest.approx(0.25)
    assert pairs[(0, 2)] == pytest.approx(0.10)
    assert pairs[(1, 1)] == pytest.approx(0.50)


def test_load_pair_tsv_missing_kinship_column_raises(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("sample_a\tsample_b\tphi\n0\t1\t0.25\n")
    with pytest.raises(ValueError, match="kinship"):
        driver.load_pair_tsv(p)


def test_load_pair_tsv_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="TSV not found"):
        driver.load_pair_tsv(tmp_path / "no.tsv")


# ---------------------------------------------------------------------------
# matrix_to_pair_map
# ---------------------------------------------------------------------------

def test_matrix_to_pair_map_upper_triangle_with_diag():
    sample_ids = [0, 1, 2]
    m = [
        [0.5, 0.1, -0.1],
        [0.1, 0.5, 0.2],
        [-0.1, 0.2, 0.5],
    ]
    out = driver.matrix_to_pair_map(sample_ids, m)
    assert len(out) == 6   # n(n+1)/2
    assert out[(0, 0)] == pytest.approx(0.5)
    assert out[(0, 1)] == pytest.approx(0.1)
    assert out[(0, 2)] == pytest.approx(-0.1)
    assert out[(1, 2)] == pytest.approx(0.2)


def test_matrix_to_pair_map_with_non_sequential_ids():
    """Sample IDs need not be 0..n-1 (some fixtures use haplotype IDs)."""
    sample_ids = [10, 5, 20]
    m = [
        [0.5, 0.1, -0.1],
        [0.1, 0.5, 0.2],
        [-0.1, 0.2, 0.5],
    ]
    out = driver.matrix_to_pair_map(sample_ids, m)
    # (10, 5) row→col entry should be stored under canonical (5, 10).
    assert out[(5, 10)] == pytest.approx(0.1)
    assert out[(10, 20)] == pytest.approx(-0.1)


# ---------------------------------------------------------------------------
# join_pairs + rel-err
# ---------------------------------------------------------------------------

def test_join_pairs_zero_rel_err_on_identity():
    egrm = {(0, 0): 1.0, (0, 1): 0.5, (1, 1): 1.0}
    rows = driver.join_pairs(egrm, dict(egrm))
    assert len(rows) == 3
    for r in rows:
        assert r.abs_diff == 0.0
        assert r.rel_err == 0.0


def test_join_pairs_known_perturbation_scales_linearly():
    """Perturbing graphpop by +eps relative to a unit-magnitude
    reference produces rel-err ≈ eps."""
    egrm = {(0, 0): 1.0, (0, 1): -2.0, (1, 1): 0.5}
    eps = 1e-5
    gp = {k: v * (1 + eps) for k, v in egrm.items()}
    rows = driver.join_pairs(egrm, gp)
    for r in rows:
        assert r.rel_err == pytest.approx(eps, rel=1e-6)


def test_join_pairs_missing_key_raises():
    egrm = {(0, 0): 1.0, (0, 1): 0.5}
    gp = {(0, 0): 1.0}    # missing (0, 1)
    with pytest.raises(KeyError, match="missing from graphpop"):
        driver.join_pairs(egrm, gp)


def test_join_pairs_handles_near_zero_reference():
    """The 1e-12 floor on the denominator prevents division blowups
    when the reference entry is near zero."""
    egrm = {(0, 0): 0.0, (0, 1): 1e-15}
    gp = {(0, 0): 1e-13, (0, 1): 2e-15}
    rows = driver.join_pairs(egrm, gp)
    # No NaN / inf produced.
    for r in rows:
        assert r.rel_err == r.rel_err          # not NaN
        assert r.rel_err < float("inf")


def test_join_pairs_carries_h4_optional_field():
    egrm = {(0, 0): 1.0, (0, 1): 0.5}
    gp = dict(egrm)
    h4 = {(0, 0): 0.999999, (0, 1): 0.500001}
    rows = driver.join_pairs(egrm, gp, egrm_h4=h4)
    assert rows[0].egrm_h4 == pytest.approx(0.999999)


# ---------------------------------------------------------------------------
# write_panel_csv + summarise
# ---------------------------------------------------------------------------

def test_write_panel_csv_schema(tmp_path):
    egrm = {(0, 0): 1.0, (0, 1): 0.5}
    rows = driver.join_pairs(egrm, dict(egrm))
    out = tmp_path / "panel.csv"
    driver.write_panel_csv(rows, out)
    with open(out) as fh:
        reader = csv.reader(fh)
        header = next(reader)
        body = list(reader)
    assert header == [
        "sample_a", "sample_b",
        "egrm_ref", "graphpop_branch_grm",
        "abs_diff", "rel_err", "egrm_h4",
    ]
    assert len(body) == 2
    # egrm_h4 column is empty (None) on the identity join.
    assert body[0][-1] == ""


def test_summarise_reports_max_and_mean():
    egrm = {(0, 0): 1.0, (0, 1): 1.0}
    gp = {(0, 0): 1.0 + 1e-7, (0, 1): 1.0 + 1e-9}
    rows = driver.join_pairs(egrm, gp)
    s = driver.summarise(rows)
    assert s["n_pairs"] == 2
    assert s["max_rel_err"] == pytest.approx(1e-7, rel=1e-3)
    assert s["mean_rel_err"] > s["max_rel_err"] / 2   # both rel-errs counted
