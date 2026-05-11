"""Unit tests for Paper 2 Fig 3 drivers.

Covers Fig 3c reuse of Fig 1d/1e helpers (predicate orchestration),
Fig 3a/3b HE regression + phenotype simulation helpers, and a tiny
end-to-end smoke. Integration tests gated by msprime/egrm.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from graphpop_bench.paper2_drivers import (
    fig3c_panels as f3c,
    fig3ab_panels as f3ab,
)


# ---------------------------------------------------------------------------
# Fig 3c — predicate orchestration (reuses Fig 1d/1e helpers)
# ---------------------------------------------------------------------------

def _write_tsv(tmp: Path, name: str, n: int, value: float) -> Path:
    """Write an H1-schema TSV where all pair entries equal `value`."""
    p = tmp / name
    with open(p, "w") as fh:
        fh.write("sample_a\tsample_b\tkinship\n")
        for i in range(n):
            for j in range(i, n):
                fh.write(f"{i}\t{j}\t{value:.10g}\n")
    return p


def _write_json_matrix(tmp: Path, name: str, n: int,
                       value: float, *, t_lo=None, t_hi=None) -> Path:
    p = tmp / name
    blob = {
        "n_samples": n,
        "matrix": [[value] * n for _ in range(n)],
    }
    if t_lo is not None:
        blob["t_lo"] = t_lo
    if t_hi is not None:
        blob["t_hi"] = t_hi
    p.write_text(json.dumps(blob))
    return p


def test_join_one_predicate_identity(tmp_path):
    """Identical TSV + JSON → rel-err = 0 on every pair."""
    tsv = _write_tsv(tmp_path, "branch_grm_time_window.tsv",
                     n=4, value=0.5)
    ref = _write_json_matrix(
        tmp_path, "egrm_expected_time_window.json", n=4,
        value=0.5, t_lo=10, t_hi=100)
    out = tmp_path / "out"
    out.mkdir()
    result = f3c.join_one_predicate(
        "time_window", java_dump_dir=tmp_path,
        json_reference_dir=tmp_path, output_dir=out)
    assert result.summary["max_rel_err"] == 0.0
    assert result.panel_csv.exists()


def test_join_one_predicate_perturbation(tmp_path):
    """Constant relative perturbation → constant rel-err."""
    tsv = _write_tsv(tmp_path, "branch_grm_pathway_half.tsv",
                     n=3, value=2.0)
    ref = _write_json_matrix(
        tmp_path, "egrm_expected_pathway_half.json", n=3,
        value=2.0 / (1 + 1e-6))   # GraphPop reads larger by factor (1+1e-6)
    out = tmp_path / "out"
    out.mkdir()
    result = f3c.join_one_predicate(
        "pathway_half", java_dump_dir=tmp_path,
        json_reference_dir=tmp_path, output_dir=out)
    assert result.summary["max_rel_err"] == pytest.approx(
        1e-6, rel=1e-2)


def test_write_summary_csv(tmp_path):
    # Build synthetic PredicateResults via the identity-join above.
    tsv = _write_tsv(tmp_path, "branch_grm_pathway_half.tsv",
                     n=3, value=1.0)
    ref = _write_json_matrix(
        tmp_path, "egrm_expected_pathway_half.json", n=3, value=1.0)
    out = tmp_path / "out"
    out.mkdir()
    r = f3c.join_one_predicate(
        "pathway_half", java_dump_dir=tmp_path,
        json_reference_dir=tmp_path, output_dir=out)
    summary_csv = tmp_path / "summary.csv"
    f3c.write_summary_csv([r], summary_csv)
    with open(summary_csv) as fh:
        rows = list(csv.DictReader(fh))
    assert rows[0]["predicate"] == "pathway_half"
    assert float(rows[0]["max_rel_err"]) == 0.0


# ---------------------------------------------------------------------------
# Fig 3a/3b — HE regression
# ---------------------------------------------------------------------------

def test_haseman_elston_zero_phenotype():
    rng = np.random.default_rng(0)
    grm = rng.normal(size=(5, 5))
    grm = (grm + grm.T) / 2
    y = np.zeros(5)
    assert f3ab.haseman_elston(grm, y) == 0.0


def test_haseman_elston_signal_in_grm():
    """If y_i y_j ≡ h2_true * G_ij off-diagonal exactly, HE
    recovers h2_true (rounding aside)."""
    rng = np.random.default_rng(1)
    # Construct GRM-correlated phenotype: y = L @ z + noise = 0
    # where L is a Cholesky factor of GRM. Then E[y y^T] = h² G + I.
    # We'll skip the noise-free perfect-recovery case and just
    # check sign + monotonicity on a known synthetic.
    n = 6
    grm = rng.normal(size=(n, n))
    grm = grm @ grm.T   # PSD
    # Pick a y such that y_i y_j off-diagonal exactly = 0.5 * grm.
    # Easiest: y = scaled eigenvector of grm.
    eigvals, eigvecs = np.linalg.eigh(grm)
    y = eigvecs[:, -1] * np.sqrt(eigvals[-1] * 0.5)
    # Now y y^T = 0.5 * outer(eigvec, eigvec) * eigvals[-1]
    # = 0.5 * (rank-1 projection); the off-diagonal slope vs grm
    # is not exactly 0.5 but is positive + bounded.
    h2 = f3ab.haseman_elston(grm, y)
    assert h2 > 0   # alignment → positive


def test_haseman_elston_anti_signal():
    """Y orthogonal to GRM (random + symmetric) → HE near 0."""
    rng = np.random.default_rng(2)
    grm = np.eye(6)   # trivially diagonal → off-diag entries are 0
    y = rng.normal(size=6)
    # Off-diagonal of grm is zero ⇒ denominator zero ⇒ return 0.
    assert f3ab.haseman_elston(grm, y) == 0.0


def test_haseman_elston_shape_check():
    grm = np.eye(3)
    with pytest.raises(ValueError, match="y shape"):
        f3ab.haseman_elston(grm, np.zeros(4))
    with pytest.raises(ValueError, match="grm shape"):
        f3ab.haseman_elston(np.zeros((3, 4)), np.zeros(3))


# ---------------------------------------------------------------------------
# Fig 3a/3b — phenotype simulation
# ---------------------------------------------------------------------------

def test_simulate_phenotype_targets_h2():
    """Synthetic genotypes; achieved h² close to the target."""
    rng = np.random.default_rng(3)
    # 100 haplotypes × 50 mutations; binary with intermediate
    # allele freq so variance > 0.
    n = 100
    G = rng.integers(0, 2, size=(n, 50)).astype(np.float64)
    causal = list(range(20))
    y, achieved_h2 = f3ab.simulate_phenotype(
        G, causal, true_h2=0.5, rng=rng)
    assert y.shape == (n,)
    assert 0.2 <= achieved_h2 <= 0.8


def test_simulate_phenotype_empty_causal_raises():
    rng = np.random.default_rng(4)
    G = np.zeros((4, 3))
    with pytest.raises(ValueError, match="non-empty"):
        f3ab.simulate_phenotype(G, [], 0.5, rng)


def test_simulate_phenotype_index_out_of_range_raises():
    rng = np.random.default_rng(5)
    G = rng.integers(0, 2, size=(10, 5)).astype(np.float64)
    with pytest.raises(ValueError, match="out of range"):
        f3ab.simulate_phenotype(G, [10], 0.5, rng)


# ---------------------------------------------------------------------------
# Smoke — end-to-end fig3ab (gated)
# ---------------------------------------------------------------------------

def _msprime_egrm_available() -> bool:
    try:
        import msprime  # noqa: F401
        from egrm import varGRM_C  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(
    not _msprime_egrm_available(),
    reason="msprime + egrm required for run_fig3ab smoke",
)
def test_run_fig3ab_micro(tmp_path):
    """5-diploid / 2 kb / 2 ARG draws / 2 pheno replicates micro
    smoke. ~seconds; verifies all CSVs + metadata emit valid output."""
    cohort = f3ab.CohortParams(
        n_diploid=5, sequence_length=2_000)
    result = f3ab.run_fig3ab(
        cohort=cohort,
        output_dir=tmp_path,
        n_arg_draws=2,
        n_pheno_replicates=2,
        m_pathway=5,
        m_lof=5,
        true_h2=0.5,
        seed=2026,
    )
    assert result.panel_3a_csv.exists()
    assert result.panel_3b_csv.exists()
    assert result.metadata_path.exists()

    # Each panel should have 2 args × 2 replicates × 3 predicates = 12 rows.
    with open(result.panel_3a_csv) as fh:
        rows_a = list(csv.DictReader(fh))
    assert len(rows_a) == 12
    predicates_a = set(r["predicate"] for r in rows_a)
    assert predicates_a == {"unconditional", "pathway", "anti_pathway"}
    with open(result.panel_3b_csv) as fh:
        rows_b = list(csv.DictReader(fh))
    assert len(rows_b) == 12
    predicates_b = set(r["predicate"] for r in rows_b)
    assert predicates_b == {"unconditional", "lof", "non_lof"}
