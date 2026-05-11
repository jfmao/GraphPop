"""Unit tests for Paper 2 Fig 4 driver helpers.

Covers aggregation, log-log fit, tool registry detection, and a
tiny end-to-end smoke (gated by msprime + tskit/egrm).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from graphpop_bench.paper2_drivers import fig4_panels as panels


# ---------------------------------------------------------------------------
# CellResult + aggregate_by_tool
# ---------------------------------------------------------------------------

def _mk(n_dip: int, tool: str, rep: int,
        wall: float, rss: float) -> panels.CellResult:
    return panels.CellResult(
        n_diploid=n_dip, n_haploid=n_dip * 2,
        sequence_length=30_000, tool=tool, replicate=rep,
        wall_clock_s=wall, rss_peak_mb=rss,
        user_cpu_s=0.0, system_cpu_s=0.0,
        exit_code=0, backend="rusage",
        receipt_path="/dev/null")


def test_aggregate_by_tool_means_and_n():
    rows = [
        _mk(25, "tskit", 0, 1.0, 10.0),
        _mk(25, "tskit", 1, 2.0, 12.0),
        _mk(25, "tskit", 2, 3.0, 14.0),
        _mk(25, "egrm", 0, 5.0, 30.0),
    ]
    agg = panels.aggregate_by_tool(rows)
    tskit_key = (50, "tskit")
    egrm_key = (50, "egrm")
    assert agg[tskit_key]["n"] == 3
    assert agg[tskit_key]["wall_mean"] == pytest.approx(2.0)
    assert agg[tskit_key]["rss_mean"] == pytest.approx(12.0)
    assert agg[tskit_key]["wall_std"] == pytest.approx(1.0)
    assert agg[egrm_key]["n"] == 1
    # Single-sample std is 0 by convention.
    assert agg[egrm_key]["wall_std"] == 0.0


def test_aggregate_empty_returns_empty():
    assert panels.aggregate_by_tool([]) == {}


# ---------------------------------------------------------------------------
# log_log_fit
# ---------------------------------------------------------------------------

def test_log_log_fit_exact_power_law():
    """y = 2 · x^1.5 → slope ≈ 1.5, intercept ≈ log 2."""
    import math
    xs = [1.0, 10.0, 100.0, 1000.0]
    ys = [2.0 * x ** 1.5 for x in xs]
    slope, intercept = panels.log_log_fit(xs, ys)
    assert slope == pytest.approx(1.5, rel=1e-6)
    assert intercept == pytest.approx(math.log(2.0), rel=1e-6)


def test_log_log_fit_quadratic():
    xs = [1.0, 2.0, 4.0, 8.0]
    ys = [x ** 2 for x in xs]
    slope, _ = panels.log_log_fit(xs, ys)
    assert slope == pytest.approx(2.0, rel=1e-9)


def test_log_log_fit_too_few_points_raises():
    with pytest.raises(ValueError, match="≥ 2"):
        panels.log_log_fit([1.0], [2.0])


def test_log_log_fit_collinear_x_raises():
    """All xs identical → zero denominator."""
    with pytest.raises(ValueError, match="collinear"):
        panels.log_log_fit([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])


def test_log_log_fit_length_mismatch_raises():
    with pytest.raises(ValueError, match="length mismatch"):
        panels.log_log_fit([1.0, 2.0], [1.0])


# ---------------------------------------------------------------------------
# Tool registry
# ---------------------------------------------------------------------------

def test_is_tool_available_known_tools_return_bool():
    for tool in panels.ALL_TOOLS:
        assert isinstance(panels.is_tool_available(tool), bool)


def test_is_tool_available_unknown_raises():
    with pytest.raises(ValueError, match="unknown tool"):
        panels.is_tool_available("invented_tool_99")


# ---------------------------------------------------------------------------
# CSV emission
# ---------------------------------------------------------------------------

def test_write_panel_csv_schema(tmp_path):
    rows = [_mk(25, "tskit", 0, 1.0, 10.0)]
    p = tmp_path / "panel.csv"
    panels.write_panel_csv(rows, p)
    with open(p) as fh:
        reader = csv.DictReader(fh)
        out = list(reader)
    assert len(out) == 1
    assert set(out[0].keys()) >= {
        "n_diploid", "n_haploid", "tool", "replicate",
        "wall_clock_s", "rss_peak_mb", "receipt_path",
    }
    assert int(out[0]["n_haploid"]) == 50


# ---------------------------------------------------------------------------
# Smoke — end-to-end (gated by msprime + at least one tool)
# ---------------------------------------------------------------------------

def _msprime_and_tskit_available() -> bool:
    try:
        import msprime  # noqa: F401
    except ImportError:
        return False
    return panels.is_tool_available("tskit")


@pytest.mark.skipif(
    not _msprime_and_tskit_available(),
    reason="msprime + tskit required for run_fig4 smoke",
)
def test_run_fig4_micro_tskit_only(tmp_path):
    """5-diploid / 2 kb / 1 replicate / tskit-only smoke."""
    config = panels.Fig4SweepConfig(
        n_diploid_values=[5, 10],
        n_replicates=1,
        sequence_length=2_000,
        seed=2026,
    )
    result = panels.run_fig4(
        config=config, output_dir=tmp_path,
        tools=("tskit",),
    )
    assert result.panel_csv.exists()
    assert result.metadata_path.exists()
    with open(result.panel_csv) as fh:
        rows = list(csv.DictReader(fh))
    # 2 N values × 1 replicate × 1 tool = 2 rows.
    assert len(rows) == 2
    for r in rows:
        assert r["tool"] == "tskit"
        assert int(r["exit_code"]) == 0
        assert float(r["wall_clock_s"]) > 0


@pytest.mark.skipif(
    not _msprime_and_tskit_available(),
    reason="msprime + tskit required for run_fig4 skip test",
)
def test_run_fig4_skips_missing_tool(tmp_path):
    """Asking for an unavailable tool writes a skip entry.

    Uses `king` for the negative case — KING is not installed in
    the local dev env. If KING ever IS installed, this test
    becomes a no-op (assertion guarded by availability check).
    """
    config = panels.Fig4SweepConfig(
        n_diploid_values=[5],
        n_replicates=1,
        sequence_length=2_000,
    )
    result = panels.run_fig4(
        config=config, output_dir=tmp_path,
        tools=("tskit", "king"),
    )
    skipped = json.loads(result.skip_path.read_text())
    if not panels.is_tool_available("king"):
        assert "king" in skipped["skipped"]
        assert "tskit" in skipped["ran"]
