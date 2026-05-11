"""Unit tests for the P2.5-bis SINGER driver + wrapper.

Wrapper-level tests run without SINGER installed; integration
test runs the full SINGER MCMC and is gated on the binary
being available.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from graphpop_bench.competitors import SingerRunner
from graphpop_bench.competitors.singer import (
    _build_singer_cmd, _find_singer_binary, _find_convert_binary,
)


def test_is_available_returns_bool():
    assert isinstance(SingerRunner.is_available(), bool)


def test_build_singer_cmd_required_args():
    cmd = _build_singer_cmd(
        "singer_master", Path("cohort"), Path("/tmp/out/singer"),
        Ne=10_000, mutation_rate=1.25e-8,
        start=0, end=1_000_000,
        n_samples=20, thin=10, ratio=1.0, polar=0.5, ploidy=2,
        seed=42,
    )
    assert cmd[0] == "singer_master"
    assert "-Ne" in cmd and "10000" in cmd
    assert "-m" in cmd and "1.25e-08" in cmd
    assert "-vcf" in cmd and "cohort" in cmd
    assert "-output" in cmd and "/tmp/out/singer" in cmd
    assert "-start" in cmd and "0" in cmd
    assert "-end" in cmd and "1000000" in cmd
    assert "-n" in cmd and "20" in cmd
    assert "-thin" in cmd and "10" in cmd
    assert "-seed" in cmd and "42" in cmd


def test_build_singer_cmd_seed_optional():
    cmd = _build_singer_cmd(
        "singer_master", Path("cohort"), Path("out"),
        Ne=10_000, mutation_rate=1e-8,
        start=0, end=1_000,
        n_samples=10, thin=5, ratio=1.0, polar=0.5, ploidy=2,
        seed=None,
    )
    assert "-seed" not in cmd


def test_run_raises_when_singer_missing(tmp_path):
    """Force-construct a runner with both binaries None."""
    runner = SingerRunner.__new__(SingerRunner)
    runner.singer_binary = None
    runner.convert_binary = "/bin/true"
    runner.python_binary = "/usr/bin/python3"
    with pytest.raises(FileNotFoundError, match="singer_master"):
        runner.run(
            tmp_path / "x.vcf", tmp_path / "out",
            Ne=10_000, mutation_rate=1e-8,
            start=0, end=1_000)


def test_run_raises_when_convert_missing(tmp_path):
    runner = SingerRunner.__new__(SingerRunner)
    runner.singer_binary = "/bin/true"
    runner.convert_binary = None
    runner.python_binary = "/usr/bin/python3"
    with pytest.raises(FileNotFoundError, match="convert_to_tskit"):
        runner.run(
            tmp_path / "x.vcf", tmp_path / "out",
            Ne=10_000, mutation_rate=1e-8,
            start=0, end=1_000)


def test_run_raises_when_vcf_missing(tmp_path):
    """If singer_master + convert are present but VCF file is not."""
    runner = SingerRunner.__new__(SingerRunner)
    runner.singer_binary = "/bin/true"
    runner.convert_binary = "/bin/true"
    runner.python_binary = "/usr/bin/python3"
    with pytest.raises(FileNotFoundError, match=r"SINGER expects"):
        runner.run(
            tmp_path / "missing.vcf", tmp_path / "out",
            Ne=10_000, mutation_rate=1e-8,
            start=0, end=1_000)


# ---------------------------------------------------------------------------
# Integration: gated on SINGER installation
# ---------------------------------------------------------------------------

def _singer_and_egrm_available() -> bool:
    if not SingerRunner.is_available():
        return False
    try:
        import msprime  # noqa: F401
        from egrm import varGRM_C  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(
    not _singer_and_egrm_available(),
    reason="SINGER + msprime + egrm required for fig2de-singer smoke",
)
def test_run_fig2de_singer_micro(tmp_path):
    """Tiny end-to-end smoke: 3+3 diploid / 5 kb / N_posterior=3."""
    from graphpop_bench.paper2_drivers import (
        fig2de_panels as f2de_v1,
        fig2de_singer_panels as f2de_singer,
    )
    params = f2de_v1.TwoPopParams(
        n_diploid_per_pop={"AFR": 3, "EUR": 3},
        sequence_length=5_000)
    result = f2de_singer.run_fig2de_singer(
        params=params,
        n_posterior=3,
        thin=5,
        Ne=10_000,
        seed=2026,
        output_dir=tmp_path,
    )
    assert result.fig2d_csv.exists()
    assert result.fig2e_csv.exists()
    assert result.metadata_path.exists()
    assert result.n_posterior_converted >= 2
    with open(result.fig2e_csv) as fh:
        rows = list(csv.DictReader(fh))
    labels = {r["pop_pair_label"] for r in rows}
    assert labels == {"AFR-AFR", "AFR-EUR", "EUR-EUR"}
