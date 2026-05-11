"""Tests for the profiling harness."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from graphpop_bench import ProfilingResult, profile_command, write_receipt
from graphpop_bench.profiling import (
    _parse_elapsed,
    _parse_gnu_time_verbose,
    host_fingerprint,
)


def test_profile_simple_command_returns_populated_result(tmp_path):
    r = profile_command(["true"], tmp_path)
    assert isinstance(r, ProfilingResult)
    assert r.exit_code == 0
    assert r.wall_clock_s >= 0.0
    assert r.rss_peak_mb >= 0.0
    assert r.backend in {"gnu_time", "rusage"}
    assert r.cmd == ["true"]


def test_profile_captures_exit_code():
    r = profile_command(["false"], "/tmp")
    assert r.exit_code != 0


def test_profile_captures_stdout(tmp_path):
    r = profile_command(
        [sys.executable, "-c", "print('hello')"], tmp_path)
    assert r.exit_code == 0
    assert "hello" in r.stdout


def test_profile_wall_clock_meaningful(tmp_path):
    # A python sleep of 0.2s should produce wall_clock_s >= ~0.15 to
    # leave room for clock noise; in CI this is plenty.
    r = profile_command(
        [sys.executable, "-c", "import time; time.sleep(0.2)"],
        tmp_path,
    )
    assert r.wall_clock_s >= 0.15


def test_write_receipt_round_trip(tmp_path):
    r = profile_command(["true"], tmp_path)
    receipt = r.as_receipt(
        tool="testtool", tool_version="0.0.1", seed=42,
        graphpop_commit="abcdef1",
    )
    path = write_receipt(receipt, tmp_path)
    assert path.exists()
    loaded = json.loads(path.read_text())
    assert loaded["tool"] == "testtool"
    assert loaded["tool_version"] == "0.0.1"
    assert loaded["seed"] == 42
    assert loaded["graphpop_commit"] == "abcdef1"
    assert loaded["exit_code"] == 0
    assert "host_fingerprint" in loaded
    assert "timestamp_iso" in loaded


def test_parse_elapsed_handles_hms_and_ms_forms():
    assert _parse_elapsed("0:00.50") == pytest.approx(0.50)
    assert _parse_elapsed("0:12") == pytest.approx(12.0)
    assert _parse_elapsed("1:00:00") == pytest.approx(3600.0)


def test_parse_gnu_time_verbose_minimal():
    raw = (
        "\tUser time (seconds): 1.23\n"
        "\tSystem time (seconds): 0.45\n"
        "\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:01.68\n"
        "\tMaximum resident set size (kbytes): 1024000\n"
    )
    m = _parse_gnu_time_verbose(raw)
    assert m["user_cpu_s"] == pytest.approx(1.23)
    assert m["system_cpu_s"] == pytest.approx(0.45)
    assert m["wall_clock_s"] == pytest.approx(1.68)
    assert m["rss_peak_mb"] == pytest.approx(1000.0)


def test_host_fingerprint_format():
    fp = host_fingerprint()
    # Should look like "linux-x86_64-64gb" or similar; at minimum
    # contains the OS family.
    assert "-" in fp
    assert any(seg in fp for seg in ("linux", "darwin", "windows"))


def test_receipt_includes_cmd_and_metrics(tmp_path):
    r = profile_command(["true"], tmp_path)
    receipt = r.as_receipt(tool="x")
    assert receipt["cmd"] == ["true"]
    assert "wall_clock_s" in receipt
    assert "rss_peak_mb" in receipt
    assert "backend" in receipt


def test_profile_creates_output_dir(tmp_path):
    out = tmp_path / "nested" / "dir"
    assert not out.exists()
    profile_command(["true"], out)
    assert out.is_dir()
