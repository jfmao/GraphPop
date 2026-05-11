"""Profile a subprocess: wall-clock + RSS peak + CPU times + receipt.

Two profiling backends, selected at runtime:

- **Linux** — wraps the command under ``/usr/bin/time -v`` and
  parses the verbose output. Captures Max RSS, wall clock,
  user/system CPU, exit code.
- **Fallback** — uses Python's ``resource.getrusage(RUSAGE_CHILDREN)``
  delta + ``time.monotonic()``. Less detailed but portable; tests
  run on this path.

Each profiled run emits a JSON receipt alongside the user's
output directory containing the GraphPop commit SHA, tool name +
version, seed, hardware fingerprint, and the captured metrics.
This is the reproducibility receipt required by
``paper/paper2_kinship_arg/benchmark_plan.md`` § 8.
"""
from __future__ import annotations

import json
import os
import platform
import re
import resource
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

# Path to GNU time on Linux; macOS BSD time has a different flag set.
_GNU_TIME = "/usr/bin/time"


@dataclass
class ProfilingResult:
    """Outcome of a profiled command.

    All fields are populated on Linux; on fallback only the
    coarse-grain fields are filled and the verbose-time fields
    remain ``None``.
    """

    cmd: list[str]
    exit_code: int
    wall_clock_s: float
    rss_peak_mb: float
    user_cpu_s: float
    system_cpu_s: float
    stdout: str
    stderr: str
    backend: str  # "gnu_time" or "rusage"
    raw_time_output: str | None = None

    def as_receipt(
        self,
        tool: str,
        tool_version: str | None = None,
        seed: int | None = None,
        graphpop_commit: str | None = None,
        extra: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Build the JSON-serialisable receipt dict.

        Reproducibility receipt schema (per benchmark_plan.md § 8):
        graphpop_commit, tool, tool_version, seed, metrics,
        hardware fingerprint, timestamp.
        """
        return {
            "tool": tool,
            "tool_version": tool_version,
            "seed": seed,
            "graphpop_commit": graphpop_commit,
            "cmd": list(self.cmd),
            "exit_code": self.exit_code,
            "wall_clock_s": self.wall_clock_s,
            "rss_peak_mb": self.rss_peak_mb,
            "user_cpu_s": self.user_cpu_s,
            "system_cpu_s": self.system_cpu_s,
            "backend": self.backend,
            "host_fingerprint": host_fingerprint(),
            "timestamp_iso": datetime.now(timezone.utc).isoformat(),
            **(dict(extra) if extra else {}),
        }


def profile_command(
    cmd: Sequence[str],
    output_dir: str | Path,
    *,
    cwd: str | Path | None = None,
    env: Mapping[str, str] | None = None,
    timeout: float | None = None,
) -> ProfilingResult:
    """Run ``cmd`` under a profiling backend; return metrics.

    cmd        : argv list to execute.
    output_dir : directory created if missing; the caller can drop
                 receipt + output files alongside.
    cwd, env, timeout : passed to ``subprocess.run``.

    Returns a ``ProfilingResult``. Does NOT write the receipt JSON
    automatically — callers do that via ``write_receipt``.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd_list = list(cmd)

    if _gnu_time_available():
        return _profile_with_gnu_time(
            cmd_list, output_dir, cwd=cwd, env=env, timeout=timeout)
    return _profile_with_rusage(
        cmd_list, cwd=cwd, env=env, timeout=timeout)


def write_receipt(
    receipt: Mapping[str, Any],
    output_dir: str | Path,
    filename: str = "receipt.json",
) -> Path:
    """Write the receipt dict next to the benchmark output.

    Returns the path written.
    """
    path = Path(output_dir) / filename
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(receipt, indent=2, sort_keys=True))
    return path


def host_fingerprint() -> str:
    """Compact host fingerprint for receipts."""
    return "-".join(
        [
            platform.system().lower(),
            platform.machine(),
            f"{_total_ram_gb()}gb",
        ]
    )


def _total_ram_gb() -> int:
    """Approximate total RAM in GB via /proc/meminfo on Linux."""
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    return max(1, round(kb / (1024 * 1024)))
    except (FileNotFoundError, ValueError, OSError):
        pass
    return 0


def _gnu_time_available() -> bool:
    """Whether ``/usr/bin/time -v`` exists and supports verbose."""
    if not os.path.exists(_GNU_TIME):
        return False
    # macOS BSD time exists at the same path but doesn't accept -v.
    # Probe it cheaply.
    try:
        r = subprocess.run(
            [_GNU_TIME, "-v", "true"],
            capture_output=True, text=True, timeout=5)
        return r.returncode == 0
    except (OSError, subprocess.SubprocessError):
        return False


def _profile_with_gnu_time(
    cmd: list[str],
    output_dir: Path,
    *,
    cwd: str | Path | None = None,
    env: Mapping[str, str] | None = None,
    timeout: float | None = None,
) -> ProfilingResult:
    """Linux path: wrap with /usr/bin/time -v -o tmp; parse output."""
    time_log = output_dir / ".graphpop_bench_time.log"
    full_cmd = [_GNU_TIME, "-v", "-o", str(time_log)] + cmd
    result = subprocess.run(
        full_cmd,
        capture_output=True, text=True,
        cwd=cwd, env=dict(env) if env else None, timeout=timeout,
    )
    raw = time_log.read_text() if time_log.exists() else ""
    metrics = _parse_gnu_time_verbose(raw)
    try:
        time_log.unlink()
    except OSError:
        pass
    return ProfilingResult(
        cmd=cmd,
        exit_code=result.returncode,
        wall_clock_s=metrics.get("wall_clock_s", 0.0),
        rss_peak_mb=metrics.get("rss_peak_mb", 0.0),
        user_cpu_s=metrics.get("user_cpu_s", 0.0),
        system_cpu_s=metrics.get("system_cpu_s", 0.0),
        stdout=result.stdout or "",
        stderr=result.stderr or "",
        backend="gnu_time",
        raw_time_output=raw,
    )


_TIME_KEY_PATTERNS = {
    "wall_clock_s": re.compile(
        r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (.+)"),
    "rss_peak_kb": re.compile(
        r"Maximum resident set size \(kbytes\): (\d+)"),
    "user_cpu_s": re.compile(r"User time \(seconds\): ([\d.]+)"),
    "system_cpu_s": re.compile(r"System time \(seconds\): ([\d.]+)"),
}


def _parse_gnu_time_verbose(raw: str) -> dict[str, float]:
    """Parse the multi-line ``time -v`` output into a metric dict."""
    out: dict[str, float] = {}
    for key, pat in _TIME_KEY_PATTERNS.items():
        m = pat.search(raw)
        if not m:
            continue
        val = m.group(1).strip()
        if key == "wall_clock_s":
            out["wall_clock_s"] = _parse_elapsed(val)
        elif key == "rss_peak_kb":
            out["rss_peak_mb"] = int(val) / 1024.0
        else:
            try:
                out[key] = float(val)
            except ValueError:
                pass
    return out


def _parse_elapsed(s: str) -> float:
    """Convert ``h:mm:ss`` or ``m:ss[.f]`` to seconds."""
    parts = s.split(":")
    if len(parts) == 3:
        h, m, sec = parts
        return int(h) * 3600 + int(m) * 60 + float(sec)
    if len(parts) == 2:
        m, sec = parts
        return int(m) * 60 + float(sec)
    return float(parts[0])


def _profile_with_rusage(
    cmd: list[str],
    *,
    cwd: str | Path | None = None,
    env: Mapping[str, str] | None = None,
    timeout: float | None = None,
) -> ProfilingResult:
    """Fallback path: time.monotonic() + getrusage(CHILDREN) delta."""
    pre = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.monotonic()
    result = subprocess.run(
        cmd,
        capture_output=True, text=True,
        cwd=cwd, env=dict(env) if env else None, timeout=timeout,
    )
    wall = time.monotonic() - t0
    post = resource.getrusage(resource.RUSAGE_CHILDREN)
    user_cpu = max(0.0, post.ru_utime - pre.ru_utime)
    system_cpu = max(0.0, post.ru_stime - pre.ru_stime)
    # ru_maxrss is KB on Linux, bytes on macOS — use sensible
    # default (KB → MB) here.
    rss_mb = max(0.0, post.ru_maxrss / 1024.0)
    return ProfilingResult(
        cmd=cmd,
        exit_code=result.returncode,
        wall_clock_s=wall,
        rss_peak_mb=rss_mb,
        user_cpu_s=user_cpu,
        system_cpu_s=system_cpu,
        stdout=result.stdout or "",
        stderr=result.stderr or "",
        backend="rusage",
        raw_time_output=None,
    )
