"""PLINK 2.0 GRM wrapper for cross-paper benchmarks.

Shells out to ``plink2 --make-grm-bin`` (or ``plink`` 1.9 as
fallback), parses the binary GRM (`.grm.bin` lower-triangle
float32) + sample IDs (`.grm.id`), and emits a normalised TSV
with the GraphPop kinship-schema columns
(``sample_a, sample_b, kinship``).

See ``PLAN_plink_grm.md`` in this folder for the design + scope.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence

import numpy as np

from ..profiling import ProfilingResult, profile_command, write_receipt


@dataclass
class PlinkGrmResult:
    """Outcome of a PLINK GRM run."""

    grm: np.ndarray              # (n, n) symmetric
    sample_ids: List[str]
    output_dir: Path
    normalised_tsv: Path
    profiling: ProfilingResult


class PlinkGrmRunner:
    """Wrapper for ``plink2 --make-grm-bin``."""

    def __init__(self, plink_binary: str | None = None):
        self.plink_binary = plink_binary or _find_plink_binary()

    @staticmethod
    def is_available() -> bool:
        """True if ``plink2`` (preferred) or ``plink`` is on PATH."""
        return _find_plink_binary() is not None

    @staticmethod
    def detect_version(binary: str) -> str | None:
        """Return the first line of ``<binary> --version`` output."""
        try:
            r = subprocess.run(
                [binary, "--version"],
                capture_output=True, text=True, timeout=10,
            )
            line = (r.stdout or r.stderr).strip().splitlines()
            return line[0] if line else None
        except (OSError, subprocess.SubprocessError):
            return None

    def run(
        self,
        input_path: str | Path,
        output_dir: str | Path,
        *,
        input_kind: str = "auto",
        extra_args: Sequence[str] | None = None,
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> PlinkGrmResult:
        """Run PLINK GRM on ``input_path``; emit TSV + receipt.

        input_path : a `.vcf(.gz)` file or a BED-prefix (no extension).
        input_kind : "vcf" | "bfile" | "auto" (default: auto-detect).
        extra_args : additional PLINK args (e.g. `["--maf", "0.05"]`).
        """
        if self.plink_binary is None:
            raise FileNotFoundError(
                "PLINK binary not on PATH; install plink2 (preferred) "
                "or plink to run this wrapper")

        input_path = Path(input_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = output_dir / "plink_grm"

        cmd = _build_plink_cmd(
            self.plink_binary, input_path, input_kind,
            out_prefix, extra_args or [],
        )
        prof = profile_command(cmd, output_dir)
        if prof.exit_code != 0:
            raise RuntimeError(
                f"PLINK exited with code {prof.exit_code}; "
                f"stderr tail:\n{prof.stderr[-2000:]}")

        grm, sample_ids = parse_grm_bin(
            grm_bin=out_prefix.with_suffix(".grm.bin"),
            grm_id=out_prefix.with_suffix(".grm.id"),
        )
        tsv_path = output_dir / "plink_grm.tsv"
        _build_normalised_tsv(grm, sample_ids, tsv_path)

        receipt = prof.as_receipt(
            tool="plink_grm",
            tool_version=self.detect_version(self.plink_binary),
            seed=seed,
            graphpop_commit=graphpop_commit,
            extra={
                "input_path": str(input_path),
                "input_kind": _resolve_input_kind(input_path, input_kind),
                "normalised_tsv": str(tsv_path),
                "n_samples": len(sample_ids),
            },
        )
        write_receipt(receipt, output_dir)

        return PlinkGrmResult(
            grm=grm,
            sample_ids=sample_ids,
            output_dir=output_dir,
            normalised_tsv=tsv_path,
            profiling=prof,
        )


def parse_grm_bin(
    grm_bin: str | Path, grm_id: str | Path,
) -> tuple[np.ndarray, List[str]]:
    """Parse PLINK's `.grm.bin` + `.grm.id` into ``(grm, sample_ids)``.

    grm_bin : path to the binary GRM (``float32`` lower-triangle).
    grm_id  : path to the FID/IID TSV.

    Returns:
        grm : ``(n, n)`` symmetric numpy ``float64`` array.
        sample_ids : length-``n`` list of IIDs (column 2 of `.grm.id`).
    """
    grm_id = Path(grm_id)
    grm_bin = Path(grm_bin)
    if not grm_id.exists():
        raise FileNotFoundError(
            f"PLINK .grm.id not found: {grm_id}")
    if not grm_bin.exists():
        raise FileNotFoundError(
            f"PLINK .grm.bin not found: {grm_bin}")

    sample_ids: List[str] = []
    for line in grm_id.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split()
        # PLINK writes "FID IID"; we keep IID. Fall back to the
        # single column if the file isn't 2-column.
        iid = parts[1] if len(parts) >= 2 else parts[0]
        sample_ids.append(iid)
    n = len(sample_ids)
    if n == 0:
        raise ValueError(f"PLINK .grm.id had no sample IDs: {grm_id}")

    expected_entries = n * (n + 1) // 2
    flat = np.fromfile(str(grm_bin), dtype=np.float32)
    if flat.size != expected_entries:
        raise ValueError(
            f"PLINK .grm.bin size mismatch: got {flat.size} float32 "
            f"entries; expected {expected_entries} for {n} samples")

    grm = np.zeros((n, n), dtype=np.float64)
    cursor = 0
    for i in range(n):
        row_len = i + 1
        row = flat[cursor:cursor + row_len].astype(np.float64)
        grm[i, : i + 1] = row
        grm[: i + 1, i] = row  # symmetric fill
        cursor += row_len
    return grm, sample_ids


def _build_normalised_tsv(
    grm: np.ndarray, sample_ids: Sequence[str], path: Path,
) -> None:
    """Emit ``sample_a, sample_b, kinship`` TSV — one row per
    unordered pair (i ≤ j), including the diagonal i == j.

    The output respects the GraphPop kinship-export schema so
    figure-generation code can join GraphPop and PLINK outputs by
    `(sample_a, sample_b)` directly.
    """
    n = len(sample_ids)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("sample_a\tsample_b\tkinship\n")
        for i in range(n):
            for j in range(i, n):
                fh.write(
                    f"{sample_ids[i]}\t{sample_ids[j]}\t"
                    f"{grm[i, j]:.8g}\n"
                )


def _find_plink_binary() -> str | None:
    """Prefer plink2 (PLINK 2.0+) over plink (1.9)."""
    for name in ("plink2", "plink"):
        path = shutil.which(name)
        if path is not None:
            return path
    return None


def _resolve_input_kind(input_path: Path, kind: str) -> str:
    if kind != "auto":
        if kind not in {"vcf", "bfile"}:
            raise ValueError(f"unknown input_kind: {kind!r}")
        return kind
    s = str(input_path)
    if s.endswith(".vcf") or s.endswith(".vcf.gz"):
        return "vcf"
    # BED-prefix: no extension; user passed `prefix` not `prefix.bed`.
    return "bfile"


def _build_plink_cmd(
    binary: str,
    input_path: Path,
    input_kind: str,
    out_prefix: Path,
    extra_args: Iterable[str],
) -> list[str]:
    kind = _resolve_input_kind(input_path, input_kind)
    cmd: list[str] = [binary]
    if kind == "vcf":
        cmd += ["--vcf", str(input_path)]
    else:
        cmd += ["--bfile", str(input_path)]
    cmd += ["--make-grm-bin", "--out", str(out_prefix)]
    cmd += list(extra_args)
    return cmd
