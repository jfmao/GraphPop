"""SINGER wrapper — Deng/Nielsen/Song Bayesian ARG posterior sampler.

Wraps the pre-compiled `singer_master` binary (github.com/
popgenmethods/SINGER, release 0.1.9-beta tested) + the
`convert_to_tskit` companion script. Produces N posterior `.trees`
files from an input VCF; downstream consumers can run any
ARG-statistic on each draw and aggregate.

Used by the P2.5-bis Fig 2d/2e driver to replace the
independent-draw proxy with a real Bayesian ARG posterior.

PATH detection: looks for `singer_master` + `convert_to_tskit`.
The convert helper is a Python script that needs `tskit` +
`numpy` (already in the bench's `compare` extra).
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence

from ..profiling import ProfilingResult, profile_command, write_receipt


@dataclass
class SingerResult:
    """Outcome of a SINGER posterior run."""

    trees_paths: List[Path]
    output_dir: Path
    profiling: ProfilingResult
    n_samples_requested: int
    n_samples_converted: int


class SingerRunner:
    """Wrapper for `singer_master` + `convert_to_tskit`."""

    def __init__(
        self,
        singer_binary: str | None = None,
        convert_binary: str | None = None,
        python_binary: str | None = None,
    ):
        self.singer_binary = singer_binary or _find_singer_binary()
        self.convert_binary = (
            convert_binary or _find_convert_binary())
        self.python_binary = python_binary or sys.executable

    @staticmethod
    def is_available() -> bool:
        """True iff both `singer_master` + `convert_to_tskit` are on PATH."""
        return (
            _find_singer_binary() is not None
            and _find_convert_binary() is not None
        )

    @staticmethod
    def detect_version(binary: str) -> str | None:
        """Probe `singer_master --help` for the version banner.

        SINGER doesn't print a version string by default; fall back
        to the directory name of the install path (e.g.,
        `singer-0.1.9-beta-linux-x86_64`).
        """
        try:
            r = subprocess.run(
                [binary, "--help"],
                capture_output=True, text=True, timeout=10,
            )
            text = (r.stdout or "") + (r.stderr or "")
            # No explicit version flag — return the install-dir name
            # so the receipt still has a stable identifier.
            real_path = Path(binary).resolve()
            return real_path.parent.name
        except (OSError, subprocess.SubprocessError):
            return None

    def run(
        self,
        vcf_path: str | Path,
        output_dir: str | Path,
        *,
        Ne: int,
        mutation_rate: float,
        start: float,
        end: float,
        n_samples: int = 20,
        thin: int = 10,
        ratio: float = 1.0,
        polar: float = 0.5,
        ploidy: int = 2,
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> SingerResult:
        """Run SINGER on `vcf_path`; convert output → `.trees` files.

        Notes on the SINGER CLI:
        - `-vcf <prefix>` — pass WITHOUT the `.vcf` extension.
        - `-output <prefix>` — also a prefix; SINGER appends
          `_branches_{i}.txt` etc.
        - `-Ne` is diploid Ne (haploid = 2 Ne).
        - `-start` + `-end` are base-pair coordinates within the
          VCF (1-indexed, both inclusive of the VCF range).

        Wall-clock + RSS are captured via `profile_command` so the
        receipt is comparable to other competitor wrappers.
        """
        if self.singer_binary is None:
            raise FileNotFoundError(
                "singer_master not on PATH; install SINGER 0.1.9+ "
                "(github.com/popgenmethods/SINGER)")
        if self.convert_binary is None:
            raise FileNotFoundError(
                "convert_to_tskit not on PATH; install SINGER's "
                "companion helper")

        vcf_path = Path(vcf_path)
        if vcf_path.suffix == ".vcf":
            vcf_prefix = vcf_path.with_suffix("")
        else:
            vcf_prefix = vcf_path
        # SINGER appends `.vcf` to the prefix; verify the file exists.
        expected_vcf = vcf_prefix.with_suffix(".vcf")
        if not expected_vcf.exists():
            raise FileNotFoundError(
                f"SINGER expects {expected_vcf} (vcf_path with "
                f".vcf extension); not found")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        singer_out_prefix = output_dir / "singer"

        cmd_singer = _build_singer_cmd(
            self.singer_binary, vcf_prefix, singer_out_prefix,
            Ne=Ne, mutation_rate=mutation_rate,
            start=start, end=end, n_samples=n_samples,
            thin=thin, ratio=ratio, polar=polar, ploidy=ploidy,
            seed=seed,
        )
        prof = profile_command(cmd_singer, output_dir)
        if prof.exit_code != 0:
            raise RuntimeError(
                f"singer_master exited with code {prof.exit_code}; "
                f"stderr tail:\n{prof.stderr[-2000:]}")

        trees_out_prefix = output_dir / "ts"
        cmd_convert = [
            self.python_binary, self.convert_binary,
            "-input", str(singer_out_prefix),
            "-output", str(trees_out_prefix),
            "-start", "0",
            "-end", str(n_samples),
            "-step", "1",
        ]
        r = subprocess.run(
            cmd_convert, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            raise RuntimeError(
                f"convert_to_tskit exited with code {r.returncode}; "
                f"stderr tail:\n{r.stderr[-2000:]}")

        # Collect the converted .trees files.
        trees_paths = sorted(
            output_dir.glob("ts_*.trees"),
            key=lambda p: int(p.stem.split("_")[-1]),
        )

        receipt = prof.as_receipt(
            tool="singer",
            tool_version=self.detect_version(self.singer_binary),
            seed=seed,
            graphpop_commit=graphpop_commit,
            extra={
                "vcf_path": str(expected_vcf),
                "Ne": Ne, "mutation_rate": mutation_rate,
                "start": start, "end": end,
                "n_samples_requested": n_samples,
                "n_samples_converted": len(trees_paths),
                "thin": thin, "ratio": ratio, "polar": polar,
                "ploidy": ploidy,
            },
        )
        write_receipt(receipt, output_dir)

        return SingerResult(
            trees_paths=list(trees_paths),
            output_dir=output_dir,
            profiling=prof,
            n_samples_requested=n_samples,
            n_samples_converted=len(trees_paths),
        )


# ---------------------------------------------------------------------------
# Binary detection + cmd construction
# ---------------------------------------------------------------------------

def _find_singer_binary() -> str | None:
    return shutil.which("singer_master")


def _find_convert_binary() -> str | None:
    return shutil.which("convert_to_tskit")


def _build_singer_cmd(
    binary: str,
    vcf_prefix: Path,
    output_prefix: Path,
    *,
    Ne: int,
    mutation_rate: float,
    start: float,
    end: float,
    n_samples: int,
    thin: int,
    ratio: float,
    polar: float,
    ploidy: int,
    seed: int | None,
) -> List[str]:
    cmd = [
        binary,
        "-Ne", str(int(Ne)),
        "-m", str(mutation_rate),
        "-vcf", str(vcf_prefix),
        "-output", str(output_prefix),
        "-start", str(int(start)),
        "-end", str(int(end)),
        "-n", str(int(n_samples)),
        "-thin", str(int(thin)),
        "-ratio", str(ratio),
        "-polar", str(polar),
        "-ploidy", str(int(ploidy)),
    ]
    if seed is not None:
        cmd += ["-seed", str(int(seed))]
    return cmd
