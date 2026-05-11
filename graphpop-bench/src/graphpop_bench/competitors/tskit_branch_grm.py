"""tskit branch_grm wrapper for cross-paper benchmarks.

Wraps ``TreeSequence.genetic_relatedness_matrix(mode='branch')``
(Tang & Chiang 2025 *Genetics*) so it's a first-class competitor
in the benchmark harness alongside the H1/H2 tools.

Key difference from PLINK / KING: tskit is a Python library, not
a binary. To keep RSS measurements isolated from the calling
process (pytest / click / numpy baseline contaminates in-process
profiling), the wrapper spawns a Python subprocess that does the
actual tskit work and dumps `.npy` + `.ids` files. The same
``profile_command`` harness used by H1/H2 captures wall-clock +
RSS of the subprocess.

The module is *both* importable AND executable as a script:

    python -m graphpop_bench.competitors.tskit_branch_grm \\
        <input.trees> <output_dir> [<mode>]

See ``PLAN_tskit_branch_grm.md`` for the design + scope.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence

import numpy as np

from ..profiling import ProfilingResult, profile_command, write_receipt


@dataclass
class TskitBranchGrmResult:
    """Outcome of a tskit branch GRM run."""

    grm: np.ndarray              # (n, n) symmetric
    sample_ids: List[str]
    output_dir: Path
    normalised_tsv: Path
    profiling: ProfilingResult


class TskitBranchGrmRunner:
    """Wrapper around ``ts.genetic_relatedness_matrix(mode='branch')``."""

    def __init__(self, python_binary: str | None = None):
        self.python_binary = python_binary or sys.executable

    @staticmethod
    def is_available() -> bool:
        """True if tskit is importable in the current env."""
        try:
            import tskit  # noqa: F401
            return True
        except ImportError:
            return False

    @staticmethod
    def detect_version() -> str | None:
        """``tskit.__version__`` if importable; else ``None``."""
        try:
            import tskit
            return tskit.__version__
        except ImportError:
            return None

    def run(
        self,
        input_path: str | Path,
        output_dir: str | Path,
        *,
        mode: str = "branch",
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> TskitBranchGrmResult:
        """Compute the tskit GRM on ``input_path``; emit TSV + receipt.

        input_path : path to a `.trees` file.
        mode       : "branch" (default; the Paper-2 statistic) or
                     "site". v1 only exercises "branch" in tests.
        """
        if mode not in {"branch", "site"}:
            raise ValueError(
                f"unknown mode {mode!r}; expected 'branch' or 'site'")

        input_path = Path(input_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.python_binary, "-m",
            "graphpop_bench.competitors.tskit_branch_grm",
            str(input_path), str(output_dir), mode,
        ]
        prof = profile_command(cmd, output_dir)
        if prof.exit_code != 0:
            raise RuntimeError(
                f"tskit subprocess exited with code {prof.exit_code}; "
                f"stderr tail:\n{prof.stderr[-2000:]}")

        grm, sample_ids = _load_grm_files(
            grm_npy=output_dir / "grm.npy",
            grm_ids=output_dir / "grm.ids",
        )
        tsv_path = output_dir / "tskit_branch_grm.tsv"
        _build_normalised_tsv(grm, sample_ids, tsv_path)

        receipt = prof.as_receipt(
            tool="tskit_branch_grm",
            tool_version=self.detect_version(),
            seed=seed,
            graphpop_commit=graphpop_commit,
            extra={
                "input_path": str(input_path),
                "mode": mode,
                "normalised_tsv": str(tsv_path),
                "n_samples": len(sample_ids),
            },
        )
        write_receipt(receipt, output_dir)

        return TskitBranchGrmResult(
            grm=grm,
            sample_ids=sample_ids,
            output_dir=output_dir,
            normalised_tsv=tsv_path,
            profiling=prof,
        )


def _inner_compute_and_dump(
    input_path: str | Path,
    output_dir: str | Path,
    mode: str = "branch",
) -> None:
    """Subprocess entry: load `.trees`, compute GRM, dump artefacts.

    Writes:
      <output_dir>/grm.npy — float64 (n, n) matrix
      <output_dir>/grm.ids — newline-separated sample node ids
    """
    import tskit

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ts = tskit.load(str(input_path))
    samples = list(ts.samples())
    sample_sets = [[int(s)] for s in samples]
    # genetic_relatedness_matrix returns shape (n_windows, n_sets,
    # n_sets); with windows=None we get a single window, squeeze.
    arr = ts.genetic_relatedness_matrix(
        sample_sets=sample_sets, mode=mode)
    if arr.ndim == 3:
        arr = arr[0]
    grm = np.asarray(arr, dtype=np.float64)
    np.save(str(output_dir / "grm.npy"), grm)
    (output_dir / "grm.ids").write_text(
        "\n".join(str(s) for s in samples) + "\n")


def _load_grm_files(
    grm_npy: Path, grm_ids: Path,
) -> tuple[np.ndarray, List[str]]:
    """Read the artefacts dumped by ``_inner_compute_and_dump``."""
    if not grm_npy.exists():
        raise FileNotFoundError(f"tskit GRM file not found: {grm_npy}")
    if not grm_ids.exists():
        raise FileNotFoundError(f"tskit IDs file not found: {grm_ids}")
    grm = np.load(str(grm_npy))
    sample_ids = [
        line for line in grm_ids.read_text().splitlines() if line.strip()
    ]
    if grm.ndim != 2 or grm.shape[0] != grm.shape[1]:
        raise ValueError(
            f"tskit GRM file has bad shape: {grm.shape}")
    if grm.shape[0] != len(sample_ids):
        raise ValueError(
            f"tskit GRM size {grm.shape} mismatches "
            f"{len(sample_ids)} sample IDs")
    return grm, sample_ids


def _build_normalised_tsv(
    grm: np.ndarray, sample_ids: Sequence[str], path: Path,
) -> None:
    """Emit ``sample_a, sample_b, kinship`` TSV.

    One row per unordered pair (i ≤ j) **including diagonal** —
    matches the H1 (PLINK GRM) schema for direct join.
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


def main() -> int:
    """Module-as-script entry: load .trees + dump artefacts.

    Usage:
      python -m graphpop_bench.competitors.tskit_branch_grm \\
          <input.trees> <output_dir> [<mode>]
    """
    if len(sys.argv) < 3:
        sys.stderr.write(
            "usage: python -m graphpop_bench.competitors.tskit_branch_grm "
            "<input.trees> <output_dir> [<mode>]\n")
        return 2
    input_path = sys.argv[1]
    output_dir = sys.argv[2]
    mode = sys.argv[3] if len(sys.argv) > 3 else "branch"
    _inner_compute_and_dump(input_path, output_dir, mode)
    return 0


if __name__ == "__main__":
    sys.exit(main())
