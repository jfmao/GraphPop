"""egrm reference wrapper for cross-paper benchmarks.

Wraps the Fan, Mancuso & Chiang 2022 *PLOS Genetics* `egrm`
package (`pip install egrm`, currently 0.1). The package
exposes ``varGRM_C(trees, ...) -> (egrm, vargrm, total_mu)``
where `egrm` is the N×N double-centred eGRM matrix used as
the reference ground-truth for GraphPop's M4.1 branch-GRM
procedure.

For Paper 2 this wrapper produces:
  - Fig 1d/1e reference matrix (rel-err vs GraphPop's output).
  - Fig 4a/b head-to-head timing competitor alongside PLINK,
    KING, tskit.

Like H3 (tskit_branch_grm), egrm is a Python library, so the
wrapper spawns a Python subprocess to keep RSS measurements
isolated from the calling process.

See ``PLAN_egrm.md`` for the design + scope.
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence

import numpy as np

from ..profiling import ProfilingResult, profile_command, write_receipt


@dataclass
class EgrmResult:
    """Outcome of an egrm run."""

    egrm: np.ndarray                 # (N, N) double-centred
    vargrm: np.ndarray | None        # (N, N) or None
    total_mu: float
    sample_ids: List[str]
    output_dir: Path
    normalised_tsv: Path
    profiling: ProfilingResult


class EgrmRunner:
    """Wrapper around ``egrm.varGRM_C`` (Fan et al. 2022)."""

    def __init__(self, python_binary: str | None = None):
        self.python_binary = python_binary or sys.executable

    @staticmethod
    def is_available() -> bool:
        """True if ``egrm`` is importable in the current env."""
        try:
            import egrm  # noqa: F401
            return True
        except ImportError:
            return False

    @staticmethod
    def detect_version() -> str | None:
        """``egrm.__version__`` if importable; else ``None``."""
        try:
            import egrm
            return getattr(egrm, "__version__", None)
        except ImportError:
            return None

    def run(
        self,
        input_path: str | Path,
        output_dir: str | Path,
        *,
        compute_var: bool = True,
        rlim: float = 0.0,
        alim: float | None = None,
        left: float = 0.0,
        right: float | None = None,
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> EgrmResult:
        """Compute the egrm on ``input_path``; emit TSV + receipt.

        compute_var : if True (default), also compute varGRM and
                      write `vargrm.npy`.
        rlim, alim  : most-recent / most-ancient time bounds.
        left, right : leftmost / rightmost base-pair bounds.
        """
        if rlim < 0:
            raise ValueError(f"rlim must be ≥ 0, got {rlim!r}")
        if alim is not None and alim <= rlim:
            raise ValueError(
                f"alim ({alim!r}) must exceed rlim ({rlim!r})")
        if left < 0:
            raise ValueError(f"left must be ≥ 0, got {left!r}")
        if right is not None and right <= left:
            raise ValueError(
                f"right ({right!r}) must exceed left ({left!r})")

        input_path = Path(input_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            self.python_binary, "-m",
            "graphpop_bench.competitors.egrm",
            str(input_path), str(output_dir),
            "var" if compute_var else "novar",
            str(rlim),
            "inf" if alim is None else str(alim),
            str(left),
            "inf" if right is None else str(right),
        ]
        prof = profile_command(cmd, output_dir)
        if prof.exit_code != 0:
            raise RuntimeError(
                f"egrm subprocess exited with code {prof.exit_code}; "
                f"stderr tail:\n{prof.stderr[-2000:]}")

        egrm, vargrm, total_mu, sample_ids = _load_egrm_files(
            egrm_npy=output_dir / "egrm.npy",
            vargrm_npy=output_dir / "vargrm.npy",
            egrm_ids=output_dir / "egrm.ids",
            meta_json=output_dir / "meta.json",
            compute_var=compute_var,
        )
        tsv_path = output_dir / "egrm.tsv"
        _build_normalised_tsv(egrm, sample_ids, tsv_path)

        receipt = prof.as_receipt(
            tool="egrm",
            tool_version=self.detect_version(),
            seed=seed,
            graphpop_commit=graphpop_commit,
            extra={
                "input_path": str(input_path),
                "compute_var": compute_var,
                "rlim": rlim,
                "alim": "inf" if alim is None else alim,
                "left": left,
                "right": "inf" if right is None else right,
                "total_mu": total_mu,
                "normalised_tsv": str(tsv_path),
                "n_samples": len(sample_ids),
            },
        )
        write_receipt(receipt, output_dir)

        return EgrmResult(
            egrm=egrm,
            vargrm=vargrm,
            total_mu=total_mu,
            sample_ids=sample_ids,
            output_dir=output_dir,
            normalised_tsv=tsv_path,
            profiling=prof,
        )


def _inner_compute_and_dump(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    compute_var: bool,
    rlim: float,
    alim: float,
    left: float,
    right: float,
) -> None:
    """Subprocess entry: load `.trees`, run varGRM_C, dump artefacts.

    Writes:
      <output_dir>/egrm.npy    — float64 (N, N) eGRM matrix
      <output_dir>/vargrm.npy  — float64 (N, N) varGRM (if compute_var)
      <output_dir>/egrm.ids    — newline-separated sample node ids
      <output_dir>/meta.json   — {"total_mu": float, "compute_var": bool}
    """
    import tskit
    from egrm import varGRM_C

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ts = tskit.load(str(input_path))
    samples = list(ts.samples())

    egrm_mat, vargrm_mat, total_mu = varGRM_C(
        ts,
        rlim=rlim,
        alim=alim,
        left=left,
        right=right,
        var=compute_var,
    )
    np.save(str(output_dir / "egrm.npy"),
            np.asarray(egrm_mat, dtype=np.float64))
    if compute_var and vargrm_mat is not None:
        np.save(str(output_dir / "vargrm.npy"),
                np.asarray(vargrm_mat, dtype=np.float64))
    (output_dir / "egrm.ids").write_text(
        "\n".join(str(s) for s in samples) + "\n")
    (output_dir / "meta.json").write_text(
        json.dumps(
            {"total_mu": float(total_mu),
             "compute_var": bool(compute_var)}))


def _load_egrm_files(
    egrm_npy: Path,
    vargrm_npy: Path,
    egrm_ids: Path,
    meta_json: Path,
    *,
    compute_var: bool,
) -> tuple[np.ndarray, np.ndarray | None, float, List[str]]:
    """Read the artefacts dumped by ``_inner_compute_and_dump``."""
    if not egrm_npy.exists():
        raise FileNotFoundError(f"eGRM file not found: {egrm_npy}")
    if not egrm_ids.exists():
        raise FileNotFoundError(f"IDs file not found: {egrm_ids}")
    if not meta_json.exists():
        raise FileNotFoundError(f"meta file not found: {meta_json}")

    egrm = np.load(str(egrm_npy))
    sample_ids = [
        line for line in egrm_ids.read_text().splitlines() if line.strip()
    ]
    meta = json.loads(meta_json.read_text())
    total_mu = float(meta["total_mu"])

    if egrm.ndim != 2 or egrm.shape[0] != egrm.shape[1]:
        raise ValueError(f"eGRM file has bad shape: {egrm.shape}")
    if egrm.shape[0] != len(sample_ids):
        raise ValueError(
            f"eGRM size {egrm.shape} mismatches "
            f"{len(sample_ids)} sample IDs")

    vargrm: np.ndarray | None = None
    if compute_var:
        if not vargrm_npy.exists():
            raise FileNotFoundError(
                f"varGRM file expected but not found: {vargrm_npy}")
        vargrm = np.load(str(vargrm_npy))
        if vargrm.shape != egrm.shape:
            raise ValueError(
                f"varGRM shape {vargrm.shape} mismatches "
                f"eGRM shape {egrm.shape}")

    return egrm, vargrm, total_mu, sample_ids


def _build_normalised_tsv(
    egrm: np.ndarray, sample_ids: Sequence[str], path: Path,
) -> None:
    """Emit ``sample_a, sample_b, kinship`` TSV.

    One row per unordered pair (i ≤ j) **including diagonal** —
    matches the H1 (PLINK GRM) + H3 (tskit) schema for direct join.
    """
    n = len(sample_ids)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("sample_a\tsample_b\tkinship\n")
        for i in range(n):
            for j in range(i, n):
                fh.write(
                    f"{sample_ids[i]}\t{sample_ids[j]}\t"
                    f"{egrm[i, j]:.8g}\n"
                )


def _parse_float_with_inf(s: str) -> float:
    """Accept ``'inf'`` (case-insensitive) alongside numeric strings."""
    if s.strip().lower() in {"inf", "infinity"}:
        return math.inf
    return float(s)


def main() -> int:
    """Module-as-script entry: load .trees + dump artefacts.

    Usage:
      python -m graphpop_bench.competitors.egrm \\
          <input.trees> <output_dir> [<var|novar>] \\
          [<rlim>] [<alim>] [<left>] [<right>]
    """
    if len(sys.argv) < 3:
        sys.stderr.write(
            "usage: python -m graphpop_bench.competitors.egrm "
            "<input.trees> <output_dir> [<var|novar>] "
            "[<rlim>] [<alim>] [<left>] [<right>]\n")
        return 2
    input_path = sys.argv[1]
    output_dir = sys.argv[2]
    var_flag = sys.argv[3] if len(sys.argv) > 3 else "var"
    rlim = _parse_float_with_inf(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    alim = _parse_float_with_inf(sys.argv[5]) if len(sys.argv) > 5 else math.inf
    left = _parse_float_with_inf(sys.argv[6]) if len(sys.argv) > 6 else 0.0
    right = _parse_float_with_inf(sys.argv[7]) if len(sys.argv) > 7 else math.inf
    if var_flag not in {"var", "novar"}:
        sys.stderr.write(
            f"var flag must be 'var' or 'novar', got {var_flag!r}\n")
        return 2
    _inner_compute_and_dump(
        input_path, output_dir,
        compute_var=(var_flag == "var"),
        rlim=rlim, alim=alim, left=left, right=right,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
