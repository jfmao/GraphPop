"""KING-robust wrapper for cross-paper benchmarks.

Shells out to ``king --kinship`` (Manichaikul et al. 2010,
DOI 10.1093/bioinformatics/btq559); parses the resulting
``.kin0`` (between-family) and ``.kin`` (within-family) tab-
separated outputs; emits a normalised TSV with the GraphPop
kinship-schema columns ``sample_a, sample_b, kinship``.

KING accepts PLINK BED input only. For VCF inputs, the wrapper
auto-converts via PLINK (requires PLINK on PATH for that path).

NOTE: KING reports the kinship coefficient phi (Manichaikul
2010), NOT a GRM entry. phi(self) = 0.5 by definition;
phi(parent-child) ≈ 0.25; phi(unrelated) ≈ 0. The normalised
TSV preserves whatever KING computes; figure-generation code
aligns the schema with PLINK GRM (which is a different
statistic) as needed.

See ``PLAN_king.md`` in this folder for the design + scope.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from ..profiling import ProfilingResult, profile_command, write_receipt
from .plink_grm import _find_plink_binary

PairKey = Tuple[str, str]


@dataclass
class KingResult:
    """Outcome of a KING run."""

    sample_ids: List[str]
    pair_kinship: Dict[PairKey, float]
    output_dir: Path
    normalised_tsv: Path
    profiling: ProfilingResult


class KingRunner:
    """Wrapper for ``king --kinship`` (or ``king --related``)."""

    def __init__(self, king_binary: str | None = None):
        self.king_binary = king_binary or _find_king_binary()

    @staticmethod
    def is_available() -> bool:
        """True if ``king`` is on PATH."""
        return _find_king_binary() is not None

    @staticmethod
    def detect_version(binary: str) -> str | None:
        """Return the first line of ``<binary>`` output (KING prints
        version info to stdout/stderr when called with no args)."""
        try:
            r = subprocess.run(
                [binary],
                capture_output=True, text=True, timeout=10,
            )
            text = (r.stdout or "") + (r.stderr or "")
            lines = text.strip().splitlines()
            return lines[0] if lines else None
        except (OSError, subprocess.SubprocessError):
            return None

    def run(
        self,
        input_path: str | Path,
        output_dir: str | Path,
        *,
        input_kind: str = "auto",
        mode: str = "kinship",
        related_degree: int = 3,
        extra_args: Sequence[str] | None = None,
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> KingResult:
        """Run KING on ``input_path``; emit normalised TSV + receipt.

        input_path     : `.vcf(.gz)` or BED-prefix.
        input_kind     : "vcf" | "bfile" | "auto".
        mode           : "kinship" (default) | "related".
        related_degree : passed only when mode="related".
        extra_args     : additional KING args (e.g. `["--prevalence", "0.01"]`).
        """
        if self.king_binary is None:
            raise FileNotFoundError(
                "KING binary not on PATH; install KING (Manichaikul "
                "et al. 2010) to run this wrapper")

        if mode not in {"kinship", "related"}:
            raise ValueError(
                f"unknown mode {mode!r}; expected 'kinship' or 'related'")

        input_path = Path(input_path)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        bed_prefix = self._ensure_bed_input(
            input_path, _resolve_input_kind(input_path, input_kind),
            output_dir)
        out_prefix = output_dir / "king_out"

        cmd = _build_king_cmd(
            self.king_binary, bed_prefix, mode, related_degree,
            out_prefix, extra_args or [])
        prof = profile_command(cmd, output_dir)
        if prof.exit_code != 0:
            raise RuntimeError(
                f"KING exited with code {prof.exit_code}; "
                f"stderr tail:\n{prof.stderr[-2000:]}")

        pair_kinship, sample_ids = parse_kin_files(
            kin0_path=_existing_or_none(out_prefix.with_suffix(".kin0")),
            kin_path=_existing_or_none(out_prefix.with_suffix(".kin")),
            fam_path=_existing_or_none(bed_prefix.with_suffix(".fam")),
        )

        tsv_path = output_dir / "king.tsv"
        _build_normalised_tsv(pair_kinship, tsv_path)

        receipt = prof.as_receipt(
            tool="king",
            tool_version=self.detect_version(self.king_binary),
            seed=seed,
            graphpop_commit=graphpop_commit,
            extra={
                "input_path": str(input_path),
                "mode": mode,
                "normalised_tsv": str(tsv_path),
                "n_samples": len(sample_ids),
                "n_pairs": len(pair_kinship),
            },
        )
        write_receipt(receipt, output_dir)

        return KingResult(
            sample_ids=sample_ids,
            pair_kinship=pair_kinship,
            output_dir=output_dir,
            normalised_tsv=tsv_path,
            profiling=prof,
        )

    def _ensure_bed_input(
        self, input_path: Path, kind: str, output_dir: Path,
    ) -> Path:
        """Return a BED prefix usable by KING.

        If ``kind == "vcf"``, run PLINK to convert; otherwise treat
        ``input_path`` as the BED prefix directly.
        """
        if kind == "bfile":
            return input_path
        plink = _find_plink_binary()
        if plink is None:
            raise FileNotFoundError(
                "VCF input requires PLINK on PATH for BED conversion; "
                "install plink2 or pass a BED prefix")
        bed_prefix = output_dir / "kin_input"
        cmd = [
            plink, "--vcf", str(input_path),
            "--make-bed", "--out", str(bed_prefix),
        ]
        r = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600)
        if r.returncode != 0:
            raise RuntimeError(
                f"PLINK VCF→BED conversion failed (exit {r.returncode}); "
                f"stderr:\n{r.stderr[-2000:]}")
        return bed_prefix


def parse_kin_files(
    kin0_path: Path | None,
    kin_path: Path | None,
    fam_path: Path | None = None,
) -> tuple[Dict[PairKey, float], List[str]]:
    """Parse KING's ``.kin0`` + ``.kin`` files into a pair-kinship map.

    Returns:
        pair_kinship : dict mapping ``(iid_a, iid_b)`` (sorted) to phi.
        sample_ids   : sorted list of all IIDs observed in the kin
                       files (and the optional ``.fam`` if provided —
                       useful for cohorts where some samples have no
                       output pairs).
    """
    if kin0_path is None and kin_path is None:
        raise FileNotFoundError(
            "no KING output files (.kin0, .kin) found")

    pair_kinship: Dict[PairKey, float] = {}
    sample_ids: set[str] = set()

    if kin0_path is not None and kin0_path.exists():
        for sa, sb, phi in _iter_kin_rows(kin0_path, between_family=True):
            key = _canonical_pair(sa, sb)
            pair_kinship[key] = phi
            sample_ids.update(key)

    if kin_path is not None and kin_path.exists():
        for sa, sb, phi in _iter_kin_rows(kin_path, between_family=False):
            key = _canonical_pair(sa, sb)
            pair_kinship[key] = phi
            sample_ids.update(key)

    if fam_path is not None and fam_path.exists():
        for line in fam_path.read_text().splitlines():
            parts = line.split()
            if len(parts) >= 2:
                sample_ids.add(parts[1])  # IID is column 2 in .fam

    return pair_kinship, sorted(sample_ids)


def _iter_kin_rows(path: Path, between_family: bool):
    """Yield ``(iid_a, iid_b, kinship)`` triples from a KING .kin/.kin0.

    KING column count varies across versions. Both files have a
    header line; we key the IID and Kinship columns by header
    name, not position.
    """
    text = path.read_text()
    if not text.strip():
        return
    lines = text.splitlines()
    header = lines[0].split()
    cols = {name: idx for idx, name in enumerate(header)}

    if "Kinship" not in cols:
        raise ValueError(
            f"KING output {path} has no 'Kinship' header column; "
            f"found columns: {header}")

    # `.kin0`: columns include FID1, IID1, FID2, IID2.
    # `.kin`:  columns include FID, ID1, ID2 (single family).
    if between_family:
        a_key = _first_present(cols, ["IID1", "ID1"])
        b_key = _first_present(cols, ["IID2", "ID2"])
    else:
        a_key = _first_present(cols, ["ID1", "IID1"])
        b_key = _first_present(cols, ["ID2", "IID2"])
    if a_key is None or b_key is None:
        raise ValueError(
            f"KING output {path} lacks IID columns; header: {header}")

    a_idx = cols[a_key]
    b_idx = cols[b_key]
    phi_idx = cols["Kinship"]

    for raw in lines[1:]:
        if not raw.strip():
            continue
        parts = raw.split()
        if len(parts) <= max(a_idx, b_idx, phi_idx):
            continue
        try:
            yield parts[a_idx], parts[b_idx], float(parts[phi_idx])
        except ValueError:
            # Skip malformed numeric entries silently.
            continue


def _first_present(cols: Dict[str, int], names: Iterable[str]) -> str | None:
    for name in names:
        if name in cols:
            return name
    return None


def _canonical_pair(a: str, b: str) -> PairKey:
    """Return the pair as a sorted tuple so ``(a,b)`` == ``(b,a)`` in the map."""
    return (a, b) if a <= b else (b, a)


def _build_normalised_tsv(
    pair_kinship: Dict[PairKey, float], path: Path,
) -> None:
    """Emit ``sample_a, sample_b, kinship`` TSV.

    Off-diagonal only; pairs sorted by (sample_a, sample_b).
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("sample_a\tsample_b\tkinship\n")
        for (a, b), phi in sorted(pair_kinship.items()):
            fh.write(f"{a}\t{b}\t{phi:.8g}\n")


def _existing_or_none(path: Path) -> Path | None:
    return path if path.exists() else None


def _find_king_binary() -> str | None:
    return shutil.which("king")


def _resolve_input_kind(input_path: Path, kind: str) -> str:
    if kind != "auto":
        if kind not in {"vcf", "bfile"}:
            raise ValueError(f"unknown input_kind: {kind!r}")
        return kind
    s = str(input_path)
    if s.endswith(".vcf") or s.endswith(".vcf.gz"):
        return "vcf"
    return "bfile"


def _build_king_cmd(
    binary: str,
    bed_prefix: Path,
    mode: str,
    related_degree: int,
    out_prefix: Path,
    extra_args: Iterable[str],
) -> list[str]:
    cmd = [binary, "-b", str(bed_prefix) + ".bed"]
    if mode == "kinship":
        cmd.append("--kinship")
    else:
        cmd += ["--related", "--degree", str(related_degree)]
    cmd += ["--prefix", str(out_prefix)]
    cmd += list(extra_args)
    return cmd
