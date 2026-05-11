"""s-LDSC partitioned-heritability wrapper.

Wraps **stratified LD Score Regression** (Finucane et al. 2015,
*Nature Genetics*; tool `ldsc.py` at github.com/bulik/ldsc) as
the annotation-stratified heritability competitor for Paper 2's
Fig 3 panels (pathway-restricted + LoF-class heritability).

Mirrors the H1 / H2 binary-detection pattern (not the H3 / H4
Python-subprocess pattern) because `ldsc.py` is a shippable
executable, not a Python library we import directly.

Output schema diverges from H1 / H3 / H4: s-LDSC produces
**partitioned variance components**, not pairwise kinship. The
normalised output here is per-category-h² + total-h² TSVs.

See ``PLAN_s_ldsc.md`` for the design + scope.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from ..profiling import ProfilingResult, profile_command, write_receipt


@dataclass
class CategoryRow:
    """One row of the s-LDSC ``.results`` table."""

    category: str
    prop_snps: float
    prop_h2: float
    prop_h2_se: float
    enrichment: float
    enrichment_se: float
    enrichment_p: float


@dataclass
class SLdscResult:
    """Outcome of an s-LDSC run."""

    categories: List[CategoryRow]
    total_h2: float | None
    total_h2_se: float | None
    intercept: float | None
    intercept_se: float | None
    lambda_gc: float | None
    n_snps: int | None
    output_dir: Path
    partition_tsv: Path
    total_tsv: Path
    profiling: ProfilingResult


class SLdscRunner:
    """Wrapper for ``ldsc.py --h2 ... --overlap-annot`` (Finucane 2015)."""

    def __init__(self, ldsc_binary: str | None = None):
        self.ldsc_binary = ldsc_binary or _find_ldsc_binary()

    @staticmethod
    def is_available() -> bool:
        """True if ``ldsc.py`` (or ``ldsc``) is on PATH."""
        return _find_ldsc_binary() is not None

    @staticmethod
    def detect_version(binary: str) -> str | None:
        """Probe ``<binary> --version`` for a version string.

        LDSC sometimes prints version info to stderr; we capture
        both streams.
        """
        try:
            r = subprocess.run(
                [binary, "--version"],
                capture_output=True, text=True, timeout=10,
            )
            text = (r.stdout or "") + (r.stderr or "")
            lines = [ln for ln in text.strip().splitlines() if ln.strip()]
            return lines[0] if lines else None
        except (OSError, subprocess.SubprocessError):
            return None

    def run(
        self,
        sumstats: str | Path,
        output_dir: str | Path,
        *,
        ref_ld_chr: str | Path,
        w_ld_chr: str | Path,
        frqfile_chr: str | Path | None = None,
        overlap_annot: bool = True,
        extra_args: Sequence[str] | None = None,
        seed: int | None = None,
        graphpop_commit: str | None = None,
    ) -> SLdscResult:
        """Run s-LDSC on ``sumstats``; emit normalised TSVs + receipt.

        sumstats     : path to LDSC `.sumstats.gz` file.
        ref_ld_chr   : per-chromosome ref-LD-score prefix (LDSC
                       expects `<prefix>{1..22}.l2.ldscore.gz`).
        w_ld_chr     : per-chromosome regression-weight prefix.
        frqfile_chr  : optional per-chromosome allele-frequency prefix.
        overlap_annot: passes `--overlap-annot` (the standard mode
                       for s-LDSC partitioned heritability).
        extra_args   : passed through to `ldsc.py` verbatim.
        """
        if self.ldsc_binary is None:
            raise FileNotFoundError(
                "ldsc.py binary not on PATH; install LDSC "
                "(Bulik-Sullivan et al. 2015) to run this wrapper")

        sumstats = Path(sumstats)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = output_dir / "ldsc"

        cmd = _build_ldsc_cmd(
            self.ldsc_binary, sumstats, ref_ld_chr, w_ld_chr,
            frqfile_chr, overlap_annot, out_prefix, extra_args or [])
        prof = profile_command(cmd, output_dir)
        if prof.exit_code != 0:
            raise RuntimeError(
                f"ldsc.py exited with code {prof.exit_code}; "
                f"stderr tail:\n{prof.stderr[-2000:]}")

        results_path = out_prefix.with_suffix(".results")
        log_path = out_prefix.with_suffix(".log")
        categories = parse_results(results_path)
        (total_h2, total_h2_se, intercept, intercept_se,
         lambda_gc, n_snps) = parse_log(log_path)

        partition_tsv = output_dir / "partition_h2.tsv"
        total_tsv = output_dir / "total_h2.tsv"
        _build_partition_tsv(categories, partition_tsv)
        _build_total_tsv(
            total_h2, total_h2_se, intercept, intercept_se,
            lambda_gc, n_snps, total_tsv)

        receipt = prof.as_receipt(
            tool="s_ldsc",
            tool_version=self.detect_version(self.ldsc_binary),
            seed=seed,
            graphpop_commit=graphpop_commit,
            extra={
                "sumstats": str(sumstats),
                "ref_ld_chr": str(ref_ld_chr),
                "w_ld_chr": str(w_ld_chr),
                "frqfile_chr": (None if frqfile_chr is None
                                else str(frqfile_chr)),
                "overlap_annot": overlap_annot,
                "partition_tsv": str(partition_tsv),
                "total_tsv": str(total_tsv),
                "n_categories": len(categories),
                "total_h2": total_h2,
                "intercept": intercept,
            },
        )
        write_receipt(receipt, output_dir)

        return SLdscResult(
            categories=categories,
            total_h2=total_h2,
            total_h2_se=total_h2_se,
            intercept=intercept,
            intercept_se=intercept_se,
            lambda_gc=lambda_gc,
            n_snps=n_snps,
            output_dir=output_dir,
            partition_tsv=partition_tsv,
            total_tsv=total_tsv,
            profiling=prof,
        )


# ---------------------------------------------------------------------------
# Parsers — .results table and .log file
# ---------------------------------------------------------------------------

_RESULTS_REQUIRED_COLS = (
    "Category",
    "Prop._SNPs",
    "Prop._h2",
    "Prop._h2_std_error",
    "Enrichment",
    "Enrichment_std_error",
    "Enrichment_p",
)


def parse_results(path: Path) -> List[CategoryRow]:
    """Parse an LDSC ``<prefix>.results`` partitioned-h² table.

    The column order has been stable across LDSC versions since
    2018, but we key columns by header name (not position) to
    survive minor reorderings.
    """
    if not path.exists():
        raise FileNotFoundError(f".results file not found: {path}")
    text = path.read_text()
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split()
    missing = [c for c in _RESULTS_REQUIRED_COLS if c not in header]
    if missing:
        raise ValueError(
            f"LDSC .results table at {path} missing required "
            f"columns: {missing}; header was: {header}")
    idx = {name: i for i, name in enumerate(header)}

    rows: List[CategoryRow] = []
    for raw in lines[1:]:
        parts = raw.split()
        if len(parts) < len(header):
            continue
        try:
            rows.append(CategoryRow(
                category=parts[idx["Category"]],
                prop_snps=float(parts[idx["Prop._SNPs"]]),
                prop_h2=float(parts[idx["Prop._h2"]]),
                prop_h2_se=float(parts[idx["Prop._h2_std_error"]]),
                enrichment=float(parts[idx["Enrichment"]]),
                enrichment_se=float(parts[idx["Enrichment_std_error"]]),
                enrichment_p=float(parts[idx["Enrichment_p"]]),
            ))
        except ValueError:
            # Skip malformed numeric rows silently — LDSC sometimes
            # emits "NA" for failed-to-estimate categories.
            continue
    return rows


_LOG_PATTERNS = {
    "total_h2": re.compile(
        r"Total Observed scale h2:\s*([-\d.eE+nan]+)\s*\(([-\d.eE+nan]+)\)"),
    "intercept": re.compile(
        r"Intercept:\s*([-\d.eE+nan]+)\s*\(([-\d.eE+nan]+)\)"),
    "lambda_gc": re.compile(
        r"Lambda GC:\s*([-\d.eE+nan]+)"),
    # Prefer the regression-SNP count (the SNPs actually used in
    # the regression); fall back to any "SNPs remain" line.
    "n_snps_regression": re.compile(
        r"After merging with regression SNP LD,\s*(\d+)\s+SNPs remain"),
    "n_snps_any": re.compile(
        r"After merging.*?,\s*(\d+)\s+SNPs remain"),
}


def parse_log(
    path: Path,
) -> tuple[float | None, float | None,
           float | None, float | None,
           float | None, int | None]:
    """Extract total h², intercept, λ_GC, and regression SNP count
    from an LDSC ``<prefix>.log``.

    Returns ``(h2, h2_se, intercept, intercept_se, lambda_gc,
    n_snps)`` with ``None`` for any pattern that doesn't match.
    """
    if not path.exists():
        return (None, None, None, None, None, None)
    text = path.read_text()

    def _to_float(s: str | None) -> float | None:
        if s is None:
            return None
        try:
            return float(s)
        except ValueError:
            return None

    m = _LOG_PATTERNS["total_h2"].search(text)
    h2 = _to_float(m.group(1)) if m else None
    h2_se = _to_float(m.group(2)) if m else None

    m = _LOG_PATTERNS["intercept"].search(text)
    intercept = _to_float(m.group(1)) if m else None
    intercept_se = _to_float(m.group(2)) if m else None

    m = _LOG_PATTERNS["lambda_gc"].search(text)
    lambda_gc = _to_float(m.group(1)) if m else None

    m = _LOG_PATTERNS["n_snps_regression"].search(text)
    if m is None:
        m = _LOG_PATTERNS["n_snps_any"].search(text)
    n_snps = int(m.group(1)) if m else None

    return (h2, h2_se, intercept, intercept_se, lambda_gc, n_snps)


# ---------------------------------------------------------------------------
# TSV emission — per-category + totals
# ---------------------------------------------------------------------------

def _build_partition_tsv(
    rows: Iterable[CategoryRow], path: Path,
) -> None:
    """Emit ``category, prop_snps, prop_h2, prop_h2_se, enrichment,
    enrichment_se, enrichment_p`` TSV.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(
            "category\tprop_snps\tprop_h2\tprop_h2_se\t"
            "enrichment\tenrichment_se\tenrichment_p\n"
        )
        for r in rows:
            fh.write(
                f"{r.category}\t{r.prop_snps:.8g}\t"
                f"{r.prop_h2:.8g}\t{r.prop_h2_se:.8g}\t"
                f"{r.enrichment:.8g}\t{r.enrichment_se:.8g}\t"
                f"{r.enrichment_p:.8g}\n"
            )


def _build_total_tsv(
    total_h2: float | None,
    total_h2_se: float | None,
    intercept: float | None,
    intercept_se: float | None,
    lambda_gc: float | None,
    n_snps: int | None,
    path: Path,
) -> None:
    """Emit a single-row TSV with ``total_h2, total_h2_se, intercept,
    intercept_se, lambda_gc, n_snps``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)

    def _fmt_f(v: float | None) -> str:
        return "" if v is None else f"{v:.8g}"

    def _fmt_i(v: int | None) -> str:
        return "" if v is None else str(v)

    with open(path, "w") as fh:
        fh.write(
            "total_h2\ttotal_h2_se\tintercept\tintercept_se\t"
            "lambda_gc\tn_snps\n"
        )
        fh.write(
            f"{_fmt_f(total_h2)}\t{_fmt_f(total_h2_se)}\t"
            f"{_fmt_f(intercept)}\t{_fmt_f(intercept_se)}\t"
            f"{_fmt_f(lambda_gc)}\t{_fmt_i(n_snps)}\n"
        )


# ---------------------------------------------------------------------------
# Binary detection + cmd construction
# ---------------------------------------------------------------------------

def _find_ldsc_binary() -> str | None:
    """Return the path to ``ldsc.py`` (preferred) or ``ldsc``.

    LDSC's canonical script name is ``ldsc.py``; some distros ship
    a wrapper called ``ldsc``.
    """
    return shutil.which("ldsc.py") or shutil.which("ldsc")


def _build_ldsc_cmd(
    binary: str,
    sumstats: Path,
    ref_ld_chr: str | Path,
    w_ld_chr: str | Path,
    frqfile_chr: str | Path | None,
    overlap_annot: bool,
    out_prefix: Path,
    extra_args: Iterable[str],
) -> list[str]:
    cmd: list[str] = [
        binary,
        "--h2", str(sumstats),
        "--ref-ld-chr", str(ref_ld_chr),
        "--w-ld-chr", str(w_ld_chr),
        "--out", str(out_prefix),
    ]
    if frqfile_chr is not None:
        cmd += ["--frqfile-chr", str(frqfile_chr)]
    if overlap_annot:
        cmd.append("--overlap-annot")
    cmd += list(extra_args)
    return cmd
