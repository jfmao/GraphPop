"""1000G sample-subset extraction helper.

R0.2 utility: parse the standard 1000G sample panel + extract
per-super-pop or per-sub-pop VCF subsets from the pre-downloaded
chrN VCFs, shelling out to `bcftools view -S samples.txt` for
fast BGZF-streaming I/O.

Used by the four Phase-2 real-data drivers (Fig 5b, 3e, 2d/e,
3f) to standardise cohort selection without bespoke filtering
in each panel driver.

Inputs (pre-existing on disk):
- 1000G phased chrN VCFs (NYGC 2022 release, 3,202 samples)
- Panel TSV: 2,504 unrelated samples × (sample, pop, super_pop,
  gender)
"""
from __future__ import annotations

import csv
import json
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence

DEFAULT_VCF_DIR = Path(
    "/mnt/data/GraphPop/data/raw/1000g/vcf")
DEFAULT_PANEL_PATH = Path(
    "/mnt/data/GraphPop/data/raw/1000g/"
    "integrated_call_samples_v3.20130502.ALL.panel")

VCF_FILENAME_TEMPLATE = (
    "1kGP_high_coverage_Illumina.chr{chr}."
    "filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

# Five super-populations in the canonical 1000G panel.
SUPER_POPS = ("AFR", "AMR", "EAS", "EUR", "SAS")


# ---------------------------------------------------------------------------
# Panel parsing
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PanelRow:
    sample: str
    pop: str        # sub-pop e.g. YRI / GBR / JPT
    super_pop: str  # AFR / AMR / EAS / EUR / SAS
    gender: str     # male / female (string in panel)


def load_panel(path: Path = DEFAULT_PANEL_PATH) -> List[PanelRow]:
    """Parse the 4-col 1000G panel TSV (tab-separated).

    Header: `sample\\tpop\\tsuper_pop\\tgender` (some releases
    add trailing tabs in the header line; we ignore those).
    """
    if not path.exists():
        raise FileNotFoundError(f"panel file not found: {path}")
    rows: List[PanelRow] = []
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        # Tolerate trailing-empty-field headers ("sample pop super_pop gender ").
        header = [c for c in header if c]
        if header[:4] != ["sample", "pop", "super_pop", "gender"]:
            raise ValueError(
                f"unexpected panel header: {header!r}")
        for row in reader:
            if not row or not row[0]:
                continue
            if len(row) < 4:
                raise ValueError(f"malformed panel row: {row!r}")
            rows.append(PanelRow(
                sample=row[0], pop=row[1],
                super_pop=row[2], gender=row[3]))
    return rows


def samples_for_super_pop(
    panel: Sequence[PanelRow], super_pop: str,
) -> List[str]:
    """Sample IDs in a single super-pop (e.g. 'EUR' → 503 IDs)."""
    return [r.sample for r in panel if r.super_pop == super_pop]


def samples_for_sub_pop(
    panel: Sequence[PanelRow], pop: str,
) -> List[str]:
    """Sample IDs in a single sub-pop (e.g. 'YRI' → 108 IDs)."""
    return [r.sample for r in panel if r.pop == pop]


def samples_for_super_pop_set(
    panel: Sequence[PanelRow], super_pops: Sequence[str],
) -> List[str]:
    """Sample IDs across a set of super-pops (e.g. ['EUR', 'AFR'])."""
    keep = set(super_pops)
    return [r.sample for r in panel if r.super_pop in keep]


# ---------------------------------------------------------------------------
# VCF extraction (bcftools shell-out)
# ---------------------------------------------------------------------------

@dataclass
class SubsetResult:
    """Outcome of a single subset extraction."""

    input_vcf: Path
    output_vcf: Path
    samples_file: Path
    receipt_path: Path
    n_samples: int
    n_variants: int
    wall_clock_s: float


def extract_subset_vcf(
    *,
    input_vcf: Path,
    output_vcf: Path,
    sample_ids: Sequence[str],
    bcftools_binary: str = "bcftools",
    region: str | None = None,
    timeout: float = 1200.0,
) -> SubsetResult:
    """Extract a sub-VCF + tabix index for the given sample IDs.

    region : optional `chr:start-end` region restriction (passed
             to `bcftools view -r`). Useful for fast micro-tests.
    """
    if not input_vcf.exists():
        raise FileNotFoundError(f"input VCF not found: {input_vcf}")
    if not sample_ids:
        raise ValueError("sample_ids must be non-empty")
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    samples_file = output_vcf.with_suffix(".samples.txt")
    samples_file.write_text("\n".join(sample_ids) + "\n")

    cmd = [bcftools_binary, "view",
           "-S", str(samples_file),
           "-Oz", "-o", str(output_vcf)]
    if region:
        cmd += ["-r", region]
    cmd.append(str(input_vcf))

    t0 = time.monotonic()
    r = subprocess.run(
        cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(
            f"bcftools view exited with code {r.returncode}; "
            f"stderr tail:\n{r.stderr[-2000:]}")

    # Build a tabix index.
    r_idx = subprocess.run(
        [bcftools_binary, "index", "--tbi", str(output_vcf)],
        capture_output=True, text=True, timeout=timeout)
    if r_idx.returncode != 0:
        raise RuntimeError(
            f"bcftools index exited with code {r_idx.returncode}; "
            f"stderr tail:\n{r_idx.stderr[-2000:]}")
    wall = time.monotonic() - t0

    n_variants = _count_variants(output_vcf, bcftools_binary)

    receipt = {
        "input_vcf": str(input_vcf),
        "output_vcf": str(output_vcf),
        "samples_file": str(samples_file),
        "region": region,
        "n_samples_requested": len(sample_ids),
        "n_variants_written": n_variants,
        "wall_clock_s": wall,
        "tool": "bcftools",
    }
    receipt_path = output_vcf.with_suffix(".subset_receipt.json")
    receipt_path.write_text(json.dumps(receipt, indent=2))

    return SubsetResult(
        input_vcf=input_vcf, output_vcf=output_vcf,
        samples_file=samples_file, receipt_path=receipt_path,
        n_samples=len(sample_ids), n_variants=n_variants,
        wall_clock_s=wall,
    )


def _count_variants(
    vcf_path: Path, bcftools_binary: str,
) -> int:
    """`bcftools view -H | wc -l` style variant count."""
    r = subprocess.run(
        [bcftools_binary, "view", "-H", str(vcf_path)],
        capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        return -1
    return r.stdout.count("\n")


# ---------------------------------------------------------------------------
# Orchestration helpers
# ---------------------------------------------------------------------------

def vcf_path_for_chr(
    chr_label: str | int,
    vcf_dir: Path = DEFAULT_VCF_DIR,
) -> Path:
    """Resolve `chr_label` (e.g. 22 or '22') → the NYGC 2022 VCF path."""
    chr_str = str(chr_label).lstrip("chr")
    return vcf_dir / VCF_FILENAME_TEMPLATE.format(chr=chr_str)


def extract_super_pop_subsets(
    *,
    panel: Sequence[PanelRow],
    super_pops: Sequence[str],
    chr_label: str | int,
    output_dir: Path,
    vcf_dir: Path = DEFAULT_VCF_DIR,
    region: str | None = None,
    bcftools_binary: str = "bcftools",
) -> dict[str, SubsetResult]:
    """One sub-VCF per super-pop for the requested chromosome."""
    input_vcf = vcf_path_for_chr(chr_label, vcf_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, SubsetResult] = {}
    for sp in super_pops:
        sample_ids = samples_for_super_pop(panel, sp)
        if not sample_ids:
            continue
        out = output_dir / f"chr{chr_label}_{sp}.vcf.gz"
        results[sp] = extract_subset_vcf(
            input_vcf=input_vcf, output_vcf=out,
            sample_ids=sample_ids, region=region,
            bcftools_binary=bcftools_binary,
        )
    return results
