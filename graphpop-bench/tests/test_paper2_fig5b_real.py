"""Unit tests for Paper 2 Fig 5b real-data driver (R1)."""
from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

import pytest

from graphpop_bench.paper2_drivers import fig5b_real_panels as f5br
from graphpop_bench.paper2_drivers.subset_1000g import (
    DEFAULT_PANEL_PATH,
    extract_subset_vcf,
    load_panel,
    samples_for_super_pop,
    vcf_path_for_chr,
)


def _real_data_available() -> bool:
    return (
        DEFAULT_PANEL_PATH.exists()
        and vcf_path_for_chr(22).exists()
        and shutil.which("bcftools") is not None
    )


def _cyvcf2_available() -> bool:
    try:
        import cyvcf2  # noqa: F401
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# load_real_cohort_haploids (gated)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not (_real_data_available() and _cyvcf2_available()),
    reason="1000G chr22 + bcftools + cyvcf2 required",
)
def test_load_real_cohort_haploids_bi_allelic_only(tmp_path):
    """Tiny region: 5 EUR samples × chr22:50_000_000-50_010_000."""
    panel = load_panel()
    samples = samples_for_super_pop(panel, "EUR")[:5]
    cohort_vcf = tmp_path / "tiny.vcf.gz"
    extract_subset_vcf(
        input_vcf=vcf_path_for_chr(22),
        output_vcf=cohort_vcf,
        sample_ids=samples,
        region="chr22:50000000-50010000",
    )
    haploids, positions, returned_samples = (
        f5br.load_real_cohort_haploids(cohort_vcf, samples))
    assert haploids.shape[0] == 10           # 5 diploids = 10 haploids
    assert haploids.shape[1] == positions.size
    # Bi-allelic + polymorphic filter: entries are 0 or 1.
    assert ((haploids == 0) | (haploids == 1)).all()
    assert returned_samples == samples


@pytest.mark.skipif(
    not _cyvcf2_available(),
    reason="cyvcf2 required",
)
def test_load_real_cohort_haploids_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError, match="VCF not found"):
        f5br.load_real_cohort_haploids(tmp_path / "nope.vcf.gz")


# ---------------------------------------------------------------------------
# Smoke run_fig5b_real (gated)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not (_real_data_available() and _cyvcf2_available()),
    reason="1000G chr22 + bcftools + cyvcf2 + PLINK required",
)
def test_run_fig5b_real_micro(tmp_path):
    """Tiny end-to-end: 20-sample EUR × 1 fraction × 1 replicate."""
    panel = load_panel()
    samples = samples_for_super_pop(panel, "EUR")[:20]
    # Pre-subset to a small region so the full driver runs in seconds.
    cohort_vcf = tmp_path / "EUR_micro.vcf.gz"
    extract_subset_vcf(
        input_vcf=vcf_path_for_chr(22),
        output_vcf=cohort_vcf,
        sample_ids=samples,
        region="chr22:50000000-50200000",
    )
    # Mock the extract step inside run_fig5b_real by setting up the
    # cohort VCF where run_fig5b_real will find it.
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    target = output_dir / "EUR_micro.vcf.gz"
    target.write_bytes(cohort_vcf.read_bytes())
    target.with_suffix(".gz.tbi").write_bytes(
        cohort_vcf.with_suffix(".gz.tbi").read_bytes())
    result = f5br.run_fig5b_real(
        cohort_label="EUR_micro",
        sample_ids=samples,
        vcf_path=vcf_path_for_chr(22),
        output_dir=output_dir,
        fractions=[0.5],
        n_pairs=3,
        n_replicates=1,
        threshold=0.0625,
        seed=2026,
    )
    assert result.panel_csv.exists()
    with open(result.panel_csv) as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1
    assert rows[0]["cohort_label"] == "EUR_micro"
    assert float(rows[0]["recall"]) >= 0.0
    # Metadata captures all the config.
    meta = json.loads(result.metadata_path.read_text())
    assert meta["n_samples_requested"] == 20
    assert meta["fractions"] == [0.5]
