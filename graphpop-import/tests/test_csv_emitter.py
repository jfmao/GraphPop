"""Tests for csv_emitter.py — CSV output logic for Neo4j admin-import."""

import csv
import math
import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from graphpop_import.csv_emitter import (
    CSVEmitter,
    CHROMOSOME_HEADER,
    IN_POPULATION_HEADER,
    NEXT_HEADER,
    ON_CHROMOSOME_HEADER,
    POPULATION_HEADER,
    SAMPLE_HEADER,
    VARIANT_HEADER,
    _fmt_float,
    _harmonic,
    _harmonic2,
)
from graphpop_import.vcf_parser import PopulationMap, VariantRecord


# ---------------------------------------------------------------------------
# Helpers: build minimal PopulationMap and VariantRecord fixtures
# ---------------------------------------------------------------------------


def _make_pop_map(
    sample_ids: list[str],
    sample_to_pop: dict[str, str],
    pop_ids: list[str] | None = None,
) -> PopulationMap:
    """Create a minimal PopulationMap for testing."""
    if pop_ids is None:
        pop_ids = sorted(set(sample_to_pop.values()))
    pop_to_indices = {}
    for pop in pop_ids:
        indices = [i for i, s in enumerate(sample_ids) if sample_to_pop.get(s) == pop]
        pop_to_indices[pop] = np.array(indices, dtype=np.int32)
    n_samples_per_pop = {pop: len(idx) for pop, idx in pop_to_indices.items()}
    packed_index = {s: i for i, s in enumerate(sample_ids)}
    return PopulationMap(
        sample_ids=sample_ids,
        pop_ids=pop_ids,
        sample_to_pop=sample_to_pop,
        pop_to_indices=pop_to_indices,
        n_samples_per_pop=n_samples_per_pop,
        sample_packed_index=packed_index,
        n_vcf_samples=len(sample_ids),
        sample_to_sex={"S1": 1, "S2": 2},
    )


def _make_variant(
    vid: str = "chr22:100:A:T",
    chrom: str = "chr22",
    pos: int = 100,
    ref: str = "A",
    alt: str = "T",
    n_pops: int = 2,
) -> VariantRecord:
    """Create a minimal VariantRecord for testing."""
    return VariantRecord(
        id=vid,
        chr=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        variant_type="SNP",
        ac=[5, 3][:n_pops],
        an=[20, 16][:n_pops],
        af=[0.25, 0.1875][:n_pops],
        het_count=[3, 1][:n_pops],
        hom_alt_count=[1, 1][:n_pops],
        het_exp=[0.375, 0.3046875][:n_pops],
        ac_total=8,
        an_total=36,
        af_total=0.22222222,
        call_rate=0.95,
        ancestral_allele="REF",
        is_polarized=True,
        gt_packed=bytes([0b01_00_10_01]),  # 4 samples packed
        phase_packed=bytes([0b0010]),
        ploidy_packed=b"",
    )


# ---------------------------------------------------------------------------
# _harmonic tests
# ---------------------------------------------------------------------------


class TestHarmonic:
    def test_harmonic_1(self):
        """H(1) = 1."""
        assert _harmonic(1) == 1.0

    def test_harmonic_5(self):
        """H(5) = 1 + 1/2 + 1/3 + 1/4 + 1/5 = 137/60."""
        expected = 1.0 + 1 / 2 + 1 / 3 + 1 / 4 + 1 / 5
        assert abs(_harmonic(5) - expected) < 1e-12
        assert abs(_harmonic(5) - 2.283333333333333) < 1e-10

    def test_harmonic_large(self):
        """H(100) should be close to ln(100) + Euler-Mascheroni."""
        euler = 0.5772156649015329
        approx = math.log(100) + euler
        assert abs(_harmonic(100) - approx) < 0.01  # rough check

    def test_harmonic2_1(self):
        """H2(1) = 1."""
        assert _harmonic2(1) == 1.0

    def test_harmonic2_5(self):
        """H2(5) = 1 + 1/4 + 1/9 + 1/16 + 1/25."""
        expected = 1.0 + 1 / 4 + 1 / 9 + 1 / 16 + 1 / 25
        assert abs(_harmonic2(5) - expected) < 1e-12

    def test_harmonic2_converges(self):
        """H2(n) should converge towards pi^2/6."""
        pi2_over_6 = math.pi ** 2 / 6
        assert _harmonic2(10000) < pi2_over_6
        assert abs(_harmonic2(10000) - pi2_over_6) < 0.001


# ---------------------------------------------------------------------------
# _fmt_float tests
# ---------------------------------------------------------------------------


class TestFmtFloat:
    def test_normal_float(self):
        """Normal float is formatted with up to 8 significant digits."""
        result = _fmt_float(0.25)
        assert result == "0.25"

    def test_zero(self):
        assert _fmt_float(0.0) == "0"

    def test_one(self):
        assert _fmt_float(1.0) == "1"

    def test_small_float(self):
        """Small floats use scientific notation via %g."""
        result = _fmt_float(1.23e-10)
        assert "1.23e-10" in result or "1.23e-010" in result

    def test_trailing_zeros_stripped(self):
        """%g format naturally strips trailing zeros."""
        result = _fmt_float(0.5)
        assert result == "0.5"


# ---------------------------------------------------------------------------
# CSVEmitter: static nodes (Sample, Population, IN_POPULATION)
# ---------------------------------------------------------------------------


class TestStaticNodes:
    @pytest.fixture()
    def pop_map(self):
        return _make_pop_map(
            sample_ids=["S1", "S2", "S3"],
            sample_to_pop={"S1": "POP_A", "S2": "POP_A", "S3": "POP_B"},
        )

    def test_sample_nodes_csv(self, tmp_path, pop_map):
        """Sample CSV has correct header and one row per sample."""
        emitter = CSVEmitter(tmp_path, pop_map)
        emitter.write_static_nodes()

        rows = list(csv.reader(open(tmp_path / "sample_nodes.csv")))
        assert rows[0] == SAMPLE_HEADER
        assert len(rows) == 4  # header + 3 samples
        # Check first sample row
        assert rows[1][0] == "S1"
        assert rows[1][1] == "Sample"
        assert rows[1][2] == "POP_A"

    def test_population_nodes_csv(self, tmp_path, pop_map):
        """Population CSV has one row per population with harmonic numbers."""
        emitter = CSVEmitter(tmp_path, pop_map)
        emitter.write_static_nodes()

        rows = list(csv.reader(open(tmp_path / "population_nodes.csv")))
        assert rows[0] == POPULATION_HEADER
        assert len(rows) == 3  # header + POP_A + POP_B

        # POP_A has 2 samples → 2n-1 = 3 → a_n = 1 + 1/2 + 1/3
        pop_a_row = [r for r in rows[1:] if r[0] == "POP_A"][0]
        n_samples = int(pop_a_row[3])
        assert n_samples == 2
        expected_a_n = _harmonic(3)
        assert abs(float(pop_a_row[4]) - expected_a_n) < 1e-6

    def test_in_population_edges_csv(self, tmp_path, pop_map):
        """IN_POPULATION edges link each sample to its population."""
        emitter = CSVEmitter(tmp_path, pop_map)
        emitter.write_static_nodes()

        rows = list(csv.reader(open(tmp_path / "in_population_edges.csv")))
        assert rows[0] == IN_POPULATION_HEADER
        assert len(rows) == 4  # header + 3 edges
        # S3 → POP_B
        s3_row = [r for r in rows[1:] if r[0] == "S3"][0]
        assert s3_row[1] == "POP_B"
        assert s3_row[2] == "IN_POPULATION"


# ---------------------------------------------------------------------------
# CSVEmitter: streaming variant/edge CSVs
# ---------------------------------------------------------------------------


class TestStreamingCSVs:
    @pytest.fixture()
    def pop_map(self):
        return _make_pop_map(
            sample_ids=["S1", "S2"],
            sample_to_pop={"S1": "POP_A", "S2": "POP_B"},
        )

    def test_variant_nodes_csv(self, tmp_path, pop_map):
        """Variant CSV has correct header and row data."""
        emitter = CSVEmitter(tmp_path, pop_map)
        variant = _make_variant()
        emitter.process_chunk([variant])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "variant_nodes.csv")))
        assert rows[0] == VARIANT_HEADER
        assert len(rows) == 2  # header + 1 variant
        row = rows[1]
        assert row[0] == "chr22:100:A:T"  # variantId
        assert row[1] == "Variant"  # :LABEL
        assert row[2] == "chr22"  # chr
        assert row[3] == "100"  # pos
        assert row[4] == "A"  # ref
        assert row[5] == "T"  # alt
        assert row[6] == "SNP"  # variant_type

    def test_variant_pop_arrays(self, tmp_path, pop_map):
        """Per-population arrays are semicolon-delimited."""
        emitter = CSVEmitter(tmp_path, pop_map)
        variant = _make_variant()
        emitter.process_chunk([variant])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "variant_nodes.csv")))
        row = rows[1]
        # pop_ids (column 7)
        assert row[7] == "POP_A;POP_B"
        # ac (column 8)
        assert row[8] == "5;3"
        # an (column 9)
        assert row[9] == "20;16"

    def test_next_edges_created(self, tmp_path, pop_map):
        """NEXT edges connect consecutive variants on the same chromosome."""
        emitter = CSVEmitter(tmp_path, pop_map)
        v1 = _make_variant(vid="chr22:100:A:T", pos=100)
        v2 = _make_variant(vid="chr22:200:G:C", pos=200, ref="G", alt="C")
        v3 = _make_variant(vid="chr22:350:A:G", pos=350, ref="A", alt="G")
        emitter.process_chunk([v1, v2, v3])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "next_edges.csv")))
        assert rows[0] == NEXT_HEADER
        assert len(rows) == 3  # header + 2 NEXT edges
        # First NEXT edge: v1 → v2, distance 100
        assert rows[1][0] == "chr22:100:A:T"
        assert rows[1][1] == "chr22:200:G:C"
        assert rows[1][2] == "NEXT"
        assert rows[1][3] == "100"
        # Second NEXT edge: v2 → v3, distance 150
        assert rows[2][3] == "150"

    def test_no_next_edge_across_chromosomes(self, tmp_path, pop_map):
        """NEXT edges do not span different chromosomes."""
        emitter = CSVEmitter(tmp_path, pop_map)
        v1 = _make_variant(vid="chr21:100:A:T", chrom="chr21", pos=100)
        v2 = _make_variant(vid="chr22:200:G:C", chrom="chr22", pos=200, ref="G", alt="C")
        emitter.process_chunk([v1, v2])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "next_edges.csv")))
        assert len(rows) == 1  # header only, no NEXT edges

    def test_on_chromosome_edges(self, tmp_path, pop_map):
        """Every variant gets an ON_CHROMOSOME edge."""
        emitter = CSVEmitter(tmp_path, pop_map)
        v1 = _make_variant(vid="chr22:100:A:T", pos=100)
        v2 = _make_variant(vid="chr22:200:G:C", pos=200, ref="G", alt="C")
        emitter.process_chunk([v1, v2])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "on_chromosome_edges.csv")))
        assert rows[0] == ON_CHROMOSOME_HEADER
        assert len(rows) == 3  # header + 2 edges
        assert rows[1] == ["chr22:100:A:T", "chr22", "ON_CHROMOSOME"]

    def test_chromosome_nodes(self, tmp_path, pop_map):
        """Chromosome nodes CSV has one entry per chromosome seen."""
        emitter = CSVEmitter(tmp_path, pop_map)
        v1 = _make_variant(vid="chr22:100:A:T", chrom="chr22", pos=100)
        v2 = _make_variant(vid="chr21:50:G:A", chrom="chr21", pos=50, ref="G", alt="A")
        emitter.process_chunk([v1, v2])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "chromosome_nodes.csv")))
        assert rows[0] == CHROMOSOME_HEADER
        assert len(rows) == 3  # header + 2 chromosomes
        chroms = {r[0] for r in rows[1:]}
        assert chroms == {"chr21", "chr22"}
        # Check GRCh38 length for chr22
        chr22_row = [r for r in rows[1:] if r[0] == "chr22"][0]
        assert chr22_row[2] == "50818468"

    def test_contig_length_override(self, tmp_path, pop_map):
        """Custom contig_lengths override GRCh38 defaults."""
        emitter = CSVEmitter(tmp_path, pop_map, contig_lengths={"chr22": 99999})
        v = _make_variant()
        emitter.process_chunk([v])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "chromosome_nodes.csv")))
        chr22_row = [r for r in rows[1:] if r[0] == "chr22"][0]
        assert chr22_row[2] == "99999"

    def test_counters(self, tmp_path, pop_map):
        """Emitter counters track variants, NEXT edges, and ON_CHROMOSOME edges."""
        emitter = CSVEmitter(tmp_path, pop_map)
        v1 = _make_variant(vid="chr22:100:A:T", pos=100)
        v2 = _make_variant(vid="chr22:200:G:C", pos=200, ref="G", alt="C")
        emitter.process_chunk([v1, v2])
        emitter.finalize()

        assert emitter.n_variants == 2
        assert emitter.n_next == 1
        assert emitter.n_on_chrom == 2
        assert emitter.chromosomes_seen == {"chr22"}

    def test_multiple_chunks(self, tmp_path, pop_map):
        """Processing multiple chunks accumulates correctly."""
        emitter = CSVEmitter(tmp_path, pop_map)
        v1 = _make_variant(vid="chr22:100:A:T", pos=100)
        v2 = _make_variant(vid="chr22:200:G:C", pos=200, ref="G", alt="C")
        emitter.process_chunk([v1])
        emitter.process_chunk([v2])
        emitter.finalize()

        assert emitter.n_variants == 2
        assert emitter.n_next == 1  # NEXT state persists across chunks

    def test_ancestral_allele_and_polarization(self, tmp_path, pop_map):
        """Ancestral allele and is_polarized fields are written correctly."""
        emitter = CSVEmitter(tmp_path, pop_map)
        v = _make_variant()
        emitter.process_chunk([v])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "variant_nodes.csv")))
        row = rows[1]
        # ancestral_allele is column 18, is_polarized is column 19
        assert row[18] == "REF"
        assert row[19] == "true"

    def test_empty_packed_arrays(self, tmp_path, pop_map):
        """Variants with empty packed arrays produce empty CSV fields."""
        emitter = CSVEmitter(tmp_path, pop_map)
        v = _make_variant()
        v.gt_packed = b""
        v.phase_packed = b""
        v.ploidy_packed = b""
        emitter.process_chunk([v])
        emitter.finalize()

        rows = list(csv.reader(open(tmp_path / "variant_nodes.csv")))
        row = rows[1]
        # gt_packed (col 20), phase_packed (col 21), ploidy_packed (col 22)
        assert row[20] == ""
        assert row[21] == ""
        assert row[22] == ""
