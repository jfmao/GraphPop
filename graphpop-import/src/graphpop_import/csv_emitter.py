"""Emit Neo4j admin-import CSV files from parsed VCF data."""

from __future__ import annotations

import csv
import logging
from io import TextIOWrapper
from pathlib import Path
from typing import TextIO

from .vcf_parser import CarriesRecord, PopulationMap, VCFParser, VariantRecord

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Chromosome lengths (GRCh38)
# ---------------------------------------------------------------------------

CHR_LENGTHS: dict[str, int] = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    # Also support without 'chr' prefix
    "1": 248956422, "2": 242193529, "3": 198295559,
    "4": 190214555, "5": 181538259, "6": 170805979,
    "7": 159345973, "8": 145138636, "9": 138394717,
    "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189,
    "16": 90338345, "17": 83257441, "18": 80373285,
    "19": 58617616, "20": 64444167, "21": 46709983,
    "22": 50818468, "X": 156040895, "Y": 57227415,
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _harmonic(n: int) -> float:
    """Compute Σ(1/i for i in 1..n)."""
    return sum(1.0 / i for i in range(1, n + 1))


def _harmonic2(n: int) -> float:
    """Compute Σ(1/i² for i in 1..n)."""
    return sum(1.0 / (i * i) for i in range(1, n + 1))


def _fmt_float(v: float) -> str:
    """Format a float, stripping trailing zeros."""
    return f"{v:.8g}"


# ---------------------------------------------------------------------------
# CSV headers (neo4j-admin import format)
# ---------------------------------------------------------------------------

VARIANT_HEADER = [
    "variantId:ID(Variant)", ":LABEL",
    "chr", "pos:long", "ref", "alt", "variant_type",
    "pop_ids:string[]",
    "ac:int[]", "an:int[]", "af:float[]",
    "het_count:int[]", "hom_alt_count:int[]",
    "ac_total:int", "an_total:int", "af_total:float", "call_rate:float",
    "het_exp:float[]",
    "ancestral_allele", "is_polarized:boolean",
]

SAMPLE_HEADER = ["sampleId:ID(Sample)", ":LABEL", "population"]

POPULATION_HEADER = [
    "populationId:ID(Population)", ":LABEL",
    "name", "n_samples:int", "a_n:float", "a_n2:float",
]

CHROMOSOME_HEADER = ["chromosomeId:ID(Chromosome)", ":LABEL", "length:long"]

CARRIES_HEADER = [":START_ID(Sample)", ":END_ID(Variant)", ":TYPE", "gt:int", "phase:int"]

NEXT_HEADER = [":START_ID(Variant)", ":END_ID(Variant)", ":TYPE", "distance_bp:long"]

ON_CHROMOSOME_HEADER = [":START_ID(Variant)", ":END_ID(Chromosome)", ":TYPE"]

IN_POPULATION_HEADER = [":START_ID(Sample)", ":END_ID(Population)", ":TYPE"]


# ---------------------------------------------------------------------------
# CSVEmitter
# ---------------------------------------------------------------------------


class CSVEmitter:
    """Write node and relationship CSVs for neo4j-admin import.

    Produces the full set of CSV files needed by ``neo4j-admin database import``
    to construct the GraphPop property graph from parsed VCF data.
    """

    def __init__(
        self,
        out_dir: Path,
        pop_map: PopulationMap,
        *,
        array_delimiter: str = ";",
    ) -> None:
        self._out_dir = Path(out_dir)
        self._pop_map = pop_map
        self._delim = array_delimiter

        # Streaming file handles (opened in _open_streaming_files)
        self._variant_fh: TextIO | None = None
        self._carries_fh: TextIO | None = None
        self._next_fh: TextIO | None = None
        self._on_chrom_fh: TextIO | None = None

        self._variant_writer: csv.writer | None = None
        self._carries_writer: csv.writer | None = None
        self._next_writer: csv.writer | None = None
        self._on_chrom_writer: csv.writer | None = None

        # NEXT-edge tracking: chr → (prev_variant_id, prev_pos)
        self._prev_variant: dict[str, tuple[str, int]] = {}

        # Chromosomes seen (for chromosome_nodes.csv)
        self._chromosomes_seen: set[str] = set()

        # Counters
        self._n_variants = 0
        self._n_carries = 0
        self._n_next = 0
        self._n_on_chrom = 0

    # -- Static nodes (written before streaming) ----------------------------

    def write_static_nodes(self) -> None:
        """Write Sample, Population, and IN_POPULATION CSVs from PopulationMap."""
        self._out_dir.mkdir(parents=True, exist_ok=True)
        self._write_sample_nodes()
        self._write_population_nodes()
        self._write_in_population_edges()
        logger.info(
            "Static CSVs written: %d samples, %d populations",
            len(self._pop_map.sample_ids),
            len(self._pop_map.pop_ids),
        )

    def _write_sample_nodes(self) -> None:
        path = self._out_dir / "sample_nodes.csv"
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(SAMPLE_HEADER)
            for sid in self._pop_map.sample_ids:
                w.writerow([sid, "Sample", self._pop_map.sample_to_pop[sid]])

    def _write_population_nodes(self) -> None:
        path = self._out_dir / "population_nodes.csv"
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(POPULATION_HEADER)
            for pop in self._pop_map.pop_ids:
                n = self._pop_map.n_samples_per_pop[pop]
                # For n diploid samples: a_n = Σ(1/i for i in 1..2n-1)
                two_n_minus_1 = 2 * n - 1
                a_n = _harmonic(two_n_minus_1) if two_n_minus_1 > 0 else 0.0
                a_n2 = _harmonic2(two_n_minus_1) if two_n_minus_1 > 0 else 0.0
                w.writerow([pop, "Population", pop, n, _fmt_float(a_n), _fmt_float(a_n2)])

    def _write_in_population_edges(self) -> None:
        path = self._out_dir / "in_population_edges.csv"
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(IN_POPULATION_HEADER)
            for sid in self._pop_map.sample_ids:
                w.writerow([sid, self._pop_map.sample_to_pop[sid], "IN_POPULATION"])

    # -- Streaming files ----------------------------------------------------

    def _open_streaming_files(self) -> None:
        """Open file handles for variant/carries/next/on_chromosome CSVs."""
        self._variant_fh = open(self._out_dir / "variant_nodes.csv", "w", newline="")
        self._carries_fh = open(self._out_dir / "carries_edges.csv", "w", newline="")
        self._next_fh = open(self._out_dir / "next_edges.csv", "w", newline="")
        self._on_chrom_fh = open(self._out_dir / "on_chromosome_edges.csv", "w", newline="")

        self._variant_writer = csv.writer(self._variant_fh)
        self._carries_writer = csv.writer(self._carries_fh)
        self._next_writer = csv.writer(self._next_fh)
        self._on_chrom_writer = csv.writer(self._on_chrom_fh)

        # Write headers
        self._variant_writer.writerow(VARIANT_HEADER)
        self._carries_writer.writerow(CARRIES_HEADER)
        self._next_writer.writerow(NEXT_HEADER)
        self._on_chrom_writer.writerow(ON_CHROMOSOME_HEADER)

    # -- Chunk processing ---------------------------------------------------

    def process_chunk(self, chunk: list[VariantRecord]) -> None:
        """Append Variant nodes, CARRIES/NEXT/ON_CHROMOSOME edges for a chunk."""
        if self._variant_writer is None:
            self._open_streaming_files()

        delim = self._delim
        pop_ids = self._pop_map.pop_ids
        pop_ids_str = delim.join(pop_ids)

        vw = self._variant_writer
        cw = self._carries_writer
        nw = self._next_writer
        ow = self._on_chrom_writer

        for rec in chunk:
            # -- Variant node --
            ancestral = getattr(rec, "ancestral_allele", None) or ""
            is_pol = getattr(rec, "is_polarized", False)
            vw.writerow([
                rec.id,
                "Variant",
                rec.chr,
                rec.pos,
                rec.ref,
                rec.alt,
                rec.variant_type,
                pop_ids_str,
                delim.join(str(x) for x in rec.ac),
                delim.join(str(x) for x in rec.an),
                delim.join(_fmt_float(x) for x in rec.af),
                delim.join(str(x) for x in rec.het_count),
                delim.join(str(x) for x in rec.hom_alt_count),
                rec.ac_total,
                rec.an_total,
                _fmt_float(rec.af_total),
                _fmt_float(rec.call_rate),
                delim.join(_fmt_float(x) for x in rec.het_exp),
                ancestral,
                str(is_pol).lower(),
            ])
            self._n_variants += 1

            # -- CARRIES edges --
            for c in rec.carries:
                cw.writerow([c.sample_id, c.variant_id, "CARRIES", c.gt, c.phase])
                self._n_carries += 1

            # -- ON_CHROMOSOME edge --
            chrom = rec.chr
            ow.writerow([rec.id, chrom, "ON_CHROMOSOME"])
            self._n_on_chrom += 1
            self._chromosomes_seen.add(chrom)

            # -- NEXT edge --
            prev = self._prev_variant.get(chrom)
            if prev is not None:
                prev_id, prev_pos = prev
                distance = rec.pos - prev_pos
                nw.writerow([prev_id, rec.id, "NEXT", distance])
                self._n_next += 1
            self._prev_variant[chrom] = (rec.id, rec.pos)

    # -- Finalize -----------------------------------------------------------

    def finalize(self) -> None:
        """Write Chromosome nodes, close streaming file handles, log summary."""
        # Write chromosome nodes
        path = self._out_dir / "chromosome_nodes.csv"
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(CHROMOSOME_HEADER)
            for chrom in sorted(self._chromosomes_seen):
                length = CHR_LENGTHS.get(chrom, 0)
                w.writerow([chrom, "Chromosome", length])

        # Close streaming handles
        for fh in (self._variant_fh, self._carries_fh, self._next_fh, self._on_chrom_fh):
            if fh is not None:
                fh.close()

        self._variant_fh = None
        self._carries_fh = None
        self._next_fh = None
        self._on_chrom_fh = None

        logger.info(
            "CSV emission complete: %d variants, %d carries edges, "
            "%d next edges, %d on_chromosome edges, %d chromosomes",
            self._n_variants,
            self._n_carries,
            self._n_next,
            self._n_on_chrom,
            len(self._chromosomes_seen),
        )

    # -- Convenience runner -------------------------------------------------

    @staticmethod
    def run(
        parser: VCFParser,
        out_dir: str | Path,
        *,
        chunk_size: int = 100_000,
    ) -> CSVEmitter:
        """Wire parser → emitter end-to-end and return the emitter.

        Usage::

            emitter = CSVEmitter.run(parser, Path("output/csv"))
        """
        out_dir = Path(out_dir)
        emitter = CSVEmitter(out_dir, parser.pop_map)

        emitter.write_static_nodes()

        for chunk in parser.iter_chunks(chunk_size):
            emitter.process_chunk(chunk)

        emitter.finalize()
        return emitter

    # -- Properties for inspection ------------------------------------------

    @property
    def n_variants(self) -> int:
        return self._n_variants

    @property
    def n_carries(self) -> int:
        return self._n_carries

    @property
    def n_next(self) -> int:
        return self._n_next

    @property
    def n_on_chrom(self) -> int:
        return self._n_on_chrom

    @property
    def chromosomes_seen(self) -> set[str]:
        return self._chromosomes_seen
