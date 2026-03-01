"""Parse VCF files and extract variant/genotype data for Neo4j import."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
from cyvcf2 import VCF

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class CarriesRecord:
    """One per non-reference genotype (HET or HOM_ALT)."""

    sample_id: str
    variant_id: str
    gt: int  # 1 = HET, 2 = HOM_ALT
    phase: int  # 0 = unphased, 1 = phased


@dataclass(slots=True)
class VariantRecord:
    """Summary for a single biallelic variant site."""

    # Identity
    id: str  # "chr22:12345:A:T"
    chr: str
    pos: int
    ref: str
    alt: str
    variant_type: str  # SNP, INDEL, SV

    # Per-population arrays (indexed parallel to PopulationMap.pop_ids)
    ac: list[int]
    an: list[int]
    af: list[float]
    het_count: list[int]
    hom_alt_count: list[int]
    het_exp: list[float]

    # Global summaries
    ac_total: int = 0
    an_total: int = 0
    af_total: float = 0.0
    call_rate: float = 0.0

    # Sparse genotype edges
    carries: list[CarriesRecord] = field(default_factory=list)


@dataclass(slots=True)
class PopulationMap:
    """Maps samples to populations; built once from the panel file."""

    sample_ids: list[str]
    pop_ids: list[str]  # sorted unique population labels
    sample_to_pop: dict[str, str]
    pop_to_indices: dict[str, np.ndarray]  # int32 index arrays into sample_ids
    n_samples_per_pop: dict[str, int]


# ---------------------------------------------------------------------------
# Panel loader
# ---------------------------------------------------------------------------


def _load_panel(panel_path: str | Path, stratify_by: str) -> dict[str, str]:
    """Load a 1000 Genomes panel/PED file and return {sample_id: pop_label}.

    Auto-detects format from the header line:
    - PED format (3202 samples): whitespace-separated, has 'SampleID' column
    - Panel format (2504 samples): tab-separated, has 'sample' column
    """
    panel_path = Path(panel_path)
    first_line = panel_path.open().readline()

    if "SampleID" in first_line or "sampleID" in first_line:
        # PED-style: whitespace-separated, columns include SampleID, Population, Superpopulation
        df = pd.read_csv(panel_path, sep=r"\s+", dtype=str)
        # Normalise column names to lowercase for matching
        col_map = {c.lower(): c for c in df.columns}
        sample_col = col_map.get("sampleid", col_map.get("sample"))
        if stratify_by == "superpopulation":
            pop_col = col_map.get("superpopulation", col_map.get("super_pop"))
        else:
            pop_col = col_map.get("population", col_map.get("pop"))
    else:
        # Panel-style: tab-separated, columns 'sample', 'pop', 'super_pop'
        df = pd.read_csv(panel_path, sep="\t", dtype=str)
        col_map = {c.lower(): c for c in df.columns}
        sample_col = col_map.get("sample")
        if stratify_by == "superpopulation":
            pop_col = col_map.get("super_pop", col_map.get("superpopulation"))
        else:
            pop_col = col_map.get("pop", col_map.get("population"))

    if sample_col is None or pop_col is None:
        raise ValueError(
            f"Cannot find sample/population columns in {panel_path}. "
            f"Columns found: {list(df.columns)}"
        )

    return dict(zip(df[sample_col], df[pop_col], strict=False))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SV_PREFIXES = ("<DEL", "<DUP", "<INV", "<INS", "<CNV", "<BND")


def _classify_variant(ref: str, alt: str) -> str:
    """Return 'SNP', 'INDEL', or 'SV'."""
    if alt.startswith("<") or any(alt.startswith(p) for p in _SV_PREFIXES):
        return "SV"
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    return "INDEL"


# ---------------------------------------------------------------------------
# Main parser
# ---------------------------------------------------------------------------


class VCFParser:
    """Streaming VCF reader that yields per-variant records with population stats.

    Uses cyvcf2 for fast, low-memory access to BCF/VCF files.
    """

    def __init__(
        self,
        vcf_path: str | Path,
        panel_path: str | Path,
        *,
        stratify_by: str = "superpopulation",
        region: str | None = None,
        include_filtered: bool = False,
    ) -> None:
        self._vcf_path = str(vcf_path)
        self._region = region
        self._include_filtered = include_filtered
        self._n_variants_processed = 0
        self._n_multiallelic_skipped = 0

        # Load panel and build population map
        sample_to_pop = _load_panel(panel_path, stratify_by)
        self._pop_map = self._build_pop_map(sample_to_pop)

        logger.info(
            "VCFParser initialised: %d samples, %d populations",
            len(self._pop_map.sample_ids),
            len(self._pop_map.pop_ids),
        )

    # -- Public API ---------------------------------------------------------

    @property
    def pop_map(self) -> PopulationMap:
        return self._pop_map

    @property
    def n_variants_processed(self) -> int:
        return self._n_variants_processed

    @property
    def n_multiallelic_skipped(self) -> int:
        return self._n_multiallelic_skipped

    def __iter__(self) -> Iterator[VariantRecord]:
        """Yield one VariantRecord per qualifying variant site."""
        vcf = VCF(self._vcf_path, lazy=True, threads=2)
        try:
            yield from self._stream(vcf)
        finally:
            vcf.close()

    def iter_chunks(
        self, chunk_size: int = 100_000
    ) -> Iterator[list[VariantRecord]]:
        """Yield lists of up to *chunk_size* VariantRecords."""
        chunk: list[VariantRecord] = []
        for rec in self:
            chunk.append(rec)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk

    # -- Internals ----------------------------------------------------------

    def _build_pop_map(self, sample_to_pop: dict[str, str]) -> PopulationMap:
        """Intersect VCF samples with panel and build index arrays."""
        vcf = VCF(self._vcf_path, lazy=True)
        vcf_samples: list[str] = list(vcf.samples)
        vcf.close()

        sample_ids: list[str] = []
        mapping: dict[str, str] = {}
        missing = 0

        for s in vcf_samples:
            if s in sample_to_pop:
                sample_ids.append(s)
                mapping[s] = sample_to_pop[s]
            else:
                missing += 1

        if missing:
            logger.warning(
                "%d samples in VCF not found in panel — excluded from pop counts",
                missing,
            )

        pop_ids = sorted(set(mapping.values()))

        # Build index arrays (positions in the VCF sample order)
        vcf_sample_index = {s: i for i, s in enumerate(vcf_samples)}
        pop_to_indices: dict[str, np.ndarray] = {}
        n_samples_per_pop: dict[str, int] = {}

        for pop in pop_ids:
            indices = np.array(
                [vcf_sample_index[s] for s in sample_ids if mapping[s] == pop],
                dtype=np.int32,
            )
            pop_to_indices[pop] = indices
            n_samples_per_pop[pop] = len(indices)

        return PopulationMap(
            sample_ids=sample_ids,
            pop_ids=pop_ids,
            sample_to_pop=mapping,
            pop_to_indices=pop_to_indices,
            n_samples_per_pop=n_samples_per_pop,
        )

    def _stream(self, vcf: VCF) -> Iterator[VariantRecord]:
        """Core streaming loop over VCF records."""
        pop_map = self._pop_map
        pop_ids = pop_map.pop_ids
        pop_indices = [pop_map.pop_to_indices[p] for p in pop_ids]
        n_pops = len(pop_ids)
        sample_ids_arr = np.array(vcf.samples)
        n_total_samples = len(sample_ids_arr)

        iterator = vcf(self._region) if self._region else vcf

        for v in iterator:
            # --- Filter checks ---
            if not self._include_filtered and v.FILTER is not None:
                continue

            alts = v.ALT
            if not alts or alts[0] in (".", "<*>", "*"):
                continue

            if len(alts) > 1:
                self._n_multiallelic_skipped += 1
                continue

            # --- Identity ---
            chrom = v.CHROM
            pos = v.POS
            ref = v.REF
            alt = alts[0]
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            variant_type = _classify_variant(ref, alt)

            # --- Genotypes ---
            # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN/MISSING, 3=HOM_ALT
            gt_types = v.gt_types  # numpy int8, shape (n_samples,)
            gt_phases = v.gt_phases  # numpy bool, shape (n_samples,)

            # --- Per-population stats ---
            ac = [0] * n_pops
            an = [0] * n_pops
            af = [0.0] * n_pops
            het_count = [0] * n_pops
            hom_alt_count = [0] * n_pops
            het_exp = [0.0] * n_pops

            for k in range(n_pops):
                idx = pop_indices[k]
                gt_k = gt_types[idx]

                n_miss = int(np.sum(gt_k == 2))
                n_het = int(np.sum(gt_k == 1))
                n_hom_alt = int(np.sum(gt_k == 3))
                n_called = len(idx) - n_miss

                ac_k = n_het + 2 * n_hom_alt
                an_k = 2 * n_called

                ac[k] = ac_k
                an[k] = an_k
                af_k = ac_k / an_k if an_k > 0 else 0.0
                af[k] = af_k
                het_count[k] = n_het
                hom_alt_count[k] = n_hom_alt
                het_exp[k] = 2.0 * af_k * (1.0 - af_k)

            ac_total = sum(ac)
            an_total = sum(an)
            af_total = ac_total / an_total if an_total > 0 else 0.0

            n_missing_total = int(np.sum(gt_types == 2))
            call_rate = (
                (n_total_samples - n_missing_total) / n_total_samples
                if n_total_samples > 0
                else 0.0
            )

            # --- CARRIES records (sparse) ---
            nonref_mask = (gt_types == 1) | (gt_types == 3)
            nonref_idx = np.flatnonzero(nonref_mask)
            carries: list[CarriesRecord] = []
            for i in nonref_idx:
                gt_val = gt_types[i]
                carries.append(
                    CarriesRecord(
                        sample_id=sample_ids_arr[i],
                        variant_id=variant_id,
                        gt=1 if gt_val == 1 else 2,
                        phase=int(gt_phases[i]),
                    )
                )

            self._n_variants_processed += 1
            if self._n_variants_processed % 100_000 == 0:
                logger.info(
                    "Processed %d variants (current: %s)",
                    self._n_variants_processed,
                    variant_id,
                )

            yield VariantRecord(
                id=variant_id,
                chr=chrom,
                pos=pos,
                ref=ref,
                alt=alt,
                variant_type=variant_type,
                ac=ac,
                an=an,
                af=af,
                het_count=het_count,
                hom_alt_count=hom_alt_count,
                het_exp=het_exp,
                ac_total=ac_total,
                an_total=an_total,
                af_total=af_total,
                call_rate=call_rate,
                carries=carries,
            )

        logger.info(
            "VCF parsing complete: %d variants processed",
            self._n_variants_processed,
        )
        if self._n_multiallelic_skipped > 0:
            logger.warning(
                "%d multi-allelic sites skipped (%.1f%% of encountered sites). "
                "To retain these variants, decompose before import with: "
                "bash scripts/prepare-vcf.sh INPUT.vcf.gz OUTPUT.vcf.gz",
                self._n_multiallelic_skipped,
                100.0
                * self._n_multiallelic_skipped
                / (self._n_variants_processed + self._n_multiallelic_skipped),
            )
