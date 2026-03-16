"""Parse VCF files and extract variant/genotype data for Neo4j import."""

from __future__ import annotations

import logging
from dataclasses import dataclass
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

    # Ancestral allele annotation (optional, from Ensembl EPO FASTA)
    ancestral_allele: str | None = None  # "REF", "ALT", or None
    is_polarized: bool = False  # True if high-confidence EPO alignment

    # Packed genotype arrays (2 bits/sample for genotype, 1 bit/sample for phase)
    # gt_packed: 00=HomRef, 01=Het, 10=HomAlt, 11=Missing
    # phase_packed: which haplotype carries ALT (0 or 1), for hets
    gt_packed: bytes = b""
    phase_packed: bytes = b""

    # Ploidy: 1 bit/sample, same layout as phase_packed. Bit=1 means haploid.
    # Empty for all-diploid chromosomes (backward compatible, zero overhead).
    ploidy_packed: bytes = b""


@dataclass(slots=True)
class PopulationMap:
    """Maps samples to populations; built once from the panel file."""

    sample_ids: list[str]
    pop_ids: list[str]  # sorted unique population labels
    sample_to_pop: dict[str, str]
    pop_to_indices: dict[str, np.ndarray]  # int32 index arrays into sample_ids
    n_samples_per_pop: dict[str, int]
    sample_packed_index: dict[str, int]  # sampleId → VCF column index (for packed arrays)
    n_vcf_samples: int = 0  # total VCF samples (including those not in panel)
    sample_to_sex: dict[str, int] = None  # sampleId → sex (1=male, 2=female, 0=unknown)

    def __post_init__(self):
        if self.sample_to_sex is None:
            self.sample_to_sex = {}


# ---------------------------------------------------------------------------
# Panel loader
# ---------------------------------------------------------------------------


def _load_panel(
    panel_path: str | Path, stratify_by: str
) -> tuple[dict[str, str], dict[str, int]]:
    """Load a panel/PED file and return ({sample_id: pop_label}, {sample_id: sex}).

    Sex encoding: 1=male, 2=female, 0=unknown.

    Auto-detects format from the header line:
    - PED format (3202 samples): whitespace-separated, has 'SampleID' column
    - Panel format (2504 samples): tab-separated, has 'sample' column
    """
    panel_path = Path(panel_path)
    first_line = panel_path.open().readline()

    if "SampleID" in first_line or "sampleID" in first_line:
        # PED-style: columns include SampleID, Population, Superpopulation
        # Use tab separator if file is tab-separated, else whitespace
        sep = "\t" if "\t" in first_line else r"\s+"
        df = pd.read_csv(panel_path, sep=sep, dtype=str)
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

    sample_to_pop = dict(zip(df[sample_col], df[pop_col], strict=False))

    # Extract sex metadata (1=male, 2=female, 0=unknown)
    sample_to_sex: dict[str, int] = {}
    sex_col = col_map.get("sex", col_map.get("gender"))
    if sex_col is not None:
        for sid, val in zip(df[sample_col], df[sex_col], strict=False):
            val_str = str(val).strip().lower()
            if val_str in ("1", "male"):
                sample_to_sex[sid] = 1
            elif val_str in ("2", "female"):
                sample_to_sex[sid] = 2
            else:
                sample_to_sex[sid] = 0

    return sample_to_pop, sample_to_sex


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_SV_PREFIXES = ("<DEL", "<DUP", "<INV", "<INS", "<CNV", "<BND")

# Remap cyvcf2 gt_types → packed 2-bit encoding
# cyvcf2: 0=HOM_REF, 1=HET, 2=MISSING, 3=HOM_ALT
# packed: 00=HomRef, 01=Het, 10=HomAlt, 11=Missing
_GT_REMAP = np.array([0, 1, 3, 2], dtype=np.uint8)


def _detect_ploidy(v) -> tuple[str, np.ndarray]:
    """Detect per-sample ploidy from a cyvcf2 Variant.

    Returns (mode, haploid_flags) where mode is
    'all_diploid', 'all_haploid', or 'mixed'.
    haploid_flags is a boolean array (True = haploid for that sample).
    """
    gts = v.genotypes  # list of [a0, a1, is_phased] or [a0, is_phased]
    n = len(gts)
    haploid_flags = np.zeros(n, dtype=bool)
    for i, g in enumerate(gts):
        if len(g) == 2:  # [allele, is_phased] → haploid
            haploid_flags[i] = True
    n_hap = int(haploid_flags.sum())
    if n_hap == 0:
        return "all_diploid", haploid_flags
    if n_hap == n:
        return "all_haploid", haploid_flags
    return "mixed", haploid_flags


def _build_ploidy_packed(haploid_flags: np.ndarray) -> bytes:
    """Pack boolean haploid flags into bytes (1 bit/sample, LSB first)."""
    return np.packbits(haploid_flags.astype(np.uint8), bitorder="little").tobytes()


def _vectorized_gt_pack(gt_types: np.ndarray) -> bytes:
    """Pack gt_types into 2 bits/sample using vectorized numpy operations."""
    gt_remapped = _GT_REMAP[gt_types]
    n_all = len(gt_remapped)
    n_padded = ((n_all + 3) // 4) * 4
    padded = np.zeros(n_padded, dtype=np.uint8)
    padded[:n_all] = gt_remapped
    groups = padded.reshape(-1, 4)
    packed = (
        groups[:, 0]
        | (groups[:, 1] << 2)
        | (groups[:, 2] << 4)
        | (groups[:, 3] << 6)
    )
    return packed.tobytes()


def _load_ancestral_fasta(fasta_path: Path) -> str:
    """Load an Ensembl ancestral allele FASTA into a single sequence string.

    Uppercase = high-confidence EPO alignment, lowercase = low-confidence.
    '.' or 'N' = unknown. Returns 0-indexed sequence.
    """
    lines = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            lines.append(line.strip())
    seq = "".join(lines)
    logger.info("Loaded ancestral FASTA: %d positions from %s", len(seq), fasta_path)
    return seq


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
        ancestral_fasta: str | Path | None = None,
        ploidy: str = "auto",
        contigs: list[str] | None = None,
    ) -> None:
        self._vcf_path = str(vcf_path)
        self._region = region
        self._include_filtered = include_filtered
        self._n_variants_processed = 0
        self._n_multiallelic_skipped = 0
        self._ploidy_mode_setting = ploidy  # "auto" or "diploid"
        self._contig_whitelist = set(contigs) if contigs else None

        # Load optional ancestral FASTA for allele polarization
        self._ancestral_seq: str | None = None
        if ancestral_fasta is not None:
            self._ancestral_seq = _load_ancestral_fasta(Path(ancestral_fasta))

        # Load panel and build population map
        sample_to_pop, sample_to_sex = _load_panel(panel_path, stratify_by)
        self._pop_map = self._build_pop_map(sample_to_pop, sample_to_sex)

        # Extract contig lengths from VCF ##contig headers (species-agnostic)
        self._contig_lengths: dict[str, int] = {}
        vcf_tmp = VCF(self._vcf_path, lazy=True)
        for name, length in zip(vcf_tmp.seqnames, vcf_tmp.seqlens):
            if length > 0:
                self._contig_lengths[name] = length
        vcf_tmp.close()

        logger.info(
            "VCFParser initialised: %d samples, %d populations, %d contigs",
            len(self._pop_map.sample_ids),
            len(self._pop_map.pop_ids),
            len(self._contig_lengths),
        )

    # -- Public API ---------------------------------------------------------

    @property
    def pop_map(self) -> PopulationMap:
        return self._pop_map

    @property
    def contig_lengths(self) -> dict[str, int]:
        """Chromosome lengths from VCF ##contig headers."""
        return self._contig_lengths

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

    def _build_pop_map(
        self,
        sample_to_pop: dict[str, str],
        sample_to_sex: dict[str, int] | None = None,
    ) -> PopulationMap:
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

        # Build packed_index: VCF column position for each sample in panel
        sample_packed_index = {s: vcf_sample_index[s] for s in sample_ids}

        return PopulationMap(
            sample_ids=sample_ids,
            pop_ids=pop_ids,
            sample_to_pop=mapping,
            pop_to_indices=pop_to_indices,
            n_samples_per_pop=n_samples_per_pop,
            sample_packed_index=sample_packed_index,
            n_vcf_samples=len(vcf_samples),
            sample_to_sex=sample_to_sex or {},
        )

    def _stream(self, vcf: VCF) -> Iterator[VariantRecord]:
        """Core streaming loop over VCF records."""
        pop_map = self._pop_map
        pop_ids = pop_map.pop_ids
        pop_indices = [pop_map.pop_to_indices[p] for p in pop_ids]
        n_pops = len(pop_ids)
        sample_ids_arr = np.array(vcf.samples)
        n_total_samples = len(sample_ids_arr)

        # Ploidy cache: chr → (mode, haploid_flags, ploidy_packed_bytes)
        chr_ploidy_cache: dict[str, tuple[str, np.ndarray, bytes]] = {}
        force_diploid = self._ploidy_mode_setting == "diploid"

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

            # Contig whitelist filter
            if self._contig_whitelist and chrom not in self._contig_whitelist:
                continue

            variant_id = f"{chrom}:{pos}:{ref}:{alt}"
            variant_type = _classify_variant(ref, alt)

            # --- Genotypes ---
            # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN/MISSING, 3=HOM_ALT
            gt_types = v.gt_types  # numpy int8, shape (n_samples,)

            # --- Ploidy detection (once per chromosome) ---
            if force_diploid:
                ploidy_mode = "all_diploid"
                haploid_flags = None
                ploidy_packed_bytes = b""
            elif chrom in chr_ploidy_cache:
                ploidy_mode, haploid_flags, ploidy_packed_bytes = chr_ploidy_cache[chrom]
            else:
                ploidy_mode, haploid_flags = _detect_ploidy(v)
                if ploidy_mode == "all_diploid":
                    ploidy_packed_bytes = b""
                    haploid_flags = None
                else:
                    ploidy_packed_bytes = _build_ploidy_packed(haploid_flags)
                chr_ploidy_cache[chrom] = (ploidy_mode, haploid_flags, ploidy_packed_bytes)
                if ploidy_mode != "all_diploid":
                    logger.info(
                        "Chromosome %s: ploidy mode=%s (%d haploid samples)",
                        chrom, ploidy_mode,
                        int(haploid_flags.sum()) if haploid_flags is not None else 0,
                    )

            # --- Per-population stats ---
            ac = [0] * n_pops
            an = [0] * n_pops
            af = [0.0] * n_pops
            het_count = [0] * n_pops
            hom_alt_count = [0] * n_pops
            het_exp = [0.0] * n_pops

            if ploidy_mode == "all_diploid":
                # FAST PATH: standard diploid counting (unchanged)
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

            elif ploidy_mode == "all_haploid":
                # All samples haploid (MT, chloroplast)
                for k in range(n_pops):
                    idx = pop_indices[k]
                    gt_k = gt_types[idx]

                    n_miss = int(np.sum(gt_k == 2))
                    n_hom_alt = int(np.sum(gt_k == 3))
                    n_called = len(idx) - n_miss

                    ac_k = n_hom_alt  # 1 allele per ALT call
                    an_k = n_called   # 1 allele per called sample

                    ac[k] = ac_k
                    an[k] = an_k
                    af_k = ac_k / an_k if an_k > 0 else 0.0
                    af[k] = af_k
                    het_count[k] = 0  # impossible for haploid
                    hom_alt_count[k] = n_hom_alt
                    het_exp[k] = 2.0 * af_k * (1.0 - af_k)

            else:
                # MIXED: per-sample diploid/haploid split (chrX, chrY)
                for k in range(n_pops):
                    idx = pop_indices[k]
                    gt_k = gt_types[idx]
                    hap_k = haploid_flags[idx]

                    # Diploid subset
                    dip_mask = ~hap_k
                    gt_dip = gt_k[dip_mask]
                    n_miss_dip = int(np.sum(gt_dip == 2))
                    n_het_dip = int(np.sum(gt_dip == 1))
                    n_hom_alt_dip = int(np.sum(gt_dip == 3))
                    n_called_dip = len(gt_dip) - n_miss_dip
                    ac_dip = n_het_dip + 2 * n_hom_alt_dip
                    an_dip = 2 * n_called_dip

                    # Haploid subset
                    gt_hap = gt_k[hap_k]
                    n_miss_hap = int(np.sum(gt_hap == 2))
                    n_hom_alt_hap = int(np.sum(gt_hap == 3))
                    n_called_hap = len(gt_hap) - n_miss_hap
                    ac_hap = n_hom_alt_hap
                    an_hap = n_called_hap

                    ac_k = ac_dip + ac_hap
                    an_k = an_dip + an_hap

                    ac[k] = ac_k
                    an[k] = an_k
                    af_k = ac_k / an_k if an_k > 0 else 0.0
                    af[k] = af_k
                    het_count[k] = n_het_dip  # only diploid samples can be het
                    hom_alt_count[k] = n_hom_alt_dip + n_hom_alt_hap
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

            # --- Packed genotype arrays (vectorized) ---
            gt_packed_data = _vectorized_gt_pack(gt_types)

            # --- Phase packing (for het sites) ---
            n_all = len(gt_types)
            phase_packed_len = (n_all + 7) >> 3
            phase_packed = bytearray(phase_packed_len)

            het_idx = np.flatnonzero(gt_types == 1)
            if len(het_idx) > 0:
                genotypes = v.genotypes
                for i in het_idx:
                    # genotypes[i] = [a0, a1, is_phased]
                    # a1=1 means ALT on second haplotype → phase=1
                    g = genotypes[i]
                    if len(g) >= 3 and g[1] == 1:
                        phase_packed[i >> 3] |= 1 << (i & 7)

            # --- Ancestral allele annotation (optional) ---
            ancestral_allele: str | None = None
            is_polarized = False
            if self._ancestral_seq is not None:
                fasta_idx = pos - 1  # 0-indexed
                if 0 <= fasta_idx < len(self._ancestral_seq):
                    anc_char = self._ancestral_seq[fasta_idx]
                    anc_upper = anc_char.upper()
                    if anc_upper not in (".", "N", "-"):
                        if anc_upper == ref[0].upper():
                            ancestral_allele = "REF"
                            is_polarized = anc_char.isupper()
                        elif anc_upper == alt[0].upper():
                            ancestral_allele = "ALT"
                            is_polarized = anc_char.isupper()

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
                ancestral_allele=ancestral_allele,
                is_polarized=is_polarized,
                gt_packed=gt_packed_data,
                phase_packed=bytes(phase_packed),
                ploidy_packed=ploidy_packed_bytes,
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
