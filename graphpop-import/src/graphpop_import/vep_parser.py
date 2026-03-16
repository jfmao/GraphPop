"""Parse VEP CSQ / SnpEff ANN annotations from a VCF and emit Gene/HAS_CONSEQUENCE CSVs."""

from __future__ import annotations

import csv
import logging
import re
from pathlib import Path

from cyvcf2 import VCF

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# CSV headers (LOAD CSV compatible, neo4j-admin naming convention)
# ---------------------------------------------------------------------------

GENE_HEADER = [
    "geneId:ID(Gene)", ":LABEL", "symbol", "biotype",
]

HAS_CONSEQUENCE_HEADER = [
    ":START_ID(Variant)", ":END_ID(Gene)", ":TYPE",
    "consequence", "impact", "feature", "feature_type",
    "sift_score:float", "sift_pred",
    "polyphen_score:float", "polyphen_pred",
    "cadd_phred:float", "revel:float",
]

# Regex for SIFT/PolyPhen format: "prediction(score)"
_PRED_SCORE_RE = re.compile(r"^(.+?)\(([0-9.]+)\)$")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse_pred_score(value: str) -> tuple[str, str]:
    """Parse 'deleterious(0.02)' → ('deleterious', '0.02').

    Returns ('', '') if value is empty or unparseable.
    """
    if not value:
        return ("", "")
    m = _PRED_SCORE_RE.match(value)
    if m:
        return (m.group(1), m.group(2))
    return (value, "")


def _load_variant_id_set(variant_csv_path: str | Path) -> set[str]:
    """Read the first column of variant_nodes.csv to build a variant ID set."""
    ids: set[str] = set()
    path = Path(variant_csv_path)
    with open(path) as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            if row:
                ids.add(row[0])
    logger.info("Loaded %d variant IDs from %s", len(ids), path)
    return ids


# ---------------------------------------------------------------------------
# Annotation format descriptors
# ---------------------------------------------------------------------------

# Unified field names used internally (mapped from CSQ or ANN)
_FIELD_ALLELE = "allele"
_FIELD_CONSEQUENCE = "consequence"
_FIELD_IMPACT = "impact"
_FIELD_SYMBOL = "symbol"
_FIELD_GENE = "gene"
_FIELD_FEATURE_TYPE = "feature_type"
_FIELD_FEATURE = "feature"
_FIELD_BIOTYPE = "biotype"
_FIELD_SIFT = "sift"
_FIELD_POLYPHEN = "polyphen"
_FIELD_CADD = "cadd"
_FIELD_REVEL = "revel"

# CSQ field name → internal name
_CSQ_FIELD_MAP = {
    "Allele": _FIELD_ALLELE,
    "Consequence": _FIELD_CONSEQUENCE,
    "IMPACT": _FIELD_IMPACT,
    "SYMBOL": _FIELD_SYMBOL,
    "Gene": _FIELD_GENE,
    "Feature_type": _FIELD_FEATURE_TYPE,
    "Feature": _FIELD_FEATURE,
    "BIOTYPE": _FIELD_BIOTYPE,
    "SIFT": _FIELD_SIFT,
    "PolyPhen": _FIELD_POLYPHEN,
    "CADD_PHRED": _FIELD_CADD,
    "REVEL": _FIELD_REVEL,
}

# ANN field index → internal name (SnpEff standard 16-field layout)
_ANN_INDEX_MAP = {
    0: _FIELD_ALLELE,         # Allele
    1: _FIELD_CONSEQUENCE,    # Annotation
    2: _FIELD_IMPACT,         # Annotation_Impact
    3: _FIELD_SYMBOL,         # Gene_Name
    4: _FIELD_GENE,           # Gene_ID
    5: _FIELD_FEATURE_TYPE,   # Feature_Type
    6: _FIELD_FEATURE,        # Feature_ID
    7: _FIELD_BIOTYPE,        # Transcript_BioType
    # Fields 8-15: Rank, HGVS.c, HGVS.p, cDNA, CDS, AA, Distance, Errors
}


# ---------------------------------------------------------------------------
# VEPParser
# ---------------------------------------------------------------------------


class VEPParser:
    """Parse VEP CSQ or SnpEff ANN annotations from a VCF.

    Auto-detects CSQ (VEP/gnomAD) vs ANN (SnpEff) format from VCF header.

    Emits two CSV files:
    - gene_nodes.csv: unique Gene nodes (geneId, symbol, biotype)
    - has_consequence_edges.csv: HAS_CONSEQUENCE edges from Variant→Gene

    Usage::

        # VEP (gnomAD)
        parser = VEPParser("gnomad.chr22_VEP.vcf.gz", "output/")
        parser.run()

        # SnpEff (e.g., rice 3kRG)
        parser = VEPParser(
            "snpeff_annotated.vcf.gz", "output/",
            variant_csv_path="output/variant_nodes.csv",
            chrom_map={"1": "Chr1", "2": "Chr2", ...},
        )
        parser.run()
    """

    def __init__(
        self,
        vep_vcf_path: str | Path,
        out_dir: str | Path,
        *,
        variant_id_set: set[str] | None = None,
        variant_csv_path: str | Path | None = None,
        chrom_map: dict[str, str] | None = None,
    ) -> None:
        """
        Parameters
        ----------
        vep_vcf_path
            Path to the annotation VCF (VEP or SnpEff).
        out_dir
            Output directory for gene_nodes.csv and has_consequence_edges.csv.
        variant_id_set
            Optional set of variant IDs (chr:pos:ref:alt) to filter on.
            Only consequences for these variants are emitted.
        variant_csv_path
            Path to an existing variant_nodes.csv; used to build variant_id_set
            if variant_id_set is not provided.
        chrom_map
            Optional mapping from annotation VCF chromosome names to main VCF
            chromosome names (e.g., {"1": "Chr1"} for rice SnpEff → main VCF).
        """
        self._vcf_path = str(vep_vcf_path)
        self._out_dir = Path(out_dir)
        self._out_dir.mkdir(parents=True, exist_ok=True)
        self._chrom_map = chrom_map or {}

        if variant_id_set is not None:
            self._variant_ids = variant_id_set
        elif variant_csv_path is not None:
            self._variant_ids = _load_variant_id_set(variant_csv_path)
        else:
            self._variant_ids = None

        # Stats
        self._n_variants_seen = 0
        self._n_variants_matched = 0
        self._n_variants_skipped = 0
        self._n_edges = 0
        self._n_genes = 0
        self._ann_format: str = ""  # "CSQ" or "ANN", set in run()

    # -- Public API ---------------------------------------------------------

    @property
    def n_variants_seen(self) -> int:
        return self._n_variants_seen

    @property
    def n_variants_matched(self) -> int:
        return self._n_variants_matched

    @property
    def n_edges(self) -> int:
        return self._n_edges

    @property
    def n_genes(self) -> int:
        return self._n_genes

    @property
    def ann_format(self) -> str:
        return self._ann_format

    def run(self) -> None:
        """Parse the VCF and write gene_nodes.csv + has_consequence_edges.csv."""
        vcf = VCF(self._vcf_path, lazy=True, threads=2)

        # 1. Detect annotation format and build field index map
        info_key, field_indices = self._detect_format(vcf)
        self._ann_format = info_key
        logger.info("Detected annotation format: %s", info_key)

        i_allele = field_indices[_FIELD_ALLELE]
        i_consequence = field_indices[_FIELD_CONSEQUENCE]
        i_impact = field_indices[_FIELD_IMPACT]
        i_symbol = field_indices[_FIELD_SYMBOL]
        i_gene = field_indices[_FIELD_GENE]
        i_feature_type = field_indices[_FIELD_FEATURE_TYPE]
        i_feature = field_indices[_FIELD_FEATURE]
        i_biotype = field_indices[_FIELD_BIOTYPE]
        i_sift = field_indices.get(_FIELD_SIFT)
        i_polyphen = field_indices.get(_FIELD_POLYPHEN)
        i_cadd = field_indices.get(_FIELD_CADD)
        i_revel = field_indices.get(_FIELD_REVEL)

        is_ann = info_key == "ANN"

        # 2. Open output CSVs
        genes: dict[str, tuple[str, str]] = {}  # geneId → (symbol, biotype)

        edge_path = self._out_dir / "has_consequence_edges.csv"
        edge_fh = open(edge_path, "w", newline="")
        edge_writer = csv.writer(edge_fh)
        edge_writer.writerow(HAS_CONSEQUENCE_HEADER)

        try:
            # 3. Stream VCF records
            vcf_iter = iter(vcf)
            while True:
                try:
                    v = next(vcf_iter)
                except StopIteration:
                    break
                except Exception as exc:
                    if "bcf_read" in str(exc) or "htslib" in str(exc):
                        logger.warning(
                            "VCF read error after %d records (truncated file?): %s",
                            self._n_variants_seen, exc,
                        )
                        break
                    raise

                self._n_variants_seen += 1

                chrom = v.CHROM
                pos = v.POS
                ref = v.REF
                alts = v.ALT

                if not alts:
                    continue

                # Get annotation string from INFO
                ann_raw = v.INFO.get(info_key)
                if ann_raw is None:
                    continue

                # Map chrom for variant ID construction
                vid_chrom = self._chrom_map.get(chrom, chrom)

                # Build variant IDs for each ALT allele
                alt_to_vid: dict[str, str] = {}
                for alt in alts:
                    vid = f"{vid_chrom}:{pos}:{ref}:{alt}"
                    alt_to_vid[alt] = vid

                # Filter: if we have a variant_id_set, check if any ALT matches
                if self._variant_ids is not None:
                    matched_alts = {
                        alt for alt, vid in alt_to_vid.items()
                        if vid in self._variant_ids
                    }
                    if not matched_alts:
                        self._n_variants_skipped += 1
                        continue
                else:
                    matched_alts = set(alts)

                self._n_variants_matched += 1

                # 4. Parse annotation entries (comma-separated, pipe-delimited)
                for entry in ann_raw.split(","):
                    fields = entry.split("|")

                    allele = fields[i_allele] if i_allele < len(fields) else ""
                    gene_id = fields[i_gene] if i_gene < len(fields) else ""

                    # Skip entries with no gene ID
                    if not gene_id:
                        continue

                    # Match allele to ALT
                    if is_ann:
                        # SnpEff ANN: Allele = VCF ALT directly
                        variant_id = self._resolve_ann_variant_id(
                            allele, alt_to_vid, matched_alts,
                        )
                    else:
                        # VEP CSQ: Allele uses simplified representation
                        variant_id = self._resolve_variant_id(
                            allele, ref, alts, alt_to_vid, matched_alts,
                        )
                    if variant_id is None:
                        continue

                    # Collect gene info (deduplicate)
                    symbol = fields[i_symbol] if i_symbol < len(fields) else ""
                    biotype = fields[i_biotype] if i_biotype < len(fields) else ""
                    if gene_id not in genes:
                        genes[gene_id] = (symbol, biotype)

                    # Parse optional scores (only present in VEP/CSQ)
                    sift_pred, sift_score = _parse_pred_score(
                        fields[i_sift] if i_sift is not None and i_sift < len(fields) else ""
                    )
                    pp_pred, pp_score = _parse_pred_score(
                        fields[i_polyphen] if i_polyphen is not None and i_polyphen < len(fields) else ""
                    )
                    cadd = fields[i_cadd] if i_cadd is not None and i_cadd < len(fields) else ""
                    revel = fields[i_revel] if i_revel is not None and i_revel < len(fields) else ""

                    consequence = fields[i_consequence] if i_consequence < len(fields) else ""
                    impact = fields[i_impact] if i_impact < len(fields) else ""
                    feature = fields[i_feature] if i_feature < len(fields) else ""
                    feature_type = fields[i_feature_type] if i_feature_type < len(fields) else ""

                    edge_writer.writerow([
                        variant_id,
                        gene_id,
                        "HAS_CONSEQUENCE",
                        consequence,
                        impact,
                        feature,
                        feature_type,
                        sift_score,
                        sift_pred,
                        pp_score,
                        pp_pred,
                        cadd,
                        revel,
                    ])
                    self._n_edges += 1

                # Progress logging
                if self._n_variants_seen % 500_000 == 0:
                    logger.info(
                        "Processed %d records (%d matched, %d edges so far)",
                        self._n_variants_seen,
                        self._n_variants_matched,
                        self._n_edges,
                    )

        finally:
            edge_fh.close()
            vcf.close()

        # 5. Write gene_nodes.csv (deduplicated)
        gene_path = self._out_dir / "gene_nodes.csv"
        with open(gene_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(GENE_HEADER)
            for gene_id, (symbol, biotype) in sorted(genes.items()):
                w.writerow([gene_id, "Gene", symbol, biotype])

        self._n_genes = len(genes)

        logger.info(
            "%s parsing complete: %d VCF records seen, %d matched, "
            "%d skipped, %d genes, %d HAS_CONSEQUENCE edges",
            info_key,
            self._n_variants_seen,
            self._n_variants_matched,
            self._n_variants_skipped,
            self._n_genes,
            self._n_edges,
        )

    # -- Internals ----------------------------------------------------------

    @staticmethod
    def _detect_format(vcf: VCF) -> tuple[str, dict[str, int]]:
        """Auto-detect CSQ (VEP) or ANN (SnpEff) from VCF header.

        Returns (info_key, field_indices) where field_indices maps internal
        field names to their position in the pipe-delimited annotation string.
        """
        csq_desc = None
        ann_desc = None

        for header_item in vcf.header_iter():
            info = header_item.info()
            hid = info.get("ID")
            if hid == "CSQ":
                csq_desc = info.get("Description", "")
            elif hid == "ANN":
                ann_desc = info.get("Description", "")

        # Prefer CSQ if both present
        if csq_desc is not None:
            return ("CSQ", VEPParser._parse_csq_format_from_desc(csq_desc))

        if ann_desc is not None:
            return ("ANN", VEPParser._parse_ann_format(ann_desc))

        raise ValueError("No CSQ or ANN INFO field found in VCF header")

    @staticmethod
    def _parse_csq_format_from_desc(desc: str) -> dict[str, int]:
        """Parse CSQ Description → internal field index map."""
        idx = desc.find("Format: ")
        if idx == -1:
            raise ValueError(
                f"CSQ header found but no 'Format: ' in Description: {desc}"
            )
        fmt_str = desc[idx + len("Format: "):].rstrip('"').rstrip("'")
        csq_fields = fmt_str.split("|")

        field_idx: dict[str, int] = {}
        name_to_pos = {name: i for i, name in enumerate(csq_fields)}

        for csq_name, internal_name in _CSQ_FIELD_MAP.items():
            if csq_name in name_to_pos:
                field_idx[internal_name] = name_to_pos[csq_name]

        # Verify required fields
        for req in (_FIELD_ALLELE, _FIELD_CONSEQUENCE, _FIELD_GENE):
            if req not in field_idx:
                raise ValueError(f"CSQ header missing required field: {req}")

        return field_idx

    @staticmethod
    def _parse_ann_format(desc: str) -> dict[str, int]:
        """Parse ANN Description → internal field index map.

        SnpEff ANN has a fixed 16-field layout. We map the known positions
        to internal field names. SIFT/PolyPhen/CADD/REVEL are absent.
        """
        # The ANN_INDEX_MAP gives us fixed positions
        field_idx: dict[str, int] = {}
        for pos, internal_name in _ANN_INDEX_MAP.items():
            field_idx[internal_name] = pos
        return field_idx

    @staticmethod
    def _parse_csq_format(vcf: VCF) -> list[str]:
        """Extract CSQ field names from VCF header (legacy, kept for compat).

        Looks for a line like:
          ##INFO=<ID=CSQ,...,Description="...Format: Allele|Consequence|IMPACT|...">
        """
        for header_item in vcf.header_iter():
            info = header_item.info()
            if info.get("ID") == "CSQ":
                desc = info.get("Description", "")
                idx = desc.find("Format: ")
                if idx == -1:
                    raise ValueError(
                        f"CSQ header found but no 'Format: ' in Description: {desc}"
                    )
                fmt_str = desc[idx + len("Format: "):].rstrip('"')
                return fmt_str.split("|")
        raise ValueError("No CSQ INFO field found in VCF header")

    @staticmethod
    def _resolve_ann_variant_id(
        ann_allele: str,
        alt_to_vid: dict[str, str],
        matched_alts: set[str],
    ) -> str | None:
        """Map SnpEff ANN Allele field to a variant ID.

        SnpEff uses the VCF ALT allele directly in the Allele field,
        so matching is simpler than VEP.
        """
        if ann_allele in alt_to_vid and ann_allele in matched_alts:
            return alt_to_vid[ann_allele]
        return None

    @staticmethod
    def _resolve_variant_id(
        vep_allele: str,
        ref: str,
        alts: list[str],
        alt_to_vid: dict[str, str],
        matched_alts: set[str],
    ) -> str | None:
        """Map VEP's Allele field back to a variant ID.

        VEP represents alleles differently from VCF:
        - SNP: ALT base (same as VCF)
        - Deletion: '-'
        - Insertion: inserted bases only (no ref padding)

        Returns the variant_id string or None if not matched.
        """
        # Fast path: direct match (works for SNPs and many cases)
        if vep_allele in alt_to_vid and vep_allele in matched_alts:
            return alt_to_vid[vep_allele]

        # For indels, try to match VEP representation to VCF ALT
        for alt in matched_alts:
            if alt not in alt_to_vid:
                continue
            if len(ref) > len(alt):
                # Deletion: VEP uses '-' or the remaining bases
                if vep_allele == "-":
                    return alt_to_vid[alt]
            elif len(alt) > len(ref):
                # Insertion: VEP uses the inserted bases (without ref padding)
                inserted = alt[len(ref):]
                if vep_allele == inserted:
                    return alt_to_vid[alt]

        return None
