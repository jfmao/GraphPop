"""Parse VEP CSQ annotations from a gnomAD VCF and emit Gene/HAS_CONSEQUENCE CSVs."""

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
# VEPParser
# ---------------------------------------------------------------------------


class VEPParser:
    """Parse VEP CSQ annotations from a gnomAD sites-only VCF.

    Emits two CSV files:
    - gene_nodes.csv: unique Gene nodes (geneId, symbol, biotype)
    - has_consequence_edges.csv: HAS_CONSEQUENCE edges from Variant→Gene

    Usage::

        parser = VEPParser("gnomad.chr22_VEP.vcf.gz", "output/")
        parser.run()
    """

    def __init__(
        self,
        vep_vcf_path: str | Path,
        out_dir: str | Path,
        *,
        variant_id_set: set[str] | None = None,
        variant_csv_path: str | Path | None = None,
    ) -> None:
        """
        Parameters
        ----------
        vep_vcf_path
            Path to the gnomAD VEP-annotated VCF (sites-only).
        out_dir
            Output directory for gene_nodes.csv and has_consequence_edges.csv.
        variant_id_set
            Optional set of variant IDs (chr:pos:ref:alt) to filter on.
            Only consequences for these variants are emitted.
        variant_csv_path
            Path to an existing variant_nodes.csv; used to build variant_id_set
            if variant_id_set is not provided.
        """
        self._vcf_path = str(vep_vcf_path)
        self._out_dir = Path(out_dir)
        self._out_dir.mkdir(parents=True, exist_ok=True)

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

    def run(self) -> None:
        """Parse the VEP VCF and write gene_nodes.csv + has_consequence_edges.csv."""
        vcf = VCF(self._vcf_path, lazy=True, threads=2)

        # 1. Parse CSQ format from header
        csq_fields = self._parse_csq_format(vcf)
        logger.info("CSQ format: %d fields", len(csq_fields))

        # Build field name → index map
        field_idx = {name: i for i, name in enumerate(csq_fields)}

        # Required field indices
        i_allele = field_idx["Allele"]
        i_consequence = field_idx["Consequence"]
        i_impact = field_idx["IMPACT"]
        i_symbol = field_idx["SYMBOL"]
        i_gene = field_idx["Gene"]
        i_feature_type = field_idx["Feature_type"]
        i_feature = field_idx["Feature"]
        i_biotype = field_idx["BIOTYPE"]
        i_sift = field_idx.get("SIFT")
        i_polyphen = field_idx.get("PolyPhen")
        i_cadd = field_idx.get("CADD_PHRED")
        i_revel = field_idx.get("REVEL")

        n_fields = len(csq_fields)

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

                # Get CSQ string from INFO
                csq_raw = v.INFO.get("CSQ")
                if csq_raw is None:
                    continue

                # Build variant IDs for each ALT allele
                alt_to_vid: dict[str, str] = {}
                for alt in alts:
                    vid = f"{chrom}:{pos}:{ref}:{alt}"
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

                # 4. Parse CSQ entries (comma-separated, pipe-delimited)
                for entry in csq_raw.split(","):
                    fields = entry.split("|")
                    if len(fields) < n_fields:
                        # Pad with empty strings if entry is short
                        fields.extend([""] * (n_fields - len(fields)))

                    allele = fields[i_allele]
                    gene_id = fields[i_gene]

                    # Skip entries with no gene ID
                    if not gene_id:
                        continue

                    # Match allele to ALT — VEP Allele field uses the ALT
                    # For SNPs this is the single base; for indels it can differ.
                    # We match by checking if the allele corresponds to a matched ALT.
                    # VEP uses a simplified allele representation: for SNPs it's
                    # the ALT base, for deletions it's '-', for insertions it's
                    # the inserted sequence. We need to map back to VCF ALT.
                    variant_id = self._resolve_variant_id(
                        allele, ref, alts, alt_to_vid, matched_alts,
                    )
                    if variant_id is None:
                        continue

                    # Collect gene info (deduplicate)
                    symbol = fields[i_symbol]
                    biotype = fields[i_biotype]
                    if gene_id not in genes:
                        genes[gene_id] = (symbol, biotype)

                    # Parse optional scores
                    sift_pred, sift_score = _parse_pred_score(
                        fields[i_sift] if i_sift is not None else ""
                    )
                    pp_pred, pp_score = _parse_pred_score(
                        fields[i_polyphen] if i_polyphen is not None else ""
                    )
                    cadd = fields[i_cadd] if i_cadd is not None else ""
                    revel = fields[i_revel] if i_revel is not None else ""

                    edge_writer.writerow([
                        variant_id,
                        gene_id,
                        "HAS_CONSEQUENCE",
                        fields[i_consequence],
                        fields[i_impact],
                        fields[i_feature],
                        fields[i_feature_type],
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
                        "Processed %d VEP records (%d matched, %d edges so far)",
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
            "VEP parsing complete: %d VCF records seen, %d matched, "
            "%d skipped, %d genes, %d HAS_CONSEQUENCE edges",
            self._n_variants_seen,
            self._n_variants_matched,
            self._n_variants_skipped,
            self._n_genes,
            self._n_edges,
        )

    # -- Internals ----------------------------------------------------------

    @staticmethod
    def _parse_csq_format(vcf: VCF) -> list[str]:
        """Extract CSQ field names from VCF header.

        Looks for a line like:
          ##INFO=<ID=CSQ,...,Description="...Format: Allele|Consequence|IMPACT|...">
        """
        for header_item in vcf.header_iter():
            info = header_item.info()
            if info.get("ID") == "CSQ":
                desc = info.get("Description", "")
                # Extract format string after "Format: "
                idx = desc.find("Format: ")
                if idx == -1:
                    raise ValueError(
                        f"CSQ header found but no 'Format: ' in Description: {desc}"
                    )
                fmt_str = desc[idx + len("Format: "):].rstrip('"')
                return fmt_str.split("|")
        raise ValueError("No CSQ INFO field found in VCF header")

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
