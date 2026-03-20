#!/usr/bin/env python3
"""Annotate human 1000G VCFs with functional consequences via SnpEff.

Generates gene_nodes.csv and has_consequence_edges.csv for all 22 human
chromosomes to enable annotation-conditioned analyses in GraphPop.

Pipeline per chromosome:
  bcftools view -G (strip genotypes) → snpEff ann → parse ANN → accumulate

Usage:
    # Prerequisites: conda install -c bioconda snpeff && snpEff download GRCh38.p14

    # Annotate all chromosomes (~2-3 hours):
    python scripts/annotate_human_1000g.py

    # Annotate specific chromosomes:
    python scripts/annotate_human_1000g.py --chromosomes 22

    # Then load into Neo4j:
    bash scripts/load-annotations.sh data/raw/1000g/full_genome_annotations/
"""

import argparse
import logging
import os
import subprocess
import sys
import time

log = logging.getLogger("annotate_human")

# ── Constants ─────────────────────────────────────────────────────────

VCF_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "1000g", "vcf")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "1000g", "full_genome_annotations")

CHROMOSOMES = [str(i) for i in range(1, 23)]
VCF_PATTERN = "1kGP_high_coverage_Illumina.chr{chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

SNPEFF_GENOME = "GRCh38.p14"
SNPEFF_JAVA_OPTS = "-Xmx8g"

# Variant ID format in our database: chr22:pos:ref:alt
# 1000G VCF uses chr22 chromosome names, so we match directly


def annotate_chromosome(chrom, vcf_dir, genes, edges, progress_file):
    """Annotate one chromosome: bcftools view -G | snpEff ann | parse ANN."""
    vcf = os.path.join(vcf_dir, VCF_PATTERN.format(chr=chrom))
    if not os.path.exists(vcf):
        log.warning(f"VCF not found: {vcf}")
        return 0, 0

    log.info(f"--- Chromosome {chrom} ---")
    t0 = time.time()

    # Pipeline: strip genotypes → annotate → parse
    env = os.environ.copy()
    env["_JAVA_OPTIONS"] = SNPEFF_JAVA_OPTS

    # bcftools strips genotype columns (VCF goes from ~3KB/line to ~200B/line)
    bcftools = subprocess.Popen(
        ["bcftools", "view", "-G", vcf],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )

    # SnpEff annotates (only needs CHROM/POS/REF/ALT/INFO)
    snpeff = subprocess.Popen(
        ["snpEff", "ann", "-noStats", "-noLog", SNPEFF_GENOME],
        stdin=bcftools.stdout, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, env=env,
    )
    bcftools.stdout.close()  # allow bcftools to get SIGPIPE

    n_variants = 0
    n_annotated = 0
    n_edges_before = len(edges)

    # Parse streaming output
    for raw_line in snpeff.stdout:
        line = raw_line.decode("utf-8", errors="replace")
        if line.startswith("#"):
            continue

        fields = line.split("\t", 8)
        if len(fields) < 8:
            continue

        chrom_name = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt_field = fields[4]
        info = fields[7]

        # Handle multi-allelic: take first ALT only
        alt = alt_field.split(",")[0]

        # Build variant ID matching our database format
        variant_id = f"{chrom_name}:{pos}:{ref}:{alt}"

        n_variants += 1

        # Extract ANN field
        ann_start = info.find("ANN=")
        if ann_start == -1:
            continue

        ann_end = info.find(";", ann_start)
        ann_str = info[ann_start + 4:ann_end] if ann_end != -1 else info[ann_start + 4:].rstrip("\n")

        n_annotated += 1
        seen = set()

        for entry in ann_str.split(","):
            parts = entry.split("|")
            if len(parts) < 8:
                continue

            # ANN: ALT|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|...
            consequence = parts[1]
            impact = parts[2]
            gene_symbol = parts[3]
            gene_id = parts[4]
            feature_type = parts[5]
            feature_id = parts[6]
            biotype = parts[7] if len(parts) > 7 else ""

            if not gene_id or not consequence:
                continue

            # Deduplicate per variant
            edge_key = (variant_id, gene_id, consequence)
            if edge_key in seen:
                continue
            seen.add(edge_key)

            # Register gene
            if gene_id not in genes:
                genes[gene_id] = {
                    "symbol": gene_symbol,
                    "biotype": biotype or "unknown",
                }

            edges.append(
                f"{variant_id},{gene_id},HAS_CONSEQUENCE,"
                f"{consequence},{impact},{feature_id},{feature_type},"
                f",,,,,\n"
            )

    snpeff.wait()
    bcftools.wait()

    elapsed = time.time() - t0
    n_new_edges = len(edges) - n_edges_before
    log.info(
        f"  chr{chrom}: {n_variants:,} variants, {n_annotated:,} annotated, "
        f"{n_new_edges:,} edges, {len(genes):,} genes total ({elapsed:.0f}s)"
    )

    # Save progress
    with open(progress_file, "a") as f:
        f.write(f"chr{chrom}\t{n_variants}\t{n_annotated}\t{n_new_edges}\t{elapsed:.1f}\n")

    return n_variants, n_annotated


def write_csvs(genes, edges, output_dir):
    """Write gene_nodes.csv and has_consequence_edges.csv."""
    os.makedirs(output_dir, exist_ok=True)

    # Gene nodes
    gene_path = os.path.join(output_dir, "gene_nodes.csv")
    with open(gene_path, "w") as f:
        f.write("geneId:ID(Gene),:LABEL,symbol,biotype\n")
        for gene_id in sorted(genes.keys()):
            info = genes[gene_id]
            # Escape commas in gene symbols
            symbol = info["symbol"].replace(",", "")
            f.write(f"{gene_id},Gene,{symbol},{info['biotype']}\n")
    log.info(f"Wrote {len(genes):,} genes to {gene_path}")

    # HAS_CONSEQUENCE edges (already pre-formatted as CSV strings)
    edge_path = os.path.join(output_dir, "has_consequence_edges.csv")
    with open(edge_path, "w") as f:
        f.write(
            ":START_ID(Variant),:END_ID(Gene),:TYPE,consequence,impact,"
            "feature,feature_type,sift_score:float,sift_pred,"
            "polyphen_score:float,polyphen_pred,cadd_phred:float,revel:float\n"
        )
        for edge_line in edges:
            f.write(edge_line)
    log.info(f"Wrote {len(edges):,} edges to {edge_path}")


def load_completed_chromosomes(progress_file):
    """Load already-completed chromosomes from progress file."""
    completed = set()
    if os.path.exists(progress_file):
        with open(progress_file) as f:
            for line in f:
                parts = line.strip().split("\t")
                if parts:
                    completed.add(parts[0])
    return completed


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(description="Annotate human 1000G VCFs with SnpEff")
    parser.add_argument("--chromosomes", nargs="+", default=CHROMOSOMES,
                        help="Chromosomes to annotate (default: 1-22)")
    parser.add_argument("--vcf-dir", default=VCF_DIR,
                        help="Directory with 1000G VCFs")
    parser.add_argument("--output-dir", default=OUTPUT_DIR,
                        help="Output directory for annotation CSVs")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    progress_file = os.path.join(args.output_dir, "progress.tsv")

    # Check prerequisites
    for cmd in ["bcftools", "snpEff"]:
        try:
            subprocess.run([cmd, "--help"], capture_output=True, timeout=5)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            if cmd == "snpEff":
                # SnpEff exits with error for --help, that's fine
                pass
            else:
                log.error(f"Required tool not found: {cmd}")
                sys.exit(1)

    # Resume support: check which chromosomes are already done
    completed = load_completed_chromosomes(progress_file)
    if completed:
        log.info(f"Resuming: {len(completed)} chromosomes already completed: {sorted(completed)}")

    genes = {}
    edges = []
    total_variants = 0
    total_annotated = 0
    t_start = time.time()

    for chrom in args.chromosomes:
        chrom_key = f"chr{chrom}"
        if chrom_key in completed:
            log.info(f"Skipping {chrom_key} (already completed)")
            continue

        nv, na = annotate_chromosome(chrom, args.vcf_dir, genes, edges, progress_file)
        total_variants += nv
        total_annotated += na

    elapsed = time.time() - t_start

    if not edges and not genes:
        log.info("No new annotations generated (all chromosomes already completed?)")
        # Try to load existing CSVs
        gene_path = os.path.join(args.output_dir, "gene_nodes.csv")
        edge_path = os.path.join(args.output_dir, "has_consequence_edges.csv")
        if os.path.exists(gene_path) and os.path.exists(edge_path):
            log.info(f"Existing CSVs found at {args.output_dir}")
        return

    # Write CSVs
    write_csvs(genes, edges, args.output_dir)

    log.info(f"\n=== Summary ===")
    log.info(f"Total variants:  {total_variants:,}")
    log.info(f"Total annotated: {total_annotated:,}")
    log.info(f"Total genes:     {len(genes):,}")
    log.info(f"Total edges:     {len(edges):,}")
    log.info(f"Total time:      {elapsed:.0f}s ({elapsed/60:.1f} min)")
    log.info(f"Output:          {args.output_dir}")
    log.info(f"\nNext steps:")
    log.info(f"  bash scripts/load-annotations.sh {args.output_dir}")


if __name__ == "__main__":
    main()
