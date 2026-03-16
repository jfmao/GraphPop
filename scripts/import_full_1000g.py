#!/usr/bin/env python3
"""Full 1000 Genomes Import — Per-chromosome CSV generation + merge.

Processes chr1-22 sequentially (each chromosome's VCF is ~1-3 GB),
generating CSVs into per-chromosome temp directories, then merges
into a single output directory for neo4j-admin import.

Usage:
    conda run -n graphevo python scripts/import_full_1000g.py
    conda run -n graphevo python scripts/import_full_1000g.py --workers 4
    conda run -n graphevo python scripts/import_full_1000g.py --start-chr 5  # resume from chr5
    conda run -n graphevo python scripts/import_full_1000g.py --merge-only   # skip CSV gen, just merge
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import shutil
import sys
import time
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "graphpop-import" / "src"))

from graphpop_import.csv_emitter import CSVEmitter
from graphpop_import.vcf_parser import VCFParser

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(process)d] %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

VCF_DIR = Path("data/raw/1000g/vcf")
VCF_PATTERN = "1kGP_high_coverage_Illumina.chr{chr_num}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
PANEL_PATH = "data/raw/1000g/20130606_g1k_3202_samples_ped_population.txt"
ANCESTRAL_DIR = Path("data/raw/ensembl_ancestor/homo_sapiens_ancestor_GRCh38")
OUT_DIR = Path("/tmp/graphpop_csv")
CHUNK_SIZE = 100_000

CHROMOSOMES = list(range(1, 23))  # chr1-22


def get_vcf_path(chr_num: int) -> Path:
    return VCF_DIR / VCF_PATTERN.format(chr_num=chr_num)


def get_ancestral_fasta(chr_num: int) -> str | None:
    fasta = ANCESTRAL_DIR / f"homo_sapiens_ancestor_{chr_num}.fa"
    if fasta.exists():
        return str(fasta)
    logger.warning("Ancestral FASTA not found: %s", fasta)
    return None


# ---------------------------------------------------------------------------
# Per-chromosome CSV generation
# ---------------------------------------------------------------------------


def process_chromosome(chr_num: int) -> dict:
    """Process a single chromosome and return stats."""
    chrom = f"chr{chr_num}"
    vcf_path = get_vcf_path(chr_num)
    if not vcf_path.exists():
        logger.error("VCF not found: %s", vcf_path)
        return {"chr_num": chr_num, "error": f"VCF not found: {vcf_path}"}

    chrom_dir = OUT_DIR / f"_chr_{chrom}"
    chrom_dir.mkdir(parents=True, exist_ok=True)

    ancestral_fasta = get_ancestral_fasta(chr_num)

    t0 = time.time()
    try:
        parser = VCFParser(
            str(vcf_path),
            PANEL_PATH,
            stratify_by="population",
            ancestral_fasta=ancestral_fasta,
        )
        emitter = CSVEmitter(
            chrom_dir,
            parser.pop_map,
            contig_lengths=parser.contig_lengths,
        )
        # Don't write static nodes per-chromosome (we'll write once in merge)
        for chunk in parser.iter_chunks(CHUNK_SIZE):
            emitter.process_chunk(chunk)
        emitter.finalize()

        elapsed = time.time() - t0
        stats = {
            "chr_num": chr_num,
            "n_variants": emitter.n_variants,
            "n_next": emitter.n_next,
            "elapsed": elapsed,
        }
        logger.info(
            "%s: %d variants, %d NEXT edges in %.1fs (%.1f min)",
            chrom, emitter.n_variants, emitter.n_next, elapsed, elapsed / 60,
        )
        return stats
    except Exception as e:
        elapsed = time.time() - t0
        logger.error("%s: FAILED after %.1fs — %s", chrom, elapsed, e)
        return {"chr_num": chr_num, "error": str(e), "elapsed": elapsed}


# ---------------------------------------------------------------------------
# Merge per-chromosome CSVs
# ---------------------------------------------------------------------------


def merge_csvs() -> None:
    """Merge per-chromosome CSVs into the final output directory."""
    logger.info("Merging per-chromosome CSVs...")

    # Write static nodes once from chr1's parser
    vcf_path = get_vcf_path(1)
    parser = VCFParser(str(vcf_path), PANEL_PATH, stratify_by="population")
    emitter = CSVEmitter(OUT_DIR, parser.pop_map, contig_lengths=parser.contig_lengths)
    emitter.write_static_nodes()
    logger.info("Static nodes written (sample_nodes, population_nodes, in_population_edges)")

    # Merge streaming CSVs: variant_nodes, next_edges, on_chromosome_edges
    for csv_name in ("variant_nodes.csv", "next_edges.csv", "on_chromosome_edges.csv"):
        out_path = OUT_DIR / csv_name
        total_rows = 0
        with open(out_path, "w") as out_f:
            header_written = False
            for chr_num in CHROMOSOMES:
                chrom_dir = OUT_DIR / f"_chr_chr{chr_num}"
                chrom_path = chrom_dir / csv_name
                if not chrom_path.exists():
                    logger.warning("Missing %s for chr%d", csv_name, chr_num)
                    continue
                with open(chrom_path) as in_f:
                    header = in_f.readline()
                    if not header_written:
                        out_f.write(header)
                        header_written = True
                    for line in in_f:
                        out_f.write(line)
                        total_rows += 1
        logger.info("Merged %s: %d rows", csv_name, total_rows)

    # Merge chromosome_nodes.csv (deduplicate)
    chroms_seen: dict[str, list] = {}
    for chr_num in CHROMOSOMES:
        chrom_dir = OUT_DIR / f"_chr_chr{chr_num}"
        chrom_path = chrom_dir / "chromosome_nodes.csv"
        if not chrom_path.exists():
            continue
        with open(chrom_path) as f:
            reader = csv.reader(f)
            header = next(reader)
            for row in reader:
                if row and row[0] not in chroms_seen:
                    chroms_seen[row[0]] = row

    chrom_out = OUT_DIR / "chromosome_nodes.csv"
    with open(chrom_out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chromosomeId:ID(Chromosome)", ":LABEL", "length:long"])
        for cid in sorted(chroms_seen.keys()):
            w.writerow(chroms_seen[cid])
    logger.info("Wrote chromosome_nodes.csv with %d chromosomes", len(chroms_seen))

    # Clean up per-chromosome directories
    for chr_num in CHROMOSOMES:
        chrom_dir = OUT_DIR / f"_chr_chr{chr_num}"
        if chrom_dir.exists():
            shutil.rmtree(chrom_dir)
    logger.info("Cleaned up per-chromosome directories")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Full 1000 Genomes Import — CSV generation + merge"
    )
    ap.add_argument(
        "--workers", type=int, default=2,
        help="Number of parallel workers (default: 2, each chr ~2-6 GB RAM)",
    )
    ap.add_argument(
        "--start-chr", type=int, default=1,
        help="Start from this chromosome number (for resume)",
    )
    ap.add_argument(
        "--end-chr", type=int, default=22,
        help="End at this chromosome number (inclusive)",
    )
    ap.add_argument(
        "--merge-only", action="store_true",
        help="Skip CSV generation, just merge existing per-chr CSVs",
    )
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if not args.merge_only:
        # Validate VCF files exist
        chroms_to_process = list(range(args.start_chr, args.end_chr + 1))
        missing = [c for c in chroms_to_process if not get_vcf_path(c).exists()]
        if missing:
            logger.error(
                "Missing VCFs for chromosomes: %s",
                ", ".join(f"chr{c}" for c in missing),
            )
            logger.error("Run scripts/download_1000g_full.sh first")
            sys.exit(1)

        logger.info(
            "Starting CSV generation: chr%d-chr%d (%d chromosomes, %d workers)",
            args.start_chr, args.end_chr, len(chroms_to_process), args.workers,
        )

        t0 = time.time()

        if args.workers == 1:
            results = [process_chromosome(c) for c in chroms_to_process]
        else:
            with Pool(args.workers) as pool:
                results = pool.map(process_chromosome, chroms_to_process)

        # Report
        elapsed = time.time() - t0
        successes = [r for r in results if "error" not in r]
        failures = [r for r in results if "error" in r]
        total_variants = sum(r.get("n_variants", 0) for r in successes)
        total_next = sum(r.get("n_next", 0) for r in successes)

        logger.info(
            "CSV generation done: %d chr OK, %d failed, %d variants, %.1f min",
            len(successes), len(failures), total_variants, elapsed / 60,
        )

        if failures:
            for f in failures:
                logger.error("  chr%d: %s", f["chr_num"], f["error"])
            sys.exit(1)

        # Print per-chromosome stats
        print("\nChromosome | Variants    | NEXT edges  | Time (min)")
        print("-----------|-------------|-------------|----------")
        for r in results:
            if "error" not in r:
                print(
                    f"chr{r['chr_num']:<7d} | {r['n_variants']:>11,} | "
                    f"{r['n_next']:>11,} | {r['elapsed']/60:.1f}"
                )
        print(f"{'TOTAL':<10s} | {total_variants:>11,} | {total_next:>11,} | {elapsed/60:.1f}")

    # Merge
    t_merge = time.time()
    merge_csvs()
    logger.info("Merge done in %.1fs", time.time() - t_merge)

    # Final summary
    print(f"\nOutput directory: {OUT_DIR}")
    print("\nNext steps:")
    print(f"  sudo bash scripts/neo4j-import.sh {OUT_DIR} neo4j")


if __name__ == "__main__":
    main()
