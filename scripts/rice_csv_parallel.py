#!/usr/bin/env python3
"""Parallel CSV generation for rice 3K — one process per chromosome.

Spawns 12 parallel workers (one per rice chromosome), each running
VCFParser + CSVEmitter on a single chromosome region. Then merges
the per-chromosome CSVs into a single output directory.

Usage:
    conda run -n graphevo python scripts/rice_csv_parallel.py [--workers 6]
"""

from __future__ import annotations

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

VCF_PATH = "data/raw/3kRG_data/NB_final_snp.vcf.gz"
PANEL_PATH = "data/raw/3kRG_data/rice_3k_panel.txt"
OUT_DIR = Path("/tmp/graphpop_rice_csv")
CHUNK_SIZE = 100_000

CHROMOSOMES = [
    "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6",
    "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12",
]


def process_chromosome(chrom: str) -> dict:
    """Process a single chromosome and return stats."""
    chrom_dir = OUT_DIR / f"_chr_{chrom}"
    chrom_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    parser = VCFParser(
        VCF_PATH, PANEL_PATH,
        stratify_by="population",
        region=chrom,
    )
    emitter = CSVEmitter(
        chrom_dir, parser.pop_map,
        contig_lengths=parser.contig_lengths,
    )
    # Don't write static nodes per-chromosome (we'll write once in merge)
    for chunk in parser.iter_chunks(CHUNK_SIZE):
        emitter.process_chunk(chunk)
    emitter.finalize()

    elapsed = time.time() - t0
    stats = {
        "chrom": chrom,
        "n_variants": emitter.n_variants,
        "n_next": emitter.n_next,
        "elapsed": elapsed,
    }
    logger.info(
        "%s: %d variants, %d NEXT edges in %.1fs",
        chrom, emitter.n_variants, emitter.n_next, elapsed,
    )
    return stats


def merge_csvs() -> None:
    """Merge per-chromosome CSVs into the final output directory."""
    logger.info("Merging per-chromosome CSVs...")

    # Write static nodes once from any chromosome's parser
    parser = VCFParser(VCF_PATH, PANEL_PATH, stratify_by="population", region="Chr1")
    emitter = CSVEmitter(OUT_DIR, parser.pop_map, contig_lengths=parser.contig_lengths)
    emitter.write_static_nodes()

    # Merge streaming CSVs: variant_nodes, next_edges, on_chromosome_edges
    for csv_name in ("variant_nodes.csv", "next_edges.csv", "on_chromosome_edges.csv"):
        out_path = OUT_DIR / csv_name
        with open(out_path, "w") as out_f:
            header_written = False
            for chrom in CHROMOSOMES:
                chrom_path = OUT_DIR / f"_chr_{chrom}" / csv_name
                if not chrom_path.exists():
                    continue
                with open(chrom_path) as in_f:
                    header = in_f.readline()
                    if not header_written:
                        out_f.write(header)
                        header_written = True
                    # Stream the rest
                    for line in in_f:
                        out_f.write(line)
        logger.info("Merged %s", csv_name)

    # Merge chromosome_nodes.csv (deduplicate)
    chroms_seen: dict[str, str] = {}  # chromosomeId → full row
    for chrom in CHROMOSOMES:
        chrom_path = OUT_DIR / f"_chr_{chrom}" / "chromosome_nodes.csv"
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
    for chrom in CHROMOSOMES:
        chrom_dir = OUT_DIR / f"_chr_{chrom}"
        if chrom_dir.exists():
            shutil.rmtree(chrom_dir)
    logger.info("Cleaned up per-chromosome directories")


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6,
                    help="Number of parallel workers (default: 6)")
    args = ap.parse_args()

    # Clean output directory
    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    OUT_DIR.mkdir(parents=True)

    t0 = time.time()

    logger.info("Starting parallel CSV generation: %d chromosomes, %d workers",
                len(CHROMOSOMES), args.workers)

    with Pool(args.workers) as pool:
        results = pool.map(process_chromosome, CHROMOSOMES)

    total_variants = sum(r["n_variants"] for r in results)
    total_next = sum(r["n_next"] for r in results)
    parallel_time = time.time() - t0
    logger.info(
        "Parallel phase done: %d variants, %d NEXT edges in %.1fs (%.1f min)",
        total_variants, total_next, parallel_time, parallel_time / 60,
    )

    merge_csvs()

    total_time = time.time() - t0
    logger.info("=== Complete: %d variants in %.1fs (%.1f min) ===",
                total_variants, total_time, total_time / 60)

    # Print per-chromosome stats
    print("\nChromosome | Variants | NEXT edges | Time (s)")
    print("-----------|----------|------------|--------")
    for r in results:
        print(f"{r['chrom']:10s} | {r['n_variants']:>8,} | {r['n_next']:>10,} | {r['elapsed']:.1f}")
    print(f"{'TOTAL':10s} | {total_variants:>8,} | {total_next:>10,} | {total_time:.1f}")


if __name__ == "__main__":
    main()
