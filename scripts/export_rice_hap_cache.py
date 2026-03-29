#!/usr/bin/env python3
"""
One-time export of phased haplotype data from Neo4j to per-chromosome numpy .npz cache
for the Rice 3K dataset.

Encoding mirrors PackedGenotypeReader.java:
  gt_packed   : 2 bits/sample (4/byte). Sample i -> byte i//4, bits (i%4)*2 and (i%4)*2+1.
                0=HOM_REF, 1=HET, 2=HOM_ALT, 3=MISSING
  phase_packed: 1 bit/sample (8/byte). Sample i -> byte i//8, bit i%8.
                0 or 1: which haplotype carries ALT (for HET sites only)

Output per chromosome  (data/rice_hap_cache/{chr}.npz):
  pos          int64   (M,)               variant positions (sorted)
  hap          uint8   (M, ceil(2N/8))    bit-packed haplotype matrix
                                           column 2i   = hap0 for sample i
                                           column 2i+1 = hap1 for sample i
  n_samples    scalar  -- total sample count (N)
  n_haplotypes scalar  -- 2*N

Also writes data/rice_hap_cache/sample_meta.json:
  {n_samples, sample_ids (VCF/packed_index order), pop_ids (population per sample)}

Usage:
    python scripts/export_rice_hap_cache.py                  # all 12 chromosomes
    python scripts/export_rice_hap_cache.py --chr Chr9       # single chromosome (test)
    python scripts/export_rice_hap_cache.py --workers 2      # parallel export
    python scripts/export_rice_hap_cache.py --force           # re-export existing
"""

import argparse
import json
import logging
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from neo4j import GraphDatabase

# -- Constants -----------------------------------------------------------------

NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB   = "neo4j"

CHROMOSOMES = [f"Chr{i}" for i in range(1, 13)]
CACHE_DIR   = "data/rice_hap_cache"

log = logging.getLogger("export_rice_hap")


# -- Packed-array decoding (mirrors PackedGenotypeReader.java) -----------------

def decode_haplotype_row(gt_bytes: bytes, phase_bytes: bytes,
                         idx: np.ndarray) -> np.ndarray:
    """
    Decode one variant's packed genotype + phase arrays into a haplotype row.

    idx: pre-built arange(n_samples, dtype=int32)

    Returns uint8 array of length 2*n_samples:
      [hap0_s0, hap1_s0, hap0_s1, hap1_s1, ...]
    """
    gt_arr = np.frombuffer(gt_bytes, dtype=np.uint8)
    gt = (gt_arr[idx >> 2].astype(np.int32) >> ((idx & 3) << 1)) & 3

    ph_arr = np.frombuffer(phase_bytes, dtype=np.uint8)
    phase  = (ph_arr[idx >> 3].astype(np.int32) >> (idx & 7)) & 1

    hom_alt = gt == 2
    het     = gt == 1

    # hap0=1 when HOM_ALT, or (HET and phase==0)
    # hap1=1 when HOM_ALT, or (HET and phase==1)
    hap0 = (hom_alt | (het & (phase == 0))).astype(np.uint8)
    hap1 = (hom_alt | (het & (phase == 1))).astype(np.uint8)

    row = np.empty(2 * len(idx), dtype=np.uint8)
    row[0::2] = hap0
    row[1::2] = hap1
    return row


# -- Sample metadata -----------------------------------------------------------

def build_sample_meta(driver, cache_dir: str, force: bool = False) -> dict:
    """
    Build sample_meta.json from Neo4j Sample nodes, ordered by packed_index.
    {n_samples, sample_ids, pop_ids}
    """
    meta_path = os.path.join(cache_dir, "sample_meta.json")
    if os.path.exists(meta_path) and not force:
        log.info("sample_meta.json exists, skipping")
        with open(meta_path) as f:
            return json.load(f)

    with driver.session(database=NEO4J_DB) as session:
        result = session.run(
            "MATCH (s:Sample) "
            "RETURN s.sampleId AS sid, s.packed_index AS idx, s.population AS pop "
            "ORDER BY s.packed_index"
        )
        records = list(result)

    sample_ids = [rec["sid"] for rec in records]
    pop_ids    = [rec["pop"] or "UNKNOWN" for rec in records]

    # Verify contiguous packed_index
    indices = [rec["idx"] for rec in records]
    if indices != list(range(len(indices))):
        log.warning("packed_index is not contiguous 0..N-1; gaps may exist")

    meta = {
        "n_samples":  len(sample_ids),
        "sample_ids": sample_ids,
        "pop_ids":    pop_ids,
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f)

    n_pops = len(set(pop_ids) - {"UNKNOWN"})
    log.info("sample_meta.json: %d samples, %d populations", len(sample_ids), n_pops)
    return meta


# -- Per-chromosome export -----------------------------------------------------

def export_chromosome(driver, chrom: str, sample_meta: dict,
                      cache_dir: str, force: bool = False) -> None:
    out_path = os.path.join(cache_dir, f"{chrom}.npz")
    if os.path.exists(out_path) and not force:
        log.info("  %s: cache exists, skipping", chrom)
        return

    n_samples   = sample_meta["n_samples"]
    n_hap_cols  = 2 * n_samples
    n_hap_bytes = (n_hap_cols + 7) // 8
    idx = np.arange(n_samples, dtype=np.int32)

    # Collect into lists (no pre-known variant count for rice)
    log.info("  %s: starting export (%d samples, %d haplotypes)...",
             chrom, n_samples, n_hap_cols)

    pos_list  = []
    rows_list = []
    t0 = time.time()
    v  = 0

    with driver.session(database=NEO4J_DB) as s:
        result = s.run(
            "MATCH (v:Variant) WHERE v.chr = $chr "
            "RETURN v.pos AS pos, v.gt_packed AS gt, v.phase_packed AS ph "
            "ORDER BY v.pos",
            chr=chrom,
        )
        for rec in result:
            pos      = rec["pos"]
            gt_bytes = bytes(rec["gt"])
            ph_bytes = bytes(rec["ph"])

            row = decode_haplotype_row(gt_bytes, ph_bytes, idx)

            pos_list.append(pos)
            rows_list.append(np.packbits(row, bitorder="big"))

            v += 1
            if v % 100_000 == 0:
                log.info("  %s: %d variants decoded (%.0fs)", chrom, v, time.time() - t0)

    if v == 0:
        log.warning("  %s: no variants found, skipping", chrom)
        return

    pos_arr    = np.array(pos_list, dtype=np.int64)
    hap_packed = np.stack(rows_list, axis=0)

    t_decode = time.time() - t0
    log.info("  %s: decoded %d variants in %.1fs, saving...", chrom, v, t_decode)

    np.savez_compressed(
        out_path,
        pos          = pos_arr,
        hap          = hap_packed,
        n_samples    = np.int32(n_samples),
        n_haplotypes = np.int32(n_hap_cols),
    )
    size_mb = os.path.getsize(out_path) / 1024 / 1024
    log.info("  %s: done in %.1fs -> %s (%.0f MB)", chrom, time.time() - t0, out_path, size_mb)


# -- Main ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Export rice haplotype data to numpy cache")
    ap.add_argument("--chr", nargs="+", default=CHROMOSOMES)
    ap.add_argument("--workers", type=int, default=1,
                    help="Parallel workers (default 1; Neo4j I/O bound)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--cache-dir", default=CACHE_DIR)
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("export_rice_hap_cache.log", mode="a"),
        ],
    )

    os.makedirs(args.cache_dir, exist_ok=True)
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    sample_meta = build_sample_meta(driver, args.cache_dir, force=args.force)
    log.info("Exporting %d chromosome(s) with %d worker(s)...",
             len(args.chr), args.workers)

    errors = []
    if args.workers == 1:
        # Sequential — simpler and avoids Neo4j connection pool issues
        for c in args.chr:
            try:
                export_chromosome(driver, c, sample_meta, args.cache_dir, args.force)
            except Exception as e:
                log.error("  %s: FAILED -- %s", c, e)
                errors.append((c, str(e)))
    else:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(export_chromosome, driver, c,
                            sample_meta, args.cache_dir, args.force): c
                for c in args.chr
            }
            for fut in as_completed(futures):
                chrom = futures[fut]
                try:
                    fut.result()
                except Exception as e:
                    log.error("  %s: FAILED -- %s", chrom, e)
                    errors.append((chrom, str(e)))

    driver.close()
    if errors:
        log.error("Export completed with %d error(s):", len(errors))
        for chrom, msg in errors:
            log.error("  %s: %s", chrom, msg)
        sys.exit(1)
    else:
        log.info("Export complete. Cache at: %s", args.cache_dir)


if __name__ == "__main__":
    main()
