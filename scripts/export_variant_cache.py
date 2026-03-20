#!/usr/bin/env python3
"""
One-time export of variant data from Neo4j to per-chromosome numpy .npz cache.

For each chromosome, exports:
  - pos, ac, an, af, het_count, hom_alt_count (full 26-population arrays)
  - ancestral_allele and variant_type as int8 codes
  - has_missense / has_synonymous boolean masks (from HAS_CONSEQUENCE graph join)
  - pop_ids ordering

Also writes data/variants_cache/pop_meta.json with a_n, a_n2, n_samples per population.

Usage:
    python scripts/export_variant_cache.py               # all 22 chromosomes
    python scripts/export_variant_cache.py --chr chr1 chr22
    python scripts/export_variant_cache.py --workers 4   # parallel chromosome export
    python scripts/export_variant_cache.py --force       # re-export even if cache exists
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

# ── Constants ──────────────────────────────────────────────────────────────────

NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB   = "neo4j"

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]
CACHE_DIR   = "data/variants_cache"

CONSEQUENCE_TYPES = ["missense_variant", "synonymous_variant"]

# Encoding maps
ANC_MAP  = {None: 0, "REF": 1, "ALT": 2}   # int8: 0=unknown,1=REF,2=ALT
VTYPE_MAP = {None: 0, "SNP": 1, "INDEL": 2} # int8: 0=unknown,1=SNP,2=INDEL

log = logging.getLogger("export_cache")


# ── Population metadata ────────────────────────────────────────────────────────

def export_pop_meta(driver, cache_dir):
    """Fetch a_n, a_n2, n_samples from Population nodes; also resolve pop_ids order."""
    meta_path = os.path.join(cache_dir, "pop_meta.json")

    with driver.session(database=NEO4J_DB) as s:
        pop_recs = s.run(
            "MATCH (p:Population) "
            "RETURN p.populationId AS pid, p.a_n AS a_n, "
            "       p.a_n2 AS a_n2, p.n_samples AS n_samples "
            "ORDER BY p.populationId"
        ).data()

    if not pop_recs:
        raise RuntimeError("No Population nodes found in database.")

    meta = {r["pid"]: {"a_n": float(r["a_n"]),
                        "a_n2": float(r["a_n2"]),
                        "n_samples": int(r["n_samples"])}
            for r in pop_recs}

    # Resolve the pop_ids ordering from an actual Variant node (first on chr1)
    with driver.session(database=NEO4J_DB) as s:
        row = s.run(
            "MATCH (v:Variant) WHERE v.chr='chr1' "
            "RETURN v.pop_ids AS pop_ids LIMIT 1"
        ).single()
    if row is None:
        raise RuntimeError("No Variant nodes found on chr1.")
    meta["pop_ids"] = list(row["pop_ids"])

    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    log.info("Saved pop_meta.json  (%d populations, order: %s...)",
             len(meta["pop_ids"]), meta["pop_ids"][:4])
    return meta


# ── Per-chromosome export ──────────────────────────────────────────────────────

def export_chromosome(driver, chrom, cache_dir, force=False):
    out_path = os.path.join(cache_dir, f"{chrom}.npz")
    if os.path.exists(out_path) and not force:
        log.info("  %s: cache exists, skipping (use --force to re-export)", chrom)
        return

    t0 = time.time()

    # ── Phase A: Stream all variants sorted by position into numpy arrays ────────
    # Use cursor iteration (not .data()) to avoid building millions of Python dicts
    log.info("  %s: streaming variants...", chrom)
    pos_list = []
    ac_list  = []
    an_list  = []
    af_list  = []
    het_list = []
    hom_list = []
    anc_list = []
    vt_list  = []

    with driver.session(database=NEO4J_DB) as s:
        result = s.run(
            "MATCH (v:Variant) WHERE v.chr = $chr "
            "RETURN v.pos AS pos, "
            "       v.ac AS ac, v.an AS an, v.af AS af, "
            "       v.het_count AS het, v.hom_alt_count AS hom_alt, "
            "       v.ancestral_allele AS anc, v.variant_type AS vtype "
            "ORDER BY v.pos",
            chr=chrom
        )
        for rec in result:          # stream one record at a time — avoids .data() bulk
            pos_list.append(rec["pos"])
            ac_list.append(rec["ac"])
            an_list.append(rec["an"])
            af_list.append(rec["af"])
            het_list.append(rec["het"])
            hom_list.append(rec["hom_alt"])
            anc_list.append(ANC_MAP.get(rec["anc"], 0))
            vt_list.append(VTYPE_MAP.get(rec["vtype"], 0))

    M = len(pos_list)
    if M == 0:
        log.warning("  %s: no variants found, skipping", chrom)
        return

    t_load = time.time() - t0
    log.info("  %s: %d variants streamed in %.1fs, building arrays...", chrom, M, t_load)

    # Convert to numpy — one array at a time so peak RAM stays low
    pos       = np.array(pos_list,  dtype=np.int64);   del pos_list
    ac        = np.array(ac_list,   dtype=np.int32);   del ac_list   # (M,26)
    an        = np.array(an_list,   dtype=np.int32);   del an_list
    af        = np.array(af_list,   dtype=np.float32); del af_list   # float32: ~600 MB saved vs float64
    het       = np.array(het_list,  dtype=np.int32);   del het_list
    hom_alt   = np.array(hom_list,  dtype=np.int32);   del hom_list
    anc_code  = np.array(anc_list,  dtype=np.int8);    del anc_list
    vtype_code= np.array(vt_list,   dtype=np.int8);    del vt_list

    # ── Phase B: Consequence masks (graph join) ────────────────────────────────
    log.info("  %s: querying consequence masks...", chrom)
    with driver.session(database=NEO4J_DB) as s:
        csq_recs = s.run(
            "MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene) "
            "WHERE v.chr = $chr AND c.consequence IN $csqs "
            "RETURN DISTINCT v.pos AS pos, c.consequence AS csq",
            chr=chrom, csqs=CONSEQUENCE_TYPES
        ).data()

    pos_to_idx = {int(p): i for i, p in enumerate(pos.tolist())}
    has_missense   = np.zeros(M, dtype=bool)
    has_synonymous = np.zeros(M, dtype=bool)
    for r in csq_recs:
        idx = pos_to_idx.get(int(r["pos"]))
        if idx is not None:
            if r["csq"] == "missense_variant":
                has_missense[idx] = True
            elif r["csq"] == "synonymous_variant":
                has_synonymous[idx] = True
    del csq_recs, pos_to_idx

    n_mis = int(has_missense.sum())
    n_syn = int(has_synonymous.sum())
    log.info("  %s: %d missense, %d synonymous sites", chrom, n_mis, n_syn)

    # ── Phase C: Save compressed .npz ─────────────────────────────────────────
    np.savez_compressed(
        out_path,
        pos            = pos,
        ac             = ac,
        an             = an,
        af             = af,
        het_count      = het,
        hom_alt_count  = hom_alt,
        ancestral_allele = anc_code,
        variant_type   = vtype_code,
        has_missense   = has_missense,
        has_synonymous = has_synonymous,
    )

    size_mb = os.path.getsize(out_path) / 1024 / 1024
    log.info("  %s: done in %.1fs → %s (%.0f MB, %d variants)",
             chrom, time.time() - t0, out_path, size_mb, M)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Export variant data to numpy cache")
    ap.add_argument("--chr", nargs="+", default=CHROMOSOMES,
                    help="Chromosomes to export (default: all 22)")
    ap.add_argument("--workers", type=int, default=4,
                    help="Parallel chromosome export workers (default: 4)")
    ap.add_argument("--force", action="store_true",
                    help="Re-export even if .npz already exists")
    ap.add_argument("--cache-dir", default=CACHE_DIR,
                    help=f"Output directory (default: {CACHE_DIR})")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler("export_variant_cache.log", mode="a")],
    )

    os.makedirs(args.cache_dir, exist_ok=True)
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    # Population metadata (fast, single pass)
    meta_path = os.path.join(args.cache_dir, "pop_meta.json")
    if not os.path.exists(meta_path) or args.force:
        export_pop_meta(driver, args.cache_dir)
    else:
        log.info("pop_meta.json exists, skipping (use --force to re-export)")

    # Chromosome export — parallel across chromosomes (all READ-only queries)
    log.info("Exporting %d chromosome(s) with %d worker(s)...",
             len(args.chr), args.workers)

    errors = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(export_chromosome, driver, c, args.cache_dir, args.force): c
            for c in args.chr
        }
        for fut in as_completed(futures):
            chrom = futures[fut]
            try:
                fut.result()
            except Exception as e:
                log.error("  %s: FAILED — %s", chrom, e)
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
