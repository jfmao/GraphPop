#!/usr/bin/env python3
"""Investigation 1: Pathway-Level Population Differentiation via Hudson FST.

Question: Which biological pathways are most divergently selected across human populations?

Method:
  Graph traversal: Pathway <-[IN_PATHWAY]- Gene <-[HAS_CONSEQUENCE]- Variant
  Per variant: compute Hudson FST (Bhatia et al. 2013) for key population pairs
  Per pathway:  mean(numerators) / mean(denominators)  (unbiased ratio-of-averages)

Why novel: No VCF tool computes FST at pathway level natively.
  One Cypher traversal replaces: VCF filtering + gene annotation + pathway join + FST tool.

Output: data/results/pathway_fst.tsv — pathways ranked by FST

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_01_pathway_fst.py
"""

import csv
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase
import os

# ── Config ──────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"

OUT_DIR    = Path("data/results")
OUT_TSV    = OUT_DIR / "pathway_fst.tsv"
OUT_JSON   = OUT_DIR / "pathway_fst.json"

# Population pairs to analyse (contrasting continental groups)
# Indices into pop_ids array (ACB=0 ASW=1 BEB=2 CDX=3 CEU=4 CHB=5 CHS=6
# CLM=7 ESN=8 FIN=9 GBR=10 GIH=11 GWD=12 IBS=13 ITU=14 JPT=15 KHV=16
# LWK=17 MSL=18 MXL=19 PEL=20 PJL=21 PUR=22 STU=23 TSI=24 YRI=25)
FOCAL_PAIRS = [
    ("YRI", "CEU"),   # Africa vs Europe
    ("YRI", "CHB"),   # Africa vs East Asia
    ("CEU", "CHB"),   # Europe vs East Asia
    ("YRI", "JPT"),   # Africa vs East Asia (alt)
    ("CEU", "JPT"),   # Europe vs East Asia (alt)
]

MIN_VARIANTS = 5      # minimum variants per pathway to report
MIN_AN       = 10     # minimum allele count total to include a variant


# ── Cypher query ─────────────────────────────────────────────────────────────
# Traverse graph: Pathway → Gene → Variant, pull allele count arrays
CYPHER = """
MATCH (p:Pathway)<-[:IN_PATHWAY]-(g:Gene)<-[:HAS_CONSEQUENCE]-(v:Variant)
WHERE v.ac IS NOT NULL AND v.an IS NOT NULL
RETURN p.pathwayId      AS pathway_id,
       p.name           AS pathway_name,
       v.variantId      AS variant_id,
       v.ac             AS ac,
       v.an             AS an,
       v.pop_ids        AS pop_ids
"""


# ── Hudson FST (Bhatia et al. 2013) ─────────────────────────────────────────

def hudson_fst_components(p1, n1, p2, n2):
    """Return (numerator, denominator) for Hudson FST estimator.

    p1, p2: allele frequencies (float)
    n1, n2: diploid sample sizes (int)

    Numerator   = (p1 - p2)^2  - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)
    Denominator = p1(1-p2) + p2(1-p1)

    Genome-wide FST = sum(num) / sum(denom)   (Bhatia et al. 2013 eq. 10)
    """
    if n1 < 2 or n2 < 2:
        return None, None
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    if den <= 0:
        return None, None
    return num, den


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 1: Pathway-Level Population Differentiation ===\n", flush=True)

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()

    # ── Step 1: Stream Pathway → Gene → Variant data ─────────────────────────
    print("[1/3] Streaming Pathway→Gene→Variant from Neo4j ...", flush=True)
    t0 = time.time()

    # pathway_data[pathway_id] = {name, pair_tag: {nums:[], dens:[]}, n_variants}
    pathway_meta = {}       # pathway_id → name
    # For each pathway and each focal pair: accumulate FST components
    # Structure: {pathway_id: {pair_tag: {'num': [], 'den': []}}}
    pathway_fst = defaultdict(lambda: defaultdict(lambda: {"num": [], "den": []}))
    n_variants_per_pathway = defaultdict(set)  # pathway_id → set of variant_ids

    pop_ids_global = None   # loaded once from first row

    n_rows = 0
    with driver.session(database=NEO4J_DB) as session:
        result = session.run(CYPHER)
        for record in result:
            n_rows += 1
            if n_rows % 50_000 == 0:
                print(f"  {n_rows:,} rows streamed ...", flush=True)

            pid   = record["pathway_id"]
            pname = record["pathway_name"]
            vid   = record["variant_id"]
            ac    = record["ac"]
            an    = record["an"]
            pop_ids = record["pop_ids"]

            if pop_ids_global is None:
                pop_ids_global = list(pop_ids)
                pop_idx = {p: i for i, p in enumerate(pop_ids_global)}
                print(f"  pop_ids: {pop_ids_global}", flush=True)

            pathway_meta[pid] = pname
            n_variants_per_pathway[pid].add(vid)

            # Compute FST components for each focal pair
            for pop_a, pop_b in FOCAL_PAIRS:
                i, j = pop_idx[pop_a], pop_idx[pop_b]
                an_i, an_j = an[i], an[j]
                if an_i < MIN_AN or an_j < MIN_AN:
                    continue
                p_i = ac[i] / an_i
                p_j = ac[j] / an_j
                n_i = an_i // 2   # diploid sample size
                n_j = an_j // 2
                num, den = hudson_fst_components(p_i, n_i, p_j, n_j)
                if num is not None:
                    tag = f"{pop_a}_vs_{pop_b}"
                    pathway_fst[pid][tag]["num"].append(num)
                    pathway_fst[pid][tag]["den"].append(den)

    print(f"  Streamed {n_rows:,} rows in {time.time()-t0:.1f}s", flush=True)
    print(f"  {len(pathway_meta):,} pathways seen", flush=True)

    # ── Step 2: Compute per-pathway FST ──────────────────────────────────────
    print("\n[2/3] Computing per-pathway Hudson FST ...", flush=True)

    results = []
    for pid, pname in pathway_meta.items():
        n_var = len(n_variants_per_pathway[pid])
        if n_var < MIN_VARIANTS:
            continue

        row = {
            "pathway_id":   pid,
            "pathway_name": pname,
            "n_variants":   n_var,
        }
        any_fst = False
        for pop_a, pop_b in FOCAL_PAIRS:
            tag = f"{pop_a}_vs_{pop_b}"
            comps = pathway_fst[pid].get(tag, {})
            nums = comps.get("num", [])
            dens = comps.get("den", [])
            if len(nums) >= MIN_VARIANTS:
                fst = np.sum(nums) / np.sum(dens)
                fst = max(0.0, fst)   # clamp negative (sampling noise)
                row[f"fst_{tag}"]  = round(float(fst), 6)
                row[f"n_{tag}"]    = len(nums)
                any_fst = True
            else:
                row[f"fst_{tag}"] = None
                row[f"n_{tag}"]   = len(nums)

        # Mean FST across all focal pairs (where available)
        fst_vals = [row[f"fst_{pop_a}_vs_{pop_b}"]
                    for pop_a, pop_b in FOCAL_PAIRS
                    if row.get(f"fst_{pop_a}_vs_{pop_b}") is not None]
        row["mean_fst"] = round(float(np.mean(fst_vals)), 6) if fst_vals else None

        if any_fst:
            results.append(row)

    # Sort by mean FST descending
    results.sort(key=lambda r: r["mean_fst"] if r["mean_fst"] is not None else -1,
                 reverse=True)

    # ── Step 3: Write outputs ─────────────────────────────────────────────────
    print(f"\n[3/3] Writing {len(results):,} pathways to {OUT_TSV} ...", flush=True)

    pair_tags = [f"{a}_vs_{b}" for a, b in FOCAL_PAIRS]
    fieldnames = (["pathway_id", "pathway_name", "n_variants", "mean_fst"] +
                  [f"fst_{t}" for t in pair_tags] +
                  [f"n_{t}" for t in pair_tags])

    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)

    driver.close()

    # Print top 20
    print("\n=== Top 20 Pathways by Mean FST ===")
    print(f"{'Rank':<5} {'FST':>7}  {'N_var':>6}  Pathway")
    print("-" * 80)
    for i, r in enumerate(results[:20], 1):
        fst_str = f"{r['mean_fst']:.4f}" if r["mean_fst"] is not None else "   N/A"
        name = r["pathway_name"][:55]
        print(f"  {i:<4} {fst_str}  {r['n_variants']:>6}  {name}")

    print(f"\nFull results: {OUT_TSV}")
    print(f"JSON:         {OUT_JSON}")


if __name__ == "__main__":
    main()
