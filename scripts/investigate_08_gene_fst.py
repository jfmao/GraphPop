#!/usr/bin/env python3
"""Investigation 8: Gene-Level Population Differentiation via Hudson FST.

Question: Which individual genes show the strongest population differentiation?
Pathway-level FST (Inv.1) averages across many genes; here we resolve to single genes.

Method:
  Stream Variant→Gene (HAS_CONSEQUENCE) with allele count arrays.
  Compute Hudson FST per gene per focal population pair (ratio-of-averages).
  Rank genes by mean FST across all pairs.
  Cross-annotate with convergent sweep signal (Inv.2) and GO terms.

Why novel: Gene-level FST natively from a graph traversal. No VCF tool does this
  without a separate annotation step + custom script per gene.

Output: data/results/gene_fst.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_08_gene_fst.py
"""

import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase
import os

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"

OUT_DIR = Path("data/results")
OUT_TSV = OUT_DIR / "gene_fst.tsv"
OUT_JSON = OUT_DIR / "gene_fst.json"

FOCAL_PAIRS = [
    ("YRI", "CEU"),
    ("YRI", "CHB"),
    ("CEU", "CHB"),
    ("YRI", "GIH"),
    ("CEU", "GIH"),
]
MIN_VARIANTS = 3     # minimum variants per gene per pair
MIN_AN       = 10    # minimum allele count to use variant

# ── Hudson FST ────────────────────────────────────────────────────────────────

def hudson_fst_components(p1, n1, p2, n2):
    if n1 < 2 or n2 < 2:
        return None, None
    num = (p1 - p2)**2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    den = p1*(1-p2) + p2*(1-p1)
    if den <= 0:
        return None, None
    return num, den


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 8: Gene-Level Hudson FST ===\n", flush=True)

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()

    # ── Load cross-reference data ─────────────────────────────────────────────
    conv_genes = {}
    conv_json = OUT_DIR / "convergent_sweeps.json"
    if conv_json.exists():
        with open(conv_json) as f:
            conv_genes = json.load(f)["convergent_genes"]

    go_enriched_genes = set()
    go_json = OUT_DIR / "go_enrichment_summary.json"
    if go_json.exists():
        with open(go_json) as f:
            go_data = json.load(f)
        for pop_data in go_data["populations"].values():
            for entry in pop_data.get("top5", []):
                go_enriched_genes.add(entry.get("go_id", ""))  # approximate

    # ── Stream Variant→Gene from Neo4j ────────────────────────────────────────
    print("[1/3] Streaming Variant→Gene pairs from Neo4j ...", flush=True)
    t0 = time.time()

    # gene → pair_tag → {num:[], den:[]}
    gene_fst  = defaultdict(lambda: defaultdict(lambda: {"num": [], "den": []}))
    gene_meta = {}   # gene_id → {n_variants, chr}
    pop_ids_global = None
    pop_idx        = {}

    n_rows = 0
    with driver.session(database=NEO4J_DB) as session:
        result = session.run("""
            MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
            WHERE v.ac IS NOT NULL AND v.an IS NOT NULL
            RETURN g.geneId AS gene_id,
                   v.variantId AS vid,
                   v.ac AS ac, v.an AS an, v.pop_ids AS pop_ids
        """)
        for rec in result:
            n_rows += 1
            if n_rows % 200_000 == 0:
                print(f"  {n_rows:,} rows ...", flush=True)

            gene_id = rec["gene_id"]
            ac      = rec["ac"]
            an      = rec["an"]
            pids    = rec["pop_ids"]

            if pop_ids_global is None and pids:
                pop_ids_global = list(pids)
                pop_idx = {p: i for i, p in enumerate(pop_ids_global)}
                print(f"  pop_ids: {pop_ids_global}", flush=True)

            if gene_id not in gene_meta:
                gene_meta[gene_id] = {"n_variants": 0}
            gene_meta[gene_id]["n_variants"] += 1

            for pop_a, pop_b in FOCAL_PAIRS:
                i = pop_idx.get(pop_a)
                j = pop_idx.get(pop_b)
                if i is None or j is None:
                    continue
                if an[i] < MIN_AN or an[j] < MIN_AN:
                    continue
                p_i = ac[i] / an[i]
                p_j = ac[j] / an[j]
                n_i, n_j = an[i] // 2, an[j] // 2
                num, den = hudson_fst_components(p_i, n_i, p_j, n_j)
                if num is not None:
                    tag = f"{pop_a}_vs_{pop_b}"
                    gene_fst[gene_id][tag]["num"].append(num)
                    gene_fst[gene_id][tag]["den"].append(den)

    print(f"  {n_rows:,} rows streamed in {time.time()-t0:.1f}s", flush=True)
    print(f"  {len(gene_meta):,} genes", flush=True)

    # ── Compute per-gene FST ──────────────────────────────────────────────────
    print("\n[2/3] Computing per-gene Hudson FST ...", flush=True)

    results = []
    for gene_id, meta in gene_meta.items():
        row = {"gene_id": gene_id, "n_variants": meta["n_variants"]}
        fst_vals = []
        for pop_a, pop_b in FOCAL_PAIRS:
            tag   = f"{pop_a}_vs_{pop_b}"
            comps = gene_fst[gene_id].get(tag, {})
            nums  = comps.get("num", [])
            dens  = comps.get("den", [])
            if len(nums) >= MIN_VARIANTS:
                fst = max(0.0, float(np.sum(nums) / np.sum(dens)))
                row[f"fst_{tag}"] = round(fst, 6)
                fst_vals.append(fst)
            else:
                row[f"fst_{tag}"] = None

        if not fst_vals:
            continue
        row["mean_fst"] = round(float(np.mean(fst_vals)), 6)
        row["max_fst"]  = round(float(np.max(fst_vals)),  6)

        # Cross-annotations
        row["in_convergent_sweep"] = gene_id in conv_genes
        row["n_conv_groups"]       = conv_genes[gene_id]["n_groups"] if gene_id in conv_genes else 0
        row["max_h12"]             = conv_genes[gene_id]["max_h12"]  if gene_id in conv_genes else 0.0

        results.append(row)

    results.sort(key=lambda r: r["mean_fst"], reverse=True)
    print(f"  {len(results):,} genes with FST", flush=True)

    # ── Write outputs ─────────────────────────────────────────────────────────
    print(f"\n[3/3] Writing outputs ...", flush=True)
    pair_tags = [f"{a}_vs_{b}" for a, b in FOCAL_PAIRS]
    with open(OUT_TSV, "w") as f:
        header = "gene_id\tn_variants\tmean_fst\tmax_fst\tin_sweep\tn_conv_groups\tmax_h12\t"
        header += "\t".join(f"fst_{t}" for t in pair_tags) + "\n"
        f.write(header)
        for r in results:
            f.write(
                f"{r['gene_id']}\t{r['n_variants']}\t{r['mean_fst']}\t{r['max_fst']}\t"
                f"{r['in_convergent_sweep']}\t{r['n_conv_groups']}\t{r['max_h12']}\t"
                + "\t".join(str(r.get(f"fst_{t}", "")) for t in pair_tags) + "\n"
            )

    with open(OUT_JSON, "w") as f:
        json.dump({
            "n_genes":      len(results),
            "focal_pairs":  [f"{a}_vs_{b}" for a, b in FOCAL_PAIRS],
            "min_variants": MIN_VARIANTS,
            "top_genes":    results[:100],
        }, f, indent=2)

    driver.close()

    # ── Print top results ─────────────────────────────────────────────────────
    print("\n=== Top 30 Genes by Mean FST ===")
    print(f"{'Rank':<5} {'Gene':<14} {'mean_FST':>9} {'max_FST':>9} {'N_var':>6} "
          f"{'Sweep':>6}  YRI/CEU  YRI/CHB  CEU/CHB")
    print("-" * 85)
    for i, r in enumerate(results[:30], 1):
        sweep = f"★{r['n_conv_groups']}grp" if r["in_convergent_sweep"] else ""
        f1 = f"{r.get('fst_YRI_vs_CEU', 0) or 0:.4f}"
        f2 = f"{r.get('fst_YRI_vs_CHB', 0) or 0:.4f}"
        f3 = f"{r.get('fst_CEU_vs_CHB', 0) or 0:.4f}"
        print(f"  {i:<4} {r['gene_id']:<14} {r['mean_fst']:>9.5f} {r['max_fst']:>9.5f} "
              f"{r['n_variants']:>6} {sweep:>6}  {f1}   {f2}   {f3}")

    print(f"\nOutput: {OUT_TSV}")
    print(f"        {OUT_JSON}")


if __name__ == "__main__":
    main()
