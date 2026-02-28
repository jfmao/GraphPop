#!/usr/bin/env python3
"""Benchmark: Conditioned Tajima's D — GraphPop vs classical pipeline.

Computes Tajima's D for missense variants in a specific Reactome pathway,
comparing the graph-native single-query approach against the classical
multi-step pipeline (filter VCF → extract sites → recompute stats).

Usage:
    conda run -n graphevo python scripts/benchmark-conditioned.py
"""

import csv
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import allel

CSV_DIR = Path("data/raw/1000g/csv_out")
CHR = "chr22"
POP = "AFR"
PATHWAY = "Cell surface interactions at the vascular wall"
CONSEQUENCE = "missense_variant"

NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

# ---------------------------------------------------------------------------
# Classical pipeline: multi-step filtering
# ---------------------------------------------------------------------------

def classical_pipeline(csv_dir, chr_val, pop, consequence, pathway):
    """Simulate the classical pipeline:
    1. Read VCF/CSV (all variants)
    2. Read consequence annotations → filter to consequence type
    3. Read pathway annotations → filter to pathway genes
    4. Intersect: only variants that are missense AND in pathway
    5. Extract allele counts for target population
    6. Compute Tajima's D
    """
    t0 = time.perf_counter()

    # Step 1: Load all variant positions + allele counts for this chr
    variant_csv = csv_dir / "variant_nodes.csv"
    pop_index = None
    all_variants = {}  # variantId -> (ac, an)

    t_read_start = time.perf_counter()
    with open(variant_csv) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row[2] != chr_val:
                continue
            if pop_index is None:
                pop_ids = row[7].split(";")
                pop_index = pop_ids.index(pop)

            vid = row[0]
            ac = [int(x) for x in row[8].split(";")][pop_index]
            an = [int(x) for x in row[9].split(";")][pop_index]
            all_variants[vid] = (ac, an)
    t_read_variants = time.perf_counter() - t_read_start

    # Step 2: Load consequence annotations → filter missense
    consequence_csv = csv_dir / "has_consequence_edges.csv"
    missense_variant_genes = {}  # variantId -> set of geneIds

    t_read_start = time.perf_counter()
    with open(consequence_csv) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if len(row) < 4:
                continue
            vid = row[0]
            gene_id = row[1]
            cons = row[3] if len(row) > 3 else ""
            if cons == consequence and vid in all_variants:
                if vid not in missense_variant_genes:
                    missense_variant_genes[vid] = set()
                missense_variant_genes[vid].add(gene_id)
    t_read_consequences = time.perf_counter() - t_read_start

    # Step 3: Load pathway annotations → find genes in target pathway
    in_pathway_csv = csv_dir / "in_pathway_edges.csv"
    pathway_csv = csv_dir / "pathway_nodes.csv"

    # First find the pathway ID
    pathway_id = None
    t_read_start = time.perf_counter()
    with open(pathway_csv) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row[2] == pathway:
                pathway_id = row[0]
                break

    # Then find genes in that pathway
    pathway_genes = set()
    with open(in_pathway_csv) as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row[1] == pathway_id:
                pathway_genes.add(row[0])
    t_read_pathway = time.perf_counter() - t_read_start

    # Step 4: Intersect — missense variants whose gene is in the pathway
    t_filter_start = time.perf_counter()
    filtered_variants = []
    for vid, gene_set in missense_variant_genes.items():
        if gene_set & pathway_genes:
            ac, an = all_variants[vid]
            if an >= 2:
                filtered_variants.append((ac, an))
    t_filter = time.perf_counter() - t_filter_start

    # Step 5: Compute Tajima's D
    t_compute_start = time.perf_counter()
    ac_alt = np.array([v[0] for v in filtered_variants])
    an_arr = np.array([v[1] for v in filtered_variants])
    ac_ref = an_arr - ac_alt
    ac_2d = np.column_stack([ac_ref, ac_alt])
    tajima_d = allel.tajima_d(allel.AlleleCountsArray(ac_2d))
    t_compute = time.perf_counter() - t_compute_start

    t_total = time.perf_counter() - t0

    return {
        "tajima_d": tajima_d,
        "n_variants": len(filtered_variants),
        "n_pathway_genes": len(pathway_genes),
        "t_read_variants": t_read_variants,
        "t_read_consequences": t_read_consequences,
        "t_read_pathway": t_read_pathway,
        "t_filter": t_filter,
        "t_compute": t_compute,
        "t_total": t_total,
    }


# ---------------------------------------------------------------------------
# GraphPop pipeline: single Cypher call
# ---------------------------------------------------------------------------

def graphpop_pipeline(chr_val, pop, consequence, pathway):
    """Single GraphPop procedure call with conditioned query."""
    t0 = time.perf_counter()

    # Use the whole chromosome since we're filtering by consequence + pathway
    cmd = [
        "cypher-shell",
        "-u", NEO4J_USER, "-p", NEO4J_PASS,
        "--format", "plain",
        f"CALL graphpop.diversity('{chr_val}', 0, 999999999, '{pop}', "
        f"{{consequence: '{consequence}', pathway: '{pathway}'}}) "
        "YIELD tajima_d, n_variants, n_segregating "
        "RETURN tajima_d, n_variants, n_segregating;",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    t_total = time.perf_counter() - t0

    if result.returncode != 0:
        print(f"ERROR: {result.stderr}", file=sys.stderr)
        return None

    lines = result.stdout.strip().split("\n")
    if len(lines) < 2:
        print(f"ERROR: unexpected output: {result.stdout}", file=sys.stderr)
        return None

    vals = lines[1].split(",")
    return {
        "tajima_d": float(vals[0].strip()),
        "n_variants": int(float(vals[1].strip())),
        "n_segregating": int(float(vals[2].strip())),
        "t_total": t_total,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 78)
    print("BENCHMARK: Conditioned Tajima's D — GraphPop vs Classical Pipeline")
    print(f"  Chromosome:  {CHR}")
    print(f"  Population:  {POP}")
    print(f"  Consequence: {CONSEQUENCE}")
    print(f"  Pathway:     {PATHWAY}")
    print("=" * 78)
    print()

    # Classical
    print("Running CLASSICAL pipeline (read CSVs → filter → scikit-allel)...")
    classical = classical_pipeline(CSV_DIR, CHR, POP, CONSEQUENCE, PATHWAY)
    print(f"  Tajima's D   = {classical['tajima_d']:.10f}")
    print(f"  Variants     = {classical['n_variants']}")
    print(f"  Pathway genes= {classical['n_pathway_genes']}")
    print(f"  Timing breakdown:")
    print(f"    Read variants:    {classical['t_read_variants']:.3f}s")
    print(f"    Read consequences:{classical['t_read_consequences']:.3f}s")
    print(f"    Read pathways:    {classical['t_read_pathway']:.3f}s")
    print(f"    Filter/intersect: {classical['t_filter']:.3f}s")
    print(f"    Compute stats:    {classical['t_compute']:.3f}s")
    print(f"    TOTAL:            {classical['t_total']:.3f}s")
    print()

    # GraphPop
    print("Running GRAPHPOP pipeline (single Cypher CALL)...")
    graphpop = graphpop_pipeline(CHR, POP, CONSEQUENCE, PATHWAY)
    if graphpop is None:
        print("  FAILED — see error above")
        return 1

    print(f"  Tajima's D   = {graphpop['tajima_d']:.10f}")
    print(f"  Variants     = {graphpop['n_variants']}")
    print(f"    TOTAL:            {graphpop['t_total']:.3f}s")
    print()

    # Comparison
    print("-" * 78)
    print("COMPARISON")
    print("-" * 78)

    if classical['tajima_d'] != 0:
        rel_err = abs(classical['tajima_d'] - graphpop['tajima_d']) / abs(classical['tajima_d']) * 100
    else:
        rel_err = abs(graphpop['tajima_d']) * 100

    speedup = classical['t_total'] / graphpop['t_total'] if graphpop['t_total'] > 0 else float('inf')

    print(f"  Tajima's D relative error: {rel_err:.6f}%")
    print(f"  Classical time:  {classical['t_total']:.3f}s")
    print(f"  GraphPop time:   {graphpop['t_total']:.3f}s")
    print(f"  Speedup:         {speedup:.1f}x")
    print()

    # Variant count note
    if classical['n_variants'] != graphpop['n_variants']:
        print(f"  Note: variant counts differ "
              f"(classical={classical['n_variants']}, graphpop={graphpop['n_variants']})")
        print(f"    This is expected — GraphPop's conditioned query uses a graph traversal")
        print(f"    that may include different variant-gene-pathway paths.")

    print()
    print("=" * 78)
    if rel_err < 1.0:
        print("BENCHMARK PASSED")
    else:
        print("BENCHMARK: results diverge — investigate")
    print("=" * 78)

    return 0


if __name__ == "__main__":
    sys.exit(main())
