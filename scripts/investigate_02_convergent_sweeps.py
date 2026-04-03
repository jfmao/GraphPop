#!/usr/bin/env python3
"""Investigation 2: Convergent Positive Selection Across Populations.

Question: Which genes and pathways show hard sweep signals in multiple independent
population groups simultaneously — evidence of convergent adaptation?

Method:
  1. Pull all GenomicWindow nodes with sweep_type='hard' → {chr, start, end, population}
  2. Pull all Variant→Gene pairs via HAS_CONSEQUENCE → {variantId, chr, pos, geneId}
  3. Python interval join: for each gene variant, find overlapping hard-sweep windows
  4. Classify populations into continental groups; require sweeps in 2+ groups
  5. Traverse Gene → IN_PATHWAY → Pathway → aggregate convergent genes per pathway

Why novel: Convergent adaptation across independent lineages detected as a multi-layer
graph pattern (GenomicWindow → position → Variant → Gene → Pathway) in one analysis.
No VCF tool integrates sweep detection + gene annotation + pathway enrichment natively.

Output: data/results/convergent_sweeps.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_02_convergent_sweeps.py
"""

import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase

# ── Config ───────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"

OUT_DIR = Path("data/results")
OUT_TSV = OUT_DIR / "convergent_sweeps.tsv"
OUT_JSON = OUT_DIR / "convergent_sweeps.json"

H12_THRESHOLD  = 0.3   # lower threshold to capture more sweep candidates
MIN_CONV_GROUPS = 2     # minimum continental groups showing a sweep in a gene

# Continental group assignments (1000 Genomes sub-populations)
CONTINENTAL = {
    "African":     {"ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"},
    "European":    {"CEU", "FIN", "GBR", "IBS", "TSI"},
    "East_Asian":  {"CDX", "CHB", "CHS", "JPT", "KHV"},
    "South_Asian": {"BEB", "GIH", "ITU", "PJL", "STU"},
    "American":    {"CLM", "MXL", "PEL", "PUR"},
}
POP_TO_GROUP = {pop: grp for grp, pops in CONTINENTAL.items() for pop in pops}


# ── Step 1: Load hard sweep windows ──────────────────────────────────────────

SWEEP_CYPHER = """
MATCH (w:GenomicWindow)
WHERE w.h12 >= $threshold AND w.sweep_type = 'hard'
RETURN w.chr AS chr, w.start AS start, w.end AS end,
       w.population AS population, w.h12 AS h12
"""

def load_sweep_windows(session, threshold):
    """Return dict: chr → sorted list of (start, end, population, h12)."""
    print(f"  Loading sweep windows (h12 ≥ {threshold}, sweep_type='hard') ...", flush=True)
    windows = defaultdict(list)
    n = 0
    result = session.run(SWEEP_CYPHER, threshold=threshold)
    for rec in result:
        windows[rec["chr"]].append((rec["start"], rec["end"],
                                    rec["population"], rec["h12"]))
        n += 1
    # Sort by start for binary search
    for chrom in windows:
        windows[chrom].sort()
    pops_seen = set()
    for wins in windows.values():
        for _, _, pop, _ in wins:
            pops_seen.add(pop)
    print(f"  {n:,} hard sweep windows across {len(windows)} chromosomes, "
          f"{len(pops_seen)} populations", flush=True)
    return dict(windows)


# ── Step 2: Load Variant → Gene pairs ────────────────────────────────────────

VARIANT_GENE_CYPHER = """
MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
WHERE v.chr IS NOT NULL AND v.pos IS NOT NULL
RETURN v.variantId AS vid, v.chr AS chr, v.pos AS pos, g.geneId AS gene_id
"""

def load_variant_gene(session):
    """Return list of (chr, pos, geneId)."""
    print("  Loading Variant→Gene pairs (HAS_CONSEQUENCE) ...", flush=True)
    pairs = []
    n = 0
    result = session.run(VARIANT_GENE_CYPHER)
    for rec in result:
        pairs.append((rec["chr"], rec["pos"], rec["gene_id"]))
        n += 1
        if n % 100_000 == 0:
            print(f"    {n:,} pairs loaded ...", flush=True)
    print(f"  {n:,} Variant→Gene pairs loaded", flush=True)
    return pairs


# ── Step 3: Load Gene → Pathway ───────────────────────────────────────────────

GENE_PATHWAY_CYPHER = """
MATCH (g:Gene)-[:IN_PATHWAY]->(p:Pathway)
RETURN g.geneId AS gene_id, p.pathwayId AS pathway_id, p.name AS pathway_name
"""

def load_gene_pathway(session):
    """Return dict: gene_id → list of (pathway_id, name)."""
    print("  Loading Gene→Pathway pairs (IN_PATHWAY) ...", flush=True)
    mapping = defaultdict(list)
    result = session.run(GENE_PATHWAY_CYPHER)
    n = 0
    for rec in result:
        mapping[rec["gene_id"]].append((rec["pathway_id"], rec["pathway_name"]))
        n += 1
    print(f"  {n:,} Gene→Pathway pairs for {len(mapping):,} genes", flush=True)
    return dict(mapping)


# ── Step 4: Interval overlap (binary search) ──────────────────────────────────

def find_overlapping_windows(chrom_windows, chrom, pos):
    """Return list of (population, h12) for windows overlapping position."""
    wins = chrom_windows.get(chrom, [])
    if not wins:
        return []
    # Binary search for windows that could overlap pos
    # Window overlaps pos if start <= pos <= end
    starts = [w[0] for w in wins]
    # Find rightmost window with start <= pos
    lo, hi = 0, len(wins)
    idx = 0
    while lo < hi:
        mid = (lo + hi) // 2
        if starts[mid] <= pos:
            idx = mid + 1
            lo = mid + 1
        else:
            hi = mid
    # Check all windows up to idx that might overlap
    hits = []
    for i in range(idx - 1, -1, -1):
        start, end, pop, h12 = wins[i]
        if end < pos:
            break
        if start <= pos <= end:
            hits.append((pop, h12))
    return hits


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 2: Convergent Positive Selection ===\n", flush=True)

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()

    with driver.session(database=NEO4J_DB) as session:
        # Load all data from Neo4j
        print("[1/4] Loading sweep windows from Neo4j ...", flush=True)
        sweep_wins = load_sweep_windows(session, H12_THRESHOLD)

        print("\n[2/4] Loading Variant→Gene pairs from Neo4j ...", flush=True)
        t0 = time.time()
        vg_pairs = load_variant_gene(session)
        print(f"  Done in {time.time()-t0:.1f}s", flush=True)

        print("\n[3/4] Loading Gene→Pathway mapping from Neo4j ...", flush=True)
        gene_pathway = load_gene_pathway(session)

    driver.close()

    # ── Step 4: Interval join ─────────────────────────────────────────────────
    print("\n[4/4] Computing convergent sweep genes ...", flush=True)
    t0 = time.time()

    # gene_id → set of continental groups showing a sweep for any variant in that gene
    gene_groups = defaultdict(set)
    # gene_id → set of populations
    gene_pops = defaultdict(set)
    # gene_id → max h12 across all overlapping windows
    gene_max_h12 = defaultdict(float)

    n_overlaps = 0
    for i, (chrom, pos, gene_id) in enumerate(vg_pairs):
        if i % 100_000 == 0 and i > 0:
            print(f"  {i:,} / {len(vg_pairs):,} variants processed ...", flush=True)
        hits = find_overlapping_windows(sweep_wins, chrom, pos)
        for pop, h12 in hits:
            grp = POP_TO_GROUP.get(pop)
            if grp:
                gene_groups[gene_id].add(grp)
                gene_pops[gene_id].add(pop)
                gene_max_h12[gene_id] = max(gene_max_h12[gene_id], h12)
                n_overlaps += 1

    print(f"  {n_overlaps:,} variant-window overlaps found", flush=True)
    print(f"  {len(gene_groups):,} genes with ≥1 sweep overlap", flush=True)

    # Filter: genes with sweeps in 2+ continental groups
    convergent_genes = {
        gid: {"groups": sorted(grps), "pops": sorted(gene_pops[gid]),
              "max_h12": round(gene_max_h12[gid], 4),
              "n_groups": len(grps)}
        for gid, grps in gene_groups.items()
        if len(grps) >= MIN_CONV_GROUPS
    }
    print(f"  {len(convergent_genes):,} genes with sweeps in ≥{MIN_CONV_GROUPS} continental groups", flush=True)
    print(f"  Done in {time.time()-t0:.1f}s", flush=True)

    # ── Pathway enrichment ────────────────────────────────────────────────────
    # For each pathway: count convergent genes, total genes, list populations
    pathway_conv  = defaultdict(list)    # pathway_id → list of convergent gene_ids
    pathway_total = defaultdict(set)     # pathway_id → all gene_ids (for background)

    all_genes_in_pathways = set()
    for gid, pairs in gene_pathway.items():
        for pid, pname in pairs:
            pathway_total[pid].add(gid)
            all_genes_in_pathways.add(gid)
            if gid in convergent_genes:
                pathway_conv[pid].append(gid)

    # Build result rows
    n_conv_total = len(set(convergent_genes.keys()) & all_genes_in_pathways)
    n_total      = len(all_genes_in_pathways)

    from scipy.stats import fisher_exact
import os
    pathway_results = []
    for pid, conv_genes in pathway_conv.items():
        n_pathway       = len(pathway_total[pid])
        n_conv_pathway  = len(conv_genes)
        if n_conv_pathway == 0:
            continue
        # Fisher's exact test (one-sided)
        #           In pathway | Not in pathway
        # Conv         a              b
        # Not conv     c              d
        a = n_conv_pathway
        b = n_conv_total - n_conv_pathway
        c = n_pathway - n_conv_pathway
        d = n_total - n_pathway - b
        _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        or_val = (a * d) / (b * c) if b > 0 and c > 0 else float("inf")

        # Representative pathway name from gene_pathway
        pname = gene_pathway[conv_genes[0]][0][1] if conv_genes else pid

        # Populations and groups in this pathway's convergent genes
        all_grps = set()
        all_pops = set()
        for gid in conv_genes:
            all_grps.update(convergent_genes[gid]["groups"])
            all_pops.update(convergent_genes[gid]["pops"])

        pathway_results.append({
            "pathway_id":         pid,
            "pathway_name":       pname,
            "n_convergent_genes": n_conv_pathway,
            "n_total_genes":      n_pathway,
            "pct_convergent":     round(100 * n_conv_pathway / n_pathway, 1),
            "continental_groups": sorted(all_grps),
            "n_groups":           len(all_grps),
            "populations":        sorted(all_pops),
            "fisher_pval":        round(pval, 6),
            "odds_ratio":         round(or_val, 3),
            "convergent_genes":   sorted(conv_genes),
        })

    # Sort by number of convergent genes then p-value
    pathway_results.sort(key=lambda r: (-r["n_groups"], -r["n_convergent_genes"], r["fisher_pval"]))

    # ── Write outputs ─────────────────────────────────────────────────────────
    print(f"\nWriting {len(pathway_results):,} pathways to {OUT_TSV} ...", flush=True)
    with open(OUT_TSV, "w") as f:
        header = ("pathway_id\tpathway_name\tn_convergent_genes\tn_total_genes\t"
                  "pct_convergent\tn_groups\tcontinental_groups\t"
                  "fisher_pval\todds_ratio\tconvergent_genes")
        f.write(header + "\n")
        for r in pathway_results:
            f.write(
                f"{r['pathway_id']}\t{r['pathway_name']}\t{r['n_convergent_genes']}\t"
                f"{r['n_total_genes']}\t{r['pct_convergent']}\t{r['n_groups']}\t"
                f"{','.join(r['continental_groups'])}\t{r['fisher_pval']}\t"
                f"{r['odds_ratio']}\t{','.join(r['convergent_genes'])}\n"
            )

    with open(OUT_JSON, "w") as f:
        json.dump({
            "summary": {
                "h12_threshold":       H12_THRESHOLD,
                "min_conv_groups":     MIN_CONV_GROUPS,
                "n_convergent_genes":  len(convergent_genes),
                "n_pathways":          len(pathway_results),
            },
            "convergent_genes": convergent_genes,
            "pathways":         pathway_results,
        }, f, indent=2)

    # Print top results
    print("\n=== Top 20 Pathways with Convergent Sweeps ===")
    print(f"{'Rank':<5} {'Conv':>5} {'Total':>6} {'Groups':<30}  Pathway")
    print("-" * 90)
    for i, r in enumerate(pathway_results[:20], 1):
        grps = ",".join(g[:3] for g in r["continental_groups"])
        name = r["pathway_name"][:45]
        print(f"  {i:<4} {r['n_convergent_genes']:>5} {r['n_total_genes']:>6}  {grps:<30}  {name}")

    print(f"\nFull results: {OUT_TSV}")
    print(f"  Convergent genes summary: {OUT_JSON}")


if __name__ == "__main__":
    main()
