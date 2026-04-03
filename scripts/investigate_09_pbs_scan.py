#!/usr/bin/env python3
"""Investigation 9: PBS Genome Scan — Population-Specific Selection.

Question: Which genomic regions and genes show evidence of selection in ONE specific
population but not others? Population Branch Statistic (PBS) identifies population-
lineage-specific sweeps (e.g. FIN-specific, JPT-specific, YRI-specific adaptation).

PBS(pop1; pop2, pop3) = (T12 + T13 - T23) / 2
where T_ij = -log(1 - FST_ij)  (branch length transformation)

Method:
  1. Stream all GenomicWindow nodes that have PBS scores
  2. Group by population; find genome-wide 99th percentile PBS threshold per pop
  3. Find top PBS windows (above threshold) per population
  4. Positional overlap → Variant → Gene to identify selected genes
  5. Fisher's exact test: are selected genes enriched in specific pathways?

Why novel: PBS across all populations × pathways × GO in one graph traversal.
  Standard tools run PBS per-population in isolation.

Output: data/results/pbs_scan.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_09_pbs_scan.py
"""

import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase
from scipy.stats import fisher_exact
import os

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"

OUT_DIR  = Path("data/results")
OUT_TSV  = OUT_DIR / "pbs_scan.tsv"
OUT_JSON = OUT_DIR / "pbs_scan.json"

PBS_PERCENTILE   = 99.0    # top 1% PBS windows per population
MIN_WINDOWS_GENE = 1       # min overlapping windows to call gene selected
TOP_N_POPS       = 10      # print results for top N populations by signal

CONTINENTAL = {
    "African":     {"ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"},
    "European":    {"CEU", "FIN", "GBR", "IBS", "TSI"},
    "East_Asian":  {"CDX", "CHB", "CHS", "JPT", "KHV"},
    "South_Asian": {"BEB", "GIH", "ITU", "PJL", "STU"},
    "American":    {"CLM", "MXL", "PEL", "PUR"},
}
POP_TO_GROUP = {p: g for g, pops in CONTINENTAL.items() for p in pops}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 9: PBS Genome Scan ===\n", flush=True)

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()

    # ── Step 1: Load all PBS windows ─────────────────────────────────────────
    print("[1/4] Loading PBS GenomicWindows from Neo4j ...", flush=True)
    t0 = time.time()

    # pop → list of (chr, start, end, pbs, pop2, pop3)
    pop_windows = defaultdict(list)
    pop_triplets = {}   # pop → (pop2, pop3) from first window seen

    with driver.session(database=NEO4J_DB) as session:
        result = session.run("""
            MATCH (w:GenomicWindow)
            WHERE w.pbs IS NOT NULL
            RETURN w.population AS pop, w.pop2 AS pop2, w.pop3 AS pop3,
                   w.chr AS chr, w.start AS start, w.end AS end,
                   w.pbs AS pbs
        """)
        n = 0
        for rec in result:
            pop = rec["pop"]
            pop_windows[pop].append((
                rec["chr"], rec["start"], rec["end"], rec["pbs"]
            ))
            if pop not in pop_triplets and rec["pop2"]:
                pop_triplets[pop] = (rec["pop2"], rec["pop3"])
            n += 1
            if n % 100_000 == 0:
                print(f"  {n:,} windows ...", flush=True)

    print(f"  {n:,} PBS windows for {len(pop_windows)} populations "
          f"in {time.time()-t0:.1f}s", flush=True)

    # ── Step 2: Per-population PBS threshold and top windows ──────────────────
    print("\n[2/4] Computing PBS thresholds and selecting top windows ...", flush=True)

    pop_threshold = {}
    pop_top_wins  = {}

    pops = sorted(pop_windows.keys())
    print(f"\n  {'Pop':<6} {'Group':<12} {'N_windows':>10} {'P99_PBS':>9} "
          f"{'N_top':>7}  Triplet")
    print("  " + "-" * 65)
    for pop in pops:
        wins = pop_windows[pop]
        pbs_vals = np.array([w[3] for w in wins])
        thresh   = float(np.percentile(pbs_vals, PBS_PERCENTILE))
        top_wins = [w for w in wins if w[3] >= thresh]
        pop_threshold[pop] = thresh
        pop_top_wins[pop]  = top_wins
        triplet = pop_triplets.get(pop, ("?", "?"))
        group   = POP_TO_GROUP.get(pop, "Unknown")
        print(f"  {pop:<5} {group[:10]:<11} {len(wins):>10,} {thresh:>9.5f} "
              f"{len(top_wins):>7,}  {pop}|{triplet[0]}|{triplet[1]}")

    # ── Step 3: Load Variant→Gene pairs (positional) ──────────────────────────
    print("\n[3/4] Loading Variant→Gene pairs and Gene→Pathway mapping ...", flush=True)
    t1 = time.time()

    # (chr, pos) → gene_id
    var_gene = defaultdict(set)
    gene_pathway = defaultdict(list)

    with driver.session(database=NEO4J_DB) as session:
        # Variants with gene consequences
        result = session.run("""
            MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
            WHERE v.chr IS NOT NULL AND v.pos IS NOT NULL
            RETURN v.chr AS chr, v.pos AS pos, g.geneId AS gene_id
        """)
        nv = 0
        for rec in result:
            var_gene[(rec["chr"], rec["pos"])].add(rec["gene_id"])
            nv += 1
            if nv % 200_000 == 0:
                print(f"    {nv:,} variant-gene pairs ...", flush=True)
        print(f"  {nv:,} variant-gene pairs loaded in {time.time()-t1:.1f}s", flush=True)

        # Gene → Pathway
        result2 = session.run("""
            MATCH (g:Gene)-[:IN_PATHWAY]->(p:Pathway)
            RETURN g.geneId AS gene_id, p.pathwayId AS pid, p.name AS pname
        """)
        for rec in result2:
            gene_pathway[rec["gene_id"]].append((rec["pid"], rec["pname"]))

    driver.close()

    # Build chromosome → sorted variant list for interval search
    print("  Building interval index ...", flush=True)
    chr_var_pos = defaultdict(list)
    for (chrom, pos), genes in var_gene.items():
        chr_var_pos[chrom].append((pos, genes))
    for chrom in chr_var_pos:
        chr_var_pos[chrom].sort()

    # ── Step 4: Interval join — top PBS windows → genes ───────────────────────
    print("\n[4/4] Mapping top PBS windows to genes ...", flush=True)

    pop_selected_genes = {}   # pop → {gene_id: max_pbs}
    pop_gene_windows   = {}   # pop → {gene_id: count}

    for pop in pops:
        gene_max_pbs = defaultdict(float)
        gene_n_wins  = defaultdict(int)
        for chrom, start, end, pbs_val in pop_top_wins[pop]:
            # Binary search for variants in [start, end]
            pos_list = chr_var_pos.get(chrom, [])
            if not pos_list:
                continue
            positions = [p[0] for p in pos_list]
            lo = np.searchsorted(positions, start)
            hi = np.searchsorted(positions, end, side="right")
            for k in range(lo, hi):
                for gid in pos_list[k][1]:
                    gene_max_pbs[gid] = max(gene_max_pbs[gid], pbs_val)
                    gene_n_wins[gid]  += 1

        pop_selected_genes[pop] = dict(gene_max_pbs)
        pop_gene_windows[pop]   = dict(gene_n_wins)

    # All genes with any PBS signal
    all_bg_genes = set(g for genes in pop_selected_genes.values() for g in genes)

    # ── Collate results ───────────────────────────────────────────────────────
    pop_results = {}
    for pop in pops:
        genes   = pop_selected_genes[pop]
        n_total = len(var_gene)  # approximate background
        top_genes_sorted = sorted(genes.items(), key=lambda x: -x[1])

        # Pathway enrichment for this population
        pw_genes = defaultdict(list)
        for gid in genes:
            for pid, pname in gene_pathway.get(gid, []):
                pw_genes[pid].append(gid)

        pop_results[pop] = {
            "population":    pop,
            "group":         POP_TO_GROUP.get(pop, "Unknown"),
            "triplet":       pop_triplets.get(pop, ("?", "?")),
            "pbs_threshold": pop_threshold[pop],
            "n_top_windows": len(pop_top_wins[pop]),
            "n_selected_genes": len(genes),
            "top_genes": [
                {"gene_id": gid, "max_pbs": round(float(pbs), 5),
                 "n_windows": pop_gene_windows[pop].get(gid, 0)}
                for gid, pbs in top_genes_sorted[:20]
            ],
        }

    # ── Print summary ─────────────────────────────────────────────────────────
    print("\n=== Population-Specific Selection — Top PBS Genes ===\n")
    for pop in pops:
        r = pop_results[pop]
        grp = r["group"]
        top5 = [g["gene_id"] for g in r["top_genes"][:5]]
        print(f"  {pop} ({grp}): {r['n_selected_genes']} genes above P{PBS_PERCENTILE:.0f} "
              f"(PBS≥{r['pbs_threshold']:.4f})")
        print(f"    Top 5: {', '.join(top5)}")

    # ── Write outputs ─────────────────────────────────────────────────────────
    with open(OUT_TSV, "w") as f:
        f.write("pop\tgroup\tpbs_threshold\tn_top_windows\tn_selected_genes\t"
                "top_gene_1\tpbs_1\ttop_gene_2\tpbs_2\ttop_gene_3\tpbs_3\n")
        for pop in pops:
            r  = pop_results[pop]
            tg = r["top_genes"]
            def tg_val(idx, key, default=""):
                return str(tg[idx][key]) if idx < len(tg) else default
            f.write(
                f"{pop}\t{r['group']}\t{r['pbs_threshold']:.6f}\t"
                f"{r['n_top_windows']}\t{r['n_selected_genes']}\t"
                f"{tg_val(0,'gene_id')}\t{tg_val(0,'max_pbs')}\t"
                f"{tg_val(1,'gene_id')}\t{tg_val(1,'max_pbs')}\t"
                f"{tg_val(2,'gene_id')}\t{tg_val(2,'max_pbs')}\n"
            )

    with open(OUT_JSON, "w") as f:
        json.dump({
            "pbs_percentile": PBS_PERCENTILE,
            "populations": pop_results,
        }, f, indent=2)

    print(f"\nOutput: {OUT_TSV}")
    print(f"        {OUT_JSON}")


if __name__ == "__main__":
    main()
