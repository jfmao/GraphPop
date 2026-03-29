#!/usr/bin/env python3
"""Rice Inv.R13: Gene-Level FST (Neo4j-dependent).

Compute per-gene FST across rice subpopulations by traversing
HAS_CONSEQUENCE edges to get gene-associated variants.

Data source: Neo4j (Variant→Gene via HAS_CONSEQUENCE)
Output:      data/results/rice/rice_inv13_gene_fst.{tsv,json}
Figure:      data/results/rice/figures/rice_fig13_gene_fst.png
"""

import json
import time
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from neo4j import GraphDatabase

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

# Key population pairs
PAIRS = [
    ("GJ-tmp", "XI-1A"),
    ("GJ-tmp", "GJ-trp"),
    ("XI-1A", "cA-Aus"),
]

POP_GROUPS = {
    "XI-1A": "Indica", "XI-1B": "Indica", "XI-2": "Indica",
    "XI-3": "Indica", "XI-adm": "Indica",
    "GJ-tmp": "Japonica", "GJ-trp": "Japonica",
    "GJ-sbtrp": "Japonica", "GJ-adm": "Japonica",
    "cA-Aus": "Aus/Basmati", "cB-Bas": "Aus/Basmati",
    "admix": "Admixed",
}


def compute_gene_fst(session, pop1, pop2, limit=500):
    """Compute Hudson FST per gene for a population pair."""
    t0 = time.time()
    query = """
    MATCH (v:Variant)-[r:HAS_CONSEQUENCE]->(g:Gene)
    WHERE v.pop_ids IS NOT NULL
    WITH g, collect(v) AS variants
    WHERE size(variants) >= 5
    WITH g, variants,
         [p IN $pops | apoc.coll.indexOf($pops, p)] AS pop_indices
    UNWIND variants AS v
    WITH g, v,
         v.pop_ids AS pids,
         v.ac AS ac, v.an AS an
    WITH g, v, pids, ac, an,
         apoc.coll.indexOf(pids, $pop1) AS idx1,
         apoc.coll.indexOf(pids, $pop2) AS idx2
    WHERE idx1 >= 0 AND idx2 >= 0
    WITH g,
         ac[idx1] AS ac1, an[idx1] AS an1,
         ac[idx2] AS ac2, an[idx2] AS an2
    WHERE an1 > 0 AND an2 > 0
    WITH g,
         collect({ac1: ac1, an1: an1, ac2: ac2, an2: an2}) AS alleles
    WHERE size(alleles) >= 5
    RETURN g.geneId AS geneId, g.symbol AS symbol,
           size(alleles) AS n_variants,
           alleles
    """

    # Simpler approach: use graphpop.diversity conditioned on gene
    # But that requires per-gene calls. Let's use a batch query instead.
    # Actually, compute FST from the allele count arrays directly.

    query_simple = """
    MATCH (v:Variant)-[hc:HAS_CONSEQUENCE]->(g:Gene)
    WHERE v.pop_ids IS NOT NULL
    WITH g.geneId AS geneId, g.symbol AS symbol,
         collect({
           ac: v.ac, an: v.an, pids: v.pop_ids,
           impact: hc.impact
         }) AS variants
    WHERE size(variants) >= 5
    RETURN geneId, symbol, size(variants) AS n_variants, variants
    LIMIT $limit
    """

    results = []
    records = session.run(query_simple, limit=limit)

    for rec in records:
        gene_id = rec["geneId"]
        symbol = rec["symbol"]
        variants = rec["variants"]
        n_var = rec["n_variants"]

        # Compute Hudson FST
        num_sum, den_sum = 0.0, 0.0
        n_used = 0
        for var in variants:
            pids = var["pids"]
            ac = var["ac"]
            an = var["an"]
            idx1 = pids.index(pop1) if pop1 in pids else -1
            idx2 = pids.index(pop2) if pop2 in pids else -1
            if idx1 < 0 or idx2 < 0:
                continue
            a1, n1 = ac[idx1], an[idx1]
            a2, n2 = ac[idx2], an[idx2]
            if n1 < 2 or n2 < 2:
                continue

            p1 = a1 / n1
            p2 = a2 / n2
            # Hudson FST numerator and denominator
            num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
            den = p1 * (1 - p2) + p2 * (1 - p1)
            num_sum += num
            den_sum += den
            n_used += 1

        fst = num_sum / den_sum if den_sum > 0 else 0.0
        fst = max(0.0, fst)  # floor at 0

        results.append({
            "geneId": gene_id,
            "symbol": symbol,
            "n_variants": n_var,
            "n_used": n_used,
            "fst": fst,
        })

    elapsed = time.time() - t0
    print(f"  {pop1} vs {pop2}: {len(results)} genes in {elapsed:.1f}s")
    sys.stdout.flush()
    return results


def main():
    t0 = time.time()
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    all_results = {}
    all_rows = []

    for pop1, pop2 in PAIRS:
        pair_label = f"{pop1}_vs_{pop2}"
        print(f"Computing gene-level FST: {pair_label}...")
        sys.stdout.flush()

        with driver.session() as session:
            results = compute_gene_fst(session, pop1, pop2, limit=5000)

        all_results[pair_label] = {
            "n_genes": len(results),
            "mean_fst": float(np.mean([r["fst"] for r in results])) if results else 0,
            "max_fst": float(np.max([r["fst"] for r in results])) if results else 0,
            "top_10": sorted(results, key=lambda x: x["fst"], reverse=True)[:10],
        }

        for r in results:
            r["pair"] = pair_label
            all_rows.append(r)

    driver.close()
    total_time = time.time() - t0

    # Save TSV
    with open(RESULTS_DIR / "rice_inv13_gene_fst.tsv", "w") as f:
        f.write("pair\tgeneId\tsymbol\tn_variants\tn_used\tfst\n")
        for r in sorted(all_rows, key=lambda x: (-x["fst"])):
            f.write(f"{r['pair']}\t{r['geneId']}\t{r['symbol']}\t{r['n_variants']}\t{r['n_used']}\t{r['fst']:.6f}\n")

    summary = {
        "pairs": all_results,
        "total_genes": len(set(r["geneId"] for r in all_rows)),
        "total_time_s": total_time,
    }
    with open(RESULTS_DIR / "rice_inv13_gene_fst.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # ── Figure: top genes per pair ──
    n_pairs = len(PAIRS)
    fig, axes = plt.subplots(1, n_pairs, figsize=(6 * n_pairs, 6))
    if n_pairs == 1:
        axes = [axes]

    for idx, (pop1, pop2) in enumerate(PAIRS):
        ax = axes[idx]
        pair_label = f"{pop1}_vs_{pop2}"
        pair_data = all_results[pair_label]
        top = pair_data["top_10"]

        genes = [t["symbol"][:20] for t in top]
        fst_vals = [t["fst"] for t in top]

        ax.barh(range(len(genes)), fst_vals, color="#4393c3", edgecolor="black", linewidth=0.5)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=7)
        ax.set_xlabel("Hudson FST")
        ax.set_title(f"{pop1} vs {pop2}\n({pair_data['n_genes']} genes, mean={pair_data['mean_fst']:.3f})",
                     fontweight="bold", fontsize=9)
        ax.invert_yaxis()

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig13_gene_fst.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"\n=== Rice Inv.R13: Gene-Level FST ===")
    print(f"Total time: {total_time:.1f}s")
    for pair, data in all_results.items():
        print(f"  {pair}: {data['n_genes']} genes, mean_fst={data['mean_fst']:.4f}, max_fst={data['max_fst']:.4f}")
        top3 = ', '.join(t['symbol'] + f"({t['fst']:.3f})" for t in data['top_10'][:3])
        print(f"    Top 3: {top3}")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
