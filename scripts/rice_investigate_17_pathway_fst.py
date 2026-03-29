#!/usr/bin/env python3
"""Rice Inv.R17: Pathway-Level FST (Neo4j-dependent).

Rank Plant Reactome pathways by mean population differentiation (Hudson FST),
analogous to human Inv.01.

Data source: Neo4j (Variant → HAS_CONSEQUENCE → Gene → IN_PATHWAY → Pathway)
Output:      data/results/rice/rice_inv17_pathway_fst.{tsv,json}
Figure:      data/results/rice/figures/rice_fig17_pathway_fst.png
"""

import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from neo4j import GraphDatabase

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

# Key population pairs
PAIRS = [
    ("GJ-tmp", "XI-1A"),
    ("GJ-tmp", "GJ-trp"),
    ("XI-1A", "cA-Aus"),
]


def compute_pathway_fst(session, pop1, pop2):
    """Compute mean Hudson FST per pathway via graph traversal."""
    t0 = time.time()

    query = """
    MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(pw:Pathway)
    WHERE v.pop_ids IS NOT NULL
    RETURN pw.pathwayId AS pathwayId, pw.name AS pathwayName,
           count(DISTINCT g) AS n_genes,
           v.pop_ids AS pids, v.ac AS ac, v.an AS an
    """
    # This returns one row per variant-pathway pair. Aggregate in Python.
    pathway_data = {}  # pathwayId -> {name, genes, num_sum, den_sum, n_variants}

    for rec in session.run(query):
        pid = rec["pathwayId"]
        pname = rec["pathwayName"]
        pids = rec["pids"]
        ac = rec["ac"]
        an = rec["an"]

        if pop1 not in pids or pop2 not in pids:
            continue
        idx1 = pids.index(pop1)
        idx2 = pids.index(pop2)
        a1, n1 = ac[idx1], an[idx1]
        a2, n2 = ac[idx2], an[idx2]
        if n1 < 2 or n2 < 2:
            continue

        p1, p2 = a1 / n1, a2 / n2
        num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
        den = p1 * (1 - p2) + p2 * (1 - p1)

        if pid not in pathway_data:
            pathway_data[pid] = {"name": pname, "n_genes": rec["n_genes"],
                                 "num_sum": 0, "den_sum": 0, "n_variants": 0}
        pathway_data[pid]["num_sum"] += num
        pathway_data[pid]["den_sum"] += den
        pathway_data[pid]["n_variants"] += 1

    results = []
    for pid, pd in pathway_data.items():
        fst = max(0, pd["num_sum"] / pd["den_sum"]) if pd["den_sum"] > 0 else 0
        results.append({
            "pathwayId": pid, "name": pd["name"],
            "n_genes": pd["n_genes"], "n_variants": pd["n_variants"],
            "fst": fst,
        })

    elapsed = time.time() - t0
    print(f"  {pop1} vs {pop2}: {len(results)} pathways in {elapsed:.1f}s")
    sys.stdout.flush()
    return sorted(results, key=lambda x: -x["fst"])


def main():
    t0 = time.time()
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    all_results = {}
    all_rows = []

    for pop1, pop2 in PAIRS:
        pair = f"{pop1}_vs_{pop2}"
        print(f"Computing pathway-level FST: {pair}...")
        sys.stdout.flush()
        with driver.session() as session:
            results = compute_pathway_fst(session, pop1, pop2)
        all_results[pair] = {
            "n_pathways": len(results),
            "mean_fst": float(np.mean([r["fst"] for r in results])) if results else 0,
            "top_10": results[:10],
        }
        for r in results:
            r["pair"] = pair
            all_rows.append(r)

    driver.close()
    total_time = time.time() - t0

    # Save TSV
    with open(RESULTS_DIR / "rice_inv17_pathway_fst.tsv", "w") as f:
        f.write("pair\tpathwayId\tname\tn_genes\tn_variants\tfst\n")
        for r in sorted(all_rows, key=lambda x: -x["fst"]):
            f.write(f"{r['pair']}\t{r['pathwayId']}\t{r['name']}\t{r['n_genes']}\t{r['n_variants']}\t{r['fst']:.6f}\n")

    with open(RESULTS_DIR / "rice_inv17_pathway_fst.json", "w") as f:
        json.dump({"pairs": all_results, "total_time_s": total_time}, f, indent=2, default=str)

    # ── Figure: top pathways per pair ──
    fig, axes = plt.subplots(1, len(PAIRS), figsize=(6 * len(PAIRS), 7))

    for idx, (pop1, pop2) in enumerate(PAIRS):
        ax = axes[idx]
        pair = f"{pop1}_vs_{pop2}"
        top = all_results[pair]["top_10"]

        names = [t["name"][:35] for t in top]
        fst_vals = [t["fst"] for t in top]

        ax.barh(range(len(names)), fst_vals, color="#4393c3", edgecolor="black", linewidth=0.5)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=6)
        ax.set_xlabel("Pathway mean FST")
        n_pw = all_results[pair]["n_pathways"]
        mean_fst = all_results[pair]["mean_fst"]
        ax.set_title(f"{pop1} vs {pop2}\n({n_pw} pathways, mean={mean_fst:.3f})",
                     fontweight="bold", fontsize=9)
        ax.invert_yaxis()

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig17_pathway_fst.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"\n=== Rice Inv.R17: Pathway-Level FST ===")
    print(f"Total time: {total_time:.1f}s")
    for pair, data in all_results.items():
        print(f"  {pair}: {data['n_pathways']} pathways, mean_fst={data['mean_fst']:.4f}")
        if data["top_10"]:
            top1 = data["top_10"][0]
            print(f"    Top: {top1['name']} (Fst={top1['fst']:.3f}, {top1['n_variants']} variants)")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
