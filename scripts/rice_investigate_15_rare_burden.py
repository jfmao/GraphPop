#!/usr/bin/env python3
"""Rice Inv.R15: Population-Specific Rare Variant Burden (Neo4j-dependent).

Count rare functional variants per population, stratified by impact level.

Data source: Neo4j (Variant nodes with AC arrays + HAS_CONSEQUENCE)
Output:      data/results/rice/rice_inv15_rare_burden.{tsv,json}
Figure:      data/results/rice/figures/rice_fig15_rare_burden.png
"""

import json
import time
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from neo4j import GraphDatabase
import os

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")

POPS = ["XI-1A", "XI-1B", "XI-2", "XI-3", "XI-adm",
        "GJ-tmp", "GJ-trp", "GJ-sbtrp", "GJ-adm",
        "cA-Aus", "cB-Bas", "admix"]

POP_GROUPS = {
    "XI-1A": "Indica", "XI-1B": "Indica", "XI-2": "Indica",
    "XI-3": "Indica", "XI-adm": "Indica",
    "GJ-tmp": "Japonica", "GJ-trp": "Japonica",
    "GJ-sbtrp": "Japonica", "GJ-adm": "Japonica",
    "cA-Aus": "Aus/Basmati", "cB-Bas": "Aus/Basmati",
    "admix": "Admixed",
}

GROUP_COLORS = {
    "Indica": "#e41a1c", "Japonica": "#377eb8",
    "Aus/Basmati": "#4daf4a", "Admixed": "#984ea3",
}


def count_rare_burden(session, pop, af_threshold=0.05):
    """Count rare functional variants (AF < threshold) per impact."""
    t0 = time.time()
    query = """
    MATCH (v:Variant)-[r:HAS_CONSEQUENCE]->(g:Gene)
    WHERE v.pop_ids IS NOT NULL
    RETURN r.impact AS impact, v.pop_ids AS pids, v.ac AS ac, v.an AS an
    """
    result = {"HIGH": {"n_rare": 0, "af_sum": 0.0},
              "MODERATE": {"n_rare": 0, "af_sum": 0.0},
              "LOW": {"n_rare": 0, "af_sum": 0.0}}

    for rec in session.run(query):
        impact = rec["impact"]
        if impact not in result:
            continue
        pids = rec["pids"]
        if pop not in pids:
            continue
        idx = pids.index(pop)
        ac = rec["ac"]
        an = rec["an"]
        if an[idx] <= 0 or ac[idx] <= 0:
            continue
        af = ac[idx] / an[idx]
        if af < af_threshold:
            result[impact]["n_rare"] += 1
            result[impact]["af_sum"] += af

    final = {}
    for imp, d in result.items():
        if d["n_rare"] > 0:
            final[imp] = {
                "n_rare": d["n_rare"],
                "mean_af": d["af_sum"] / d["n_rare"],
            }

    elapsed = time.time() - t0
    total = sum(d.get("n_rare", 0) for d in final.values())
    print(f"  {pop}: {total:,} rare variants in {elapsed:.1f}s")
    sys.stdout.flush()
    return final


def main():
    t0 = time.time()
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    all_results = {}
    rows = []

    print("Computing rare variant burden per population...")
    sys.stdout.flush()

    with driver.session() as session:
        for pop in POPS:
            burden = count_rare_burden(session, pop)
            all_results[pop] = burden

            total = sum(v["n_rare"] for v in burden.values())
            for impact, data in burden.items():
                rows.append({
                    "population": pop,
                    "group": POP_GROUPS.get(pop, "?"),
                    "impact": impact,
                    "n_rare": data["n_rare"],
                    "mean_af": data["mean_af"],
                    "total_rare": total,
                })

    driver.close()
    total_time = time.time() - t0

    # Save TSV
    with open(RESULTS_DIR / "rice_inv15_rare_burden.tsv", "w") as f:
        f.write("population\tgroup\timpact\tn_rare\tmean_af\ttotal_rare\n")
        for r in rows:
            f.write(f"{r['population']}\t{r['group']}\t{r['impact']}\t"
                    f"{r['n_rare']}\t{r['mean_af']:.6f}\t{r['total_rare']}\n")

    summary = {
        "populations": {
            pop: {
                "total_rare": sum(v["n_rare"] for v in data.values()),
                "per_impact": data,
            }
            for pop, data in all_results.items()
        },
        "total_time_s": total_time,
    }
    with open(RESULTS_DIR / "rice_inv15_rare_burden.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 2 panels ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel (a): Total rare burden per population
    ax = axes[0]
    pop_sorted = sorted(POPS, key=lambda p: sum(v["n_rare"] for v in all_results[p].values()),
                        reverse=True)
    totals = [sum(v["n_rare"] for v in all_results[p].values()) for p in pop_sorted]
    colors = [GROUP_COLORS.get(POP_GROUPS.get(p, "?"), "gray") for p in pop_sorted]
    bars = ax.barh(range(len(pop_sorted)), totals, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(pop_sorted)))
    ax.set_yticklabels(pop_sorted, fontsize=8)
    ax.set_xlabel("Number of rare functional variants (AF < 5%)")
    ax.set_title("a  Rare variant burden per population", fontweight="bold", loc="left")
    ax.invert_yaxis()

    for bar, v in zip(bars, totals):
        ax.text(v + 50, bar.get_y() + bar.get_height() / 2, f"{v:,}", va="center", fontsize=7)

    # Panel (b): Stacked by impact
    ax = axes[1]
    impacts = ["HIGH", "MODERATE", "LOW"]
    impact_colors = {"HIGH": "#d73027", "MODERATE": "#fc8d59", "LOW": "#4575b4"}
    bottom = np.zeros(len(pop_sorted))
    for imp in impacts:
        vals = [all_results[p].get(imp, {}).get("n_rare", 0) for p in pop_sorted]
        ax.barh(range(len(pop_sorted)), vals, left=bottom, color=impact_colors[imp],
                label=imp, edgecolor="black", linewidth=0.3)
        bottom += np.array(vals)

    ax.set_yticks(range(len(pop_sorted)))
    ax.set_yticklabels(pop_sorted, fontsize=8)
    ax.set_xlabel("Rare variants by impact level")
    ax.set_title("b  Impact stratification", fontweight="bold", loc="left")
    ax.legend(fontsize=7)
    ax.invert_yaxis()

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig15_rare_burden.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"\n=== Rice Inv.R15: Rare Variant Burden ===")
    print(f"Total time: {total_time:.1f}s")
    for pop in pop_sorted[:5]:
        total = sum(v["n_rare"] for v in all_results[pop].values())
        print(f"  {pop}: {total:,} rare variants")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
