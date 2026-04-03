#!/usr/bin/env python3
"""Rice Inv.R14: Consequence-Level FST Bias (Neo4j-dependent).

Compare FST distributions across impact levels (HIGH, MODERATE, LOW)
to test whether purifying selection constrains differentiation at
functional sites.

Data source: Neo4j (HAS_CONSEQUENCE edges with impact levels)
Output:      data/results/rice/rice_inv14_consequence_bias.{tsv,json}
Figure:      data/results/rice/figures/rice_fig14_consequence_bias.png
"""

import json
import time
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from neo4j import GraphDatabase
import os

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")

POP1 = "GJ-tmp"
POP2 = "XI-1A"
IMPACTS = ["HIGH", "MODERATE", "LOW"]


def query_fst_by_impact(session, pop1, pop2, impact, limit=50000):
    """Get per-variant FST for variants with a given impact level."""
    t0 = time.time()
    query = """
    MATCH (v:Variant)-[r:HAS_CONSEQUENCE]->(g:Gene)
    WHERE r.impact = $impact AND v.pop_ids IS NOT NULL
    WITH DISTINCT v
    RETURN v.pop_ids AS pids, v.ac AS ac, v.an AS an
    LIMIT $limit
    """
    records = list(session.run(query, impact=impact, limit=limit))

    fst_vals = []
    for rec in records:
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
        if den > 0:
            fst_vals.append(max(0, num / den))

    elapsed = time.time() - t0
    print(f"  {impact}: {len(fst_vals)} variants ({len(records)} queried) in {elapsed:.1f}s")
    sys.stdout.flush()
    return fst_vals


def main():
    t0 = time.time()
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    impact_data = {}
    rows = []

    print(f"Computing consequence-level FST: {POP1} vs {POP2}")
    sys.stdout.flush()

    with driver.session() as session:
        for impact in IMPACTS:
            fst_vals = query_fst_by_impact(session, POP1, POP2, impact)
            impact_data[impact] = {
                "n_variants": len(fst_vals),
                "mean_fst": float(np.mean(fst_vals)) if fst_vals else 0,
                "median_fst": float(np.median(fst_vals)) if fst_vals else 0,
                "std_fst": float(np.std(fst_vals)) if fst_vals else 0,
                "fst_values": fst_vals,
            }

    driver.close()

    # Mann-Whitney tests between impact levels
    mw_results = {}
    for i in range(len(IMPACTS)):
        for j in range(i + 1, len(IMPACTS)):
            imp1, imp2 = IMPACTS[i], IMPACTS[j]
            if impact_data[imp1]["fst_values"] and impact_data[imp2]["fst_values"]:
                stat, pval = stats.mannwhitneyu(
                    impact_data[imp1]["fst_values"],
                    impact_data[imp2]["fst_values"],
                    alternative="two-sided"
                )
                mw_results[f"{imp1}_vs_{imp2}"] = {"U": float(stat), "p": float(pval)}

    total_time = time.time() - t0

    # Save TSV
    with open(RESULTS_DIR / "rice_inv14_consequence_bias.tsv", "w") as f:
        f.write("impact\tn_variants\tmean_fst\tmedian_fst\tstd_fst\n")
        for imp in IMPACTS:
            d = impact_data[imp]
            f.write(f"{imp}\t{d['n_variants']}\t{d['mean_fst']:.6f}\t{d['median_fst']:.6f}\t{d['std_fst']:.6f}\n")

    summary = {
        "pair": f"{POP1}_vs_{POP2}",
        "impacts": {k: {kk: vv for kk, vv in v.items() if kk != "fst_values"}
                    for k, v in impact_data.items()},
        "mann_whitney": mw_results,
        "total_time_s": total_time,
    }
    with open(RESULTS_DIR / "rice_inv14_consequence_bias.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 2 panels ──
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel (a): FST distribution by impact
    ax = axes[0]
    colors = {"HIGH": "#d73027", "MODERATE": "#fc8d59", "LOW": "#4575b4"}
    for imp in IMPACTS:
        vals = impact_data[imp]["fst_values"]
        if vals:
            ax.hist(vals, bins=50, alpha=0.5, density=True,
                    label=f"{imp} (n={len(vals)}, mean={np.mean(vals):.4f})",
                    color=colors[imp])
    ax.set_xlabel("Per-variant FST")
    ax.set_ylabel("Density")
    ax.set_title(f"a  FST by impact ({POP1} vs {POP2})", fontweight="bold", loc="left")
    ax.legend(fontsize=7)

    # Panel (b): Box plot
    ax = axes[1]
    box_data = [impact_data[imp]["fst_values"] for imp in IMPACTS if impact_data[imp]["fst_values"]]
    box_labels = [imp for imp in IMPACTS if impact_data[imp]["fst_values"]]
    bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True, showfliers=False)
    for patch, imp in zip(bp["boxes"], box_labels):
        patch.set_facecolor(colors[imp])
    ax.set_ylabel("FST")
    ax.set_title("b  FST by consequence impact", fontweight="bold", loc="left")

    # Add significance annotations
    y_max = max(np.percentile(d, 75) for d in box_data if len(d) > 0) * 1.1
    for key, mw in mw_results.items():
        if mw["p"] < 0.001:
            ax.text(0.5, 0.95, f"Mann-Whitney: p < 0.001", transform=ax.transAxes,
                    fontsize=7, ha="center", va="top")
            break

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig14_consequence_bias.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"\n=== Rice Inv.R14: Consequence-Level FST Bias ===")
    print(f"Total time: {total_time:.1f}s")
    for imp in IMPACTS:
        d = impact_data[imp]
        print(f"  {imp:10s}: n={d['n_variants']:,}, mean_fst={d['mean_fst']:.4f}")
    for key, mw in mw_results.items():
        print(f"  {key}: U={mw['U']:.0f}, p={mw['p']:.2e}")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
