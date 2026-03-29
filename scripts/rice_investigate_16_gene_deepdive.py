#!/usr/bin/env python3
"""Rice Inv.R16: Domestication Gene Deep-Dive (GW5 + Hd1).

Detailed multi-statistic profile at two key domestication loci,
analogous to human Inv.04 (KCNE1 deep-dive).

Data source: results/rice/rice_deep_integration_results.json → sweeps
             results/rice/rice_interpretation_results.json → annotate
Output:      data/results/rice/rice_inv16_gene_deepdive.{tsv,json}
Figure:      data/results/rice/figures/rice_fig16_gene_deepdive.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

TARGET_GENES = ["GW5", "Hd1", "Wx", "PROG1", "OsC1", "DRO1"]

POP_ORDER = [
    "XI-1A", "XI-1B", "XI-2", "XI-3", "XI-adm",
    "cA-Aus", "cB-Bas",
    "GJ-trp", "GJ-sbtrp", "GJ-tmp", "GJ-adm",
    "admix",
]


def main():
    with open(ROOT / "results" / "rice" / "rice_deep_integration_results.json") as f:
        di = json.load(f)

    sweeps = di["sweeps"]
    rows = []

    for gene in TARGET_GENES:
        if gene not in sweeps:
            print(f"  {gene}: not found in sweep data")
            continue

        gd = sweeps[gene]
        for pop in POP_ORDER:
            if pop not in gd["populations"]:
                continue
            pd = gd["populations"][pop]
            rows.append({
                "gene": gene,
                "function": gd.get("function", ""),
                "chr": gd.get("chr", ""),
                "pos": gd.get("pos", 0),
                "population": pop,
                "h1": pd.get("h1", 0),
                "h12": pd.get("h12", 0),
                "h2_h1": pd.get("h2_h1", 0),
                "hap_diversity": pd.get("hap_diversity", 0),
                "sweep_type": pd.get("sweep_type", "none"),
            })

    # Save TSV
    header = ["gene", "function", "chr", "pos", "population",
              "h1", "h12", "h2_h1", "hap_diversity", "sweep_type"]
    with open(RESULTS_DIR / "rice_inv16_gene_deepdive.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")

    # Gene-level summary
    gene_summary = {}
    for gene in TARGET_GENES:
        if gene not in sweeps:
            continue
        gd = sweeps[gene]
        pops_data = gd["populations"]
        h12_vals = [pops_data[p]["h12"] for p in POP_ORDER if p in pops_data]
        n_swept = sum(1 for p in pops_data.values() if p.get("sweep_type") != "none")
        gene_summary[gene] = {
            "function": gd.get("function", ""),
            "chr": gd.get("chr", ""),
            "pos": gd.get("pos", 0),
            "mean_h12": float(np.mean(h12_vals)) if h12_vals else 0,
            "max_h12": float(np.max(h12_vals)) if h12_vals else 0,
            "n_swept_pops": n_swept,
            "n_total_pops": len(pops_data),
        }

    with open(RESULTS_DIR / "rice_inv16_gene_deepdive.json", "w") as f:
        json.dump(gene_summary, f, indent=2)

    # ── Figure: one panel per gene ──
    n_genes = len([g for g in TARGET_GENES if g in sweeps])
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()

    for idx, gene in enumerate([g for g in TARGET_GENES if g in sweeps]):
        ax = axes[idx]
        gd = sweeps[gene]
        pops_present = [p for p in POP_ORDER if p in gd["populations"]]
        h12_vals = [gd["populations"][p]["h12"] for p in pops_present]
        sweep_types = [gd["populations"][p].get("sweep_type", "none") for p in pops_present]

        colors = []
        for st in sweep_types:
            if st == "hard":
                colors.append("#d73027")
            elif st == "soft":
                colors.append("#fc8d59")
            else:
                colors.append("#4575b4")

        bars = ax.bar(range(len(pops_present)), h12_vals, color=colors,
                      edgecolor="black", linewidth=0.5)
        ax.set_xticks(range(len(pops_present)))
        ax.set_xticklabels(pops_present, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("H12")
        ax.set_title(f"{gene} — {gd.get('function', '')}", fontweight="bold", fontsize=9)
        ax.axhline(0.05, color="gray", linestyle="--", linewidth=0.5)

    # Hide unused
    for idx in range(n_genes, len(axes)):
        axes[idx].set_visible(False)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#d73027", label="Hard sweep"),
        Patch(facecolor="#fc8d59", label="Soft sweep"),
        Patch(facecolor="#4575b4", label="No sweep"),
    ]
    fig.legend(handles=legend_elements, loc="lower right", fontsize=8)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig16_gene_deepdive.png", dpi=200, bbox_inches="tight")
    plt.close()

    print("=== Rice Inv.R16: Domestication Gene Deep-Dive ===")
    for gene, gs in gene_summary.items():
        print(f"  {gene:8s} ({gs['function'][:25]:25s}): mean_H12={gs['mean_h12']:.4f}, "
              f"max_H12={gs['max_h12']:.4f}, swept={gs['n_swept_pops']}/{gs['n_total_pops']}")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
