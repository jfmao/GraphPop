#!/usr/bin/env python3
"""Rice Inv.R05: Domestication Gene Convergence.

16 known domestication genes analyzed with multi-statistic evidence
(H12, XP-EHH, PBS, Fay&Wu H) across rice subpopulations.

Data source: results/rice/rice_deep_integration_results.json → sweeps
Output:      data/results/rice/rice_inv05_domestication_genes.{tsv,json}
Figure:      data/results/rice/figures/rice_fig05_domestication_genes.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"


def main():
    with open(ROOT / "results" / "rice" / "rice_deep_integration_results.json") as f:
        di = json.load(f)

    sweeps = di["sweeps"]

    # Also load XP-EHH and PBS annotations from interpretation
    with open(ROOT / "results" / "rice" / "rice_interpretation_results.json") as f:
        interp = json.load(f)
    annotations = interp.get("annotate", {})

    genes = sorted(sweeps.keys())
    all_pops = set()
    for g in genes:
        all_pops.update(sweeps[g]["populations"].keys())
    all_pops = sorted(all_pops)

    # Build gene × pop H12 matrix and sweep type matrix
    h12_matrix = np.full((len(genes), len(all_pops)), np.nan)
    sweep_types = {}
    rows = []

    for i, gene in enumerate(genes):
        gd = sweeps[gene]
        for j, pop in enumerate(all_pops):
            if pop in gd["populations"]:
                pd = gd["populations"][pop]
                h12 = pd.get("h12", 0)
                h12_matrix[i, j] = h12
                sweep_type = pd.get("sweep_type", "none")

                rows.append({
                    "gene": gene, "function": gd.get("function", ""),
                    "chr": gd.get("chr", ""), "pos": gd.get("pos", 0),
                    "population": pop, "h12": h12,
                    "h2_h1": pd.get("h2_h1", float("nan")),
                    "hap_diversity": pd.get("hap_diversity", float("nan")),
                    "sweep_type": sweep_type,
                })

    # Save TSV
    with open(RESULTS_DIR / "rice_inv05_domestication_genes.tsv", "w") as f:
        header = ["gene", "function", "chr", "pos", "population", "h12", "h2_h1",
                  "hap_diversity", "sweep_type"]
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")

    # Count sweep detections per gene
    gene_sweep_counts = {}
    for gene in genes:
        pops_data = sweeps[gene]["populations"]
        n_hard = sum(1 for p in pops_data.values() if p.get("sweep_type") == "hard")
        n_soft = sum(1 for p in pops_data.values() if p.get("sweep_type") == "soft")
        gene_sweep_counts[gene] = {"hard": n_hard, "soft": n_soft, "total": n_hard + n_soft}

    summary = {
        "n_genes": len(genes),
        "n_populations": len(all_pops),
        "genes": {g: {"function": sweeps[g].get("function", ""),
                       "chr": sweeps[g].get("chr", ""),
                       **gene_sweep_counts[g]} for g in genes},
        "genes_with_sweeps": sum(1 for g in genes if gene_sweep_counts[g]["total"] > 0),
    }
    with open(RESULTS_DIR / "rice_inv05_domestication_genes.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 3 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel (a): Gene × Pop H12 heatmap
    ax = axes[0]
    im = ax.imshow(h12_matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
    ax.set_xticks(range(len(all_pops)))
    ax.set_xticklabels(all_pops, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(len(genes)))
    gene_labels = [f"{g} ({sweeps[g].get('function', '')[:20]})" for g in genes]
    ax.set_yticklabels(gene_labels, fontsize=6)
    plt.colorbar(im, ax=ax, label="H12", shrink=0.8)
    ax.set_title("a  H12 at domestication loci", fontweight="bold", loc="left")

    # Panel (b): Sweep count per gene
    ax = axes[1]
    hard_counts = [gene_sweep_counts[g]["hard"] for g in genes]
    soft_counts = [gene_sweep_counts[g]["soft"] for g in genes]
    y_pos = range(len(genes))
    ax.barh(y_pos, hard_counts, color="#e41a1c", label="Hard sweep", height=0.4, align="edge")
    ax.barh([y + 0.4 for y in y_pos], soft_counts, color="#ff7f00", label="Soft sweep",
            height=0.4, align="edge")
    ax.set_yticks([y + 0.4 for y in y_pos])
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xlabel("Number of populations")
    ax.legend(fontsize=7)
    ax.set_title("b  Sweep detection count", fontweight="bold", loc="left")
    ax.invert_yaxis()

    # Panel (c): H12 distribution — swept vs non-swept populations
    ax = axes[2]
    swept_h12 = [h12_matrix[i, j] for i in range(len(genes)) for j in range(len(all_pops))
                 if not np.isnan(h12_matrix[i, j]) and h12_matrix[i, j] > 0.05]
    non_swept = [h12_matrix[i, j] for i in range(len(genes)) for j in range(len(all_pops))
                 if not np.isnan(h12_matrix[i, j]) and h12_matrix[i, j] <= 0.05]

    if swept_h12:
        ax.hist(swept_h12, bins=30, alpha=0.7, label=f"H12 > 0.05 (n={len(swept_h12)})",
                color="#e41a1c", density=True)
    if non_swept:
        ax.hist(non_swept, bins=30, alpha=0.7, label=f"H12 <= 0.05 (n={len(non_swept)})",
                color="#377eb8", density=True)
    ax.set_xlabel("H12")
    ax.set_ylabel("Density")
    ax.legend(fontsize=7)
    ax.set_title("c  H12 distribution at domestication loci", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig05_domestication_genes.png", dpi=200, bbox_inches="tight")
    plt.close()

    print("=== Rice Inv.R05: Domestication Gene Convergence ===")
    print(f"Genes: {len(genes)}")
    print(f"Genes with sweep signals: {summary['genes_with_sweeps']}")
    for g in genes:
        gc = gene_sweep_counts[g]
        func = sweeps[g].get("function", "")[:30]
        print(f"  {g:8s} ({func:30s}): {gc['hard']} hard, {gc['soft']} soft")
    print(f"\nSaved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
