#!/usr/bin/env python3
"""Rice Inv.R03: Pairwise Fst Landscape + UPGMA Tree.

66 pairwise Fst values reveal the indica-japonica split and subpopulation
structure. UPGMA dendrogram recovers Wang et al. 2018 classification.

Data source: results/rice/rice_interpretation_results.json → fst_matrix, dxy_matrix, tree
Output:      data/results/rice/rice_inv03_fst_landscape.{tsv,json}
Figure:      data/results/rice/figures/rice_fig03_fst_landscape.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

POP_ORDER = [
    "XI-1A", "XI-1B", "XI-2", "XI-3", "XI-adm",
    "cA-Aus", "cB-Bas",
    "GJ-trp", "GJ-sbtrp", "GJ-tmp", "GJ-adm",
    "admix",
]

POP_GROUPS = {
    "XI-1A": "Indica", "XI-1B": "Indica", "XI-2": "Indica",
    "XI-3": "Indica", "XI-adm": "Indica",
    "GJ-tmp": "Japonica", "GJ-trp": "Japonica",
    "GJ-sbtrp": "Japonica", "GJ-adm": "Japonica",
    "cA-Aus": "Aus/Basmati", "cB-Bas": "Aus/Basmati",
    "admix": "Admixed",
}


def main():
    with open(ROOT / "results" / "rice" / "rice_interpretation_results.json") as f:
        interp = json.load(f)

    fst_matrix = interp["fst_matrix"]
    dxy_matrix = interp["dxy_matrix"]
    tree = interp["tree"]

    pops = POP_ORDER
    n = len(pops)

    # Build symmetric Fst and Dxy matrices
    fst_mat = np.zeros((n, n))
    dxy_mat = np.zeros((n, n))
    for i, p1 in enumerate(pops):
        for j, p2 in enumerate(pops):
            if i == j:
                continue
            key = f"{p1}_vs_{p2}"
            key_rev = f"{p2}_vs_{p1}"
            fst_entry = fst_matrix.get(key, fst_matrix.get(key_rev, None))
            if fst_entry is not None:
                if isinstance(fst_entry, dict):
                    fst_mat[i, j] = fst_entry.get("mean_fst", 0)
                else:
                    fst_mat[i, j] = float(fst_entry)

            dxy_entry = dxy_matrix.get(key, dxy_matrix.get(key_rev, None))
            if dxy_entry is not None:
                if isinstance(dxy_entry, dict):
                    dxy_mat[i, j] = dxy_entry.get("mean_dxy", 0)
                else:
                    dxy_mat[i, j] = float(dxy_entry)

    # Symmetrize
    fst_mat = (fst_mat + fst_mat.T) / 2
    dxy_mat = (dxy_mat + dxy_mat.T) / 2

    # Fst vs Dxy correlation (upper triangle)
    upper = np.triu_indices(n, k=1)
    fst_upper = fst_mat[upper]
    dxy_upper = dxy_mat[upper]
    valid = (fst_upper > 0) & (dxy_upper > 0)
    if valid.sum() > 2:
        rho, pval = stats.spearmanr(fst_upper[valid], dxy_upper[valid])
    else:
        rho, pval = float("nan"), float("nan")

    # Most/least divergent pairs
    pair_labels = [(pops[i], pops[j]) for i, j in zip(upper[0], upper[1])]
    nonzero = fst_upper > 0
    max_idx = np.argmax(fst_upper)
    min_idx = np.argmin(np.where(nonzero, fst_upper, np.inf))

    # Save TSV
    with open(RESULTS_DIR / "rice_inv03_fst_landscape.tsv", "w") as f:
        f.write("pop1\tpop2\tfst\tdxy\n")
        for idx, (p1, p2) in enumerate(pair_labels):
            f.write(f"{p1}\t{p2}\t{fst_upper[idx]:.6f}\t{dxy_upper[idx]:.6f}\n")

    summary = {
        "n_pairs": len(pair_labels),
        "most_divergent": {"pair": list(pair_labels[max_idx]), "fst": float(fst_upper[max_idx])},
        "least_divergent": {"pair": list(pair_labels[min_idx]), "fst": float(fst_upper[fst_upper > 0][min_idx])},
        "mean_fst": float(np.mean(fst_upper)),
        "fst_vs_dxy_rho": float(rho),
        "fst_vs_dxy_pval": float(pval),
        "newick": tree.get("newick", ""),
    }
    with open(RESULTS_DIR / "rice_inv03_fst_landscape.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 3 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # Panel (a): Fst heatmap
    ax = axes[0]
    mask = np.zeros_like(fst_mat, dtype=bool)
    np.fill_diagonal(mask, True)
    masked = np.ma.array(fst_mat, mask=mask)
    im = ax.imshow(masked, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(n))
    ax.set_xticklabels(pops, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels(pops, fontsize=7)
    plt.colorbar(im, ax=ax, label="Pairwise Fst (W&C)", shrink=0.8)
    ax.set_title("a  Pairwise Fst heatmap", fontweight="bold", loc="left")

    # Annotate cells
    for i in range(n):
        for j in range(n):
            if i != j:
                val = fst_mat[i, j]
                color = "white" if val > 0.45 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=5, color=color)

    # Panel (b): UPGMA dendrogram
    ax = axes[1]
    # Use Fst as distance for hierarchical clustering
    condensed = squareform(fst_mat, checks=False)
    Z = linkage(condensed, method="average")  # UPGMA = average linkage
    dn = dendrogram(Z, labels=pops, ax=ax, orientation="right",
                    leaf_font_size=8, color_threshold=0.3)
    ax.set_xlabel("Fst distance (UPGMA)")
    ax.set_title("b  UPGMA dendrogram", fontweight="bold", loc="left")

    # Panel (c): Fst vs Dxy scatter
    ax = axes[2]
    ax.scatter(dxy_upper, fst_upper, s=30, alpha=0.7, edgecolors="black", linewidth=0.3)
    ax.set_xlabel("Dxy (absolute divergence)")
    ax.set_ylabel("Fst (relative divergence)")
    ax.set_title(f"c  Fst vs Dxy (rho={rho:.2f})", fontweight="bold", loc="left")

    # Label extremes
    ax.annotate(f"{pair_labels[max_idx][0]}-{pair_labels[max_idx][1]}",
                (dxy_upper[max_idx], fst_upper[max_idx]),
                fontsize=7, xytext=(5, 5), textcoords="offset points")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig03_fst_landscape.png", dpi=200, bbox_inches="tight")
    plt.close()

    print("=== Rice Inv.R03: Pairwise Fst Landscape ===")
    print(f"Pairs: {len(pair_labels)}")
    print(f"Most divergent: {pair_labels[max_idx][0]} vs {pair_labels[max_idx][1]} (Fst={fst_upper[max_idx]:.3f})")
    print(f"Mean Fst: {np.mean(fst_upper):.3f}")
    print(f"Fst vs Dxy: rho={rho:.3f}")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
