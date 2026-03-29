#!/usr/bin/env python3
"""Rice Inv.R01: Evolutionary Fingerprint PCA.

12 populations profiled by 9 evolutionary features, projected onto
PCA space to reveal the indica-japonica axis and bottleneck gradient.

Data source: results/rice/rice_deep_integration_results.json → fingerprints
Output:      data/results/rice/rice_inv01_evo_fingerprint.{tsv,json}
Figure:      data/results/rice/figures/rice_fig01_evo_fingerprint.png
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Population group colors
GROUP_COLORS = {
    "Indica": "#e41a1c",
    "Japonica": "#377eb8",
    "Aus/Basmati": "#4daf4a",
    "Admixed": "#984ea3",
}

POP_GROUPS = {
    "XI-1A": "Indica", "XI-1B": "Indica", "XI-2": "Indica",
    "XI-3": "Indica", "XI-adm": "Indica",
    "GJ-tmp": "Japonica", "GJ-trp": "Japonica",
    "GJ-sbtrp": "Japonica", "GJ-adm": "Japonica",
    "cA-Aus": "Aus/Basmati", "cB-Bas": "Aus/Basmati",
    "admix": "Admixed",
}


def main():
    with open(ROOT / "results" / "rice" / "rice_deep_integration_results.json") as f:
        di = json.load(f)

    fp = di["fingerprints"]
    features = fp["feature_names"]
    profiles = fp["profiles"]
    pca = fp["pca"]

    pops = sorted(profiles.keys())
    n_features = len(features)

    # Build feature matrix for TSV
    rows = []
    for pop in pops:
        p = profiles[pop]
        coord = pca["coordinates"][pop]
        row = {
            "population": pop,
            "group": POP_GROUPS.get(pop, "Unknown"),
            "PC1": coord["PC1"],
            "PC2": coord["PC2"],
        }
        for feat in features:
            row[feat] = p.get(feat, float("nan"))
        rows.append(row)

    # Save TSV
    header = ["population", "group", "PC1", "PC2"] + features
    with open(RESULTS_DIR / "rice_inv01_evo_fingerprint.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")

    # Save JSON summary
    var_exp = pca["variance_explained"]
    loadings = pca["loadings"]
    summary = {
        "n_populations": len(pops),
        "n_features": n_features,
        "features": features,
        "variance_explained": var_exp,
        "pc1_top_loadings": sorted(
            [(feat, loadings[feat][0]) for feat in loadings],
            key=lambda x: abs(x[1]), reverse=True
        )[:5],
        "pc2_top_loadings": sorted(
            [(feat, loadings[feat][1]) for feat in loadings],
            key=lambda x: abs(x[1]), reverse=True
        )[:5],
        "indica_centroid_pc1": np.mean([
            pca["coordinates"][p]["PC1"] for p in pops if POP_GROUPS.get(p) == "Indica"
        ]),
        "japonica_centroid_pc1": np.mean([
            pca["coordinates"][p]["PC1"] for p in pops if POP_GROUPS.get(p) == "Japonica"
        ]),
    }
    with open(RESULTS_DIR / "rice_inv01_evo_fingerprint.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 3 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # Panel (a): PCA biplot
    ax = axes[0]
    for pop in pops:
        coord = pca["coordinates"][pop]
        grp = POP_GROUPS.get(pop, "Unknown")
        color = GROUP_COLORS.get(grp, "gray")
        ax.scatter(coord["PC1"], coord["PC2"], c=color, s=100, edgecolors="black",
                   linewidth=0.5, zorder=3)
        ax.annotate(pop, (coord["PC1"], coord["PC2"]),
                    fontsize=7, ha="center", va="bottom",
                    xytext=(0, 6), textcoords="offset points")

    # Loading arrows
    scale = 3.0
    for feat in features:
        lx, ly = loadings[feat][0] * scale, loadings[feat][1] * scale
        ax.annotate("", xy=(lx, ly), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=1))
        ax.text(lx * 1.15, ly * 1.15, feat, fontsize=6, color="gray",
                ha="center", va="center")

    ax.set_xlabel(f"PC1 ({var_exp[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_exp[1]*100:.1f}%)")
    ax.axhline(0, color="lightgray", lw=0.5)
    ax.axvline(0, color="lightgray", lw=0.5)
    ax.set_title("a  Evolutionary fingerprint PCA", fontweight="bold", loc="left")

    # Legend
    for grp, color in GROUP_COLORS.items():
        ax.scatter([], [], c=color, s=60, label=grp, edgecolors="black", linewidth=0.5)
    ax.legend(fontsize=7, loc="lower left")

    # Panel (b): Feature heatmap (z-scored)
    ax = axes[1]
    feat_matrix = np.array([[profiles[pop].get(feat, 0) for feat in features] for pop in pops])
    z_matrix = (feat_matrix - feat_matrix.mean(axis=0)) / (feat_matrix.std(axis=0) + 1e-10)

    im = ax.imshow(z_matrix, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(n_features))
    ax.set_xticklabels(features, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(pops)))
    ax.set_yticklabels(pops, fontsize=7)
    plt.colorbar(im, ax=ax, label="Z-score", shrink=0.8)
    ax.set_title("b  Feature z-scores", fontweight="bold", loc="left")

    # Panel (c): Dendrogram
    ax = axes[2]
    Z = linkage(z_matrix, method="ward")
    dn = dendrogram(Z, labels=pops, ax=ax, orientation="right",
                    leaf_font_size=8, color_threshold=0)
    ax.set_xlabel("Ward distance")
    ax.set_title("c  Population dendrogram", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig01_evo_fingerprint.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Print summary
    print("=== Rice Inv.R01: Evolutionary Fingerprint PCA ===")
    print(f"Populations: {len(pops)}")
    print(f"Features: {n_features} ({', '.join(features)})")
    print(f"PC1: {var_exp[0]*100:.1f}%, PC2: {var_exp[1]*100:.1f}%")
    print(f"\nPC1 top loadings (indica-japonica axis):")
    for feat, val in summary["pc1_top_loadings"]:
        print(f"  {feat:25s} {val:+.3f}")
    print(f"\nIndica centroid PC1: {summary['indica_centroid_pc1']:.2f}")
    print(f"Japonica centroid PC1: {summary['japonica_centroid_pc1']:.2f}")
    print(f"\nSaved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
