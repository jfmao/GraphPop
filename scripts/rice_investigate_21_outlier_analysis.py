#!/usr/bin/env python3
"""Rice Investigation 21: Outlier Detection in Individual Process Space.

Analogous to human Inv.17 Analysis 4. Identifies rice accessions that deviate
most from their population centroid in the statistical PCA space computed by
Rice Inv.20.

For each individual:
  - Compute z-score distance from population centroid in PC1/PC2 space
  - Flag outliers at z > 3 and z > 5
  - Report top outliers with their population and feature values

This reveals individual-level heterogeneity within rice subgroups, potentially
identifying misclassified accessions, introgression events, or unusual breeding
history.

Data sources:
  - data/results/rice_inv20_individual_trajectory.tsv (from Rice Inv.20)

Output:
  data/results/rice_inv21_outlier_analysis.tsv / .json
  data/results/figures/rice_fig21_outlier_analysis.png

Usage:
  conda run -n graphevo python -u scripts/rice_investigate_21_outlier_analysis.py
"""

import csv
import json
import time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

RESULTS = Path("data/results")
FIG_DIR = RESULTS / "figures"
INPUT_TSV = RESULTS / "rice_inv20_individual_trajectory.tsv"
OUT_TSV   = RESULTS / "rice_inv21_outlier_analysis.tsv"
OUT_JSON  = RESULTS / "rice_inv21_outlier_analysis.json"

GROUP_COLORS = {
    "Japonica":  "#5C8AE0",
    "Indica":    "#E05C5C",
    "Aus":       "#5CC45C",
    "Basmati":   "#E0A85C",
    "Admixed":   "#A85CE0",
    "Unknown":   "#888888",
}


# ── Load Rice Inv.20 individual trajectory ────────────────────────────────────

def load_individual_trajectory():
    """Load the TSV output from Rice Inv.20."""
    rows = list(csv.DictReader(open(INPUT_TSV), delimiter="\t"))
    for r in rows:
        for col in ["het_rate", "hom_alt_rate", "rare_burden", "private_burden",
                     "froh_genome", "stat_pc1", "stat_pc2", "stat_pc3", "umap1", "umap2"]:
            r[col] = float(r[col])
        r["n_roh_genome"] = int(r["n_roh_genome"])
    print(f"  Loaded {len(rows)} samples from {INPUT_TSV}", flush=True)
    return rows


# ── Outlier Analysis ──────────────────────────────────────────────────────────

def run_outlier_analysis(rows: list) -> dict:
    """Identify individuals deviating most from their population centroid in PC space."""
    print("=== Rice Outlier Analysis: Individuals in Process Space ===\n", flush=True)

    # Group by population
    pop_members = defaultdict(list)
    for i, r in enumerate(rows):
        pop_members[r["population"]].append(i)

    # Compute population centroids in PC1/PC2 space
    pop_centroid = {}
    pop_std = {}
    for pop, idx in pop_members.items():
        pc1 = np.array([rows[i]["stat_pc1"] for i in idx])
        pc2 = np.array([rows[i]["stat_pc2"] for i in idx])
        pop_centroid[pop] = (pc1.mean(), pc2.mean())
        pop_std[pop] = (pc1.std() + 1e-9, pc2.std() + 1e-9)

    # Compute per-individual deviation from population centroid (z-score distance)
    outlier_records = []
    for i, r in enumerate(rows):
        pop = r["population"]
        if pop not in pop_centroid:
            continue
        cx, cy = pop_centroid[pop]
        sx, sy = pop_std[pop]
        z_pc1 = (r["stat_pc1"] - cx) / sx
        z_pc2 = (r["stat_pc2"] - cy) / sy
        z_dist = np.sqrt(z_pc1**2 + z_pc2**2)
        outlier_records.append({
            "sampleId":       r["sampleId"],
            "population":     pop,
            "group":          r["group"],
            "z_distance":     round(float(z_dist), 4),
            "z_pc1":          round(float(z_pc1), 4),
            "z_pc2":          round(float(z_pc2), 4),
            "het_rate":       r["het_rate"],
            "hom_alt_rate":   r["hom_alt_rate"],
            "rare_burden":    r["rare_burden"],
            "private_burden": r["private_burden"],
            "froh_genome":    r["froh_genome"],
            "stat_pc1":       r["stat_pc1"],
            "stat_pc2":       r["stat_pc2"],
            "umap1":          r["umap1"],
            "umap2":          r["umap2"],
        })

    outlier_records.sort(key=lambda x: -x["z_distance"])

    # Top outliers
    top_outliers = outlier_records[:30]
    print("  Top 15 outlier accessions (z-distance from population centroid):")
    for o in top_outliers[:15]:
        print(f"    {o['sampleId']:20s} ({o['population']:8s}/{o['group'][:3]})  "
              f"z={o['z_distance']:.2f}  het={o['het_rate']:.6f}  "
              f"hom_alt={o['hom_alt_rate']:.6f}  rare={o['rare_burden']:.0f}")

    # Count outliers at thresholds
    n_z3 = sum(1 for o in outlier_records if o["z_distance"] > 3)
    n_z5 = sum(1 for o in outlier_records if o["z_distance"] > 5)
    print(f"\n  Outliers z > 3: {n_z3} ({n_z3/len(outlier_records)*100:.1f}%)")
    print(f"  Outliers z > 5: {n_z5} ({n_z5/len(outlier_records)*100:.1f}%)")

    # Distribution of z-distances by group
    grp_zdist = defaultdict(list)
    for o in outlier_records:
        grp_zdist[o["group"]].append(o["z_distance"])
    print("\n  Mean z-distance by group:")
    for grp, vals in sorted(grp_zdist.items()):
        print(f"    {grp:12s}: mean={np.mean(vals):.3f}  std={np.std(vals):.3f}  "
              f"max={np.max(vals):.3f}  n_z3={sum(1 for v in vals if v > 3)}")

    # Distribution by population
    pop_zdist = defaultdict(list)
    for o in outlier_records:
        pop_zdist[o["population"]].append(o["z_distance"])
    print("\n  Mean z-distance by population:")
    for pop, vals in sorted(pop_zdist.items()):
        print(f"    {pop:10s}: mean={np.mean(vals):.3f}  std={np.std(vals):.3f}  "
              f"max={np.max(vals):.3f}  n={len(vals)}  n_z3={sum(1 for v in vals if v > 3)}")

    # What drives outliers? Correlation of z_distance with features
    z_dists = np.array([o["z_distance"] for o in outlier_records])
    print("\n  z_distance correlations with features:")
    for feat in ["het_rate", "hom_alt_rate", "rare_burden", "private_burden", "froh_genome"]:
        vals = np.array([o[feat] for o in outlier_records])
        r, p = spearmanr(z_dists, vals)
        print(f"    z_distance vs {feat:20s}: rho={r:.3f}, p={p:.4f}")

    return {
        "outlier_records": outlier_records,
        "pop_centroid":    {p: list(v) for p, v in pop_centroid.items()},
        "pop_std":         {p: list(v) for p, v in pop_std.items()},
        "group_zdist":     {g: {"mean": float(np.mean(v)), "std": float(np.std(v)),
                                "max": float(np.max(v)),
                                "n_z3": sum(1 for x in v if x > 3),
                                "n_z5": sum(1 for x in v if x > 5)}
                            for g, v in grp_zdist.items()},
        "n_outliers_z3":   n_z3,
        "n_outliers_z5":   n_z5,
    }


# ── Visualization ─────────────────────────────────────────────────────────────

def make_figure(outlier_records: list, pop_centroid: dict, pop_std: dict):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(22, 12))
    fig.patch.set_facecolor("#0D0D0D")
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.35,
                          left=0.07, right=0.97, top=0.91, bottom=0.07)
    axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]
    for ax in axes:
        ax.set_facecolor("#161616")
        for spine in ax.spines.values():
            spine.set_color("#444")

    z_dists = np.array([o["z_distance"] for o in outlier_records])
    groups = [o["group"] for o in outlier_records]
    colors = [GROUP_COLORS.get(g, "#888") for g in groups]
    pc1 = np.array([o["stat_pc1"] for o in outlier_records])
    pc2 = np.array([o["stat_pc2"] for o in outlier_records])
    umap1 = np.array([o["umap1"] for o in outlier_records])
    umap2 = np.array([o["umap2"] for o in outlier_records])
    het = np.array([o["het_rate"] for o in outlier_records])
    hom_alt = np.array([o["hom_alt_rate"] for o in outlier_records])
    rare = np.array([o["rare_burden"] for o in outlier_records])

    # ── A: PC1 vs PC2, color by z-distance ───────────────────────────────
    ax = axes[0]
    ax.set_title("A  PC Space, colored by z-distance from centroid", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    sc = ax.scatter(pc1, pc2, c=z_dists, cmap="hot", s=6, alpha=0.6,
                    linewidths=0, vmin=0, vmax=min(8, z_dists.max()))
    # Mark z>5 outliers
    mask_z5 = z_dists > 5
    if mask_z5.any():
        ax.scatter(pc1[mask_z5], pc2[mask_z5], facecolors="none",
                   edgecolors="cyan", s=40, linewidths=1.2, zorder=3, label=f"z>5 (n={mask_z5.sum()})")
        ax.legend(fontsize=7, facecolor="#222", labelcolor="#CCC", framealpha=0.8)
    plt.colorbar(sc, ax=ax, label="z-distance", pad=0.02)
    ax.set_xlabel("PC1", color="#AAA", fontsize=8)
    ax.set_ylabel("PC2", color="#AAA", fontsize=8)
    ax.tick_params(colors="#888", labelsize=7)

    # ── B: UMAP, color by z-distance ─────────────────────────────────────
    ax = axes[1]
    ax.set_title("B  UMAP Space, colored by z-distance", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    sc = ax.scatter(umap1, umap2, c=z_dists, cmap="hot", s=6, alpha=0.6,
                    linewidths=0, vmin=0, vmax=min(8, z_dists.max()))
    if mask_z5.any():
        ax.scatter(umap1[mask_z5], umap2[mask_z5], facecolors="none",
                   edgecolors="cyan", s=40, linewidths=1.2, zorder=3)
    plt.colorbar(sc, ax=ax, label="z-distance", pad=0.02)
    ax.set_xlabel("UMAP 1", color="#AAA", fontsize=8)
    ax.set_ylabel("UMAP 2", color="#AAA", fontsize=8)
    ax.tick_params(colors="#888", labelsize=7)

    # ── C: z-distance distribution by group ──────────────────────────────
    ax = axes[2]
    ax.set_title("C  z-distance Distribution by Group", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    group_order = ["Japonica", "Indica", "Aus", "Basmati", "Admixed", "Unknown"]
    grp_data = defaultdict(list)
    for o in outlier_records:
        grp_data[o["group"]].append(o["z_distance"])
    plot_data = []
    plot_labels = []
    plot_colors = []
    for grp in group_order:
        if grp in grp_data and len(grp_data[grp]) > 3:
            plot_data.append(grp_data[grp])
            plot_labels.append(grp)
            plot_colors.append(GROUP_COLORS.get(grp, "#888"))

    if plot_data:
        vp = ax.violinplot(plot_data, positions=range(len(plot_data)),
                           showmedians=True, showextrema=False)
        for i, body in enumerate(vp["bodies"]):
            body.set_facecolor(plot_colors[i])
            body.set_alpha(0.5)
            body.set_edgecolor("none")
        vp["cmedians"].set_color("#EEE")
        vp["cmedians"].set_linewidth(2)
        ax.set_xticks(range(len(plot_labels)))
        ax.set_xticklabels(plot_labels, color="#CCC", fontsize=7, rotation=15)
    ax.axhline(3, color="#FF6666", ls="--", lw=0.8, alpha=0.7, label="z=3")
    ax.axhline(5, color="#FF3333", ls="--", lw=0.8, alpha=0.7, label="z=5")
    ax.legend(fontsize=7, facecolor="#222", labelcolor="#CCC", framealpha=0.8)
    ax.set_ylabel("z-distance", color="#AAA", fontsize=8)
    ax.tick_params(axis="y", colors="#888", labelsize=7)

    # ── D: z-distance vs het_rate ─────────────────────────────────────────
    ax = axes[3]
    ax.set_title("D  z-distance vs Het Rate", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    ax.scatter(het, z_dists, c=colors, s=6, alpha=0.5, linewidths=0)
    r, p = spearmanr(het, z_dists)
    ax.set_xlabel("Het rate (Chr1)", color="#AAA", fontsize=8)
    ax.set_ylabel("z-distance", color="#AAA", fontsize=8)
    ax.axhline(3, color="#FF6666", ls="--", lw=0.8, alpha=0.5)
    ax.text(0.02, 0.95, f"rho={r:.3f}, p={p:.2e}", transform=ax.transAxes,
            color="#CCC", fontsize=7, va="top")
    ax.tick_params(colors="#888", labelsize=7)

    # ── E: z-distance vs rare_burden ──────────────────────────────────────
    ax = axes[4]
    ax.set_title("E  z-distance vs Rare Burden", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    ax.scatter(rare, z_dists, c=colors, s=6, alpha=0.5, linewidths=0)
    r2, p2 = spearmanr(rare, z_dists)
    ax.set_xlabel("Rare burden (Chr1)", color="#AAA", fontsize=8)
    ax.set_ylabel("z-distance", color="#AAA", fontsize=8)
    ax.axhline(3, color="#FF6666", ls="--", lw=0.8, alpha=0.5)
    ax.text(0.02, 0.95, f"rho={r2:.3f}, p={p2:.2e}", transform=ax.transAxes,
            color="#CCC", fontsize=7, va="top")
    ax.tick_params(colors="#888", labelsize=7)

    # ── F: Top 20 outliers bar chart ──────────────────────────────────────
    ax = axes[5]
    ax.set_title("F  Top 20 Outliers by z-distance", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    top20 = outlier_records[:20]
    names = [f"{o['sampleId'][:16]} ({o['population']})" for o in top20]
    z_vals = [o["z_distance"] for o in top20]
    bar_colors = [GROUP_COLORS.get(o["group"], "#888") for o in top20]
    y_pos = range(len(top20))
    ax.barh(y_pos, z_vals[::-1], color=bar_colors[::-1], edgecolor="#333", height=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names[::-1], fontsize=6, color="#CCC")
    ax.set_xlabel("z-distance", color="#AAA", fontsize=8)
    ax.axvline(3, color="#FF6666", ls="--", lw=0.8, alpha=0.7)
    ax.axvline(5, color="#FF3333", ls="--", lw=0.8, alpha=0.7)
    ax.tick_params(axis="x", colors="#888", labelsize=7)

    fig.suptitle(
        "Rice Outlier Detection -- 3,024 Accessions, z-distance from Population Centroid in PC Space",
        color="#EEE", fontsize=13, y=0.96,
    )

    # Legend for groups
    legend_handles = [mpatches.Patch(color=col, label=grp) for grp, col in GROUP_COLORS.items()
                      if grp != "Unknown"]
    fig.legend(handles=legend_handles, loc="lower center", ncol=5,
               fontsize=8, facecolor="#222", labelcolor="#CCC", framealpha=0.8,
               bbox_to_anchor=(0.5, 0.01))

    out_path = FIG_DIR / "rice_fig21_outlier_analysis.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()
    print("=== Rice Investigation 21: Outlier Detection in Process Space ===\n")

    # 1. Load trajectory data from Inv.20
    print("Loading trajectory data...", flush=True)
    rows = load_individual_trajectory()

    # 2. Run outlier analysis
    results = run_outlier_analysis(rows)

    # 3. Figure
    print("\nGenerating figure...", flush=True)
    make_figure(
        results["outlier_records"],
        results["pop_centroid"],
        results["pop_std"],
    )

    # 4. Write TSV
    print("\nWriting outputs...", flush=True)
    fields = [
        "sampleId", "population", "group", "z_distance", "z_pc1", "z_pc2",
        "het_rate", "hom_alt_rate", "rare_burden", "private_burden",
        "froh_genome", "stat_pc1", "stat_pc2", "umap1", "umap2",
    ]
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(results["outlier_records"])

    # 5. Write JSON summary
    summary = {
        "n_samples":       len(results["outlier_records"]),
        "n_outliers_z3":   results["n_outliers_z3"],
        "n_outliers_z5":   results["n_outliers_z5"],
        "group_zdist":     results["group_zdist"],
        "pop_centroid":    results["pop_centroid"],
        "top_30_outliers": results["outlier_records"][:30],
    }
    with open(OUT_JSON, "w") as f:
        json.dump(summary, f, indent=2)

    elapsed = time.time() - t_start
    print(f"\nOutputs: {OUT_TSV}")
    print(f"         {OUT_JSON}")
    print(f"Total time: {elapsed:.1f}s")

    # Quick verification
    print(f"\nVerification:")
    print(f"  Total samples: {len(results['outlier_records'])}")
    print(f"  Outliers z>3:  {results['n_outliers_z3']}")
    print(f"  Outliers z>5:  {results['n_outliers_z5']}")
    if results["outlier_records"]:
        top = results["outlier_records"][0]
        print(f"  Top outlier:   {top['sampleId']} ({top['population']}) z={top['z_distance']:.2f}")


if __name__ == "__main__":
    main()
