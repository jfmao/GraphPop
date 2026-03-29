#!/usr/bin/env python3
"""Rice Inv.R09: Multi-Statistic Correlations.

Correlation heatmap and key scatter plots from 8 evolutionary statistics
across 12 rice subpopulations.

Data source: results/rice/rice_deep_integration_results.json → correlations
Output:      data/results/rice/rice_inv09_stat_correlations.{tsv,json}
Figure:      data/results/rice/figures/rice_fig09_stat_correlations.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

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

# Friendly stat labels
STAT_LABELS = {
    "pi": "pi", "theta_w": "theta_W", "tajima_d": "Tajima's D",
    "fay_wu_h": "Fay&Wu H", "fis": "FIS", "froh": "FROH",
    "n_sweeps": "N sweeps", "pinsps": "piN/piS",
}


def main():
    with open(ROOT / "results" / "rice" / "rice_deep_integration_results.json") as f:
        di = json.load(f)

    corr = di["correlations"]
    stat_names = corr["stat_names"]
    corr_matrix_dict = corr["correlation_matrix"]
    pop_stats = corr["population_stats"]
    notable = corr["notable_correlations"]
    expected = corr["expected_patterns"]

    pops = sorted(pop_stats.keys())
    n_stats = len(stat_names)

    # Build numeric correlation matrix
    corr_matrix = np.zeros((n_stats, n_stats))
    for i, s1 in enumerate(stat_names):
        for j, s2 in enumerate(stat_names):
            corr_matrix[i, j] = corr_matrix_dict[s1][s2]

    # Save TSV: population × stats table
    header = ["population", "group"] + stat_names
    rows = []
    for pop in pops:
        row = {"population": pop, "group": POP_GROUPS.get(pop, "Unknown")}
        for stat in stat_names:
            row[stat] = pop_stats[pop].get(stat, float("nan"))
        rows.append(row)

    with open(RESULTS_DIR / "rice_inv09_stat_correlations.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")

    # Save JSON summary
    summary = {
        "n_populations": len(pops),
        "n_statistics": n_stats,
        "stat_names": stat_names,
        "notable_correlations": notable,
        "expected_patterns": expected,
        "strongest_positive": max(notable, key=lambda x: x["spearman_r"]),
        "strongest_negative": min(notable, key=lambda x: x["spearman_r"]),
    }
    with open(RESULTS_DIR / "rice_inv09_stat_correlations.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 2 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a): Correlation heatmap
    ax = axes[0]
    labels = [STAT_LABELS.get(s, s) for s in stat_names]
    im = ax.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
    ax.set_xticks(range(n_stats))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n_stats))
    ax.set_yticklabels(labels, fontsize=7)
    # Annotate cells
    for i in range(n_stats):
        for j in range(n_stats):
            val = corr_matrix[i, j]
            color = "white" if abs(val) > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                    fontsize=6, color=color)
    plt.colorbar(im, ax=ax, label="Spearman r", shrink=0.8)
    ax.set_title("a  Stat-stat correlation (Spearman)", fontweight="bold", loc="left")

    # Panel (b): Key scatter plots — pick top 2 notable correlations
    ax = axes[1]
    # Plot pi vs n_sweeps (strongest negative) and pi vs theta_w (strongest positive)
    top_pos = summary["strongest_positive"]
    top_neg = summary["strongest_negative"]

    # Use pi vs n_sweeps as the main scatter
    stat_x = top_neg["stat1"]
    stat_y = top_neg["stat2"]
    x_vals = [pop_stats[p].get(stat_x, float("nan")) for p in pops]
    y_vals = [pop_stats[p].get(stat_y, float("nan")) for p in pops]

    for k, pop in enumerate(pops):
        grp = POP_GROUPS.get(pop, "Unknown")
        color = GROUP_COLORS.get(grp, "gray")
        ax.scatter(x_vals[k], y_vals[k], c=color, s=80, edgecolors="black",
                   linewidth=0.5, zorder=3)
        ax.annotate(pop, (x_vals[k], y_vals[k]), fontsize=6, ha="center",
                    va="bottom", xytext=(0, 5), textcoords="offset points")

    ax.set_xlabel(STAT_LABELS.get(stat_x, stat_x))
    ax.set_ylabel(STAT_LABELS.get(stat_y, stat_y))
    r_val = top_neg["spearman_r"]
    ax.set_title(f"b  {STAT_LABELS.get(stat_x, stat_x)} vs "
                 f"{STAT_LABELS.get(stat_y, stat_y)} (r={r_val:.2f})",
                 fontweight="bold", loc="left")
    # Legend
    for grp, color in GROUP_COLORS.items():
        ax.scatter([], [], c=color, s=50, label=grp, edgecolors="black", linewidth=0.5)
    ax.legend(fontsize=7, loc="best")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig09_stat_correlations.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Print summary
    print("=== Rice Inv.R09: Multi-Statistic Correlations ===")
    print(f"Populations: {len(pops)}")
    print(f"Statistics: {n_stats} ({', '.join(stat_names)})")
    print(f"\nNotable correlations (|r| > 0.5):")
    for nc in notable:
        print(f"  {nc['stat1']:10s} × {nc['stat2']:10s}: "
              f"r={nc['spearman_r']:+.3f} ({nc['direction']})")
    print(f"\nExpected patterns:")
    for ep in expected:
        obs = ep["observed_r"]
        print(f"  {ep['stat1']:10s} × {ep['stat2']:10s}: "
              f"expected={ep['expected']}, observed r={obs:+.3f}")
        print(f"    Reason: {ep['reason']}")
    print(f"\nSaved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
