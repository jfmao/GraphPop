#!/usr/bin/env python3
"""Rice Inv.R12: ROH × Sweep Cross-Correlation.

Correlate inbreeding (FROH) with selection signals (n_sweeps, H12)
across rice subpopulations — analogous to human Inv.03.

Data source: results/rice/rice_deep_integration_results.json → fingerprints, roh_distribution
Output:      data/results/rice/rice_inv12_roh_sweep_correlation.{tsv,json}
Figure:      data/results/rice/figures/rice_fig12_roh_sweep_correlation.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

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


def main():
    with open(ROOT / "results" / "rice" / "rice_deep_integration_results.json") as f:
        di = json.load(f)

    profiles = di["fingerprints"]["profiles"]
    pops = sorted(profiles.keys())

    # Extract key variables per population
    data = []
    for pop in pops:
        p = profiles[pop]
        data.append({
            "population": pop,
            "group": POP_GROUPS.get(pop, "?"),
            "pi": p["pi"],
            "froh": p["froh"],
            "n_sweeps": p["n_sweeps"],
            "hard_sweep_fraction": p["hard_sweep_fraction"],
            "fis": p["fis"],
            "tajima_d": p["tajima_d"],
            "pinsps": p.get("pinsps", 0),
        })

    # Compute correlations
    pi = np.array([d["pi"] for d in data])
    froh = np.array([d["froh"] for d in data])
    n_sweeps = np.array([d["n_sweeps"] for d in data])
    hard_frac = np.array([d["hard_sweep_fraction"] for d in data])
    fis = np.array([d["fis"] for d in data])
    tajd = np.array([d["tajima_d"] for d in data])

    correlations = {}
    pairs = [
        ("pi", "n_sweeps", pi, n_sweeps),
        ("pi", "froh", pi, froh),
        ("froh", "n_sweeps", froh, n_sweeps),
        ("froh", "fis", froh, fis),
        ("tajima_d", "n_sweeps", tajd, n_sweeps),
        ("pi", "hard_sweep_fraction", pi, hard_frac),
    ]

    for name1, name2, x, y in pairs:
        rho, pval = stats.spearmanr(x, y)
        correlations[f"{name1}_vs_{name2}"] = {"rho": float(rho), "p": float(pval)}

    # Save TSV
    header = ["population", "group", "pi", "froh", "n_sweeps", "hard_sweep_fraction",
              "fis", "tajima_d", "pinsps"]
    with open(RESULTS_DIR / "rice_inv12_roh_sweep_correlation.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for d in data:
            f.write("\t".join(str(d[h]) for h in header) + "\n")

    summary = {
        "n_populations": len(pops),
        "correlations": correlations,
        "key_finding": "pi vs n_sweeps" if abs(correlations["pi_vs_n_sweeps"]["rho"]) > 0.5 else "weak",
    }
    with open(RESULTS_DIR / "rice_inv12_roh_sweep_correlation.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 4 panels ──
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    plot_pairs = [
        (axes[0, 0], pi, n_sweeps, "pi", "n_sweeps", "a"),
        (axes[0, 1], froh, n_sweeps, "FROH", "n_sweeps", "b"),
        (axes[1, 0], pi, froh, "pi", "FROH", "c"),
        (axes[1, 1], tajd, n_sweeps, "Tajima's D", "n_sweeps", "d"),
    ]

    for ax, x, y, xlabel, ylabel, panel in plot_pairs:
        for d in data:
            grp = d["group"]
            color = GROUP_COLORS.get(grp, "gray")
            xi = d[xlabel.lower().replace("'s ", "_").replace(" ", "_")]
            yi = d[ylabel.lower()]
            ax.scatter(xi, yi, c=color, s=80, edgecolors="black", linewidth=0.5, zorder=3)
            ax.annotate(d["population"], (xi, yi), fontsize=6, xytext=(3, 3),
                        textcoords="offset points")

        rho, pval = stats.spearmanr(x, y)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{panel}  {xlabel} vs {ylabel} (rho={rho:.2f}, p={pval:.3f})",
                     fontweight="bold", loc="left", fontsize=9)

    # Legend
    for grp, color in GROUP_COLORS.items():
        axes[0, 0].scatter([], [], c=color, s=60, label=grp, edgecolors="black", linewidth=0.5)
    axes[0, 0].legend(fontsize=7, loc="upper right")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig12_roh_sweep_correlation.png", dpi=200, bbox_inches="tight")
    plt.close()

    print("=== Rice Inv.R12: ROH × Sweep Cross-Correlation ===")
    for key, c in correlations.items():
        print(f"  {key}: rho={c['rho']:.3f}, p={c['p']:.4f}")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
