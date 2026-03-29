#!/usr/bin/env python3
"""Rice Inv.R02: piN/piS Deep Analysis.

Annotation-conditioned diversity reveals the cost of domestication:
piN/piS > 1.0 in all 12 subpopulations.

Data source: results/rice/rice_interpretation_results.json → pinsps, pinsps_ratios
             results/rice/rice_deep_integration_results.json → fingerprints (for FROH)
Output:      data/results/rice/rice_inv02_pinsps.{tsv,json}
Figure:      data/results/rice/figures/rice_fig02_pinsps.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

CHROMOSOMES = [f"Chr{i}" for i in range(1, 13)]

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
    with open(ROOT / "results" / "rice" / "rice_interpretation_results.json") as f:
        interp = json.load(f)
    with open(ROOT / "results" / "rice" / "rice_deep_integration_results.json") as f:
        di = json.load(f)

    pinsps_ratios = interp["pinsps_ratios"]
    pinsps_raw = interp["pinsps"]
    profiles = di["fingerprints"]["profiles"]

    pops = sorted(pinsps_ratios.keys())

    # Build per-population per-chromosome piN/piS matrix
    chr_matrix = {}  # pop -> [piN/piS per chr]
    for pop in pops:
        vals = []
        per_chr = pinsps_ratios[pop].get("per_chr", {})
        for chrom in CHROMOSOMES:
            if chrom in per_chr:
                vals.append(per_chr[chrom].get("pi_n_pi_s", float("nan")))
            else:
                vals.append(float("nan"))
        chr_matrix[pop] = vals

    def get_pinsps(pop):
        return pinsps_ratios[pop].get("mean_pi_n_pi_s", float("nan"))

    # Sort by genome-wide piN/piS (descending)
    pop_sorted = sorted(pops, key=lambda p: get_pinsps(p), reverse=True)

    # Save TSV
    header = ["population", "group", "pinsps_genome_wide", "froh"] + CHROMOSOMES
    with open(RESULTS_DIR / "rice_inv02_pinsps.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for pop in pop_sorted:
            froh = profiles.get(pop, {}).get("froh", float("nan"))
            genome_wide = get_pinsps(pop)
            vals = [pop, POP_GROUPS.get(pop, "?"), f"{genome_wide:.4f}", f"{froh:.4f}"]
            vals += [f"{v:.4f}" if not np.isnan(v) else "NA" for v in chr_matrix[pop]]
            f.write("\t".join(vals) + "\n")

    # Correlation: FROH vs piN/piS
    froh_vals = [profiles[p]["froh"] for p in pops]
    pinsps_vals = [get_pinsps(p) for p in pops]
    rho, pval = stats.spearmanr(froh_vals, pinsps_vals)

    # Indica vs Japonica means
    indica_mean = np.mean([get_pinsps(p) for p in pops if POP_GROUPS.get(p) == "Indica"])
    japonica_mean = np.mean([get_pinsps(p) for p in pops if POP_GROUPS.get(p) == "Japonica"])

    summary = {
        "all_above_1": all(get_pinsps(p) > 1.0 for p in pops),
        "highest": {"pop": pop_sorted[0], "pinsps": get_pinsps(pop_sorted[0])},
        "lowest": {"pop": pop_sorted[-1], "pinsps": get_pinsps(pop_sorted[-1])},
        "indica_mean": indica_mean,
        "japonica_mean": japonica_mean,
        "froh_vs_pinsps_rho": rho,
        "froh_vs_pinsps_pval": pval,
        "per_chromosome_range": {
            pop: {"min": float(np.nanmin(chr_matrix[pop])), "max": float(np.nanmax(chr_matrix[pop]))}
            for pop in pops
        },
    }
    with open(RESULTS_DIR / "rice_inv02_pinsps.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 3 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # Panel (a): piN/piS bar chart per population (sorted)
    ax = axes[0]
    y_pos = range(len(pop_sorted))
    colors = [GROUP_COLORS.get(POP_GROUPS.get(p, "?"), "gray") for p in pop_sorted]
    values = [get_pinsps(p) for p in pop_sorted]
    bars = ax.barh(y_pos, values, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pop_sorted, fontsize=8)
    ax.axvline(1.0, color="black", linestyle="--", linewidth=1, label="Neutral (piN/piS=1)")
    ax.set_xlabel("piN/piS (genome-wide)")
    ax.set_title("a  piN/piS per subpopulation", fontweight="bold", loc="left")
    ax.invert_yaxis()

    for bar, v in zip(bars, values):
        ax.text(v + 0.002, bar.get_y() + bar.get_height() / 2,
                f"{v:.3f}", va="center", fontsize=7)

    # Panel (b): Per-chromosome heatmap
    ax = axes[1]
    mat = np.array([chr_matrix[pop] for pop in pop_sorted])
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0.9, vmax=1.25)
    ax.set_xticks(range(12))
    ax.set_xticklabels([f"Chr{i+1}" for i in range(12)], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(pop_sorted)))
    ax.set_yticklabels(pop_sorted, fontsize=7)
    plt.colorbar(im, ax=ax, label="piN/piS", shrink=0.8)
    ax.set_title("b  piN/piS per chromosome", fontweight="bold", loc="left")

    # Panel (c): FROH vs piN/piS scatter
    ax = axes[2]
    for pop in pops:
        grp = POP_GROUPS.get(pop, "Unknown")
        color = GROUP_COLORS.get(grp, "gray")
        ax.scatter(profiles[pop]["froh"], get_pinsps(pop),
                   c=color, s=80, edgecolors="black", linewidth=0.5, zorder=3)
        ax.annotate(pop, (profiles[pop]["froh"], get_pinsps(pop)),
                    fontsize=6, ha="left", va="bottom",
                    xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel("FROH (inbreeding)")
    ax.set_ylabel("piN/piS")
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_title(f"c  FROH vs piN/piS (rho={rho:.2f}, p={pval:.3f})",
                 fontweight="bold", loc="left", fontsize=9)

    for grp, color in GROUP_COLORS.items():
        ax.scatter([], [], c=color, s=60, label=grp, edgecolors="black", linewidth=0.5)
    ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig02_pinsps.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Print summary
    print("=== Rice Inv.R02: piN/piS Deep Analysis ===")
    print(f"All populations piN/piS > 1.0: {summary['all_above_1']}")
    print(f"Highest: {pop_sorted[0]} ({values[0]:.4f})")
    print(f"Lowest:  {pop_sorted[-1]} ({values[-1]:.4f})")
    print(f"Indica mean:   {indica_mean:.4f}")
    print(f"Japonica mean: {japonica_mean:.4f}")
    print(f"FROH vs piN/piS: rho={rho:.3f}, p={pval:.4f}")
    print(f"\nSaved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
