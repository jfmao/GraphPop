#!/usr/bin/env python3
"""Rice Inv.R06: ROH Distribution + Inbreeding Landscape.

FROH analysis across 12 rice subpopulations reveals the selfing mating
system and bottleneck signatures.

Data source: results/rice/rice_interpretation_results.json → roh_hmm_summary
             results/rice/rice_deep_integration_results.json → roh_distribution, fingerprints
Output:      data/results/rice/rice_inv06_roh_landscape.{tsv,json}
Figure:      data/results/rice/figures/rice_fig06_roh_landscape.png
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

CHROMOSOMES = [f"Chr{i}" for i in range(1, 13)]


def main():
    with open(ROOT / "results" / "rice" / "rice_interpretation_results.json") as f:
        interp = json.load(f)
    with open(ROOT / "results" / "rice" / "rice_deep_integration_results.json") as f:
        di = json.load(f)

    roh_summary = interp["roh_hmm_summary"]
    roh_dist = di["roh_distribution"]
    profiles = di["fingerprints"]["profiles"]

    pops = sorted(roh_summary.keys())

    # Sort by FROH descending
    pop_sorted = sorted(pops, key=lambda p: profiles.get(p, {}).get("froh", 0), reverse=True)

    # Save TSV
    with open(RESULTS_DIR / "rice_inv06_roh_landscape.tsv", "w") as f:
        f.write("population\tgroup\tfroh\tmean_n_roh\tmean_total_kb\tpi\n")
        for pop in pop_sorted:
            froh = profiles.get(pop, {}).get("froh", 0)
            pi = profiles.get(pop, {}).get("pi", 0)
            rs = roh_summary.get(pop, {})
            n_roh = rs.get("mean_n_roh", rs.get("n_roh_mean", 0))
            total_kb = rs.get("mean_total_kb", rs.get("total_kb_mean", 0))
            f.write(f"{pop}\t{POP_GROUPS.get(pop, '?')}\t{froh:.6f}\t"
                    f"{n_roh}\t{total_kb}\t{pi:.6f}\n")

    # FROH vs pi correlation
    froh_vals = [profiles[p]["froh"] for p in pops]
    pi_vals = [profiles[p]["pi"] for p in pops]
    rho, pval = stats.spearmanr(froh_vals, pi_vals)

    summary = {
        "highest_froh": {"pop": pop_sorted[0], "froh": profiles[pop_sorted[0]]["froh"]},
        "lowest_froh": {"pop": pop_sorted[-1], "froh": profiles[pop_sorted[-1]]["froh"]},
        "froh_vs_pi_rho": float(rho),
        "froh_vs_pi_pval": float(pval),
        "indica_mean_froh": float(np.mean([profiles[p]["froh"] for p in pops
                                            if POP_GROUPS.get(p) == "Indica"])),
        "japonica_mean_froh": float(np.mean([profiles[p]["froh"] for p in pops
                                              if POP_GROUPS.get(p) == "Japonica"])),
    }
    with open(RESULTS_DIR / "rice_inv06_roh_landscape.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 3 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # Panel (a): FROH bar chart
    ax = axes[0]
    y_pos = range(len(pop_sorted))
    colors = [GROUP_COLORS.get(POP_GROUPS.get(p, "?"), "gray") for p in pop_sorted]
    froh_list = [profiles[p]["froh"] for p in pop_sorted]
    bars = ax.barh(y_pos, froh_list, color=colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pop_sorted, fontsize=8)
    ax.set_xlabel("FROH")
    ax.set_title("a  Inbreeding coefficient (FROH)", fontweight="bold", loc="left")
    ax.invert_yaxis()
    for bar, v in zip(bars, froh_list):
        ax.text(v + 0.001, bar.get_y() + bar.get_height() / 2,
                f"{v:.4f}", va="center", fontsize=7)

    # Panel (b): Per-chromosome ROH profile
    ax = axes[1]
    per_chr = roh_dist.get("per_chromosome", {})
    if per_chr:
        chr_pops = sorted(per_chr.keys())
        chr_mat = np.zeros((len(CHROMOSOMES), len(chr_pops)))
        for j, pop in enumerate(chr_pops):
            pop_chr = per_chr[pop]
            for i, chrom in enumerate(CHROMOSOMES):
                entry = pop_chr.get(chrom, {})
                if isinstance(entry, dict):
                    chr_mat[i, j] = entry.get("froh", entry.get("total_bp", 0))
                else:
                    chr_mat[i, j] = float(entry)

        im = ax.imshow(chr_mat.T, aspect="auto", cmap="YlOrRd")
        ax.set_xticks(range(len(CHROMOSOMES)))
        ax.set_xticklabels([f"C{i+1}" for i in range(12)], fontsize=7)
        ax.set_yticks(range(len(chr_pops)))
        ax.set_yticklabels(chr_pops, fontsize=7)
        plt.colorbar(im, ax=ax, label="Total ROH (bp)", shrink=0.8)
    else:
        ax.text(0.5, 0.5, "No per-chr ROH data", ha="center", va="center",
                transform=ax.transAxes)
    ax.set_title("b  ROH per chromosome", fontweight="bold", loc="left")

    # Panel (c): FROH vs pi scatter
    ax = axes[2]
    for pop in pops:
        grp = POP_GROUPS.get(pop, "Unknown")
        color = GROUP_COLORS.get(grp, "gray")
        ax.scatter(profiles[pop]["pi"], profiles[pop]["froh"],
                   c=color, s=80, edgecolors="black", linewidth=0.5, zorder=3)
        ax.annotate(pop, (profiles[pop]["pi"], profiles[pop]["froh"]),
                    fontsize=6, xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel("pi (nucleotide diversity)")
    ax.set_ylabel("FROH")
    ax.set_title(f"c  Diversity vs inbreeding (rho={rho:.2f})",
                 fontweight="bold", loc="left", fontsize=9)

    for grp, color in GROUP_COLORS.items():
        ax.scatter([], [], c=color, s=60, label=grp, edgecolors="black", linewidth=0.5)
    ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig06_roh_landscape.png", dpi=200, bbox_inches="tight")
    plt.close()

    print("=== Rice Inv.R06: ROH Distribution + Inbreeding Landscape ===")
    print(f"Highest FROH: {pop_sorted[0]} ({froh_list[0]:.4f})")
    print(f"Lowest FROH:  {pop_sorted[-1]} ({froh_list[-1]:.4f})")
    print(f"FROH vs pi: rho={rho:.3f}, p={pval:.4f}")
    print(f"Indica mean FROH:   {summary['indica_mean_froh']:.4f}")
    print(f"Japonica mean FROH: {summary['japonica_mean_froh']:.4f}")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
