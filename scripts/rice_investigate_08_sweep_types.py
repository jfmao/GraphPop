#!/usr/bin/env python3
"""Rice Inv.R08: Hard vs Soft Sweep Classification.

Classify sweep signals at 16 known domestication genes by H12 vs H2/H1
ratio. Hard sweeps have high H12 and low H2/H1; soft sweeps have
moderate H12 and high H2/H1.

Data source: results/rice/rice_deep_integration_results.json → sweeps
Output:      data/results/rice/rice_inv08_sweep_types.{tsv,json}
Figure:      data/results/rice/figures/rice_fig08_sweep_types.png
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

SWEEP_COLORS = {"hard": "#e41a1c", "soft": "#ff7f00", "none": "#999999"}


def main():
    with open(ROOT / "results" / "rice" / "rice_deep_integration_results.json") as f:
        di = json.load(f)

    sweeps = di["sweeps"]
    genes = sorted(sweeps.keys())

    # Build rows with H12, H2/H1, sweep_type per gene × pop
    rows = []
    for gene in genes:
        gd = sweeps[gene]
        for pop, pd in gd["populations"].items():
            h12 = pd.get("h12", 0)
            h2_h1 = pd.get("h2_h1", float("nan"))
            sweep_type = pd.get("sweep_type", "none")
            rows.append({
                "gene": gene,
                "function": gd.get("function", ""),
                "chr": gd.get("chr", ""),
                "population": pop,
                "group": POP_GROUPS.get(pop, "Unknown"),
                "h12": h12,
                "h2_h1": h2_h1,
                "hap_diversity": pd.get("hap_diversity", float("nan")),
                "n_haplotypes": pd.get("n_haplotypes", 0),
                "sweep_type": sweep_type,
            })

    # Save TSV
    header = ["gene", "function", "chr", "population", "group",
              "h12", "h2_h1", "hap_diversity", "n_haplotypes", "sweep_type"]
    with open(RESULTS_DIR / "rice_inv08_sweep_types.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")

    # Count by type and group
    type_counts = {"hard": 0, "soft": 0, "none": 0}
    group_type_counts = {}
    for r in rows:
        st = r["sweep_type"]
        type_counts[st] = type_counts.get(st, 0) + 1
        grp = r["group"]
        if grp not in group_type_counts:
            group_type_counts[grp] = {"hard": 0, "soft": 0, "none": 0}
        group_type_counts[grp][st] = group_type_counts[grp].get(st, 0) + 1

    # Per-gene summary
    gene_summary = {}
    for gene in genes:
        gene_rows = [r for r in rows if r["gene"] == gene]
        n_hard = sum(1 for r in gene_rows if r["sweep_type"] == "hard")
        n_soft = sum(1 for r in gene_rows if r["sweep_type"] == "soft")
        mean_h12 = np.nanmean([r["h12"] for r in gene_rows])
        mean_h2h1 = np.nanmean([r["h2_h1"] for r in gene_rows])
        gene_summary[gene] = {
            "function": sweeps[gene].get("function", ""),
            "n_hard": n_hard,
            "n_soft": n_soft,
            "n_none": len(gene_rows) - n_hard - n_soft,
            "mean_h12": float(mean_h12),
            "mean_h2_h1": float(mean_h2h1),
        }

    summary = {
        "n_genes": len(genes),
        "n_observations": len(rows),
        "type_counts": type_counts,
        "group_type_counts": group_type_counts,
        "gene_summary": gene_summary,
    }
    with open(RESULTS_DIR / "rice_inv08_sweep_types.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 2 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a): H12 vs H2/H1 scatter with classification zones
    ax = axes[0]
    for r in rows:
        h12 = r["h12"]
        h2h1 = r["h2_h1"]
        if not np.isfinite(h2h1):
            continue
        color = SWEEP_COLORS.get(r["sweep_type"], "gray")
        marker = "o"
        ax.scatter(h12, h2h1, c=color, s=30, alpha=0.6, edgecolors="black",
                   linewidth=0.3, marker=marker, zorder=3)

    # Classification zone shading
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(0.05, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    # Label zones
    ax.text(0.3, 0.1, "Hard sweep\nzone", fontsize=8, color="#e41a1c",
            ha="center", va="center", alpha=0.7, style="italic")
    ax.text(0.3, 0.8, "Soft sweep\nzone", fontsize=8, color="#ff7f00",
            ha="center", va="center", alpha=0.7, style="italic")
    ax.text(0.01, 0.5, "Neutral\nzone", fontsize=8, color="#999999",
            ha="center", va="center", alpha=0.7, style="italic")

    ax.set_xlabel("H12")
    ax.set_ylabel("H2/H1")
    ax.set_title("a  H12 vs H2/H1 at domestication loci", fontweight="bold", loc="left")
    # Legend
    for st, color in SWEEP_COLORS.items():
        ax.scatter([], [], c=color, s=50, label=f"{st} (n={type_counts.get(st, 0)})",
                   edgecolors="black", linewidth=0.3)
    ax.legend(fontsize=8, loc="upper right")

    # Panel (b): Sweep type counts per varietal group (stacked bar)
    ax = axes[1]
    groups = sorted(group_type_counts.keys())
    hard_vals = [group_type_counts[g].get("hard", 0) for g in groups]
    soft_vals = [group_type_counts[g].get("soft", 0) for g in groups]
    none_vals = [group_type_counts[g].get("none", 0) for g in groups]
    x = np.arange(len(groups))

    ax.bar(x, hard_vals, label="Hard sweep", color="#e41a1c", edgecolor="black", linewidth=0.5)
    ax.bar(x, soft_vals, bottom=hard_vals, label="Soft sweep", color="#ff7f00",
           edgecolor="black", linewidth=0.5)
    bottom2 = [h + s for h, s in zip(hard_vals, soft_vals)]
    ax.bar(x, none_vals, bottom=bottom2, label="No sweep", color="#999999",
           edgecolor="black", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Gene x Population observations")
    ax.legend(fontsize=8)
    ax.set_title("b  Sweep classification by varietal group", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig08_sweep_types.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Print summary
    print("=== Rice Inv.R08: Hard vs Soft Sweep Classification ===")
    print(f"Genes: {len(genes)}, Observations: {len(rows)}")
    print(f"Sweep types: hard={type_counts['hard']}, soft={type_counts['soft']}, "
          f"none={type_counts['none']}")
    print(f"\nPer varietal group:")
    for grp in groups:
        gc = group_type_counts[grp]
        print(f"  {grp:15s}: hard={gc['hard']:3d}, soft={gc['soft']:3d}, none={gc['none']:3d}")
    print(f"\nPer gene:")
    for gene in genes:
        gs = gene_summary[gene]
        print(f"  {gene:8s} ({gs['function'][:25]:25s}): "
              f"hard={gs['n_hard']}, soft={gs['n_soft']}, "
              f"mean_H12={gs['mean_h12']:.4f}, mean_H2/H1={gs['mean_h2_h1']:.3f}")
    print(f"\nSaved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
