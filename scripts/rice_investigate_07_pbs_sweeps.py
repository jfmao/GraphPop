#!/usr/bin/env python3
"""Rice Inv.R07: PBS Population-Specific Sweeps.

Top PBS windows per population trio, annotated with genes from
xpehh_annotations and sweep_annotations.

Data source: results/rice/rice_interpretation_results.json → pbs_summary, annotate
Output:      data/results/rice/rice_inv07_pbs_sweeps.{tsv,json}
Figure:      data/results/rice/figures/rice_fig07_pbs_sweeps.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

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
    with open(ROOT / "results" / "rice" / "rice_interpretation_results.json") as f:
        interp = json.load(f)

    pbs_summary = interp["pbs_summary"]
    annotations = interp.get("annotate", {})
    sweep_ann = annotations.get("sweep_annotations", {})
    xpehh_ann = annotations.get("xpehh_annotations", {})

    # Collect all PBS comparison results
    rows = []
    comparison_stats = []
    for comp_id, comp in pbs_summary.items():
        focal = comp["focal"]
        sister = comp["sister"]
        outgroup = comp["outgroup"]
        n_windows = comp["n_windows"]
        mean_pbs = comp["mean_pbs"]
        top_windows = comp["top_windows"]

        comparison_stats.append({
            "comparison": comp_id,
            "focal": focal,
            "sister": sister,
            "outgroup": outgroup,
            "n_windows": n_windows,
            "mean_pbs": mean_pbs,
            "max_pbs": max(w["pbs"] for w in top_windows) if top_windows else 0,
            "n_top_windows": len(top_windows),
        })

        for w in top_windows:
            gene_symbols = list({g["symbol"] for g in w.get("genes", [])})
            known = w.get("known_genes", [])
            rows.append({
                "comparison": comp_id,
                "focal": focal,
                "sister": sister,
                "outgroup": outgroup,
                "chr": w["chr"],
                "start": w["start"],
                "end": w["end"],
                "pbs": w["pbs"],
                "fst_wc": w.get("fst_wc", float("nan")),
                "n_variants": w.get("n_variants", 0),
                "genes": ";".join(gene_symbols),
                "n_genes": len(gene_symbols),
                "known_genes": ";".join(
                    kg["gene"] if isinstance(kg, dict) else str(kg) for kg in known
                ) if known else "",
            })

    # Save TSV
    header = ["comparison", "focal", "sister", "outgroup", "chr", "start", "end",
              "pbs", "fst_wc", "n_variants", "genes", "n_genes", "known_genes"]
    with open(RESULTS_DIR / "rice_inv07_pbs_sweeps.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")

    # Build summary JSON
    # Count unique genes across all top PBS windows
    all_genes = set()
    for r in rows:
        if r["genes"]:
            all_genes.update(r["genes"].split(";"))

    # Top PBS window overall
    top_row = max(rows, key=lambda r: r["pbs"]) if rows else {}

    summary = {
        "n_comparisons": len(pbs_summary),
        "comparisons": comparison_stats,
        "total_top_windows": len(rows),
        "unique_genes_in_top_windows": len(all_genes),
        "overall_top_window": {
            "comparison": top_row.get("comparison", ""),
            "chr": top_row.get("chr", ""),
            "start": top_row.get("start", 0),
            "end": top_row.get("end", 0),
            "pbs": top_row.get("pbs", 0),
            "genes": top_row.get("genes", ""),
        },
    }
    with open(RESULTS_DIR / "rice_inv07_pbs_sweeps.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 3 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel (a): Mean PBS + Max PBS per comparison (bar chart)
    ax = axes[0]
    comp_names = [cs["focal"] for cs in comparison_stats]
    mean_vals = [cs["mean_pbs"] for cs in comparison_stats]
    max_vals = [cs["max_pbs"] for cs in comparison_stats]
    x = np.arange(len(comp_names))
    width = 0.35
    colors_mean = [GROUP_COLORS.get(POP_GROUPS.get(c, ""), "gray") for c in comp_names]
    ax.bar(x - width / 2, mean_vals, width, label="Mean PBS", color=colors_mean, alpha=0.6,
           edgecolor="black", linewidth=0.5)
    ax.bar(x + width / 2, max_vals, width, label="Max PBS", color=colors_mean,
           edgecolor="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(comp_names, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("PBS")
    ax.legend(fontsize=7)
    ax.set_title("a  PBS per focal population", fontweight="bold", loc="left")

    # Panel (b): Top 10 PBS windows across all comparisons
    ax = axes[1]
    sorted_rows = sorted(rows, key=lambda r: r["pbs"], reverse=True)[:15]
    labels = []
    pbs_vals = []
    bar_colors = []
    for r in reversed(sorted_rows):
        gene_str = r["genes"][:25] + "..." if len(r["genes"]) > 25 else r["genes"]
        lbl = f"{r['focal']} {r['chr']}:{r['start']//1_000_000}M"
        if gene_str:
            lbl += f" ({gene_str})"
        labels.append(lbl)
        pbs_vals.append(r["pbs"])
        bar_colors.append(GROUP_COLORS.get(POP_GROUPS.get(r["focal"], ""), "gray"))
    ax.barh(range(len(labels)), pbs_vals, color=bar_colors, edgecolor="black", linewidth=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("PBS")
    ax.set_title("b  Top 15 PBS windows", fontweight="bold", loc="left")

    # Panel (c): PBS vs FST scatter
    ax = axes[2]
    for r in rows:
        focal = r["focal"]
        grp = POP_GROUPS.get(focal, "Unknown")
        color = GROUP_COLORS.get(grp, "gray")
        fst_val = r["fst_wc"]
        if np.isfinite(fst_val):
            ax.scatter(fst_val, r["pbs"], c=color, s=20, alpha=0.6, edgecolors="none")
    ax.set_xlabel("FST (Weir-Cockerham)")
    ax.set_ylabel("PBS")
    ax.set_title("c  PBS vs FST at top windows", fontweight="bold", loc="left")
    # Legend
    for grp, color in GROUP_COLORS.items():
        ax.scatter([], [], c=color, s=40, label=grp)
    ax.legend(fontsize=7, loc="upper left")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig07_pbs_sweeps.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Print summary
    print("=== Rice Inv.R07: PBS Population-Specific Sweeps ===")
    print(f"Comparisons: {len(pbs_summary)}")
    print(f"Total top windows: {len(rows)}")
    print(f"Unique genes in top windows: {len(all_genes)}")
    print(f"\nPer-comparison summary:")
    for cs in comparison_stats:
        print(f"  {cs['focal']:8s} vs {cs['sister']:8s} (out={cs['outgroup']:6s}): "
              f"mean={cs['mean_pbs']:.3f}, max={cs['max_pbs']:.3f}, n_win={cs['n_windows']}")
    print(f"\nTop PBS window: {summary['overall_top_window']}")
    print(f"\nSaved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
