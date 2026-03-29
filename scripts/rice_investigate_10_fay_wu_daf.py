#!/usr/bin/env python3
"""Rice Inv.R10: Fay & Wu's H and DAF Enrichment.

Fay & Wu's H per population (genome-wide and per-chromosome) plus
derived allele frequency enrichment in PBS and XP-EHH peak regions.

Data source: results/rice/rice_interpretation_results.json → fay_wu_summary, daf_enrichment
Output:      data/results/rice/rice_inv10_fay_wu_daf.{tsv,json}
Figure:      data/results/rice/figures/rice_fig10_fay_wu_daf.png
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


def main():
    with open(ROOT / "results" / "rice" / "rice_interpretation_results.json") as f:
        interp = json.load(f)

    fw_summary = interp["fay_wu_summary"]
    daf = interp["daf_enrichment"]

    pops = sorted(fw_summary.keys())

    # Build TSV rows: per population Fay & Wu summary
    rows = []
    for pop in pops:
        fwd = fw_summary[pop]
        rows.append({
            "population": pop,
            "group": POP_GROUPS.get(pop, "Unknown"),
            "mean_h": fwd["mean_h"],
            "mean_h_norm": fwd["mean_h_norm"],
            "n_chr": len(fwd.get("per_chr", {})),
        })

    header = ["population", "group", "mean_h", "mean_h_norm", "n_chr"]
    with open(RESULTS_DIR / "rice_inv10_fay_wu_daf.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[h]) for h in header) + "\n")

    # DAF enrichment analysis
    background = daf.get("background_daf", {})
    pbs_peak = daf.get("pbs_peak_daf", {})
    xpehh_peak = daf.get("xpehh_peak_daf", {})

    daf_rows = []
    for comp_id, pd in pbs_peak.items():
        daf_rows.append({
            "test": "PBS",
            "comparison": comp_id,
            "focal": pd.get("focal", ""),
            "mean_peak_daf": pd.get("mean_peak_daf", float("nan")),
            "background_daf": pd.get("background_daf", float("nan")),
            "enrichment": pd.get("enrichment", float("nan")),
            "n_peaks": pd.get("n_peaks", 0),
        })
    for comp_id, xd in xpehh_peak.items():
        daf_rows.append({
            "test": "XP-EHH",
            "comparison": comp_id,
            "focal": xd.get("focal", ""),
            "mean_peak_daf": xd.get("mean_peak_daf", float("nan")),
            "background_daf": xd.get("background_daf", float("nan")),
            "enrichment": xd.get("enrichment", float("nan")),
            "n_peaks": xd.get("n_peaks", 0),
        })

    # Summary JSON
    # Sort pops by mean_h (most negative = strongest directional selection)
    sorted_by_h = sorted(rows, key=lambda r: r["mean_h"])
    mean_enrichment_pbs = np.nanmean([d["enrichment"] for d in daf_rows if d["test"] == "PBS"])
    mean_enrichment_xpehh = np.nanmean([d["enrichment"] for d in daf_rows if d["test"] == "XP-EHH"])

    summary = {
        "n_populations": len(pops),
        "most_negative_h": {"pop": sorted_by_h[0]["population"],
                            "mean_h": sorted_by_h[0]["mean_h"],
                            "mean_h_norm": sorted_by_h[0]["mean_h_norm"]},
        "least_negative_h": {"pop": sorted_by_h[-1]["population"],
                             "mean_h": sorted_by_h[-1]["mean_h"],
                             "mean_h_norm": sorted_by_h[-1]["mean_h_norm"]},
        "mean_daf_enrichment_pbs": float(mean_enrichment_pbs),
        "mean_daf_enrichment_xpehh": float(mean_enrichment_xpehh),
        "n_pbs_comparisons": sum(1 for d in daf_rows if d["test"] == "PBS"),
        "n_xpehh_comparisons": sum(1 for d in daf_rows if d["test"] == "XP-EHH"),
        "daf_enrichment_details": daf_rows,
    }
    with open(RESULTS_DIR / "rice_inv10_fay_wu_daf.json", "w") as f:
        json.dump(summary, f, indent=2)

    # ── Figure: 2 panels ──────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel (a): Fay & Wu's H per population (bar chart, sorted)
    ax = axes[0]
    sorted_pops = [r["population"] for r in sorted_by_h]
    h_vals = [r["mean_h"] for r in sorted_by_h]
    bar_colors = [GROUP_COLORS.get(POP_GROUPS.get(p, ""), "gray") for p in sorted_pops]

    bars = ax.bar(range(len(sorted_pops)), h_vals, color=bar_colors,
                  edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(sorted_pops)))
    ax.set_xticklabels(sorted_pops, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Mean Fay & Wu's H")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("a  Fay & Wu's H per population", fontweight="bold", loc="left")
    # Legend
    for grp, color in GROUP_COLORS.items():
        ax.bar([], [], color=color, label=grp, edgecolor="black", linewidth=0.5)
    ax.legend(fontsize=7, loc="lower right")

    # Panel (b): DAF enrichment — peak vs background
    ax = axes[1]
    # Separate PBS and XP-EHH
    pbs_rows = [d for d in daf_rows if d["test"] == "PBS"]
    xpehh_rows = [d for d in daf_rows if d["test"] == "XP-EHH"]

    all_daf = pbs_rows + xpehh_rows
    labels = []
    peak_daf = []
    bg_daf = []
    test_colors = []
    for d in all_daf:
        lbl = d["focal"]
        if d["test"] == "XP-EHH":
            lbl = f"{d['comparison']}"
        labels.append(f"{d['test']}: {lbl}")
        peak_daf.append(d["mean_peak_daf"])
        bg_daf.append(d["background_daf"])
        test_colors.append("#e41a1c" if d["test"] == "PBS" else "#377eb8")

    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, bg_daf, width, label="Background DAF", color="#cccccc",
           edgecolor="black", linewidth=0.5)
    ax.bar(x + width / 2, peak_daf, width, label="Peak DAF", color=test_colors,
           edgecolor="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=5)
    ax.set_ylabel("Mean DAF")
    ax.legend(fontsize=7)
    ax.set_title("b  DAF enrichment at selection peaks", fontweight="bold", loc="left")

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig10_fay_wu_daf.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Print summary
    print("=== Rice Inv.R10: Fay & Wu's H + DAF Enrichment ===")
    print(f"Populations: {len(pops)}")
    print(f"\nFay & Wu's H (sorted, most negative first):")
    for r in sorted_by_h:
        print(f"  {r['population']:10s} ({r['group']:12s}): "
              f"H={r['mean_h']:.6f}, H_norm={r['mean_h_norm']:.3f}")
    print(f"\nDAF enrichment:")
    print(f"  PBS mean enrichment:    {mean_enrichment_pbs:.3f}")
    print(f"  XP-EHH mean enrichment: {mean_enrichment_xpehh:.3f}")
    for d in daf_rows:
        print(f"  {d['test']:7s} {d['comparison']:35s}: "
              f"peak={d['mean_peak_daf']:.4f}, bg={d['background_daf']:.4f}, "
              f"enrich={d['enrichment']:.3f}")
    print(f"\nSaved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
