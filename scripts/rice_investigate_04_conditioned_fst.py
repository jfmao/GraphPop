#!/usr/bin/env python3
"""Rice Inv.R04: Annotation-Conditioned Fst Genome Scan.

Missense vs synonymous Fst sliding windows reveal adaptive protein evolution.
This is a GraphPop-unique analysis: no classical tool can compute consequence-
conditioned sliding-window Fst in a single operation.

Data source: results/rice/rice_interpretation_results.json → gscan
Output:      data/results/rice/rice_inv04_conditioned_fst.{tsv,json}
Figure:      data/results/rice/figures/rice_fig04_conditioned_fst.png
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

CHROMOSOMES = [f"Chr{i}" for i in range(1, 13)]

# Key population pairs for conditioned genome scans
GSCAN_PAIRS = [
    ("GJ-tmp", "XI-1A"),
    ("GJ-tmp", "GJ-trp"),
    ("XI-1A", "cA-Aus"),
]


def load_gscan_windows(gscan, pop1, pop2, consequence):
    """Extract windows for a pop pair and consequence type."""
    windows = []
    for chrom in CHROMOSOMES:
        key = f"{pop1}|{pop2}|{chrom}|{consequence}"
        if key not in gscan:
            key = f"{pop2}|{pop1}|{chrom}|{consequence}"
        if key not in gscan:
            continue
        data = gscan[key]
        if isinstance(data, dict):
            if "result" in data:
                for w in data["result"]:
                    windows.append(w)
            elif "windows" in data:
                for w in data["windows"]:
                    w.setdefault("chr", chrom)
                    windows.append(w)
        elif isinstance(data, list):
            for w in data:
                w.setdefault("chr", chrom)
                windows.append(w)
    return windows


def main():
    with open(ROOT / "results" / "rice" / "rice_interpretation_results.json") as f:
        interp = json.load(f)

    gscan = interp["gscan"]

    all_rows = []
    pair_summaries = {}

    for pop1, pop2 in GSCAN_PAIRS:
        pair_label = f"{pop1} vs {pop2}"

        mis_windows = load_gscan_windows(gscan, pop1, pop2, "missense_variant")
        syn_windows = load_gscan_windows(gscan, pop1, pop2, "synonymous_variant")

        if not mis_windows or not syn_windows:
            print(f"WARNING: No conditioned gscan data for {pair_label}")
            continue

        # Windows have different start positions per consequence type.
        # Compare per-chromosome means and individual window distributions.
        for chrom in CHROMOSOMES:
            mis_chr = [w for w in mis_windows if w.get("chr", "") == chrom]
            syn_chr = [w for w in syn_windows if w.get("chr", "") == chrom]
            if not mis_chr or not syn_chr:
                continue
            for w in mis_chr:
                fst_val = w.get("fst_wc", w.get("fst", float("nan")))
                all_rows.append({
                    "pair": pair_label, "chr": chrom,
                    "start": w.get("start", 0),
                    "consequence": "missense", "fst": fst_val,
                })
            for w in syn_chr:
                fst_val = w.get("fst_wc", w.get("fst", float("nan")))
                all_rows.append({
                    "pair": pair_label, "chr": chrom,
                    "start": w.get("start", 0),
                    "consequence": "synonymous", "fst": fst_val,
                })

        mis_fst = [w.get("fst_wc", w.get("fst", 0)) for w in mis_windows]
        syn_fst = [w.get("fst_wc", w.get("fst", 0)) for w in syn_windows]
        mis_fst = [v for v in mis_fst if not np.isnan(v)]
        syn_fst = [v for v in syn_fst if not np.isnan(v)]

        if not mis_fst or not syn_fst:
            continue

        mean_mis = np.mean(mis_fst)
        mean_syn = np.mean(syn_fst)
        ratio = mean_mis / mean_syn if mean_syn > 0 else float("nan")

        pair_summaries[pair_label] = {
            "n_windows_missense": len(mis_fst),
            "n_windows_synonymous": len(syn_fst),
            "mean_fst_missense": float(mean_mis),
            "mean_fst_synonymous": float(mean_syn),
            "ratio": float(ratio),
        }

        print(f"{pair_label}: mis={len(mis_fst)} win (Fst={mean_mis:.4f}), "
              f"syn={len(syn_fst)} win (Fst={mean_syn:.4f}), ratio={ratio:.3f}")

    # Save TSV
    if all_rows:
        with open(RESULTS_DIR / "rice_inv04_conditioned_fst.tsv", "w") as f:
            header = ["pair", "chr", "start", "consequence", "fst"]
            f.write("\t".join(header) + "\n")
            for r in all_rows:
                f.write("\t".join(str(r[h]) for h in header) + "\n")

    with open(RESULTS_DIR / "rice_inv04_conditioned_fst.json", "w") as f:
        json.dump({"pairs": pair_summaries, "n_total_rows": len(all_rows)}, f, indent=2)

    if not pair_summaries:
        print("\nNo conditioned Fst data found. Checking gscan key format...")
        sample_keys = list(gscan.keys())[:10]
        for k in sample_keys:
            print(f"  Key: {k}")
        return

    # ── Figure: one row per pair, showing missense vs synonymous distributions ──
    n_pairs = len(pair_summaries)
    fig, axes = plt.subplots(1, n_pairs, figsize=(6 * n_pairs, 5))
    if n_pairs == 1:
        axes = [axes]

    for idx, (pair_label, ps) in enumerate(pair_summaries.items()):
        ax = axes[idx]
        pair_rows = [r for r in all_rows if r["pair"] == pair_label]

        mis_vals = [r["fst"] for r in pair_rows if r["consequence"] == "missense"]
        syn_vals = [r["fst"] for r in pair_rows if r["consequence"] == "synonymous"]

        bins = np.linspace(0, 1.0, 50)
        ax.hist(mis_vals, bins=bins, alpha=0.6, label=f"Missense (mean={ps['mean_fst_missense']:.3f})",
                color="#e41a1c", density=True)
        ax.hist(syn_vals, bins=bins, alpha=0.6, label=f"Synonymous (mean={ps['mean_fst_synonymous']:.3f})",
                color="#377eb8", density=True)

        ax.axvline(ps["mean_fst_missense"], color="#e41a1c", linestyle="--", linewidth=1.5)
        ax.axvline(ps["mean_fst_synonymous"], color="#377eb8", linestyle="--", linewidth=1.5)

        ax.set_xlabel("Fst (W&C)")
        ax.set_ylabel("Density")
        ax.set_title(f"{pair_label}\nratio={ps['ratio']:.3f}",
                     fontweight="bold", fontsize=9)
        ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig04_conditioned_fst.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"\n=== Rice Inv.R04: Conditioned Fst Genome Scan ===")
    for pair, ps in pair_summaries.items():
        print(f"  {pair}: ratio={ps['ratio']:.3f} "
              f"(mis={ps['n_windows_missense']}, syn={ps['n_windows_synonymous']})")
    print(f"Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
