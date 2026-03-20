#!/usr/bin/env python3
"""Fix FROH values in human_interpretation_results.json.

The Java ROHProcedure previously computed FROH using the variant span
(last_variant_pos - first_variant_pos) instead of the true chromosome length.
For most chromosomes this is fine (<0.3% error), but for acrocentric chromosomes
with large pericentromeric gaps the span is substantially shorter than the true
length, inflating FROH by 12–26%:

  chr13: 16% inflation  (span ~98 Mb vs true 114 Mb)
  chr14: 18% inflation  (span ~91 Mb vs true 107 Mb)
  chr15: 20% inflation  (span ~85 Mb vs true 102 Mb)
  chr21: 12% inflation  (span ~42 Mb vs true 47 Mb)
  chr22: 26% inflation  (span ~40 Mb vs true 51 Mb)

This script rewrites the froh fields using the correct hg38 chromosome lengths.
Saves atomically via tmp+rename. Prints a before/after summary.

Usage:
  conda run -n graphevo python scripts/fix_froh.py [--results PATH]
"""

import argparse
import json
import os
import sys

CHR_LENGTHS_HG38 = {
    "chr1":  248956422, "chr2":  242193529, "chr3":  198295559,
    "chr4":  190214555, "chr5":  181538259, "chr6":  170805979,
    "chr7":  159345973, "chr8":  145138636, "chr9":  138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16":  90338345, "chr17":  83257441, "chr18":  80373285,
    "chr19":  58617616, "chr20":  64444167, "chr21":  46709983,
    "chr22":  50818468,
}


def fix_roh_hmm(roh_hmm: dict) -> dict:
    """Return a corrected copy of the roh_hmm dict."""
    fixed = {}
    for key, entry in roh_hmm.items():
        pop, chrom = key.split("|")
        chr_len = CHR_LENGTHS_HG38.get(chrom)
        if chr_len is None:
            fixed[key] = entry  # leave sex chromosomes untouched
            continue

        new_entry = dict(entry)
        new_per_sample = []
        for s in entry["per_sample"]:
            ns = dict(s)
            ns["froh"] = s["total_length"] / chr_len
            new_per_sample.append(ns)

        new_entry["per_sample"] = new_per_sample
        froh_vals = [s["froh"] for s in new_per_sample]
        new_entry["mean_froh"] = sum(froh_vals) / len(froh_vals) if froh_vals else 0.0
        new_entry["max_froh"]  = max(froh_vals) if froh_vals else 0.0
        fixed[key] = new_entry
    return fixed


def fix_roh_hmm_summary(roh_hmm_summary: dict, fixed_roh_hmm: dict) -> dict:
    """Recompute mean_froh_across_chr for each population."""
    # Build pop → list of per-chr mean_froh
    pop_chr_froh: dict[str, list[float]] = {}
    for key, entry in fixed_roh_hmm.items():
        pop, chrom = key.split("|")
        if chrom not in CHR_LENGTHS_HG38:
            continue
        pop_chr_froh.setdefault(pop, []).append(entry["mean_froh"])

    fixed = {}
    for pop, entry in roh_hmm_summary.items():
        new_entry = dict(entry)
        froh_list = pop_chr_froh.get(pop, [])
        if froh_list:
            new_entry["mean_froh_across_chr"] = sum(froh_list) / len(froh_list)
        fixed[pop] = new_entry
    return fixed


def print_summary(roh_hmm_old: dict, roh_hmm_new: dict):
    """Print before/after comparison for the affected chromosomes."""
    affected = ["chr13", "chr14", "chr15", "chr21", "chr22"]
    print(f"\n{'Chr':<8} {'Pop':<6} {'Old mean_froh':>14} {'New mean_froh':>14} {'Change':>8}")
    print("-" * 56)
    for key in sorted(roh_hmm_old.keys()):
        pop, chrom = key.split("|")
        if chrom not in affected or pop != "CEU":
            continue
        old_f = roh_hmm_old[key]["mean_froh"]
        new_f = roh_hmm_new[key]["mean_froh"]
        pct = (new_f - old_f) / old_f * 100 if old_f else 0
        print(f"  {chrom:<8} {pop:<6} {old_f:>14.6f} {new_f:>14.6f} {pct:>+7.1f}%")

    # Count affected entries (> 1% change)
    n_changed = sum(
        1 for k, v in roh_hmm_new.items()
        if abs(v["mean_froh"] - roh_hmm_old[k]["mean_froh"]) /
           max(roh_hmm_old[k]["mean_froh"], 1e-12) > 0.01
    )
    print(f"\n  Entries with >1% change: {n_changed} / {len(roh_hmm_new)}")


def main():
    parser = argparse.ArgumentParser(description="Fix FROH values in results JSON")
    parser.add_argument("--results", default="human_interpretation_results.json",
                        help="Path to results file (default: human_interpretation_results.json)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print summary without writing")
    args = parser.parse_args()

    results_path = args.results
    if not os.path.exists(results_path):
        sys.exit(f"Error: {results_path} not found")

    print(f"Loading {results_path} ...")
    with open(results_path) as f:
        data = json.load(f)

    if "roh_hmm" not in data:
        sys.exit("Error: 'roh_hmm' key not found in results file")

    roh_hmm_old = data["roh_hmm"]
    print(f"  {len(roh_hmm_old)} roh_hmm entries (pop × chr)")

    fixed_roh_hmm     = fix_roh_hmm(roh_hmm_old)
    fixed_roh_summary = fix_roh_hmm_summary(data.get("roh_hmm_summary", {}), fixed_roh_hmm)

    print_summary(roh_hmm_old, fixed_roh_hmm)

    if args.dry_run:
        print("\nDry run — no file written.")
        return

    data["roh_hmm"]         = fixed_roh_hmm
    data["roh_hmm_summary"] = fixed_roh_summary

    tmp_path = results_path + ".froh_fix.tmp"
    print(f"\nWriting fixed results to {tmp_path} ...")
    with open(tmp_path, "w") as f:
        json.dump(data, f)
    os.rename(tmp_path, results_path)
    print(f"Done. {results_path} updated atomically.")


if __name__ == "__main__":
    main()
