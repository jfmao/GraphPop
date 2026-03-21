#!/usr/bin/env python3
"""Investigation 3: ROH × Sweep Population Cross-Correlation.

Question: Do populations with high runs-of-homozygosity also show more intense
positive selection sweeps? Is inbreeding correlated with sweep intensity across
the 26 1000 Genomes populations?

Method:
  For each of 26 populations, extract from human_interpretation_results.json:
    - FROH (mean across chromosomes): proxy for recent inbreeding / small N_e
    - Sweep intensity metrics: n_hard_sweeps, mean_h12, max_h12
    - π (nucleotide diversity): proxy for long-term N_e
    - Tajima's D: selection/bottleneck signal
    - Mean FST vs other continental groups: population isolation
  Compute Spearman correlations across populations.

Why novel: Multi-layer cross-population correlation (ROH × sweep × diversity)
  from a single JSON checkpoint — no separate VCF tool can compute all three
  layers simultaneously. The graph integrates individual-level (ROH) and
  population-level (sweep windows) signals.

Output: data/results/roh_sweep_correlation.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_03_roh_sweep_correlation.py
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

# ── Config ────────────────────────────────────────────────────────────────────
RESULTS_JSON = Path("human_interpretation_results.json")
OUT_DIR      = Path("data/results")
OUT_TSV      = OUT_DIR / "roh_sweep_correlation.tsv"
OUT_JSON     = OUT_DIR / "roh_sweep_correlation.json"

# Continental group assignments
CONTINENTAL = {
    "African":     {"ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"},
    "European":    {"CEU", "FIN", "GBR", "IBS", "TSI"},
    "East_Asian":  {"CDX", "CHB", "CHS", "JPT", "KHV"},
    "South_Asian": {"BEB", "GIH", "ITU", "PJL", "STU"},
    "American":    {"CLM", "MXL", "PEL", "PUR"},
}
POP_TO_GROUP = {p: g for g, pops in CONTINENTAL.items() for p in pops}

H12_HARD_THRESH = 0.3   # minimum h12 to count as a hard sweep

# ── Load data ─────────────────────────────────────────────────────────────────

def load_population_metrics(data: dict) -> dict:
    """Return dict: pop → metrics dict."""
    div_table   = data["annotate"]["diversity_table"]
    sweep_ann   = data["annotate"]["sweep_annotations"]
    roh_summary = data["roh_hmm_summary"]
    fst_matrix  = data["fst_matrix"]

    # Build FST lookup: pop → list of (partner, fst)
    pop_fst = defaultdict(list)
    for pair_key, pair_val in fst_matrix.items():
        p1, _, p2 = pair_key.partition("_vs_")
        mfst = pair_val.get("mean_fst", None)
        if mfst is not None:
            pop_fst[p1].append(mfst)
            pop_fst[p2].append(mfst)

    pops = sorted(div_table.keys())
    metrics = {}
    for pop in pops:
        div = div_table[pop]
        roh = roh_summary.get(pop, {})
        sa  = sweep_ann.get(pop, {})
        sweeps = sa.get("sweeps", [])

        # Sweep metrics
        hard = [s for s in sweeps if s.get("h12", 0) >= H12_HARD_THRESH]
        mean_h12 = float(np.mean([s["h12"] for s in hard])) if hard else 0.0
        max_h12  = float(np.max([s["h12"]  for s in hard])) if hard else 0.0

        # Mean FST against all partners
        all_fst = pop_fst.get(pop, [])
        mean_pop_fst = float(np.mean(all_fst)) if all_fst else np.nan

        metrics[pop] = {
            "pop":            pop,
            "group":          POP_TO_GROUP.get(pop, "Unknown"),
            "mean_pi":        div.get("mean_pi", np.nan),
            "mean_theta_w":   div.get("mean_theta_w", np.nan),
            "mean_tajima_d":  div.get("mean_tajima_d", np.nan),
            "mean_fis":       div.get("mean_fis", np.nan),
            "froh":           roh.get("mean_froh_across_chr", np.nan),
            "n_hard_sweeps":  len(hard),
            "mean_h12":       mean_h12,
            "max_h12":        max_h12,
            "mean_fst":       mean_pop_fst,
        }

    return metrics


def spearman_report(x_arr, y_arr, x_name, y_name, pop_names) -> dict:
    """Compute Spearman ρ and p-value; return info dict."""
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    xm, ym = x_arr[mask], y_arr[mask]
    if len(xm) < 5:
        return {"x": x_name, "y": y_name, "rho": None, "pval": None, "n": int(mask.sum())}
    rho, pval = spearmanr(xm, ym)
    return {"x": x_name, "y": y_name, "rho": round(float(rho), 4),
            "pval": round(float(pval), 6), "n": int(mask.sum())}


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 3: ROH × Sweep Cross-Correlation ===\n", flush=True)

    print("[1/3] Loading human_interpretation_results.json ...", flush=True)
    with open(RESULTS_JSON) as f:
        data = json.load(f)

    print("[2/3] Computing per-population metrics ...", flush=True)
    metrics = load_population_metrics(data)
    pops = sorted(metrics.keys())
    print(f"  {len(pops)} populations loaded", flush=True)

    # Print summary table
    print(f"\n{'Pop':<6} {'Group':<12} {'FROH':>7} {'n_sweep':>8} "
          f"{'mean_h12':>9} {'pi':>8} {'Taj_D':>7} {'Fis':>7} {'mean_FST':>9}")
    print("-" * 85)
    for pop in pops:
        m = metrics[pop]
        print(f"  {pop:<5} {m['group'][:10]:<11} {m['froh']:>7.4f} {m['n_hard_sweeps']:>8} "
              f"{m['mean_h12']:>9.4f} {m['mean_pi']:>8.5f} "
              f"{m['mean_tajima_d']:>7.3f} {m['mean_fis']:>7.4f} {m['mean_fst']:>9.5f}")

    # ── Spearman correlations ────────────────────────────────────────────────
    print("\n[3/3] Computing Spearman correlations across 26 populations ...", flush=True)

    keys = ["froh", "n_hard_sweeps", "mean_h12", "max_h12",
            "mean_pi", "mean_theta_w", "mean_tajima_d", "mean_fis", "mean_fst"]
    arrays = {k: np.array([metrics[p][k] for p in pops]) for k in keys}

    # Pairs of interest
    pairs = [
        ("froh", "n_hard_sweeps"),
        ("froh", "mean_h12"),
        ("froh", "mean_pi"),
        ("froh", "mean_tajima_d"),
        ("froh", "mean_fst"),
        ("mean_pi", "n_hard_sweeps"),
        ("mean_pi", "mean_h12"),
        ("mean_pi", "mean_tajima_d"),
        ("mean_pi", "mean_fst"),
        ("mean_tajima_d", "n_hard_sweeps"),
        ("mean_tajima_d", "mean_h12"),
        ("mean_fst", "n_hard_sweeps"),
        ("mean_fst", "mean_h12"),
        ("mean_fis", "froh"),
        ("mean_fis", "n_hard_sweeps"),
    ]

    correlations = []
    print(f"\n  {'X variable':<18} {'Y variable':<18} {'rho':>7}  {'pval':>10}  {'n':>4}")
    print("  " + "-" * 65)
    for xk, yk in pairs:
        res = spearman_report(arrays[xk], arrays[yk], xk, yk, pops)
        correlations.append(res)
        rho_s  = f"{res['rho']:>7.4f}" if res["rho"] is not None else "    N/A"
        pval_s = f"{res['pval']:>10.4f}" if res["pval"] is not None else "       N/A"
        star = ""
        if res["pval"] is not None:
            if res["pval"] < 0.001: star = "***"
            elif res["pval"] < 0.01:  star = "**"
            elif res["pval"] < 0.05:  star = "*"
        print(f"  {xk:<18} {yk:<18} {rho_s}  {pval_s}  {res['n']:>4}  {star}")

    # ── Write outputs ─────────────────────────────────────────────────────────
    # TSV of per-population metrics
    header = "pop\tgroup\tfroh\tn_hard_sweeps\tmean_h12\tmax_h12\tmean_pi\tmean_theta_w\tmean_tajima_d\tmean_fis\tmean_fst\n"
    with open(OUT_TSV, "w") as f:
        f.write(header)
        for pop in pops:
            m = metrics[pop]
            f.write(f"{pop}\t{m['group']}\t{m['froh']:.6f}\t{m['n_hard_sweeps']}\t"
                    f"{m['mean_h12']:.6f}\t{m['max_h12']:.6f}\t{m['mean_pi']:.6f}\t"
                    f"{m['mean_theta_w']:.6f}\t{m['mean_tajima_d']:.6f}\t"
                    f"{m['mean_fis']:.6f}\t{m['mean_fst']:.6f}\n")

    with open(OUT_JSON, "w") as f:
        json.dump({
            "n_populations": len(pops),
            "h12_hard_threshold": H12_HARD_THRESH,
            "per_population": metrics,
            "correlations": correlations,
        }, f, indent=2)

    # Key findings
    sig = [c for c in correlations if c["pval"] is not None and c["pval"] < 0.05]
    print(f"\n=== Summary ===")
    print(f"  {len(sig)}/{len(correlations)} correlations significant (p<0.05)")
    print(f"\nKey significant correlations:")
    for c in sorted(sig, key=lambda x: abs(x["rho"]), reverse=True)[:5]:
        direction = "↑" if c["rho"] > 0 else "↓"
        print(f"  {c['x']} vs {c['y']}: ρ={c['rho']:+.4f} {direction}  (p={c['pval']:.4f})")

    print(f"\nOutput: {OUT_TSV}")
    print(f"        {OUT_JSON}")


if __name__ == "__main__":
    main()
