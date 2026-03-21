#!/usr/bin/env python3
"""Investigation 7: Out-of-Africa Genetic Gradient Analysis.

Question: How do population genetic parameters vary along the Out-of-Africa
dispersal route? Can we reconstruct the OoA bottleneck signature and rank
populations by evolutionary distance from Africa?

Method:
  For each of 26 populations:
    - Geographic distance from Addis Ababa, Ethiopia (center of OoA origin)
    - Genetic metrics from human_interpretation_results.json:
        π (nucleotide diversity), θ_W, Tajima's D, FROH, mean FST,
        sweep fraction, mean h12
  Compute Spearman correlations with geographic distance.
  PCA of the evolutionary fingerprint across 26 populations.

Why novel: Multi-layer evolutionary gradient (ROH + diversity + sweeps)
  across populations from a single graph database checkpoint.

Output: data/results/ooa_gradient.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_07_ooa_gradient.py
"""

import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

# ── Config ────────────────────────────────────────────────────────────────────
RESULTS_JSON = Path("human_interpretation_results.json")
OUT_DIR      = Path("data/results")
OUT_TSV      = OUT_DIR / "ooa_gradient.tsv"
OUT_JSON     = OUT_DIR / "ooa_gradient.json"

H12_HARD_THRESH = 0.3

# Geographic coordinates (latitude, longitude) for 1000G populations
# Sources: 1000 Genomes sample metadata + published population genetics literature
POP_COORDS = {
    "YRI": (7.38,   3.90),    # Yoruba, Ibadan, Nigeria
    "LWK": (0.27,  36.22),    # Luhya, Webuye, Kenya
    "GWD": (13.67, -14.80),   # Gambian
    "MSL": (8.50,  -13.30),   # Mende, Sierra Leone
    "ESN": (6.00,    8.00),   # Esan, Nigeria
    "ACB": (13.17, -59.54),   # African Caribbean, Barbados
    "ASW": (35.50, -83.80),   # African American SW USA
    "CEU": (45.30,   5.80),   # Utah Caucasian (European)
    "TSI": (43.78,  11.25),   # Tuscany, Italy
    "GBR": (51.51,  -0.13),   # British, England
    "FIN": (61.00,  25.70),   # Finnish
    "IBS": (40.42,  -3.70),   # Iberian, Spain
    "CHB": (39.91, 116.39),   # Han Chinese, Beijing
    "JPT": (35.68, 139.69),   # Japanese, Tokyo
    "CHS": (25.00, 120.00),   # Han Chinese, South
    "CDX": (25.00, 101.54),   # Dai Chinese, Yunnan
    "KHV": (10.82, 106.63),   # Kinh Vietnamese, Ho Chi Minh
    "BEB": (23.68,  90.35),   # Bengali, Bangladesh
    "GIH": (22.30,  70.73),   # Gujarati Indian
    "ITU": (13.08,  80.27),   # Indian Telugu, UK
    "PJL": (31.55,  74.34),   # Punjabi, Lahore
    "STU": (7.00,   80.00),   # Sri Lankan Tamil, UK
    "CLM": (-4.00, -73.00),   # Colombian
    "MXL": (23.00, -102.00),  # Mexican Ancestry, LA
    "PEL": (-12.05, -77.03),  # Peruvian, Lima
    "PUR": (18.46,  -66.10),  # Puerto Rican
}

# OoA origin: Addis Ababa, Ethiopia
OOA_LAT, OOA_LON = 9.03, 38.74

CONTINENTAL = {
    "African":     {"ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"},
    "European":    {"CEU", "FIN", "GBR", "IBS", "TSI"},
    "East_Asian":  {"CDX", "CHB", "CHS", "JPT", "KHV"},
    "South_Asian": {"BEB", "GIH", "ITU", "PJL", "STU"},
    "American":    {"CLM", "MXL", "PEL", "PUR"},
}
POP_TO_GROUP = {p: g for g, pops in CONTINENTAL.items() for p in pops}


def haversine_km(lat1, lon1, lat2, lon2) -> float:
    """Great-circle distance in km."""
    R = 6371.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = math.sin(dphi/2)**2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam/2)**2
    return R * 2 * math.asin(math.sqrt(a))


def load_metrics(data: dict) -> dict:
    div_table   = data["annotate"]["diversity_table"]
    sweep_ann   = data["annotate"]["sweep_annotations"]
    roh_summary = data["roh_hmm_summary"]
    fst_matrix  = data["fst_matrix"]

    pop_fst = defaultdict(list)
    for pair_key, pair_val in fst_matrix.items():
        p1, _, p2 = pair_key.partition("_vs_")
        mfst = pair_val.get("mean_fst")
        if mfst is not None:
            pop_fst[p1].append(mfst)
            pop_fst[p2].append(mfst)

    metrics = {}
    for pop, (lat, lon) in POP_COORDS.items():
        div = div_table.get(pop, {})
        roh = roh_summary.get(pop, {})
        sa  = sweep_ann.get(pop, {})
        sweeps = sa.get("sweeps", [])
        hard = [s for s in sweeps if s.get("h12", 0) >= H12_HARD_THRESH]

        all_h12 = [s["h12"] for s in hard]
        mean_fst_pop = float(np.mean(pop_fst[pop])) if pop_fst[pop] else np.nan

        dist_km = haversine_km(OOA_LAT, OOA_LON, lat, lon)

        metrics[pop] = {
            "pop":           pop,
            "group":         POP_TO_GROUP.get(pop, "Unknown"),
            "lat":           lat,
            "lon":           lon,
            "dist_ooa_km":   round(dist_km, 1),
            "mean_pi":       div.get("mean_pi", np.nan),
            "mean_theta_w":  div.get("mean_theta_w", np.nan),
            "mean_tajima_d": div.get("mean_tajima_d", np.nan),
            "mean_fis":      div.get("mean_fis", np.nan),
            "froh":          roh.get("mean_froh_across_chr", np.nan),
            "n_hard_sweeps": len(hard),
            "mean_h12":      float(np.mean(all_h12)) if all_h12 else 0.0,
            "mean_fst":      mean_fst_pop,
        }

    return metrics


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 7: Out-of-Africa Genetic Gradient ===\n", flush=True)

    print("[1/3] Loading human_interpretation_results.json ...", flush=True)
    with open(RESULTS_JSON) as f:
        data = json.load(f)

    print("[2/3] Computing per-population metrics and OoA distances ...", flush=True)
    metrics = load_metrics(data)
    pops = sorted(metrics.keys(), key=lambda p: metrics[p]["dist_ooa_km"])

    print(f"\n{'Pop':<6} {'Group':<12} {'Dist_km':>8} {'π':>8} {'FROH':>7} "
          f"{'n_sweep':>8} {'mean_h12':>9} {'Taj_D':>7}")
    print("-" * 80)
    for pop in pops:
        m = metrics[pop]
        print(f"  {pop:<5} {m['group'][:10]:<11} {m['dist_ooa_km']:>8.0f} "
              f"{m['mean_pi']:>8.5f} {m['froh']:>7.4f} "
              f"{m['n_hard_sweeps']:>8} {m['mean_h12']:>9.4f} {m['mean_tajima_d']:>7.3f}")

    # ── Spearman correlations with OoA distance ───────────────────────────────
    print("\n[3/3] Spearman correlations with OoA distance ...", flush=True)

    dist = np.array([metrics[p]["dist_ooa_km"] for p in pops])
    metric_keys = ["mean_pi", "mean_theta_w", "mean_tajima_d", "mean_fis",
                   "froh", "n_hard_sweeps", "mean_h12", "mean_fst"]

    correlations = []
    print(f"\n  {'Metric':<18} {'rho':>7}  {'pval':>10}  sig")
    print("  " + "-" * 45)
    for mk in metric_keys:
        arr = np.array([metrics[p][mk] for p in pops])
        mask = np.isfinite(arr) & np.isfinite(dist)
        if mask.sum() < 5:
            continue
        rho, pval = spearmanr(dist[mask], arr[mask])
        star = ""
        if pval < 0.001: star = "***"
        elif pval < 0.01:  star = "**"
        elif pval < 0.05:  star = "*"
        print(f"  {mk:<18} {rho:>+7.4f}  {pval:>10.4f}  {star}")
        correlations.append({"metric": mk, "rho": round(float(rho), 4),
                              "pval": round(float(pval), 6), "n": int(mask.sum())})

    # ── Simple PCA of evolutionary fingerprint ────────────────────────────────
    pca_keys = ["mean_pi", "mean_tajima_d", "mean_fis", "froh", "n_hard_sweeps", "mean_h12"]
    X = np.array([[metrics[p][k] for k in pca_keys] for p in pops], dtype=float)
    # Standardize
    X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-12)
    # SVD-based PCA
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    pcs = U * S   # (26, n_features) projected scores
    var_exp = (S**2) / (S**2).sum()

    print(f"\nPCA variance explained: PC1={100*var_exp[0]:.1f}%, PC2={100*var_exp[1]:.1f}%")
    print(f"\n{'Pop':<6} {'Group':<12} {'PC1':>7} {'PC2':>7} {'Dist_km':>8}")
    print("-" * 50)
    for i, pop in enumerate(pops):
        m = metrics[pop]
        print(f"  {pop:<5} {m['group'][:10]:<11} {pcs[i,0]:>+7.3f} {pcs[i,1]:>+7.3f} {m['dist_ooa_km']:>8.0f}")

    # Correlation of PC1 with OoA distance
    rho_pc1, pval_pc1 = spearmanr(dist, pcs[:, 0])
    rho_pc2, pval_pc2 = spearmanr(dist, pcs[:, 1])
    print(f"\nPC1 vs OoA distance: ρ={rho_pc1:+.4f}  p={pval_pc1:.4f}")
    print(f"PC2 vs OoA distance: ρ={rho_pc2:+.4f}  p={pval_pc2:.4f}")

    # ── Write outputs ─────────────────────────────────────────────────────────
    with open(OUT_TSV, "w") as f:
        f.write("pop\tgroup\tlat\tlon\tdist_ooa_km\tmean_pi\tmean_theta_w\t"
                "mean_tajima_d\tmean_fis\tfroh\tn_hard_sweeps\tmean_h12\tmean_fst\tpc1\tpc2\n")
        for i, pop in enumerate(pops):
            m = metrics[pop]
            f.write(f"{pop}\t{m['group']}\t{m['lat']}\t{m['lon']}\t{m['dist_ooa_km']:.1f}\t"
                    f"{m['mean_pi']:.6f}\t{m['mean_theta_w']:.6f}\t{m['mean_tajima_d']:.6f}\t"
                    f"{m['mean_fis']:.6f}\t{m['froh']:.6f}\t{m['n_hard_sweeps']}\t"
                    f"{m['mean_h12']:.6f}\t{m['mean_fst']:.6f}\t"
                    f"{pcs[i,0]:.4f}\t{pcs[i,1]:.4f}\n")

    with open(OUT_JSON, "w") as f:
        json.dump({
            "ooa_origin": {"lat": OOA_LAT, "lon": OOA_LON, "place": "Addis Ababa, Ethiopia"},
            "h12_hard_threshold": H12_HARD_THRESH,
            "pca_variance_explained": [round(float(v), 4) for v in var_exp[:4].tolist()],
            "pca_loadings": {pca_keys[j]: round(float(Vt[0, j]), 4) for j in range(len(pca_keys))},
            "pc1_vs_ooa_dist": {"rho": round(float(rho_pc1), 4), "pval": round(float(pval_pc1), 6)},
            "pc2_vs_ooa_dist": {"rho": round(float(rho_pc2), 4), "pval": round(float(pval_pc2), 6)},
            "per_population": [
                {**metrics[pops[i]], "pc1": round(float(pcs[i, 0]), 4), "pc2": round(float(pcs[i, 1]), 4)}
                for i in range(len(pops))
            ],
            "correlations_with_ooa_dist": correlations,
        }, f, indent=2)

    print(f"\nOutput: {OUT_TSV}")
    print(f"        {OUT_JSON}")


if __name__ == "__main__":
    main()
