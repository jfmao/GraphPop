#!/usr/bin/env python3
"""Investigation 10: Temporal Stratification of Selection Signals.

Question: Can we decompose selection signals by their approximate age across
26 populations, using three complementary statistics with different time scales?

Time proxies:
  VERY RECENT  (<5 kya):   FROH — runs-of-homozygosity; recent bottleneck/inbreeding
  RECENT       (~5-30 kya): iHS / nSL — haplotype-based, detects ongoing/recent sweeps
  INTERMEDIATE (~10-100 kya): H12 (Garud's H) — haplotype homozygosity, medium-age sweeps

Method:
  1. FROH: from human_interpretation_results.json roh_hmm_summary per pop
  2. iHS: per population, compute fraction of genome with |iHS| > 2 (standardized)
         and fraction with |nSL| > 2; from Variant nodes in Neo4j
  3. H12: from GenomicWindow nodes per population, compute fraction of windows
         with h12 >= 0.3 (sweep candidates)
  4. Build a 3-axis temporal fingerprint per population
  5. Cluster populations by fingerprint; identify populations with unusually
     elevated selection at specific time depths

Why novel: Temporal decomposition of selection using three graph-native layers.
  No single VCF tool integrates ROH + iHS + haplotype homozygosity simultaneously.

Output: data/results/temporal_selection.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_10_temporal_selection.py
"""

import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase
from scipy.stats import spearmanr

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"

RESULTS_JSON = Path("human_interpretation_results.json")
OUT_DIR      = Path("data/results")
OUT_TSV      = OUT_DIR / "temporal_selection.tsv"
OUT_JSON     = OUT_DIR / "temporal_selection.json"

IHS_THRESHOLD   = 2.0     # |iHS| > 2 = genome-wide significant
H12_THRESHOLD   = 0.3     # h12 >= 0.3 = sweep candidate

# iHS populations available in Neo4j Variant nodes
IHS_POPS = [
    "CEU","FIN","GBR","IBS","TSI",          # European
    "CHB","JPT","CDX","KHV",                 # East Asian
    "BEB","GIH","ITU","PJL","STU",           # South Asian
    "CLM","MXL","PEL","PUR",                 # American
]
# Populations in h12 sweep_annotations (all 26)
CONTINENTAL = {
    "African":     {"ACB","ASW","ESN","GWD","LWK","MSL","YRI"},
    "European":    {"CEU","FIN","GBR","IBS","TSI"},
    "East_Asian":  {"CDX","CHB","CHS","JPT","KHV"},
    "South_Asian": {"BEB","GIH","ITU","PJL","STU"},
    "American":    {"CLM","MXL","PEL","PUR"},
}
POP_TO_GROUP = {p: g for g, pops in CONTINENTAL.items() for p in pops}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 10: Temporal Stratification of Selection ===\n", flush=True)

    # ── Layer 1: FROH (very recent) ────────────────────────────────────────────
    print("[1/4] Loading FROH (very recent selection proxy) ...", flush=True)
    with open(RESULTS_JSON) as f:
        interp = json.load(f)

    roh_summary = interp["roh_hmm_summary"]
    froh_per_pop = {pop: v["mean_froh_across_chr"] for pop, v in roh_summary.items()}
    print(f"  {len(froh_per_pop)} populations with FROH", flush=True)

    # ── Layer 2: H12 (intermediate-age sweeps) ─────────────────────────────────
    print("\n[2/4] Computing H12 sweep fraction (intermediate timescale) ...", flush=True)
    sweep_ann = interp["annotate"]["sweep_annotations"]
    h12_fraction = {}   # pop → fraction of (non-trivial) windows with h12 >= threshold

    for pop, pop_data in sweep_ann.items():
        sweeps = pop_data.get("sweeps", [])
        if not sweeps:
            h12_fraction[pop] = 0.0
            continue
        h12_vals = np.array([s["h12"] for s in sweeps])
        h12_fraction[pop] = float(np.mean(h12_vals >= H12_THRESHOLD))

    print(f"  {len(h12_fraction)} populations with H12 data", flush=True)

    # ── Layer 3: iHS/nSL (recent selection) ────────────────────────────────────
    print("\n[3/4] Computing iHS/nSL selection fraction from Neo4j ...", flush=True)
    t0 = time.time()

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()

    ihs_fraction = {}   # pop → fraction of variants with |iHS| > threshold
    nsl_fraction = {}   # pop → fraction of variants with |nSL| > threshold
    ihs_mean_abs = {}   # pop → mean |iHS| (strength of signal)

    with driver.session(database=NEO4J_DB) as session:
        for pop in IHS_POPS:
            print(f"  {pop} ...", flush=True)
            # Count variants with |iHS| > threshold and total with any iHS
            r = session.run(f"""
                MATCH (v:Variant)
                WHERE v.ihs_{pop} IS NOT NULL
                RETURN count(v) AS n_total,
                       sum(CASE WHEN abs(v.ihs_{pop}) > $thresh THEN 1 ELSE 0 END) AS n_sig,
                       avg(abs(v.ihs_{pop})) AS mean_abs
            """, thresh=IHS_THRESHOLD).single()

            if r and r["n_total"] > 0:
                ihs_fraction[pop] = float(r["n_sig"]) / float(r["n_total"])
                ihs_mean_abs[pop] = float(r["mean_abs"]) if r["mean_abs"] else 0.0
            else:
                ihs_fraction[pop] = 0.0
                ihs_mean_abs[pop] = 0.0

            # nSL
            r2 = session.run(f"""
                MATCH (v:Variant)
                WHERE v.nsl_{pop} IS NOT NULL
                RETURN count(v) AS n_total,
                       sum(CASE WHEN abs(v.nsl_{pop}) > $thresh THEN 1 ELSE 0 END) AS n_sig
            """, thresh=IHS_THRESHOLD).single()

            if r2 and r2["n_total"] > 0:
                nsl_fraction[pop] = float(r2["n_sig"]) / float(r2["n_total"])
            else:
                nsl_fraction[pop] = 0.0

    driver.close()
    print(f"  Done in {time.time()-t0:.1f}s", flush=True)

    # ── Assemble temporal fingerprint ─────────────────────────────────────────
    print("\n[4/4] Assembling temporal fingerprint and ranking ...", flush=True)

    # All populations that have at least one layer
    all_pops = sorted(set(froh_per_pop) | set(h12_fraction) | set(ihs_fraction))

    rows = []
    for pop in all_pops:
        froh    = froh_per_pop.get(pop, np.nan)
        h12_fr  = h12_fraction.get(pop, np.nan)
        ihs_fr  = ihs_fraction.get(pop, np.nan)
        nsl_fr  = nsl_fraction.get(pop, np.nan)
        ihs_abs = ihs_mean_abs.get(pop, np.nan)

        # Classify temporal bias: which time layer is relatively elevated?
        # Normalize each across all populations and compare
        rows.append({
            "pop":           pop,
            "group":         POP_TO_GROUP.get(pop, "Unknown"),
            "froh":          froh,          # very recent
            "h12_fraction":  h12_fr,        # intermediate
            "ihs_fraction":  ihs_fr,        # recent
            "nsl_fraction":  nsl_fr,        # recent
            "ihs_mean_abs":  ihs_abs,       # recent strength
        })

    # Z-score normalization across populations for each layer
    def zscore(vals):
        arr = np.array(vals, dtype=float)
        mask = np.isfinite(arr)
        if mask.sum() < 2:
            return arr
        arr[mask] = (arr[mask] - arr[mask].mean()) / (arr[mask].std() + 1e-12)
        return arr

    froh_z   = zscore([r["froh"]           for r in rows])
    h12_z    = zscore([r["h12_fraction"]    for r in rows])
    ihs_z    = zscore([r["ihs_fraction"]    for r in rows])

    for i, row in enumerate(rows):
        row["froh_z"]  = round(float(froh_z[i]),  3) if np.isfinite(froh_z[i])  else None
        row["h12_z"]   = round(float(h12_z[i]),   3) if np.isfinite(h12_z[i])   else None
        row["ihs_z"]   = round(float(ihs_z[i]),   3) if np.isfinite(ihs_z[i])   else None

        # Temporal bias: which layer is highest z-score?
        z_vals = {
            "very_recent(FROH)":         froh_z[i],
            "recent(iHS)":               ihs_z[i],
            "intermediate(H12)":         h12_z[i],
        }
        finite_z = {k: v for k, v in z_vals.items() if np.isfinite(v)}
        if finite_z:
            row["dominant_timescale"] = max(finite_z, key=finite_z.get)
        else:
            row["dominant_timescale"] = "unknown"

    # ── Print results ──────────────────────────────────────────────────────────
    print("\n=== Temporal Selection Fingerprint — All Populations ===")
    print(f"{'Pop':<6} {'Group':<12} {'FROH':>8} {'FROH_z':>7} "
          f"{'iHS_fr':>8} {'iHS_z':>7} {'H12_fr':>8} {'H12_z':>7}  Dominant timescale")
    print("-" * 95)
    for row in sorted(rows, key=lambda r: r["group"] + r["pop"]):
        froh_s = f"{row['froh']:.4f}"       if np.isfinite(row.get('froh', np.nan))       else "  N/A"
        ihs_s  = f"{row['ihs_fraction']:.4f}" if np.isfinite(row.get('ihs_fraction', np.nan)) else "  N/A"
        h12_s  = f"{row['h12_fraction']:.4f}" if np.isfinite(row.get('h12_fraction', np.nan)) else "  N/A"
        fz  = f"{row['froh_z']:+.2f}" if row['froh_z'] is not None else "  N/A"
        iz  = f"{row['ihs_z']:+.2f}"  if row['ihs_z']  is not None else "  N/A"
        hz  = f"{row['h12_z']:+.2f}"  if row['h12_z']  is not None else "  N/A"
        print(f"  {row['pop']:<5} {row['group'][:10]:<11} "
              f"{froh_s:>8} {fz:>7} {ihs_s:>8} {iz:>7} {h12_s:>8} {hz:>7}  "
              f"{row['dominant_timescale']}")

    # ── Cross-layer correlations ───────────────────────────────────────────────
    print("\n=== Cross-layer Spearman Correlations ===")
    have_all = [r for r in rows if all(r.get(k) is not None
                for k in ["froh_z", "ihs_z", "h12_z"])]
    if len(have_all) >= 5:
        fz_arr = np.array([r["froh_z"] for r in have_all])
        iz_arr = np.array([r["ihs_z"]  for r in have_all])
        hz_arr = np.array([r["h12_z"]  for r in have_all])
        for (x, xn), (y, yn) in [
            ((fz_arr, "FROH_z"),   (iz_arr, "iHS_z")),
            ((fz_arr, "FROH_z"),   (hz_arr, "H12_z")),
            ((iz_arr, "iHS_z"),    (hz_arr, "H12_z")),
        ]:
            rho, pval = spearmanr(x, y)
            star = "***" if pval<0.001 else ("**" if pval<0.01 else ("*" if pval<0.05 else ""))
            print(f"  {xn} vs {yn}: ρ={rho:+.4f}  p={pval:.4f}  {star}")

    # ── Summary: dominant timescale per continental group ─────────────────────
    print("\n=== Dominant Selection Timescale by Continental Group ===")
    group_ts = defaultdict(list)
    for row in rows:
        group_ts[row["group"]].append(row["dominant_timescale"])
    for grp, ts_list in sorted(group_ts.items()):
        from collections import Counter
import os
        most_common = Counter(ts_list).most_common(1)[0]
        print(f"  {grp:<12}: {most_common[0]}  ({most_common[1]}/{len(ts_list)} pops)")

    # ── Write outputs ─────────────────────────────────────────────────────────
    with open(OUT_TSV, "w") as f:
        f.write("pop\tgroup\tfroh\tfroh_z\tihs_fraction\tihs_z\t"
                "nsl_fraction\th12_fraction\th12_z\tihs_mean_abs\tdominant_timescale\n")
        for row in rows:
            f.write(
                f"{row['pop']}\t{row['group']}\t"
                f"{row.get('froh','')}\t{row.get('froh_z','')}\t"
                f"{row.get('ihs_fraction','')}\t{row.get('ihs_z','')}\t"
                f"{row.get('nsl_fraction','')}\t{row.get('h12_fraction','')}\t"
                f"{row.get('h12_z','')}\t{row.get('ihs_mean_abs','')}\t"
                f"{row.get('dominant_timescale','')}\n"
            )

    with open(OUT_JSON, "w") as f:
        json.dump({
            "config": {
                "ihs_threshold":  IHS_THRESHOLD,
                "h12_threshold":  H12_THRESHOLD,
                "ihs_pops":       IHS_POPS,
            },
            "populations": rows,
        }, f, indent=2)

    print(f"\nOutput: {OUT_TSV}")
    print(f"        {OUT_JSON}")


if __name__ == "__main__":
    main()
