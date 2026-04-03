#!/usr/bin/env python3
"""Investigation 11: Functional Consequence Bias in Selection.

Question: Does positive selection preferentially act through HIGH-impact coding changes
(missense, stop-gained) or LOW-impact regulatory/synonymous variants? Does this differ
by population?

Method:
  Stream Variant→Gene (HAS_CONSEQUENCE) with allele count arrays.
  Group variants by VEP impact: HIGH (stop_gained, frameshift), MODERATE (missense),
  LOW (synonymous, splice_region), MODIFIER (intronic/UTR/intergenic).
  Per category: compute Hudson FST for YRI/CEU, YRI/CHB, CEU/CHB; compute fraction
  of variants with |iHS| > 2 (positive selection enrichment).
  Statistical test: Mann-Whitney U between FST distributions of HIGH vs LOW impact.

Why graph-native:
  HAS_CONSEQUENCE edges carry consequence + impact co-located with allele count arrays
  on Variant nodes. One Cypher traversal joins all layers simultaneously.

Output: data/results/consequence_bias.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_11_consequence_bias.py
"""

import csv
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase
from scipy.stats import mannwhitneyu
import os

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"

OUT_DIR  = Path("data/results")
OUT_TSV  = OUT_DIR / "consequence_bias.tsv"
OUT_JSON = OUT_DIR / "consequence_bias.json"

# Population order in variant cache (and in Variant.pop_ids)
POP_IDS = ['ACB','ASW','BEB','CDX','CEU','CHB','CHS','CLM','ESN','FIN',
           'GBR','GIH','GWD','IBS','ITU','JPT','KHV','LWK','MSL','MXL',
           'PEL','PJL','PUR','STU','TSI','YRI']
IDX = {p: i for i, p in enumerate(POP_IDS)}

FOCAL_PAIRS = [
    ("YRI", "CEU"),
    ("YRI", "CHB"),
    ("CEU", "CHB"),
]

# iHS populations (those with standardized iHS stored on Variant nodes)
IHS_POPS = ["CEU", "CHB", "GIH", "FIN", "JPT"]
IHS_THRESHOLD = 2.0  # |iHS| > 2 → strong positive selection signal

MIN_AN = 10   # minimum allele number per population


# ── Hudson FST (ratio-of-averages estimator) ──────────────────────────────────
def hudson_fst_num_den(ac, an, i1, i2):
    """Return numerator and denominator for Hudson FST at a single variant."""
    n1, n2 = an[i1], an[i2]
    if n1 < 2 or n2 < 2:
        return None, None
    p1 = ac[i1] / n1
    p2 = ac[i2] / n2
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    if den <= 0:
        return None, None
    return float(num), float(den)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print("=== Investigation 11: Functional Consequence Bias in Selection ===\n", flush=True)

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()

    # Per-impact accumulators:
    #   fst_nums[impact][pair] = list of numerators
    #   fst_dens[impact][pair] = list of denominators
    #   fst_per_var[impact][pair] = list of per-variant FST values (for Mann-Whitney)
    #   ihs_counts[impact][pop] = (n_abs_gt2, n_total_with_ihs)
    fst_nums    = defaultdict(lambda: defaultdict(list))
    fst_dens    = defaultdict(lambda: defaultdict(list))
    fst_per_var = defaultdict(lambda: defaultdict(list))
    ihs_counts  = defaultdict(lambda: defaultdict(lambda: [0, 0]))  # [n_enriched, n_total]
    n_variants  = defaultdict(int)
    n_unique_v  = 0  # total unique variants with HAS_CONSEQUENCE

    CYPHER = """
    MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)
    WHERE v.ac IS NOT NULL AND v.an IS NOT NULL
    RETURN c.impact     AS impact,
           v.ac         AS ac,
           v.an         AS an,
           v.pop_ids    AS pop_ids,
           v.ihs_CEU    AS ihs_ceu,
           v.ihs_CHB    AS ihs_chb,
           v.ihs_GIH    AS ihs_gih,
           v.ihs_FIN    AS ihs_fin,
           v.ihs_JPT    AS ihs_jpt
    """

    print("Streaming HAS_CONSEQUENCE edges from Neo4j...", flush=True)
    seen_variants = defaultdict(set)  # impact → set of (ac_hash) to deduplicate per-variant FST

    with driver.session(database=NEO4J_DB) as session:
        result = session.run(CYPHER)
        for i, rec in enumerate(result):
            if i % 100_000 == 0 and i > 0:
                elapsed = time.time() - t0
                print(f"  {i:,} records | {elapsed:.0f}s elapsed", flush=True)

            impact   = rec["impact"] or "MODIFIER"
            ac_raw   = rec["ac"]
            an_raw   = rec["an"]
            pop_ids  = list(rec["pop_ids"]) if rec["pop_ids"] else []

            if not ac_raw or not an_raw or not pop_ids:
                continue

            ac = np.asarray(ac_raw, dtype=np.int32)
            an = np.asarray(an_raw, dtype=np.int32)

            # Map local pop_ids to global POP_IDS indices
            local_idx = {}
            for j, p in enumerate(pop_ids):
                if p in IDX:
                    local_idx[p] = j

            n_variants[impact] += 1
            n_unique_v += 1

            # Compute FST per focal pair
            for p1_name, p2_name in FOCAL_PAIRS:
                pair_key = f"{p1_name}/{p2_name}"
                j1 = local_idx.get(p1_name)
                j2 = local_idx.get(p2_name)
                if j1 is None or j2 is None:
                    continue
                num, den = hudson_fst_num_den(ac, an, j1, j2)
                if num is not None:
                    fst_nums[impact][pair_key].append(num)
                    fst_dens[impact][pair_key].append(den)
                    # Per-variant FST for Mann-Whitney
                    if den > 0:
                        fst_per_var[impact][pair_key].append(num / den)

            # iHS enrichment
            ihs_vals = {
                "CEU": rec["ihs_ceu"],
                "CHB": rec["ihs_chb"],
                "GIH": rec["ihs_gih"],
                "FIN": rec["ihs_fin"],
                "JPT": rec["ihs_jpt"],
            }
            for pop, val in ihs_vals.items():
                if val is not None:
                    ihs_counts[impact][pop][1] += 1
                    if abs(val) > IHS_THRESHOLD:
                        ihs_counts[impact][pop][0] += 1

    elapsed = time.time() - t0
    print(f"\nStreaming complete: {n_unique_v:,} records in {elapsed:.1f}s", flush=True)

    # ── Compute mean FST per impact category ─────────────────────────────────
    IMPACTS = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
    results = []

    for impact in IMPACTS:
        row = {"impact": impact, "n_variants": n_variants.get(impact, 0)}
        for p1_name, p2_name in FOCAL_PAIRS:
            pair_key = f"{p1_name}/{p2_name}"
            nums = fst_nums[impact][pair_key]
            dens = fst_dens[impact][pair_key]
            if dens:
                fst_val = sum(nums) / sum(dens)
                row[f"fst_{p1_name}_{p2_name}"] = round(fst_val, 6)
                row[f"n_{p1_name}_{p2_name}"] = len(nums)
            else:
                row[f"fst_{p1_name}_{p2_name}"] = None
                row[f"n_{p1_name}_{p2_name}"] = 0

        # iHS enrichment
        for pop in IHS_POPS:
            n_enr, n_tot = ihs_counts[impact][pop]
            row[f"ihs_enrich_{pop}"] = round(n_enr / n_tot, 4) if n_tot > 0 else None
            row[f"ihs_n_{pop}"] = n_tot

        results.append(row)

    # ── Mann-Whitney U: HIGH vs LOW impact FST (main focal pair) ──────────────
    mw_results = {}
    for p1_name, p2_name in FOCAL_PAIRS:
        pair_key = f"{p1_name}/{p2_name}"
        high_fsts = fst_per_var["HIGH"][pair_key]
        low_fsts  = fst_per_var["LOW"][pair_key]
        mod_fsts  = fst_per_var["MODERATE"][pair_key]
        mw_key    = pair_key

        mw_results[mw_key] = {}
        for name, grp in [("HIGH_vs_LOW", (high_fsts, low_fsts)),
                          ("HIGH_vs_MODERATE", (high_fsts, mod_fsts)),
                          ("MODERATE_vs_LOW", (mod_fsts, low_fsts))]:
            a, b = grp
            if len(a) >= 5 and len(b) >= 5:
                stat, pval = mannwhitneyu(a, b, alternative="two-sided")
                mw_results[mw_key][name] = {
                    "U": round(float(stat), 2),
                    "pval": float(f"{pval:.3e}"),
                    "n1": len(a),
                    "n2": len(b),
                    "mean1": round(float(np.mean(a)), 4),
                    "mean2": round(float(np.mean(b)), 4),
                }

    driver.close()

    # ── Biological interpretation ─────────────────────────────────────────────
    # HIGH < LOW FST → purifying selection dominates HIGH-impact variants
    # HIGH > LOW FST → some HIGH-impact variants are positively selected
    for row in results:
        fst_high = row.get("fst_YRI_CEU")
        fst_low_row = next((r for r in results if r["impact"] == "LOW"), None)
        if fst_low_row and fst_high is not None:
            fst_low_val = fst_low_row.get("fst_YRI_CEU")
            if fst_low_val:
                row["interpretation"] = (
                    "purifying_selection_dominant"
                    if fst_high < fst_low_val else "mixed_or_positive_selection"
                )

    # ── Write outputs ─────────────────────────────────────────────────────────
    fieldnames = ["impact", "n_variants"]
    for p1, p2 in FOCAL_PAIRS:
        fieldnames += [f"fst_{p1}_{p2}", f"n_{p1}_{p2}"]
    for pop in IHS_POPS:
        fieldnames += [f"ihs_enrich_{pop}", f"ihs_n_{pop}"]
    fieldnames.append("interpretation")

    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    out_data = {
        "n_total_records": n_unique_v,
        "focal_pairs": [f"{p1}/{p2}" for p1, p2 in FOCAL_PAIRS],
        "ihs_pops": IHS_POPS,
        "ihs_threshold": IHS_THRESHOLD,
        "impact_stats": results,
        "mann_whitney": mw_results,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2, default=str)

    print(f"\n=== Results ===")
    for row in results:
        print(f"  {row['impact']:10s}: n={row['n_variants']:6,d}  "
              f"FST(YRI/CEU)={row.get('fst_YRI_CEU', 'N/A')}  "
              f"iHS_CEU_enrich={row.get('ihs_enrich_CEU', 'N/A')}")

    print(f"\n  Mann-Whitney (YRI/CEU, HIGH vs LOW):")
    mw = mw_results.get("YRI/CEU", {})
    if "HIGH_vs_LOW" in mw:
        m = mw["HIGH_vs_LOW"]
        print(f"    HIGH mean_FST={m['mean1']:.4f}  LOW mean_FST={m['mean2']:.4f}  p={m['pval']}")

    print(f"\nOutputs: {OUT_TSV}  {OUT_JSON}")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
