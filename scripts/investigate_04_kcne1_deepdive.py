#!/usr/bin/env python3
"""Investigation 4: KCNE1 Convergent Sweep Deep-Dive.

Question: Is the KCNE1 hard sweep signal (observed in ALL 5 continental groups
in Investigation 2) evidence of an ancient pre-OoA sweep vs. parallel convergent
adaptation? How does FST at the KCNE1 locus compare to the genome-wide background?

KCNE1 (KCNE1 — potassium voltage-gated channel subfamily E member 1):
  - Cardiac K+ channel beta subunit (IKs channel with KCNQ1)
  - Implicated in cardiac repolarization, long QT syndrome (LQT5)
  - Sweep in ALL 5 continental groups with h12=1.0 → extraordinary signal

Method:
  1. Query Neo4j for GenomicWindow nodes overlapping KCNE1 → get sweep windows per pop
  2. Query Variant nodes with HAS_CONSEQUENCE → KCNE1 → get allele frequencies
  3. Compute Hudson FST at KCNE1 locus (variants within sweep windows)
  4. Compare to genome-wide FST baseline from fst_matrix
  5. Test hypothesis: ancient pre-OoA sweep → LOW FST; convergent → HIGH FST
  6. Retrieve all KCNE1 consequence variants and annotate with population-level context

Output: data/results/kcne1_deepdive.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_04_kcne1_deepdive.py
"""

import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB   = "neo4j"

RESULTS_JSON = Path("human_interpretation_results.json")
OUT_DIR      = Path("data/results")
OUT_TSV      = OUT_DIR / "kcne1_deepdive.tsv"
OUT_JSON     = OUT_DIR / "kcne1_deepdive.json"

TARGET_GENE  = "KCNE1"
H12_THRESH   = 0.3    # minimum h12 for sweep classification

FOCAL_PAIRS = [
    ("YRI", "CEU"),   # Africa vs Europe
    ("YRI", "CHB"),   # Africa vs East Asia
    ("CEU", "CHB"),   # Europe vs East Asia
    ("YRI", "GIH"),   # Africa vs South Asia
    ("CEU", "GIH"),   # Europe vs South Asia
]

CONTINENTAL = {
    "African":     {"ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"},
    "European":    {"CEU", "FIN", "GBR", "IBS", "TSI"},
    "East_Asian":  {"CDX", "CHB", "CHS", "JPT", "KHV"},
    "South_Asian": {"BEB", "GIH", "ITU", "PJL", "STU"},
    "American":    {"CLM", "MXL", "PEL", "PUR"},
}
POP_TO_GROUP = {p: g for g, pops in CONTINENTAL.items() for p in pops}


def hudson_fst(p1, n1, p2, n2):
    """Return (num, den) for Hudson FST estimator (Bhatia et al. 2013)."""
    if n1 < 2 or n2 < 2:
        return None, None
    num = (p1 - p2)**2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    den = p1*(1-p2) + p2*(1-p1)
    if den <= 0:
        return None, None
    return num, den


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"=== Investigation 4: {TARGET_GENE} Convergent Sweep Deep-Dive ===\n", flush=True)

    # ── Load genome-wide FST baseline ─────────────────────────────────────────
    print("[1/5] Loading genome-wide FST baseline ...", flush=True)
    with open(RESULTS_JSON) as f:
        interp = json.load(f)
    fst_matrix = interp["fst_matrix"]

    baseline_fst = {}
    for pair_key, pair_val in fst_matrix.items():
        p1, _, p2 = pair_key.partition("_vs_")
        mfst = pair_val.get("mean_fst")
        if mfst is not None:
            baseline_fst[f"{p1}_{p2}"] = mfst
            baseline_fst[f"{p2}_{p1}"] = mfst

    print(f"  {len(baseline_fst)//2} FST pairs loaded", flush=True)

    # ── Load KCNE1 sweep windows from human_interpretation_results ────────────
    print(f"\n[2/5] Loading {TARGET_GENE} sweep windows from interpretation results ...", flush=True)
    sweep_ann = interp["annotate"]["sweep_annotations"]
    kcne1_sweeps = defaultdict(list)   # pop → list of sweep dicts with chr/start/end

    # KCNE1 is on chr21 (known biology)
    # Find sweep windows that overlap KCNE1 in the sweep annotations
    # The sweep windows don't directly tag genes in all cases, so we'll use Neo4j for that
    # But we can also use the 'genes' field if populated
    for pop, pop_data in sweep_ann.items():
        for sw in pop_data.get("sweeps", []):
            genes = sw.get("genes", []) + sw.get("known_selection_genes", [])
            if TARGET_GENE in genes:
                kcne1_sweeps[pop].append(sw)

    print(f"  {TARGET_GENE} sweep windows found in {len(kcne1_sweeps)} populations:", flush=True)
    for pop, sws in sorted(kcne1_sweeps.items()):
        for sw in sws:
            print(f"    {pop}: {sw['chr']}:{sw['start']}-{sw['end']} "
                  f"h12={sw['h12']:.4f} ({sw.get('sweep_type','?')})", flush=True)

    # ── Query Neo4j for KCNE1 variants ───────────────────────────────────────
    print(f"\n[3/5] Querying Neo4j for {TARGET_GENE} variants ...", flush=True)
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()

    with driver.session(database=NEO4J_DB) as session:
        # Get all variants with HAS_CONSEQUENCE → KCNE1
        result = session.run("""
            MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene {geneId: $gene})
            WHERE v.chr IS NOT NULL AND v.pos IS NOT NULL AND v.ac IS NOT NULL
            RETURN v.variantId AS vid, v.chr AS chr, v.pos AS pos,
                   v.ac AS ac, v.an AS an, v.af AS af,
                   v.pop_ids AS pop_ids, v.variant_type AS vtype,
                   c.consequence AS consequence, c.impact AS impact
            ORDER BY v.pos
        """, gene=TARGET_GENE)

        kcne1_variants = []
        pop_ids_global = None
        for rec in result:
            if pop_ids_global is None and rec["pop_ids"] is not None:
                pop_ids_global = list(rec["pop_ids"])
                pop_idx = {p: i for i, p in enumerate(pop_ids_global)}
            kcne1_variants.append({
                "vid":         rec["vid"],
                "chr":         rec["chr"],
                "pos":         rec["pos"],
                "ac":          list(rec["ac"]) if rec["ac"] else [],
                "an":          list(rec["an"]) if rec["an"] else [],
                "af":          list(rec["af"]) if rec["af"] else [],
                "pop_ids":     list(rec["pop_ids"]) if rec["pop_ids"] else [],
                "vtype":       rec["vtype"],
                "consequence": rec["consequence"],
                "impact":      rec["impact"],
            })

        print(f"  {len(kcne1_variants)} {TARGET_GENE} variants found", flush=True)

        # Also query GenomicWindow nodes overlapping KCNE1 (from graph sweeps)
        # KCNE1 is on chr21 — find windows with h12 >= threshold for each pop
        result2 = session.run("""
            MATCH (w:GenomicWindow)
            WHERE w.h12 >= $threshold AND w.sweep_type = 'hard'
              AND w.chr = 'chr21'
            RETURN w.population AS pop, w.chr AS chr,
                   w.start AS start, w.end AS end,
                   w.h12 AS h12, w.h1 AS h1, w.h2_h1 AS h2_h1,
                   w.hap_diversity AS hap_div
            ORDER BY w.h12 DESC
        """, threshold=H12_THRESH)

        neo4j_sweeps = [dict(r) for r in result2]
        print(f"  {len(neo4j_sweeps)} hard sweep windows on chr21 in Neo4j", flush=True)

    driver.close()

    if not kcne1_variants:
        print(f"  WARNING: No {TARGET_GENE} variants found in Neo4j!", flush=True)
        return

    # Determine KCNE1 locus extent from variant positions
    positions = [v["pos"] for v in kcne1_variants]
    locus_chr   = kcne1_variants[0]["chr"]
    locus_start = min(positions)
    locus_end   = max(positions)
    print(f"\n  {TARGET_GENE} locus: {locus_chr}:{locus_start:,}-{locus_end:,} "
          f"({(locus_end-locus_start)//1000} kb, {len(kcne1_variants)} variants)", flush=True)

    # Filter Neo4j sweep windows overlapping KCNE1 locus
    kcne1_neo4j_sweeps = [
        w for w in neo4j_sweeps
        if w["start"] <= locus_end and w["end"] >= locus_start
    ]
    print(f"\n  Neo4j hard sweep windows overlapping {TARGET_GENE} locus:", flush=True)
    pops_with_sweep = set()
    for w in sorted(kcne1_neo4j_sweeps, key=lambda x: x["pop"]):
        group = POP_TO_GROUP.get(w["pop"], "?")
        print(f"    {w['pop']} ({group}): h12={w['h12']:.4f}  "
              f"hap_div={w.get('hap_div', 'N/A')}", flush=True)
        pops_with_sweep.add(w["pop"])

    groups_with_sweep = {POP_TO_GROUP.get(p) for p in pops_with_sweep if POP_TO_GROUP.get(p)}
    print(f"\n  Continental groups with sweep: {sorted(groups_with_sweep)}", flush=True)

    # ── Compute FST at KCNE1 locus ─────────────────────────────────────────────
    print(f"\n[4/5] Computing Hudson FST at {TARGET_GENE} locus vs genome-wide ...", flush=True)

    if pop_ids_global is None:
        print("  ERROR: could not determine pop_ids order", flush=True)
        return

    locus_fst = {}
    for pop_a, pop_b in FOCAL_PAIRS:
        if pop_a not in pop_idx or pop_b not in pop_idx:
            continue
        i, j = pop_idx[pop_a], pop_idx[pop_b]
        nums, dens = [], []
        for v in kcne1_variants:
            ac = v["ac"]
            an = v["an"]
            if not ac or not an or i >= len(ac) or j >= len(ac):
                continue
            if an[i] < 10 or an[j] < 10:
                continue
            p_i = ac[i] / an[i]
            p_j = ac[j] / an[j]
            n_i, n_j = an[i] // 2, an[j] // 2
            num, den = hudson_fst(p_i, n_i, p_j, n_j)
            if num is not None:
                nums.append(num)
                dens.append(den)

        if not nums:
            continue
        fst_locus = max(0.0, np.sum(nums) / np.sum(dens))
        pair_key = f"{pop_a}_{pop_b}"
        genome_fst = baseline_fst.get(pair_key, np.nan)
        ratio = fst_locus / genome_fst if genome_fst > 0 else np.nan
        locus_fst[pair_key] = {
            "pop_a": pop_a, "pop_b": pop_b,
            "n_variants": len(nums),
            "fst_locus": round(float(fst_locus), 6),
            "fst_genome": round(float(genome_fst), 6) if np.isfinite(genome_fst) else None,
            "fst_ratio": round(float(ratio), 3) if np.isfinite(ratio) else None,
        }

    print(f"\n  {'Pair':<15} {'FST_locus':>10} {'FST_genome':>11} "
          f"{'Ratio':>7}  {'Interpretation'}")
    print("  " + "-" * 70)
    for pk, vals in sorted(locus_fst.items()):
        ratio = vals["fst_ratio"]
        if ratio is None:
            interp_str = "insufficient data"
        elif ratio < 0.5:
            interp_str = "ANCIENT sweep (FST lower than genome)"
        elif ratio < 1.5:
            interp_str = "neutral sweep signal"
        else:
            interp_str = "CONVERGENT sweep (FST elevated)"
        print(f"  {pk:<15} {vals['fst_locus']:>10.5f} {vals['fst_genome'] or 0:>11.5f} "
              f"{ratio or 0:>7.2f}x  {interp_str}")

    # ── Allele frequency summary per population ────────────────────────────────
    print(f"\n[5/5] Allele frequency summary for {TARGET_GENE} variants ...", flush=True)

    # High-impact variants
    high_impact = [v for v in kcne1_variants if v.get("impact") in {"HIGH", "MODERATE"}]
    print(f"  {len(high_impact)} HIGH/MODERATE impact variants in {TARGET_GENE}", flush=True)

    if high_impact and pop_ids_global:
        print(f"\n  HIGH/MODERATE impact variant allele frequencies (top 10 by max AF):")
        print(f"  {'Variant':<30} {'Impact':<10} {'Consequence':<25}  "
              + "  ".join(f"{p[:3]}" for p in ["YRI","CEU","CHB","GIH","MXL"]))
        print("  " + "-" * 100)
        focal_pops = ["YRI", "CEU", "CHB", "GIH", "MXL"]
        focal_idx = [pop_idx.get(p) for p in focal_pops]

        # Sort by max AF across focal pops
        def max_af(v):
            af = v["af"]
            if not af:
                return 0.0
            return max((af[i] for i in focal_idx if i is not None and i < len(af)), default=0.0)

        for v in sorted(high_impact, key=max_af, reverse=True)[:10]:
            af = v["af"]
            afs = []
            for fi in focal_idx:
                if fi is not None and af and fi < len(af):
                    afs.append(f"{af[fi]:.3f}")
                else:
                    afs.append(" N/A")
            vid_short = v["vid"][-25:] if len(v["vid"]) > 25 else v["vid"]
            print(f"  {vid_short:<30} {(v['impact'] or ''):<10} "
                  f"{(v['consequence'] or '')[:24]:<25}  " + "  ".join(afs))

    # ── Write outputs ─────────────────────────────────────────────────────────
    with open(OUT_TSV, "w") as f:
        f.write("pair\tfst_locus\tfst_genome\tfst_ratio\tn_variants\n")
        for pk, vals in sorted(locus_fst.items()):
            f.write(f"{pk}\t{vals['fst_locus']}\t{vals['fst_genome']}\t"
                    f"{vals['fst_ratio']}\t{vals['n_variants']}\n")

    with open(OUT_JSON, "w") as f:
        json.dump({
            "target_gene": TARGET_GENE,
            "locus": {"chr": locus_chr, "start": locus_start, "end": locus_end,
                      "n_variants": len(kcne1_variants)},
            "h12_threshold": H12_THRESH,
            "pops_with_hard_sweep": sorted(pops_with_sweep),
            "groups_with_hard_sweep": sorted(groups_with_sweep),
            "n_groups": len(groups_with_sweep),
            "sweep_windows_neo4j": kcne1_neo4j_sweeps,
            "sweep_windows_annotate": {
                p: sws for p, sws in kcne1_sweeps.items()
            },
            "locus_fst": locus_fst,
            "n_high_impact_variants": len(high_impact),
        }, f, indent=2)

    print(f"\n=== Conclusion ===")
    low_ratio = [v for v in locus_fst.values() if v["fst_ratio"] is not None and v["fst_ratio"] < 0.5]
    high_ratio = [v for v in locus_fst.values() if v["fst_ratio"] is not None and v["fst_ratio"] > 1.5]

    if len(low_ratio) >= len(FOCAL_PAIRS) // 2:
        print(f"  ANCIENT pre-OoA sweep: FST at {TARGET_GENE} is BELOW genome-wide average")
        print(f"  → The sweep predates the Out-of-Africa dispersal (~60-70 kya)")
        print(f"  → The haplotype was fixed BEFORE population differentiation")
    elif len(high_ratio) >= len(FOCAL_PAIRS) // 2:
        print(f"  CONVERGENT sweep: FST at {TARGET_GENE} is ABOVE genome-wide average")
        print(f"  → Independent sweep events in multiple continental populations")
    else:
        print(f"  MIXED signal: FST pattern consistent with moderate differentiation")
        print(f"  → May be intermediate-age sweep or complex selective history")
    print(f"\n  Groups with h12 sweep signal: {sorted(groups_with_sweep)}")
    print(f"\nOutput: {OUT_TSV}")
    print(f"        {OUT_JSON}")


if __name__ == "__main__":
    main()
