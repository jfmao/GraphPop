#!/usr/bin/env python3
"""Investigation 14: Population-Specific Rare Variant Burden.

Question: Which genes carry the highest burden of population-specific rare variants
(AF < 0.5% in one population, absent elsewhere)? Reveals genes accumulating private
deleterious mutations not yet purged by selection — proxy for bottlenecks and
relaxed functional constraint.

Method:
  Stream Variant→Gene (HAS_CONSEQUENCE) with allele count arrays.
  For each variant: identify "private rare" alleles:
    - ac[pop] > 0
    - ac[pop] / an[pop] < 0.005 (AF < 0.5%)
    - zero copies in ALL other populations
  Aggregate per-gene per-population count of private rare variants.
  Correlate gene rare burden with FROH from temporal_selection.json.
  Contrast burden in HIGH-impact vs LOW-impact consequences per gene.

Why graph-native:
  Population-stratified allele count arrays (ac[], an[], pop_ids[]) on each Variant
  node combined with HAS_CONSEQUENCE→Gene traversal in one streaming query.
  Standard tools require separate VCF subsets per population + external gene
  annotation + custom aggregation.

Output: data/results/rare_burden.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_14_rare_burden.py
"""

import csv
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase
from scipy.stats import spearmanr

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB   = "neo4j"

OUT_DIR  = Path("data/results")
OUT_TSV  = OUT_DIR / "rare_burden.tsv"
OUT_JSON = OUT_DIR / "rare_burden.json"

RARE_AC_MAX       = 4      # ≤ 4 copies in focal population (rare by absolute count)
RARE_AF_THRESHOLD = 0.02   # < 2% in focal population (relaxed for 1000G min AC=2)
MIN_AN            = 10     # minimum allele number per population to compute AF
TOP_N_GENES       = 50     # top genes to report in detail


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print("=== Investigation 14: Population-Specific Rare Variant Burden ===\n", flush=True)

    # ── Load FROH data for correlation ────────────────────────────────────────
    froh_by_pop = {}  # pop → FROH value
    froh_z_by_pop = {}
    temporal_path = OUT_DIR / "temporal_selection.json"
    if temporal_path.exists():
        with open(temporal_path) as f:
            temp_data = json.load(f)
        for pop_data in temp_data.get("populations", []):
            pop = pop_data.get("pop")
            froh = pop_data.get("froh")
            froh_z = pop_data.get("froh_z")
            if pop and froh is not None:
                froh_by_pop[pop] = froh
                froh_z_by_pop[pop] = froh_z if froh_z is not None else 0.0
    print(f"Loaded FROH for {len(froh_by_pop)} populations", flush=True)
    print(f"  High-FROH pops: " + ", ".join(
        f"{p}(z={froh_z_by_pop.get(p, 0):.2f})"
        for p in sorted(froh_z_by_pop, key=froh_z_by_pop.get, reverse=True)[:5]
    ), flush=True)

    # ── Stream Variant→Gene edges ─────────────────────────────────────────────
    # Accumulators: gene_id → pop → {n_private_rare, n_total, n_high_impact_rare, n_low_impact_rare}
    gene_pop_counts = defaultdict(
        lambda: defaultdict(lambda: {
            "n_private_rare": 0,
            "n_total": 0,
            "n_high_rare": 0,
            "n_low_rare": 0,
        })
    )
    gene_symbol_map = {}  # gene_id → symbol

    CYPHER = """
    MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)
    WHERE v.ac IS NOT NULL AND v.an IS NOT NULL
    RETURN g.geneId   AS gene_id,
           g.symbol   AS symbol,
           c.impact   AS impact,
           v.ac       AS ac,
           v.an       AS an,
           v.pop_ids  AS pop_ids
    """

    print("Streaming Variant→Gene consequences from Neo4j...", flush=True)
    n_records = 0
    n_private = 0

    with GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS)) as driver:
        with driver.session(database=NEO4J_DB) as session:
            result = session.run(CYPHER)
            for rec in result:
                n_records += 1
                if n_records % 100_000 == 0:
                    elapsed = time.time() - t0
                    print(f"  {n_records:,} records | {n_private:,} private-rare found | "
                          f"{elapsed:.0f}s", flush=True)

                gene_id  = rec["gene_id"]
                symbol   = rec["symbol"] or gene_id
                impact   = rec["impact"] or "MODIFIER"
                ac_raw   = rec["ac"]
                an_raw   = rec["an"]
                pop_ids  = list(rec["pop_ids"]) if rec["pop_ids"] else []

                if not ac_raw or not an_raw or not pop_ids:
                    continue

                ac = np.asarray(ac_raw, dtype=np.int32)
                an = np.asarray(an_raw, dtype=np.int32)
                n_pops = len(pop_ids)

                gene_symbol_map[gene_id] = symbol

                # For each population, check if this is a private rare variant
                for focal_idx, focal_pop in enumerate(pop_ids):
                    if focal_idx >= len(ac) or focal_idx >= len(an):
                        continue
                    ac_focal = ac[focal_idx]
                    an_focal = an[focal_idx]

                    gene_pop_counts[gene_id][focal_pop]["n_total"] += 1

                    if an_focal < MIN_AN or ac_focal == 0:
                        continue

                    if ac_focal > RARE_AC_MAX:
                        continue

                    af_focal = ac_focal / an_focal
                    if af_focal >= RARE_AF_THRESHOLD:
                        continue

                    # Check if absent in ALL other populations
                    other_ac = np.sum(ac) - ac_focal
                    if other_ac > 0:
                        continue

                    # This is a private rare variant in focal_pop
                    n_private += 1
                    gene_pop_counts[gene_id][focal_pop]["n_private_rare"] += 1
                    if impact == "HIGH":
                        gene_pop_counts[gene_id][focal_pop]["n_high_rare"] += 1
                    elif impact in ("LOW", "MODIFIER"):
                        gene_pop_counts[gene_id][focal_pop]["n_low_rare"] += 1

    elapsed = time.time() - t0
    print(f"\nStreaming complete: {n_records:,} records, {n_private:,} private-rare events "
          f"in {elapsed:.1f}s", flush=True)

    # ── Aggregate per-gene total burden across all populations ────────────────
    # For each gene: sum private rare variants across all populations
    gene_total_rare = {}
    for gene_id, pop_dict in gene_pop_counts.items():
        total = sum(v["n_private_rare"] for v in pop_dict.values())
        gene_total_rare[gene_id] = total

    n_genes_with_rare = sum(1 for v in gene_total_rare.values() if v > 0)
    print(f"  Genes with ≥1 private rare variant: {n_genes_with_rare:,}", flush=True)

    # ── Build flat TSV rows ────────────────────────────────────────────────────
    rows = []
    for gene_id, pop_dict in gene_pop_counts.items():
        symbol = gene_symbol_map.get(gene_id, gene_id)
        for pop, counts in pop_dict.items():
            if counts["n_private_rare"] == 0:
                continue
            froh = froh_by_pop.get(pop, None)
            froh_z = froh_z_by_pop.get(pop, None)
            rows.append({
                "gene_id":       gene_id,
                "symbol":        symbol,
                "population":    pop,
                "n_private_rare": counts["n_private_rare"],
                "n_high_rare":   counts["n_high_rare"],
                "n_low_rare":    counts["n_low_rare"],
                "n_total_variants": counts["n_total"],
                "rare_burden_pct": round(
                    100 * counts["n_private_rare"] / max(counts["n_total"], 1), 2
                ),
                "froh":   round(froh, 4) if froh is not None else "",
                "froh_z": round(froh_z, 3) if froh_z is not None else "",
            })

    rows.sort(key=lambda r: r["n_private_rare"], reverse=True)
    print(f"  Total gene×pop rows with rare burden: {len(rows):,}", flush=True)

    # ── Correlation: per-pop mean rare burden vs FROH ──────────────────────────
    pop_mean_burden = defaultdict(list)
    for r in rows:
        pop_mean_burden[r["population"]].append(r["n_private_rare"])

    pops_with_froh  = [p for p in pop_mean_burden if p in froh_by_pop]
    mean_burdens    = [np.mean(pop_mean_burden[p]) for p in pops_with_froh]
    froh_vals       = [froh_by_pop[p] for p in pops_with_froh]

    correlation_r, correlation_p = None, None
    if len(pops_with_froh) >= 5:
        correlation_r, correlation_p = spearmanr(froh_vals, mean_burdens)
        print(f"\n  FROH × rare burden correlation (Spearman):")
        print(f"    ρ = {correlation_r:.4f}  p = {correlation_p:.3e}")

    # Print top-burden populations
    pop_burden_summary = sorted(
        [(p, np.mean(v), froh_by_pop.get(p, 0)) for p, v in pop_mean_burden.items()],
        key=lambda x: x[1], reverse=True
    )
    print("\n  Top populations by mean per-gene rare burden:")
    for pop, mean_b, froh in pop_burden_summary[:8]:
        print(f"    {pop}: mean_rare={mean_b:.2f}  FROH={froh:.4f}")

    # ── Write outputs ─────────────────────────────────────────────────────────
    fieldnames = ["gene_id", "symbol", "population", "n_private_rare",
                  "n_high_rare", "n_low_rare", "n_total_variants",
                  "rare_burden_pct", "froh", "froh_z"]
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Top genes by total rare burden across populations
    top_genes_by_total = sorted(gene_total_rare.items(), key=lambda x: x[1], reverse=True)[:TOP_N_GENES]

    out_data = {
        "config": {
            "rare_ac_max": RARE_AC_MAX,
            "rare_af_threshold": RARE_AF_THRESHOLD,
            "min_an": MIN_AN,
        },
        "summary": {
            "n_records_processed": n_records,
            "n_private_rare_events": n_private,
            "n_genes_with_rare": n_genes_with_rare,
            "n_output_rows": len(rows),
        },
        "froh_correlation": {
            "n_populations": len(pops_with_froh),
            "spearman_r": round(float(correlation_r), 4) if correlation_r is not None else None,
            "spearman_p": float(f"{correlation_p:.3e}") if correlation_p is not None else None,
        },
        "top_genes": [
            {
                "gene_id": gid,
                "symbol": gene_symbol_map.get(gid, gid),
                "total_private_rare": n,
            }
            for gid, n in top_genes_by_total
        ],
        "pop_summary": {
            p: {"mean_rare": round(np.mean(v), 2), "froh": froh_by_pop.get(p),
                "froh_z": froh_z_by_pop.get(p)}
            for p, v in pop_mean_burden.items()
        },
    }

    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2, default=str)

    print(f"\n  Top 10 genes by total private rare burden:")
    for gid, n in top_genes_by_total[:10]:
        sym = gene_symbol_map.get(gid, gid)
        print(f"    {sym:15s}  total_rare={n}")

    print(f"\nOutputs: {OUT_TSV}  {OUT_JSON}")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
