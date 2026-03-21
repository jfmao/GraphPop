#!/usr/bin/env python3
"""Investigation 15: Sweep Spatial Extent via Positional Traversal.

Question: What is the genomic extent (bp) of confirmed convergent sweep haplotypes,
and does extent correlate with sweep strength (h12), population size (FROH), or
temporal age (from Inv.10)?

Method:
  For each convergent sweep gene (9 genes from convergent_sweeps.json):
    1. Find gene chromosome and center position via HAS_CONSEQUENCE traversal in Neo4j.
    2. Find the sweep population pair (focal pair with highest FST/sweep signal).
    3. Load variant_cache numpy arrays for the chromosome (pre-computed, fast).
    4. Compute Hudson FST in 10 kb bins radiating outward from the gene center.
    5. Define sweep extent: first bin where FST < (peak_FST × 0.5) OR
       FST < (genome_mean_FST + 2σ).
    6. Correlate sweep extent with h12 peak value (sweep strength) and
       FROH of sweep population.

Why graph-native (highest novelty):
  Uses Neo4j HAS_CONSEQUENCE graph traversal to anchor the sweep locus, then
  exploits the pre-built variant position index (variant_cache) for spatial
  FST decay computation. The gene→chromosome→position linkage is purely
  graph-native; no VCF scanning needed to locate the gene.

Output: data/results/sweep_extent.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_15_sweep_extent.py
"""

import csv
import json
import time
from pathlib import Path
from collections import defaultdict

import numpy as np
from neo4j import GraphDatabase
from scipy.stats import spearmanr

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB   = "neo4j"

OUT_DIR   = Path("data/results")
OUT_TSV   = OUT_DIR / "sweep_extent.tsv"
OUT_JSON  = OUT_DIR / "sweep_extent.json"

CACHE_DIR = Path("data/variants_cache")

BIN_SIZE    = 10_000   # 10 kb bins for FST decay
MAX_DIST    = 1_000_000  # 1 Mb search radius each direction
HALF_DECAY  = 0.5       # sweep boundary: FST < peak × 0.5
MIN_AN      = 10        # min allele number per population

# Population order in variant cache (from pop_meta.json)
POP_IDS = ['ACB','ASW','BEB','CDX','CEU','CHB','CHS','CLM','ESN','FIN',
           'GBR','GIH','GWD','IBS','ITU','JPT','KHV','LWK','MSL','MXL',
           'PEL','PJL','PUR','STU','TSI','YRI']
IDX = {p: i for i, p in enumerate(POP_IDS)}

# Population group assignments for sweep pairs
GROUP_POPS = {
    "African":    ["YRI", "LWK", "MSL", "ESN", "GWD", "ACB", "ASW"],
    "European":   ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "East_Asian": ["CHB", "JPT", "CHS", "CDX", "KHV"],
    "South_Asian":["GIH", "ITU", "PJL", "BEB", "STU"],
    "American":   ["CLM", "MXL", "PEL", "PUR"],
}

# Default focal comparison: sweep pop vs African (YRI)
DEFAULT_PAIR = ("YRI",)  # will be extended with sweep pop per gene


# ── Hudson FST (scalar) ───────────────────────────────────────────────────────
def hudson_fst_arrays(ac1, an1, ac2, an2):
    """Compute Hudson FST for arrays of variants. Returns (numerator_sum, denom_sum)."""
    valid = (an1 >= 2) & (an2 >= 2) & (an1 > 0) & (an2 > 0)
    if valid.sum() == 0:
        return 0.0, 0.0, 0
    p1 = ac1[valid].astype(float) / an1[valid]
    p2 = ac2[valid].astype(float) / an2[valid]
    n1 = an1[valid].astype(float)
    n2 = an2[valid].astype(float)
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    valid2 = den > 0
    return float(num[valid2].sum()), float(den[valid2].sum()), int(valid2.sum())


# ── Sweep extent calculation ──────────────────────────────────────────────────
def compute_sweep_extent(pos, ac, an, gene_center, pop1_idx, pop2_idx,
                          bin_size=BIN_SIZE, max_dist=MAX_DIST, half_decay=HALF_DECAY):
    """Compute FST decay in bins radiating from gene_center.

    Returns dict with:
      bins: list of {dist_center, fst, n_variants}
      peak_fst: FST at gene center bin
      genome_mean_fst: mean FST across all bins
      extent_left_bp: distance (bp) to left extent boundary
      extent_right_bp: distance to right extent boundary
      total_extent_bp: total extent
    """
    ac1 = ac[:, pop1_idx]
    ac2 = ac[:, pop2_idx]
    an1 = an[:, pop1_idx]
    an2 = an[:, pop2_idx]

    gene_lo = gene_center - max_dist
    gene_hi = gene_center + max_dist
    mask = (pos >= gene_lo) & (pos <= gene_hi)
    pos_w  = pos[mask]
    ac1_w  = ac1[mask]
    ac2_w  = ac2[mask]
    an1_w  = an1[mask]
    an2_w  = an2[mask]

    if len(pos_w) == 0:
        return None

    # Bin boundaries from gene_center outward
    n_bins_each = max_dist // bin_size
    bin_edges   = np.arange(-n_bins_each, n_bins_each + 1) * bin_size

    bins = []
    for i in range(len(bin_edges) - 1):
        lo = gene_center + bin_edges[i]
        hi = gene_center + bin_edges[i + 1]
        b_mask = (pos_w >= lo) & (pos_w < hi)
        n_var  = int(b_mask.sum())
        dist_center = int((bin_edges[i] + bin_edges[i + 1]) // 2)

        if n_var < 2:
            bins.append({"dist": dist_center, "fst": np.nan, "n_var": n_var})
            continue

        num, den, n_valid = hudson_fst_arrays(
            ac1_w[b_mask], an1_w[b_mask], ac2_w[b_mask], an2_w[b_mask]
        )
        fst = float(num / den) if den > 0 else np.nan
        bins.append({"dist": dist_center, "fst": fst, "n_var": n_var})

    # FST at gene center (dist ∈ [-bin_size/2, +bin_size/2])
    center_bins = [b for b in bins if abs(b["dist"]) <= bin_size]
    fst_vals    = [b["fst"] for b in bins if not np.isnan(b["fst"])]

    if not fst_vals:
        return None

    genome_mean_fst = float(np.nanmean(fst_vals))
    genome_std_fst  = float(np.nanstd(fst_vals))

    if center_bins:
        center_fst_vals = [b["fst"] for b in center_bins if not np.isnan(b["fst"])]
        peak_fst = float(np.nanmax(center_fst_vals)) if center_fst_vals else genome_mean_fst
    else:
        # Fall back to max FST near gene
        close_bins = [b for b in bins if abs(b["dist"]) <= 50_000]
        close_fsts = [b["fst"] for b in close_bins if not np.isnan(b["fst"])]
        peak_fst = float(np.nanmax(close_fsts)) if close_fsts else genome_mean_fst

    threshold = max(peak_fst * half_decay, genome_mean_fst + 2 * genome_std_fst)

    # Determine extent: walk outward until FST drops below threshold
    # Left extent (negative distances)
    left_bins  = sorted([b for b in bins if b["dist"] < 0], key=lambda x: x["dist"], reverse=True)
    right_bins = sorted([b for b in bins if b["dist"] >= 0], key=lambda x: x["dist"])

    def find_extent_bp(ordered_bins):
        for b in ordered_bins:
            if not np.isnan(b["fst"]) and b["fst"] < threshold:
                return abs(b["dist"])
        return max_dist  # reached boundary

    left_extent  = find_extent_bp(left_bins)
    right_extent = find_extent_bp(right_bins)
    total_extent = left_extent + right_extent

    return {
        "bins":            bins,
        "peak_fst":        round(peak_fst, 4),
        "genome_mean_fst": round(genome_mean_fst, 4),
        "genome_std_fst":  round(genome_std_fst, 4),
        "threshold":       round(threshold, 4),
        "extent_left_bp":  int(left_extent),
        "extent_right_bp": int(right_extent),
        "total_extent_bp": int(total_extent),
        "total_extent_kb": round(total_extent / 1000, 1),
    }


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print("=== Investigation 15: Sweep Spatial Extent ===\n", flush=True)

    # ── Load convergent sweeps data ────────────────────────────────────────────
    conv_path = OUT_DIR / "convergent_sweeps.json"
    if not conv_path.exists():
        raise FileNotFoundError(f"Missing {conv_path} — run investigate_02 first")

    with open(conv_path) as f:
        conv_data = json.load(f)

    conv_genes = conv_data.get("convergent_genes", {})
    if isinstance(conv_genes, list):
        # convert list format
        conv_genes = {g.get("gene_id", g): g for g in conv_genes}

    print(f"Convergent sweep genes: {list(conv_genes.keys())}", flush=True)

    # ── Load FROH and temporal data ───────────────────────────────────────────
    froh_by_pop = {}
    h12_by_pop  = {}
    temporal_path = OUT_DIR / "temporal_selection.json"
    if temporal_path.exists():
        with open(temporal_path) as f:
            temp_data = json.load(f)
        for pop_data in temp_data.get("populations", []):
            pop = pop_data.get("pop")
            if pop:
                froh_by_pop[pop] = pop_data.get("froh", np.nan)
                h12_by_pop[pop]  = pop_data.get("h12_fraction", np.nan)

    # ── Neo4j: find gene locations ─────────────────────────────────────────────
    print("Querying Neo4j for gene locations...", flush=True)
    gene_loci = {}  # gene_id → {chr, min_pos, max_pos, center}

    CYPHER = """
    MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
    WHERE g.geneId IN $gene_ids
    RETURN g.geneId AS gene_id, v.chr AS chr,
           min(v.pos) AS min_pos, max(v.pos) AS max_pos,
           count(v) AS n_variants
    """
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    with driver.session(database=NEO4J_DB) as session:
        result = session.run(CYPHER, gene_ids=list(conv_genes.keys()))
        for rec in result:
            gid = rec["gene_id"]
            if gid:
                gene_loci[gid] = {
                    "chr":      rec["chr"],
                    "min_pos":  rec["min_pos"],
                    "max_pos":  rec["max_pos"],
                    "center":   (rec["min_pos"] + rec["max_pos"]) // 2,
                    "n_variants": rec["n_variants"],
                }
    driver.close()

    print(f"  Found loci for {len(gene_loci)} of {len(conv_genes)} genes", flush=True)
    for gid, locus in sorted(gene_loci.items()):
        print(f"    {gid:15s}: {locus['chr']}:{locus['min_pos']:,}-{locus['max_pos']:,} "
              f"({locus['n_variants']} coding variants)")

    # ── Load variant cache and compute FST decay ───────────────────────────────
    results = []
    cache_loaded = {}  # chr → cache data

    for gene_id, gene_data in conv_genes.items():
        locus = gene_loci.get(gene_id)
        if not locus:
            print(f"\n  [skip] {gene_id}: no locus found in Neo4j", flush=True)
            continue

        chrom   = locus["chr"]
        center  = locus["center"]
        pops    = gene_data.get("pops", [])  # sweep populations

        print(f"\nProcessing {gene_id} ({chrom}:{center:,})...", flush=True)

        # Load variant cache for this chromosome
        if chrom not in cache_loaded:
            cache_path = CACHE_DIR / f"{chrom}.npz"
            if not cache_path.exists():
                print(f"  [skip] no variant cache for {chrom}", flush=True)
                continue
            print(f"  Loading {chrom} variant cache...", flush=True)
            c = np.load(cache_path)
            # Filter to SNPs only (variant_type == 1)
            snp_mask = c["variant_type"] == 1
            cache_loaded[chrom] = {
                "pos": c["pos"][snp_mask].astype(np.int64),
                "ac":  c["ac"][snp_mask].astype(np.int32),
                "an":  c["an"][snp_mask].astype(np.int32),
            }
            print(f"  Loaded {snp_mask.sum():,} SNPs", flush=True)

        cache = cache_loaded[chrom]
        pos = cache["pos"]
        ac  = cache["ac"]
        an  = cache["an"]

        # Choose focal population pair
        # For sweep: use the non-African sweep pop vs YRI
        # If multiple pops, use the one with highest h12
        sweep_pop = None
        for p in sorted(pops, key=lambda x: h12_by_pop.get(x, 0), reverse=True):
            if p in IDX and p != "YRI":
                sweep_pop = p
                break
        if sweep_pop is None:
            # Fall back to CEU if no non-African pop
            sweep_pop = pops[0] if pops else "CEU"

        # Compare sweep_pop vs YRI (or vs ESN if sweep_pop is African)
        if sweep_pop in GROUP_POPS.get("African", []):
            ref_pop = "CEU"
        else:
            ref_pop = "YRI"

        if sweep_pop not in IDX or ref_pop not in IDX:
            print(f"  [skip] pop {sweep_pop} or {ref_pop} not in POP_IDS", flush=True)
            continue

        pop1_idx = IDX[ref_pop]
        pop2_idx = IDX[sweep_pop]

        print(f"  Sweep pair: {ref_pop} vs {sweep_pop}", flush=True)

        # Compute FST decay
        ext = compute_sweep_extent(pos, ac, an, center, pop1_idx, pop2_idx)
        if ext is None:
            print(f"  [skip] insufficient variants near {gene_id}", flush=True)
            continue

        max_h12 = gene_data.get("max_h12", np.nan)
        n_groups = gene_data.get("n_groups", len(gene_data.get("groups", [])))

        result_row = {
            "gene_id":        gene_id,
            "chr":            chrom,
            "center_pos":     int(center),
            "sweep_pop":      sweep_pop,
            "ref_pop":        ref_pop,
            "max_h12":        max_h12,
            "n_conv_groups":  n_groups,
            "froh_sweep_pop": round(froh_by_pop.get(sweep_pop, np.nan), 4),
            "peak_fst":       ext["peak_fst"],
            "genome_mean_fst": ext["genome_mean_fst"],
            "fst_threshold":  ext["threshold"],
            "extent_left_kb": round(ext["extent_left_bp"] / 1000, 1),
            "extent_right_kb": round(ext["extent_right_bp"] / 1000, 1),
            "total_extent_kb": ext["total_extent_kb"],
            "n_bins":         len(ext["bins"]),
        }
        results.append({**result_row, "decay_bins": ext["bins"]})

        print(f"  Peak FST: {ext['peak_fst']:.4f}  "
              f"Extent: {ext['total_extent_kb']:.0f} kb "
              f"(L={ext['extent_left_bp']//1000} kb, R={ext['extent_right_bp']//1000} kb)")

    # ── Correlations ─────────────────────────────────────────────────────────
    if len(results) >= 3:
        extents  = [r["total_extent_kb"] for r in results]
        h12_vals = [r["max_h12"] for r in results]
        froh_vals = [r["froh_sweep_pop"] for r in results]
        peak_fsts = [r["peak_fst"] for r in results]

        print("\n=== Correlations ===", flush=True)
        if all(not np.isnan(v) for v in h12_vals):
            r_h12, p_h12 = spearmanr(h12_vals, extents)
            print(f"  Extent vs h12:      ρ={r_h12:.3f}  p={p_h12:.3f}")
        if all(not np.isnan(v) for v in froh_vals):
            r_froh, p_froh = spearmanr(froh_vals, extents)
            print(f"  Extent vs FROH:     ρ={r_froh:.3f}  p={p_froh:.3f}")
        if all(not np.isnan(v) for v in peak_fsts):
            r_fst, p_fst = spearmanr(peak_fsts, extents)
            print(f"  Extent vs peak_FST: ρ={r_fst:.3f}  p={p_fst:.3f}")

    # ── Write outputs ─────────────────────────────────────────────────────────
    tsv_fields = ["gene_id", "chr", "center_pos", "sweep_pop", "ref_pop",
                  "max_h12", "n_conv_groups", "froh_sweep_pop",
                  "peak_fst", "genome_mean_fst", "fst_threshold",
                  "extent_left_kb", "extent_right_kb", "total_extent_kb", "n_bins"]

    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=tsv_fields, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    # Save full decay curves in JSON
    out_data = {
        "config": {
            "bin_size_bp":   BIN_SIZE,
            "max_dist_bp":   MAX_DIST,
            "half_decay":    HALF_DECAY,
        },
        "genes": [
            {k: v for k, v in r.items() if k != "decay_bins"}
            for r in results
        ],
        "decay_curves": {
            r["gene_id"]: r.get("decay_bins", [])
            for r in results
        },
    }

    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2, default=lambda x: None if np.isnan(x) else x)

    print(f"\n=== Summary ===")
    print(f"  Genes analyzed: {len(results)}")
    if results:
        extents_kb = [r["total_extent_kb"] for r in results]
        print(f"  Mean extent: {np.mean(extents_kb):.0f} kb")
        print(f"  Min extent:  {np.min(extents_kb):.0f} kb ({results[np.argmin(extents_kb)]['gene_id']})")
        print(f"  Max extent:  {np.max(extents_kb):.0f} kb ({results[np.argmax(extents_kb)]['gene_id']})")

    print(f"\nOutputs: {OUT_TSV}  {OUT_JSON}")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
