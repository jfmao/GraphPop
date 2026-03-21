#!/usr/bin/env python3
"""Cache Neo4j data for Analysis 2: ROH × Convergent Sweep Genes.

Queries the graph for:
  1. All variants in each convergent sweep gene (chr, pos, iHS values per pop)
  2. High-H12 GenomicWindow nodes overlapping each gene (sweep region boundaries)
  3. Per-individual genotype at sweep peak variants (from hap_cache)
  4. Correlation of individual FROH with HOM_ALT dosage at sweep peaks

Since Gene nodes have no chr/start/end, variant positions are obtained via
HAS_CONSEQUENCE edges. The hap_cache provides per-individual genotypes.

Cache saved to: data/results/sweep_gene_cache.json

Usage:
  conda run -n graphevo python -u scripts/cache_sweep_gene_data.py
"""

import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

RESULTS   = Path("data/results")
HAP_DIR   = Path("data/hap_cache")
VAR_DIR   = Path("data/variants_cache")
OUT_PATH  = RESULTS / "sweep_gene_cache.json"
CHUNK     = 50_000

# IHS populations stored on Variant nodes
IHS_POPS = ["CEU", "FIN", "GBR", "IBS", "TSI",
            "CHB", "JPT", "CDX", "KHV",
            "BEB", "GIH", "ITU", "PJL", "STU",
            "CLM", "MXL", "PEL", "PUR"]

POP_GROUP = {
    "ACB":"African","ASW":"African","ESN":"African","GWD":"African",
    "LWK":"African","MSL":"African","YRI":"African",
    "CEU":"European","FIN":"European","GBR":"European","IBS":"European","TSI":"European",
    "BEB":"South_Asian","GIH":"South_Asian","ITU":"South_Asian",
    "PJL":"South_Asian","STU":"South_Asian",
    "CDX":"East_Asian","CHB":"East_Asian","CHS":"East_Asian",
    "JPT":"East_Asian","KHV":"East_Asian",
    "CLM":"American","MXL":"American","PEL":"American","PUR":"American",
}


def main():
    t0 = time.time()
    print("=== Cache: Sweep Gene Data for Analysis 2 ===\n")

    # Load convergent sweep genes
    sw_data = json.load(open(RESULTS / "convergent_sweeps.json"))
    conv_genes = sw_data["convergent_genes"]
    gene_names = list(conv_genes.keys())
    print(f"Convergent sweep genes: {gene_names}\n", flush=True)

    # Load sample metadata
    meta = json.load(open(HAP_DIR / "sample_meta.json"))
    sample_ids = meta["sample_ids"]
    pop_ids    = meta["pop_ids"]
    n_samples  = len(sample_ids)

    # Load FROH per sample from Inv.16b
    import csv
    ind_rows = {r["sampleId"]: float(r["froh_genome"])
                for r in csv.DictReader(
                    open(RESULTS / "individual_trajectory.tsv"), delimiter="\t")}
    froh_arr = np.array([ind_rows.get(s, 0.0) for s in sample_ids])

    from neo4j import GraphDatabase
    driver = GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "graphpop"))

    # ── Step 1: Query variant positions in each convergent sweep gene ─────────
    print("Step 1: Querying variant positions in convergent sweep genes...", flush=True)
    ihs_props = " ".join([f"v.ihs_{p} AS ihs_{p}," for p in IHS_POPS]).rstrip(",")
    gene_variants = {}   # gene_name → list of {chr, pos, ihs_*}

    with driver.session(database="neo4j") as s:
        for gene in gene_names:
            result = s.run(
                f"MATCH (g:Gene)<-[c:HAS_CONSEQUENCE]-(v:Variant) "
                f"WHERE g.symbol = $gene "
                f"RETURN v.chr AS chr, v.pos AS pos, {ihs_props}",
                gene=gene,
            ).data()
            variants = []
            for row in result:
                entry = {"chr": row["chr"], "pos": int(row["pos"])}
                for p in IHS_POPS:
                    v = row.get(f"ihs_{p}")
                    if v is not None:
                        entry[f"ihs_{p}"] = round(float(v), 4)
                variants.append(entry)
            gene_variants[gene] = variants
            chrom = variants[0]["chr"] if variants else "?"
            print(f"  {gene}: {len(variants)} variants on {chrom}", flush=True)

    # ── Step 2: Query H12 windows overlapping each gene's region ─────────────
    print("\nStep 2: Querying H12 windows overlapping sweep genes...", flush=True)
    gene_windows = {}

    with driver.session(database="neo4j") as s:
        for gene, variants in gene_variants.items():
            if not variants:
                gene_windows[gene] = []
                continue
            chrom   = variants[0]["chr"]
            min_pos = min(v["pos"] for v in variants)
            max_pos = max(v["pos"] for v in variants)
            # Look ±500kb around gene
            result = s.run(
                "MATCH (w:GenomicWindow) "
                "WHERE w.chr = $chr "
                "  AND w.start <= $max_pos + 500000 "
                "  AND w.end   >= $min_pos - 500000 "
                "  AND w.h12 IS NOT NULL "
                "RETURN w.population AS pop, w.chr AS chr, "
                "       w.start AS start, w.end AS end, "
                "       w.h12 AS h12, w.sweep_type AS sweep_type "
                "ORDER BY w.h12 DESC LIMIT 20",
                chr=chrom, min_pos=min_pos, max_pos=max_pos,
            ).data()
            windows = [{
                "pop": r["pop"], "chr": r["chr"],
                "start": int(r["start"]), "end": int(r["end"]),
                "h12": round(float(r["h12"]), 4),
                "sweep_type": r.get("sweep_type", ""),
            } for r in result]
            gene_windows[gene] = windows
            print(f"  {gene}: {len(windows)} nearby H12 windows "
                  f"(max_h12={max((w['h12'] for w in windows), default=0):.3f})", flush=True)

    driver.close()

    # ── Step 3: Per-individual genotype at sweep peak variants ────────────────
    print("\nStep 3: Computing per-individual genotype at sweep peak variants...", flush=True)
    gene_ind_stats = {}

    for gene, variants in gene_variants.items():
        if not variants:
            continue

        chrom = variants[0]["chr"]
        hap_path = HAP_DIR / f"{chrom}.npz"
        var_path  = VAR_DIR  / f"{chrom}.npz"
        if not hap_path.exists():
            print(f"  {gene}: hap_cache missing for {chrom}, skipping", flush=True)
            continue

        hap_data = np.load(hap_path)
        var_data  = np.load(var_path)
        pos_all   = var_data["pos"]
        vtype     = var_data["variant_type"]
        ac_all    = var_data["ac"]
        an_all    = var_data["an"]
        hap_packed = hap_data["hap"]

        # Focus on SNPs only
        snp_mask  = vtype == 1
        pos_snps  = pos_all[snp_mask]
        hap_snp   = hap_packed[snp_mask]
        ac_snps   = ac_all[snp_mask].sum(axis=1)
        an_snps   = an_all[snp_mask].sum(axis=1).clip(1)
        global_af = ac_snps / an_snps

        # Target positions: all SNPs in this gene
        target_positions = set(v["pos"] for v in variants)
        target_mask_idx  = np.where(snp_mask)[0]   # indices in full array
        # Match snp positions to gene positions
        gene_snp_idx = np.array([i for i, p in enumerate(pos_snps)
                                  if int(p) in target_positions], dtype=np.int64)

        if len(gene_snp_idx) == 0:
            print(f"  {gene}: no SNPs found in hap_cache, skipping", flush=True)
            continue

        # For each gene SNP, compute per-individual dosage (0/1/2)
        # Use chunked approach for memory safety
        n_hap = 2 * n_samples
        dosage_matrix = np.zeros((len(gene_snp_idx), n_samples), dtype=np.int8)

        for chunk_i, si in enumerate(gene_snp_idx):
            chunk = np.unpackbits(hap_snp[si:si+1], axis=1, bitorder="big")[:, :n_hap]
            dosage_matrix[chunk_i] = (chunk[0, 0::2] + chunk[0, 1::2]).astype(np.int8)

        # Per-variant stats
        variant_stats = []
        for i, si in enumerate(gene_snp_idx):
            pos_val = int(pos_snps[si])
            af_val  = float(global_af[si])
            dos     = dosage_matrix[i]
            # Frequency of HOM_ALT (dosage=2) per individual
            hom_alt_rate = float((dos == 2).mean())
            het_rate     = float((dos == 1).mean())
            # Correlation of individual FROH with dosage
            valid = froh_arr > 0
            from scipy.stats import spearmanr
            rho, p = spearmanr(froh_arr[valid], dos[valid].astype(float))
            variant_stats.append({
                "pos":          pos_val,
                "global_af":    round(af_val, 5),
                "hom_alt_rate": round(hom_alt_rate, 4),
                "het_rate":     round(het_rate, 4),
                "froh_rho":     round(float(rho), 4),
                "froh_p":       round(float(p), 6),
            })

        # Find peak variant (highest max_h12 based on iHS or AF near 0.5)
        # Identify sweep peak = variant with largest sum of |iHS| across pops
        peak_ihs_sums = []
        pos_to_variant = {v["pos"]: v for v in variants}
        for si in gene_snp_idx:
            pos_val = int(pos_snps[si])
            vinfo   = pos_to_variant.get(pos_val, {})
            ihs_sum = sum(abs(vinfo.get(f"ihs_{p}", 0) or 0) for p in IHS_POPS)
            peak_ihs_sums.append(ihs_sum)

        peak_idx = int(np.argmax(peak_ihs_sums)) if peak_ihs_sums else 0
        peak_dos = dosage_matrix[peak_idx]

        # Group-level HOM_ALT rates at peak variant
        group_homalt = {}
        for grp in ["African","European","East_Asian","South_Asian","American"]:
            idx = [i for i, p in enumerate(pop_ids) if POP_GROUP.get(p) == grp]
            if idx:
                group_homalt[grp] = round(float((peak_dos[idx] == 2).mean()), 4)

        # FROH quartile analysis for peak variant
        froh_q = np.percentile(froh_arr[froh_arr > 0], [25, 50, 75])
        q_homalt = {}
        for qi in range(1, 5):
            lo = froh_q[qi-2] if qi > 1 else 0
            hi = froh_q[qi-1] if qi < 4 else 1.0
            if qi == 4:
                mask = froh_arr >= froh_q[2]
            elif qi == 1:
                mask = (froh_arr > 0) & (froh_arr < froh_q[0])
            else:
                mask = (froh_arr >= froh_q[qi-2]) & (froh_arr < froh_q[qi-1])
            if mask.sum() > 0:
                q_homalt[f"Q{qi}"] = round(float((peak_dos[mask] == 2).mean()), 4)

        gene_ind_stats[gene] = {
            "chr":            chrom,
            "n_snps":         len(gene_snp_idx),
            "peak_pos":       int(pos_snps[gene_snp_idx[peak_idx]]),
            "peak_af":        round(float(global_af[gene_snp_idx[peak_idx]]), 5),
            "peak_ihs_sum":   round(float(peak_ihs_sums[peak_idx]), 3),
            "group_homalt_at_peak": group_homalt,
            "froh_quartile_homalt": q_homalt,
            "variant_stats":  variant_stats,
        }
        print(f"  {gene} ({chrom}): {len(gene_snp_idx)} SNPs  "
              f"peak_pos={int(pos_snps[gene_snp_idx[peak_idx]])}  "
              f"HOM_ALT by group: " +
              ", ".join(f"{g[:3]}={v:.3f}" for g,v in group_homalt.items()),
              flush=True)

    # ── Write cache ───────────────────────────────────────────────────────────
    out = {
        "convergent_genes_meta": {
            g: {"groups": info["groups"], "max_h12": info["max_h12"],
                "n_groups": info["n_groups"]}
            for g, info in conv_genes.items()
        },
        "gene_variants":   gene_variants,
        "gene_windows":    gene_windows,
        "gene_ind_stats":  gene_ind_stats,
        "n_samples":       n_samples,
        "ihs_populations": IHS_POPS,
    }

    with open(OUT_PATH, "w") as f:
        json.dump(out, f, indent=2)

    size_mb = OUT_PATH.stat().st_size / 1e6
    print(f"\nCache saved: {OUT_PATH}  ({size_mb:.1f} MB)")
    print(f"Total time: {time.time()-t0:.0f}s")

    # Print summary
    print("\n=== Summary ===")
    for gene, stats in gene_ind_stats.items():
        info = conv_genes[gene]
        print(f"\n{gene}  ({stats['chr']}, {len(info['groups'])} groups, "
              f"max_H12={info['max_h12']}):")
        print(f"  Peak pos: {stats['peak_pos']}  AF={stats['peak_af']:.3f}")
        print(f"  HOM_ALT at peak: " +
              ", ".join(f"{g}={v:.3f}" for g,v in stats['group_homalt_at_peak'].items()))
        print(f"  FROH quartile HOM_ALT: " +
              ", ".join(f"{q}={v:.3f}" for q,v in stats['froh_quartile_homalt'].items()))


if __name__ == "__main__":
    main()
