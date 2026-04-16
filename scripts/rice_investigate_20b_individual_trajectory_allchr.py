#!/usr/bin/env python3
"""Rice Investigation 20b: Genome-Wide Individual Evolutionary Trajectory Map (All 12 Chromosomes).

Extends Inv.20 (Chr1-only) to aggregate individual-level statistics across all
12 rice chromosomes for a true genome-wide picture.

For each chromosome, computes per-individual:
  1. het_count / poly_count  -> het_rate (averaged across chromosomes)
  2. hom_alt_count / poly_count -> hom_alt_rate (averaged across chromosomes)
  3. rare_burden   -> summed across chromosomes (count of rare alleles, AF < 0.05)
  4. private_burden -> summed across chromosomes (count of ultra-rare alleles, AF < 0.005)

Genome-wide aggregation:
  - het_rate, hom_alt_rate: weighted average by number of polymorphic sites per chromosome
  - rare_burden, private_burden: sum across all chromosomes
  - hom_alt_rate used as inbreeding proxy (rice FROH equivalent)

Then: PCA (numpy SVD), optional UMAP, z-distance outlier detection.

Output:
  data/results/rice/rice_inv20b_individual_trajectory_allchr.tsv
  data/results/rice/rice_inv20b_individual_trajectory_allchr.json

Usage:
  conda run -n graphmana python scripts/rice_investigate_20b_individual_trajectory_allchr.py
"""

import csv
import gc
import json
import os
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

# ── Config ────────────────────────────────────────────────────────────────────
_ROOT        = Path(os.environ.get("GRAPHPOP_ROOT", Path(__file__).resolve().parents[1]))
BASE_DIR     = _ROOT
RESULTS      = BASE_DIR / "data" / "results" / "rice"
HAP_CACHE    = BASE_DIR / "data" / "rice_hap_cache"
OUT_TSV      = RESULTS / "rice_inv20b_individual_trajectory_allchr.tsv"
OUT_JSON     = RESULTS / "rice_inv20b_individual_trajectory_allchr.json"

CHUNK        = 50_000       # variants per chunk for unpacking
RARE_AF_THRESH    = 0.05
PRIVATE_AF_THRESH = 0.005

CHROMOSOMES = [f"Chr{i}" for i in range(1, 13)]

FEATURE_NAMES = ["het_rate", "hom_alt_rate", "rare_burden", "private_burden"]

# Population group mapping for rice 3K
POP_GROUP = {
    "GJ-tmp":   "Japonica",
    "GJ-trp":   "Japonica",
    "GJ-sbtrp": "Japonica",
    "GJ-adm":   "Japonica",
    "XI-1A":    "Indica",
    "XI-1B":    "Indica",
    "XI-2":     "Indica",
    "XI-3":     "Indica",
    "XI-adm":   "Indica",
    "cA-Aus":   "Aus",
    "cB-Bas":   "Basmati",
    "admix":    "Admixed",
    "na":       "Unknown",
}

GROUP_COLORS = {
    "Japonica": "#5C8AE0",
    "Indica":   "#E05C5C",
    "Aus":      "#5CC45C",
    "Basmati":  "#E0A85C",
    "Admixed":  "#A85CE0",
    "Unknown":  "#888888",
}


# ── Per-chromosome feature computation ────────────────────────────────────────

def compute_chromosome_features(chrom: str, n_samples: int) -> dict:
    """Compute per-individual features for one chromosome. Returns raw counts."""
    npz_path = HAP_CACHE / f"{chrom}.npz"
    print(f"\n  Loading {chrom}...", flush=True)
    t0 = time.time()

    data = np.load(npz_path)
    hap_packed = data["hap"]        # (n_variants, ceil(2*n_samples/8))
    n_variants = hap_packed.shape[0]
    n_haplotypes = 2 * n_samples
    print(f"    {n_variants:,} variants, {n_haplotypes:,} haplotypes", flush=True)

    # Pass 1: global allele frequencies
    print(f"    Computing allele frequencies...", flush=True)
    global_ac = np.zeros(n_variants, dtype=np.int64)
    n_chunks = (n_variants + CHUNK - 1) // CHUNK

    for cs in range(0, n_variants, CHUNK):
        ce = min(cs + CHUNK, n_variants)
        chunk_hap = np.unpackbits(hap_packed[cs:ce], axis=1, bitorder="big")[:, :n_haplotypes]
        global_ac[cs:ce] = chunk_hap.sum(axis=1)

    global_af = global_ac / n_haplotypes

    poly_mask    = (global_af > 0) & (global_af < 1.0)
    rare_mask    = (global_af > 0) & (global_af < RARE_AF_THRESH)
    private_mask = (global_af > 0) & (global_af < PRIVATE_AF_THRESH)
    n_poly = int(poly_mask.sum())

    print(f"    Polymorphic: {n_poly:,}  Rare: {rare_mask.sum():,}  "
          f"Private: {private_mask.sum():,}  ({time.time()-t0:.1f}s)", flush=True)

    # Pass 2: per-individual counts
    print(f"    Computing per-individual features...", flush=True)
    het_count     = np.zeros(n_samples, dtype=np.int64)
    hom_alt_count = np.zeros(n_samples, dtype=np.int64)
    rare_count    = np.zeros(n_samples, dtype=np.int64)
    private_count = np.zeros(n_samples, dtype=np.int64)

    for ci, cs in enumerate(range(0, n_variants, CHUNK)):
        ce = min(cs + CHUNK, n_variants)
        chunk_hap = np.unpackbits(hap_packed[cs:ce], axis=1, bitorder="big")[:, :n_haplotypes]

        h1 = chunk_hap[:, 0::2]  # (chunk_size, n_samples)
        h2 = chunk_hap[:, 1::2]

        # Het and hom-alt at polymorphic sites only
        cm_poly = poly_mask[cs:ce]
        if cm_poly.any():
            het_count     += (h1[cm_poly] != h2[cm_poly]).sum(axis=0)
            hom_alt_count += ((h1[cm_poly] == 1) & (h2[cm_poly] == 1)).sum(axis=0)

        # Rare and private burden (any alt allele carried)
        alt_dose = h1.astype(np.int16) + h2

        cm_rare = rare_mask[cs:ce]
        if cm_rare.any():
            rare_count += (alt_dose[cm_rare] > 0).sum(axis=0)

        cm_priv = private_mask[cs:ce]
        if cm_priv.any():
            private_count += (alt_dose[cm_priv] > 0).sum(axis=0)

        if (ci + 1) % 20 == 0 or ci == n_chunks - 1:
            print(f"      chunk {ci+1}/{n_chunks}  ({time.time()-t0:.1f}s)", flush=True)

    # Clean up to free memory
    del hap_packed, data, global_ac, global_af, poly_mask, rare_mask, private_mask
    gc.collect()

    elapsed = time.time() - t0
    print(f"    {chrom} done in {elapsed:.1f}s  "
          f"het_mean={het_count.mean():.0f}  rare_mean={rare_count.mean():.0f}", flush=True)

    return {
        "het_count":     het_count,
        "hom_alt_count": hom_alt_count,
        "rare_count":    rare_count,
        "private_count": private_count,
        "n_poly":        n_poly,
    }


# ── PCA via numpy SVD ────────────────────────────────────────────────────────

def pca_svd(mat_z: np.ndarray, n_components: int = 4):
    """PCA using numpy SVD (no sklearn dependency)."""
    # Center (should already be zero-mean if standardized, but be safe)
    centered = mat_z - mat_z.mean(axis=0)
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    # Explained variance
    explained_var = (S ** 2) / (mat_z.shape[0] - 1)
    total_var = explained_var.sum()
    explained_ratio = explained_var / total_var
    # PC scores
    n_comp = min(n_components, Vt.shape[0])
    scores = centered @ Vt[:n_comp].T
    return scores, explained_ratio[:n_comp], Vt[:n_comp]


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()
    print("=" * 78)
    print("Rice Investigation 20b: Genome-Wide Individual Trajectory (All 12 Chromosomes)")
    print("=" * 78)

    # 1. Sample metadata
    print("\nLoading sample metadata...", flush=True)
    meta = json.load(open(HAP_CACHE / "sample_meta.json"))
    sample_ids = meta["sample_ids"]
    pop_ids    = meta["pop_ids"]
    n_samples  = len(sample_ids)
    print(f"  {n_samples} samples, {len(set(pop_ids))} populations", flush=True)

    for grp in sorted(set(POP_GROUP.get(p, "Unknown") for p in pop_ids)):
        cnt = sum(1 for p in pop_ids if POP_GROUP.get(p, "Unknown") == grp)
        print(f"    {grp}: {cnt}")

    # 2. Compute per-chromosome features and aggregate
    print(f"\nProcessing {len(CHROMOSOMES)} chromosomes...", flush=True)

    # Accumulators: weighted average for rates, sum for burdens
    total_het_count     = np.zeros(n_samples, dtype=np.int64)
    total_hom_alt_count = np.zeros(n_samples, dtype=np.int64)
    total_rare_burden   = np.zeros(n_samples, dtype=np.int64)
    total_private_burden= np.zeros(n_samples, dtype=np.int64)
    total_poly_sites    = 0

    chr_summaries = {}
    for chrom in CHROMOSOMES:
        result = compute_chromosome_features(chrom, n_samples)

        total_het_count      += result["het_count"]
        total_hom_alt_count  += result["hom_alt_count"]
        total_rare_burden    += result["rare_count"]
        total_private_burden += result["private_count"]
        total_poly_sites     += result["n_poly"]

        chr_summaries[chrom] = {
            "n_poly": result["n_poly"],
            "mean_het": round(float(result["het_count"].mean()), 1),
            "mean_rare": round(float(result["rare_count"].mean()), 1),
        }

        # Free per-chromosome arrays
        del result
        gc.collect()

    # Compute genome-wide rates (weighted by total polymorphic sites)
    print(f"\n  Total polymorphic sites across genome: {total_poly_sites:,}", flush=True)
    het_rate     = total_het_count / max(total_poly_sites, 1)
    hom_alt_rate = total_hom_alt_count / max(total_poly_sites, 1)
    rare_burden  = total_rare_burden.astype(np.float64)
    private_burden = total_private_burden.astype(np.float64)

    print(f"  Genome-wide het_rate:     [{het_rate.min():.6f}, {het_rate.max():.6f}]  "
          f"mean={het_rate.mean():.6f}", flush=True)
    print(f"  Genome-wide hom_alt_rate: [{hom_alt_rate.min():.6f}, {hom_alt_rate.max():.6f}]  "
          f"mean={hom_alt_rate.mean():.6f}", flush=True)
    print(f"  Genome-wide rare_burden:  [{rare_burden.min():.0f}, {rare_burden.max():.0f}]  "
          f"mean={rare_burden.mean():.1f}", flush=True)
    print(f"  Genome-wide priv_burden:  [{private_burden.min():.0f}, {private_burden.max():.0f}]  "
          f"mean={private_burden.mean():.1f}", flush=True)

    # 3. Build feature matrix
    print("\nBuilding feature matrix...", flush=True)
    mat = np.column_stack([het_rate, hom_alt_rate, rare_burden, private_burden])
    assert mat.shape == (n_samples, 4)

    # Standardize (z-score)
    means = mat.mean(axis=0)
    stds  = mat.std(axis=0)
    stds[stds == 0] = 1.0  # prevent div-by-zero
    mat_z = (mat - means) / stds

    for j, fname in enumerate(FEATURE_NAMES):
        v = mat[:, j]
        print(f"  {fname:20s}  mean={v.mean():.6f}  std={v.std():.6f}  "
              f"range=[{v.min():.6f}, {v.max():.6f}]")

    # 4. PCA via numpy SVD
    print("\nComputing PCA (numpy SVD)...", flush=True)
    pca_coords, pca_var_ratio, pca_components = pca_svd(mat_z, n_components=4)
    for i in range(len(pca_var_ratio)):
        top3 = sorted(range(4), key=lambda j: abs(pca_components[i][j]), reverse=True)[:3]
        top3_str = ", ".join(f"{FEATURE_NAMES[j]}({pca_components[i][j]:+.3f})" for j in top3)
        print(f"  PC{i+1}: {pca_var_ratio[i]*100:.1f}%  |  {top3_str}")

    # 5. UMAP (optional)
    umap_coords = None
    try:
        import umap as umap_lib
        print("\nComputing UMAP (n_neighbors=15, min_dist=0.2)...", flush=True)
        reducer = umap_lib.UMAP(
            n_components=2, n_neighbors=15, min_dist=0.2,
            random_state=42, metric="euclidean",
        )
        umap_coords = reducer.fit_transform(mat_z)
        print(f"  UMAP range: x=[{umap_coords[:,0].min():.2f},{umap_coords[:,0].max():.2f}]  "
              f"y=[{umap_coords[:,1].min():.2f},{umap_coords[:,1].max():.2f}]")
    except ImportError:
        print("\nUMAP not available, skipping.", flush=True)

    # 6. Z-distance from subpopulation centroid (outlier detection)
    print("\nComputing z-distance from subpopulation centroids...", flush=True)
    groups = [POP_GROUP.get(p, "Unknown") for p in pop_ids]
    z_distance = np.zeros(n_samples, dtype=np.float64)

    # Compute centroid per subpopulation (using pop_ids, not major groups)
    pop_centroids = {}
    for pop in set(pop_ids):
        idx = [i for i, p in enumerate(pop_ids) if p == pop]
        if len(idx) >= 3:
            centroid = mat_z[idx].mean(axis=0)
            cov = np.cov(mat_z[idx].T)
            # Regularize covariance
            cov += np.eye(4) * 1e-6
            try:
                cov_inv = np.linalg.inv(cov)
            except np.linalg.LinAlgError:
                cov_inv = np.eye(4)
            pop_centroids[pop] = (centroid, cov_inv)

    for i in range(n_samples):
        pop = pop_ids[i]
        if pop in pop_centroids:
            centroid, cov_inv = pop_centroids[pop]
            diff = mat_z[i] - centroid
            # Mahalanobis distance
            z_distance[i] = np.sqrt(diff @ cov_inv @ diff)
        else:
            # No centroid available, use Euclidean distance from global mean
            z_distance[i] = np.sqrt(np.sum(mat_z[i] ** 2))

    # Find outliers (z_distance > 3)
    outlier_mask = z_distance > 3.0
    n_outliers = outlier_mask.sum()
    print(f"  Outliers (Mahalanobis > 3): {n_outliers} / {n_samples} "
          f"({100*n_outliers/n_samples:.1f}%)", flush=True)
    if n_outliers > 0:
        outlier_idx = np.where(outlier_mask)[0]
        for oi in outlier_idx[:10]:  # show first 10
            print(f"    {sample_ids[oi]:12s}  pop={pop_ids[oi]:10s}  "
                  f"z_dist={z_distance[oi]:.2f}  "
                  f"het={het_rate[oi]:.6f}  hom_alt={hom_alt_rate[oi]:.6f}")
        if n_outliers > 10:
            print(f"    ... and {n_outliers - 10} more")

    # 7. Save outputs
    print("\nSaving outputs...", flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    # TSV
    fields = [
        "sampleId", "population", "group",
        "het_rate", "hom_alt_rate", "rare_burden", "private_burden",
        "inbreeding_proxy",  # hom_alt_rate as FROH equivalent
        "z_distance", "is_outlier",
        "stat_pc1", "stat_pc2", "stat_pc3", "stat_pc4",
    ]
    if umap_coords is not None:
        fields += ["umap1", "umap2"]

    rows = []
    for i in range(n_samples):
        row = {
            "sampleId":         sample_ids[i],
            "population":       pop_ids[i],
            "group":            groups[i],
            "het_rate":         round(float(het_rate[i]), 8),
            "hom_alt_rate":     round(float(hom_alt_rate[i]), 8),
            "rare_burden":      int(rare_burden[i]),
            "private_burden":   int(private_burden[i]),
            "inbreeding_proxy": round(float(hom_alt_rate[i]), 8),
            "z_distance":       round(float(z_distance[i]), 4),
            "is_outlier":       bool(outlier_mask[i]),
            "stat_pc1":         round(float(pca_coords[i, 0]), 4),
            "stat_pc2":         round(float(pca_coords[i, 1]), 4),
            "stat_pc3":         round(float(pca_coords[i, 2]), 4) if pca_coords.shape[1] > 2 else 0.0,
            "stat_pc4":         round(float(pca_coords[i, 3]), 4) if pca_coords.shape[1] > 3 else 0.0,
        }
        if umap_coords is not None:
            row["umap1"] = round(float(umap_coords[i, 0]), 4)
            row["umap2"] = round(float(umap_coords[i, 1]), 4)
        rows.append(row)

    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  TSV: {OUT_TSV} ({len(rows)} rows)", flush=True)

    # JSON summary
    group_stats = {}
    for grp in sorted(set(groups)):
        idx = [i for i, g in enumerate(groups) if g == grp]
        if not idx:
            continue
        group_stats[grp] = {
            "n_samples": len(idx),
            "het_rate":     {"mean": round(float(het_rate[idx].mean()), 8),
                             "std":  round(float(het_rate[idx].std()), 8)},
            "hom_alt_rate": {"mean": round(float(hom_alt_rate[idx].mean()), 8),
                             "std":  round(float(hom_alt_rate[idx].std()), 8)},
            "rare_burden":  {"mean": round(float(rare_burden[idx].mean()), 1),
                             "std":  round(float(rare_burden[idx].std()), 1)},
            "private_burden": {"mean": round(float(private_burden[idx].mean()), 1),
                               "std":  round(float(private_burden[idx].std()), 1)},
            "inbreeding_proxy": {"mean": round(float(hom_alt_rate[idx].mean()), 8),
                                 "std":  round(float(hom_alt_rate[idx].std()), 8)},
            "z_distance":   {"mean": round(float(z_distance[idx].mean()), 4),
                             "std":  round(float(z_distance[idx].std()), 4)},
            "n_outliers":   int(outlier_mask[idx].sum()),
        }

    out_data = {
        "analysis":              "rice_inv20b_individual_trajectory_allchr",
        "description":           "Genome-wide individual evolutionary trajectory from all 12 chromosomes",
        "n_samples":             n_samples,
        "n_features":            4,
        "features":              FEATURE_NAMES,
        "chromosomes":           CHROMOSOMES,
        "total_polymorphic_sites": total_poly_sites,
        "rare_af_threshold":     RARE_AF_THRESH,
        "private_af_threshold":  PRIVATE_AF_THRESH,
        "pca_explained_variance": [round(float(v), 4) for v in pca_var_ratio],
        "pca_components":        [[round(float(c), 4) for c in row] for row in pca_components],
        "chromosome_summaries":  chr_summaries,
        "group_statistics":      group_stats,
        "outlier_threshold":     3.0,
        "total_outliers":        int(n_outliers),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2)
    print(f"  JSON: {OUT_JSON}", flush=True)

    # 8. Summary
    elapsed = time.time() - t_start
    print(f"\n{'=' * 78}")
    print(f"SUMMARY")
    print(f"{'=' * 78}")
    print(f"  Chromosomes processed:     {len(CHROMOSOMES)}")
    print(f"  Total polymorphic sites:   {total_poly_sites:,}")
    print(f"  Samples:                   {n_samples}")
    print(f"  Outliers (z>3):            {n_outliers}")
    print(f"  PCA variance explained:    {[f'{v*100:.1f}%' for v in pca_var_ratio]}")
    print(f"  UMAP computed:             {umap_coords is not None}")
    print(f"  Total time:                {elapsed:.1f}s")
    print()

    for grp in ["Japonica", "Indica", "Aus", "Basmati", "Admixed", "Unknown"]:
        if grp in group_stats:
            gs = group_stats[grp]
            print(f"  {grp:12s} (n={gs['n_samples']:4d})  "
                  f"het={gs['het_rate']['mean']:.6f}  "
                  f"hom_alt={gs['hom_alt_rate']['mean']:.6f}  "
                  f"rare={gs['rare_burden']['mean']:.0f}  "
                  f"private={gs['private_burden']['mean']:.0f}  "
                  f"outliers={gs['n_outliers']}")

    print(f"\nDone.")


if __name__ == "__main__":
    main()
