#!/usr/bin/env python3
"""Investigation 19: Individual-Level Pathway Trajectory Embedding (Analysis 5).

Places 3,202 human individuals in pathway-rare-burden space:
  cell[i, j] = count of rare missense variants in pathway j carried by individual i
               (genome-wide, all 22 autosomes, global AF < 0.01)

Unlike standard PCA (ancestry-based) or the process-space embedding (Inv.16b, 5 scalar
features), this is a high-dimensional functional representation: individuals are
characterised by WHICH biological systems accumulate their private rare variation,
not just HOW MUCH.

Inspired by single-cell trajectory inference:
  individuals  ↔  cells
  pathways     ↔  genes
  rare burden  ↔  expression

Method:
  1. Load variant→pathway map (variant_pathway_cache.json, no Neo4j needed)
  2. For each of 22 chromosomes: chunked hap_cache pass to accumulate
     per-individual, per-pathway rare missense carrier counts
  3. Build sparse (3202, N_pathways) matrix; filter low-variance pathways
  4. TruncatedSVD → UMAP for nonlinear embedding
  5. Analyse: which pathways drive continental separation?
     Do pathway clusters correspond to known biology?
  6. Compare to Inv.16b individual process-space (panel F)

Output:
  data/results/pathway_trajectory.tsv / .json
  data/results/figures/fig24_pathway_trajectory.png

Runtime estimate: ~25-35 minutes (22 chromosome hap passes)
No Neo4j required — runs entirely from cached numpy + JSON files.

Usage:
  conda run -n graphevo python -u scripts/investigate_19_pathway_trajectory.py
"""

import csv
import json
import time
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.sparse import lil_matrix, csr_matrix
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")

RESULTS   = Path("data/results")
FIG_DIR   = RESULTS / "figures"
HAP_DIR   = Path("data/hap_cache")
VAR_DIR   = Path("data/variants_cache")
CACHE     = RESULTS / "variant_pathway_cache.json"
OUT_TSV   = RESULTS / "pathway_trajectory.tsv"
OUT_JSON  = RESULTS / "pathway_trajectory.json"
CHUNK     = 50_000
RARE_AF   = 0.01

# Minimum fraction of individuals carrying ≥1 rare missense in a pathway
# (filters near-zero-variance pathways before SVD/UMAP)
MIN_CARRIER_FRAC = 0.02

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]

GROUP_COLORS = {
    "African":     "#E05C5C",
    "European":    "#5C8AE0",
    "South_Asian": "#E0A85C",
    "East_Asian":  "#5CC45C",
    "American":    "#A85CE0",
}
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


# ── Step 1: Load caches ───────────────────────────────────────────────────────

def load_pathway_cache() -> tuple:
    """Load variant→pathway map and build pathway index."""
    print("Loading variant_pathway_cache.json...", flush=True)
    cache = json.load(open(CACHE))
    variant_pathways = cache["variant_pathways"]   # {chr: {pos_str: [pathway_ids]}}
    pathway_meta     = cache["pathway_meta"]        # {pathway_id: {name, mean_fst}}

    # Build sorted pathway index
    all_pathway_ids = sorted({
        pid
        for chrom_data in variant_pathways.values()
        for pid_list in chrom_data.values()
        for pid in pid_list
    })
    pathway_idx = {pid: i for i, pid in enumerate(all_pathway_ids)}
    n_pathways  = len(all_pathway_ids)
    print(f"  {n_pathways:,} pathways  |  "
          f"{sum(len(v) for v in variant_pathways.values()):,} annotated positions", flush=True)
    return variant_pathways, pathway_meta, pathway_idx, all_pathway_ids


# ── Step 2: Genome-wide accumulation ─────────────────────────────────────────

def build_pathway_matrix(
    variant_pathways: dict,
    pathway_idx: dict,
    n_samples: int,
    n_pathways: int,
) -> np.ndarray:
    """
    Accumulate per-individual, per-pathway rare missense carrier counts.
    Returns dense (n_samples, n_pathways) float32 matrix.
    """
    mat = np.zeros((n_samples, n_pathways), dtype=np.float32)
    n_hap = 2 * n_samples

    for chrom in CHROMOSOMES:
        hap_path = HAP_DIR / f"{chrom}.npz"
        var_path  = VAR_DIR  / f"{chrom}.npz"
        if not hap_path.exists() or not var_path.exists():
            print(f"  SKIP {chrom} (cache missing)", flush=True)
            continue

        chrom_data = variant_pathways.get(chrom, {})
        if not chrom_data:
            print(f"  {chrom}: no annotated positions, skipping", flush=True)
            continue

        hap_data   = np.load(hap_path)
        var_data   = np.load(var_path)
        hap_packed = hap_data["hap"]
        pos_all    = var_data["pos"]
        vtype      = var_data["variant_type"]
        ac_all     = var_data["ac"]
        an_all     = var_data["an"]
        mis_all    = var_data["has_missense"]

        snp_mask  = vtype == 1
        pos_snps  = pos_all[snp_mask]
        hap_snp   = hap_packed[snp_mask]
        ac_snps   = ac_all[snp_mask].astype(np.int64)
        an_snps   = an_all[snp_mask].astype(np.int64)
        mis_snps  = mis_all[snp_mask]
        global_af = ac_snps.sum(axis=1) / an_snps.sum(axis=1).clip(1)

        # Build index: which SNP positions (in snp-filtered array) are in chrom_data?
        pos_to_pathways = {}
        for i, pos in enumerate(pos_snps):
            pos_str = str(int(pos))
            if pos_str in chrom_data:
                pids = chrom_data[pos_str]
                col_idxs = [pathway_idx[p] for p in pids if p in pathway_idx]
                if col_idxs:
                    pos_to_pathways[i] = col_idxs

        # Filter: rare missense AND annotated positions
        # Build arrays for fast chunked access
        target_snp_indices = np.array([
            i for i in pos_to_pathways
            if mis_snps[i] and global_af[i] < RARE_AF
        ], dtype=np.int64)

        n_target = len(target_snp_indices)
        t_chr = time.time()
        print(f"  {chrom}: {n_target:,} rare missense positions with pathway annotation...",
              end="", flush=True)

        if n_target == 0:
            print(" skip", flush=True)
            continue

        # Process each target position individually (they're scattered, not contiguous)
        # Chunk by groups of CHUNK positions for memory efficiency
        for batch_start in range(0, n_target, CHUNK):
            batch_end = min(batch_start + CHUNK, n_target)
            batch_idx = target_snp_indices[batch_start:batch_end]

            for snp_i in batch_idx:
                # Unpack single row
                row = np.unpackbits(
                    hap_snp[snp_i:snp_i+1], axis=1, bitorder="big"
                )[0, :n_hap]
                dosage = (row[0::2] + row[1::2]).astype(np.float32)  # (n_samples,)
                carrier = (dosage > 0).astype(np.float32)             # 0 or 1

                # Increment all pathway columns for this position
                for col in pos_to_pathways[snp_i]:
                    mat[:, col] += carrier

        print(f"  {time.time()-t_chr:.0f}s", flush=True)

    return mat


# ── Step 3: Filter, normalise, embed ─────────────────────────────────────────

def filter_matrix(mat: np.ndarray, all_pathway_ids: list,
                  pathway_meta: dict, min_frac: float = MIN_CARRIER_FRAC):
    """Remove pathways where too few individuals carry any rare missense."""
    carrier_frac = (mat > 0).mean(axis=0)
    keep = carrier_frac >= min_frac
    mat_filt      = mat[:, keep]
    pids_kept     = [all_pathway_ids[i] for i in range(len(all_pathway_ids)) if keep[i]]
    print(f"  Pathways with ≥{min_frac*100:.0f}% carrier fraction: "
          f"{keep.sum():,} / {len(all_pathway_ids):,}", flush=True)
    return mat_filt, pids_kept


def embed(mat_filt: np.ndarray, n_svd_components: int = 50):
    """TruncatedSVD then UMAP."""
    from sklearn.decomposition import TruncatedSVD
    from sklearn.preprocessing import normalize
    import umap as umap_lib

    print(f"  TruncatedSVD ({n_svd_components} components)...", flush=True)
    svd = TruncatedSVD(n_components=n_svd_components, random_state=42)
    svd_coords = svd.fit_transform(mat_filt)
    ev = svd.explained_variance_ratio_
    print(f"  SVD: cumulative variance (top-{n_svd_components}): "
          f"{ev.sum()*100:.1f}%  (PC1={ev[0]*100:.1f}%, PC2={ev[1]*100:.1f}%)", flush=True)

    print(f"  UMAP (n_neighbors=15, min_dist=0.15, metric=cosine)...", flush=True)
    reducer = umap_lib.UMAP(
        n_components=2, n_neighbors=15, min_dist=0.15,
        metric="cosine", random_state=42,
    )
    umap_coords = reducer.fit_transform(svd_coords)
    return svd_coords, svd, umap_coords


# ── Step 4: Top pathways per group ───────────────────────────────────────────

def top_pathways_per_group(mat_filt, pids_kept, pathway_meta, pop_ids, top_n=10):
    """For each continental group, find pathways with highest mean burden."""
    groups = ["African", "European", "East_Asian", "South_Asian", "American"]
    grp_top = {}
    for grp in groups:
        idx = [i for i, p in enumerate(pop_ids) if POP_GROUP.get(p) == grp]
        if not idx:
            continue
        means = mat_filt[idx].mean(axis=0)
        top_idx = np.argsort(means)[::-1][:top_n]
        grp_top[grp] = [
            {"pathway_id": pids_kept[j],
             "name":       pathway_meta.get(pids_kept[j], {}).get("name", pids_kept[j])[:60],
             "mean_burden": round(float(means[j]), 4),
             "mean_fst":   pathway_meta.get(pids_kept[j], {}).get("mean_fst")}
            for j in top_idx
        ]
    return grp_top


# ── Step 5: Figure ────────────────────────────────────────────────────────────

def make_figure(sample_ids, pop_ids, svd_coords, umap_coords, svd,
                mat_filt, pids_kept, pathway_meta,
                grp_top, ind_trajectory_path):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    groups  = [POP_GROUP.get(p, "Other") for p in pop_ids]
    colors  = [GROUP_COLORS.get(g, "#888") for g in groups]

    fig = plt.figure(figsize=(22, 16))
    fig.patch.set_facecolor("#0D0D0D")
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.35,
                          left=0.07, right=0.97, top=0.93, bottom=0.07)
    axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]
    for ax in axes:
        ax.set_facecolor("#161616")
        for spine in ax.spines.values():
            spine.set_color("#444")

    # ── A: Pathway-space UMAP ─────────────────────────────────────────────
    ax = axes[0]
    ax.set_title("A  Pathway Rare-Burden Space (UMAP, n=3,202)", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
               c=colors, s=6, alpha=0.55, linewidths=0, zorder=2)
    for grp, col in GROUP_COLORS.items():
        ax.scatter([], [], c=col, s=30, label=grp.replace("_", " "))
    ax.legend(fontsize=7, loc="best", facecolor="#222",
              labelcolor="#CCC", framealpha=0.8)
    ax.set_xlabel("Pathway-UMAP 1", color="#AAA", fontsize=8)
    ax.set_ylabel("Pathway-UMAP 2", color="#AAA", fontsize=8)
    ax.tick_params(colors="#888", labelsize=7)

    # ── B: SVD PC1 vs PC2 ─────────────────────────────────────────────────
    ax = axes[1]
    ev = svd.explained_variance_ratio_
    ax.set_title(f"B  Pathway SVD  (PC1={ev[0]*100:.1f}%, PC2={ev[1]*100:.1f}%)",
                 color="#EEE", fontsize=11, pad=6, loc="left")
    ax.scatter(svd_coords[:, 0], svd_coords[:, 1],
               c=colors, s=6, alpha=0.55, linewidths=0, zorder=2)
    ax.set_xlabel(f"SVD-PC1 ({ev[0]*100:.1f}%)", color="#AAA", fontsize=8)
    ax.set_ylabel(f"SVD-PC2 ({ev[1]*100:.1f}%)", color="#AAA", fontsize=8)
    ax.axhline(0, color="#333", lw=0.5); ax.axvline(0, color="#333", lw=0.5)
    ax.tick_params(colors="#888", labelsize=7)

    # ── C: Top SVD-PC1 pathway loadings ───────────────────────────────────
    ax = axes[2]
    ax.set_title("C  Pathway SVD-PC1 top loadings", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    loadings = svd.components_[0]
    top_pos = np.argsort(loadings)[::-1][:8]
    top_neg = np.argsort(loadings)[:8]
    top_idx = np.concatenate([top_neg, top_pos])
    load_vals = loadings[top_idx]
    names = [pathway_meta.get(pids_kept[j], {}).get("name", pids_kept[j])[:30]
             for j in top_idx]
    bar_colors = ["#E05C5C" if v > 0 else "#5C8AE0" for v in load_vals]
    ax.barh(range(len(top_idx)), load_vals[::-1],
            color=bar_colors[::-1], edgecolor="#333")
    ax.set_yticks(range(len(top_idx)))
    ax.set_yticklabels([n[:28] for n in names[::-1]], fontsize=5.5, color="#CCC")
    ax.set_xlabel("SVD-PC1 loading", color="#AAA", fontsize=8)
    ax.axvline(0, color="#555", lw=0.8)
    ax.tick_params(axis="x", colors="#888", labelsize=7)

    # ── D: Total pathway burden by group (violin) ─────────────────────────
    ax = axes[3]
    ax.set_title("D  Total pathway-annotated rare missense by group",
                 color="#EEE", fontsize=11, pad=6, loc="left")
    group_order = ["African","European","East_Asian","South_Asian","American"]
    total_burden = mat_filt.sum(axis=1)
    for xi, grp in enumerate(group_order):
        col = GROUP_COLORS.get(grp, "#888")
        idx = [i for i, p in enumerate(pop_ids) if POP_GROUP.get(p) == grp]
        vals = total_burden[idx]
        if len(vals) > 3:
            vp = ax.violinplot([vals], positions=[xi], widths=0.7,
                               showmedians=True, showextrema=False)
            for body in vp["bodies"]:
                body.set_facecolor(col); body.set_alpha(0.6); body.set_edgecolor("none")
            vp["cmedians"].set_color(col); vp["cmedians"].set_linewidth(2)
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels([g.replace("_","\n") for g in group_order],
                       color="#CCC", fontsize=7)
    ax.tick_params(axis="y", colors="#888", labelsize=7)
    ax.set_ylabel("Total rare missense burden\n(pathway-annotated, genome-wide)",
                  color="#AAA", fontsize=8)

    # ── E: Compare pathway-UMAP-1 vs process-space PC1 ────────────────────
    ax = axes[4]
    ax.set_title("E  Pathway-UMAP-1 vs Process-Space PC1 (Inv.16b)",
                 color="#EEE", fontsize=11, pad=6, loc="left")
    if ind_trajectory_path.exists():
        ind_rows = {r["sampleId"]: float(r["stat_pc1"])
                    for r in csv.DictReader(
                        open(ind_trajectory_path), delimiter="\t")}
        proc_pc1 = np.array([ind_rows.get(s, np.nan) for s in sample_ids])
        valid = ~np.isnan(proc_pc1)
        ax.scatter(proc_pc1[valid], umap_coords[valid, 0],
                   c=[colors[i] for i in np.where(valid)[0]],
                   s=5, alpha=0.45, linewidths=0, zorder=2)
        if valid.sum() > 10:
            rho, p = spearmanr(proc_pc1[valid], umap_coords[valid, 0])
            ax.set_title(f"E  Pathway-UMAP-1 vs Process-PC1  (ρ={rho:.3f}, p={p:.2e})",
                         color="#EEE", fontsize=9, pad=6, loc="left")
        ax.set_xlabel("Individual Process-Space PC1 (Inv.16b)", color="#AAA", fontsize=8)
        ax.set_ylabel("Pathway-UMAP 1", color="#AAA", fontsize=8)
        ax.tick_params(colors="#888", labelsize=7)

    # ── F: Top pathways per group heatmap ─────────────────────────────────
    ax = axes[5]
    ax.set_title("F  Top pathway burden by continental group (mean variants/individual)",
                 color="#EEE", fontsize=11, pad=6, loc="left")
    # Collect top 5 pathways from African (highest burden) for heatmap rows
    afr_top5_pids = [x["pathway_id"] for x in grp_top.get("African", [])[:5]]
    afr_top5_names = [pathway_meta.get(p, {}).get("name", p)[:35] for p in afr_top5_pids]
    group_order2 = ["African", "European", "East_Asian", "South_Asian", "American"]
    hmap = []
    for pid in afr_top5_pids:
        if pid not in {pids_kept[i]: i for i in range(len(pids_kept))}:
            continue
        col_i = pids_kept.index(pid) if pid in pids_kept else -1
        if col_i < 0:
            continue
        row = []
        for grp in group_order2:
            idx = [i for i, p in enumerate(pop_ids) if POP_GROUP.get(p) == grp]
            row.append(mat_filt[idx, col_i].mean() if idx else 0.0)
        hmap.append(row)
    if hmap:
        hmap = np.array(hmap)
        im = ax.imshow(hmap, aspect="auto", cmap="YlOrRd", vmin=0)
        ax.set_xticks(range(len(group_order2)))
        ax.set_xticklabels([g.replace("_","\n") for g in group_order2],
                           color="#CCC", fontsize=7)
        ax.set_yticks(range(len(afr_top5_names)))
        ax.set_yticklabels(afr_top5_names, fontsize=6, color="#CCC")
        plt.colorbar(im, ax=ax, fraction=0.03, pad=0.03,
                     label="Mean rare variants/individual").ax.yaxis.set_tick_params(
                         color="#888", labelcolor="#888")

    fig.suptitle(
        "Individual-Level Pathway Trajectory — 3,202 Humans in Pathway Rare-Burden Space\n"
        f"Dimensions: genome-wide rare missense burden per Reactome pathway "
        f"(all 22 autosomes, AF<0.01)",
        color="#EEE", fontsize=12, y=0.97,
    )
    out_path = FIG_DIR / "fig24_pathway_trajectory.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()
    print("=== Investigation 19: Individual-Level Pathway Trajectory Embedding ===\n")

    # Load sample metadata
    meta = json.load(open(HAP_DIR / "sample_meta.json"))
    sample_ids = meta["sample_ids"]
    pop_ids    = meta["pop_ids"]
    n_samples  = len(sample_ids)
    print(f"Samples: {n_samples}  |  Populations: {len(set(pop_ids))}\n")

    # Step 1: Load pathway cache
    variant_pathways, pathway_meta, pathway_idx, all_pathway_ids = load_pathway_cache()
    n_pathways = len(all_pathway_ids)

    # Step 2: Build genome-wide pathway burden matrix
    print(f"\nBuilding {n_samples} × {n_pathways} pathway burden matrix...", flush=True)
    print("(Processing 22 chromosomes — estimated 25-35 minutes)\n", flush=True)
    mat = build_pathway_matrix(variant_pathways, pathway_idx, n_samples, n_pathways)

    total_counts = mat.sum()
    print(f"\nMatrix complete: {total_counts:,.0f} total carrier-pathway events", flush=True)
    print(f"Mean per individual: {mat.sum(axis=1).mean():.1f}  "
          f"Range: [{mat.sum(axis=1).min():.0f}, {mat.sum(axis=1).max():.0f}]", flush=True)

    # Step 3: Filter low-variance pathways
    print("\nFiltering pathways...", flush=True)
    mat_filt, pids_kept = filter_matrix(mat, all_pathway_ids, pathway_meta)
    n_kept = mat_filt.shape[1]

    # Step 4: SVD + UMAP embedding
    print(f"\nEmbedding {n_samples} × {n_kept} matrix...", flush=True)
    n_svd = min(50, n_kept - 1)
    svd_coords, svd, umap_coords = embed(mat_filt, n_svd_components=n_svd)
    print(f"  UMAP range: x=[{umap_coords[:,0].min():.2f},{umap_coords[:,0].max():.2f}]  "
          f"y=[{umap_coords[:,1].min():.2f},{umap_coords[:,1].max():.2f}]", flush=True)

    # Step 5: Top pathways per group
    print("\nComputing top pathways per continental group...", flush=True)
    grp_top = top_pathways_per_group(mat_filt, pids_kept, pathway_meta, pop_ids)
    for grp, top in grp_top.items():
        print(f"  {grp} top 3: " +
              ", ".join(f"{t['name'][:30]}({t['mean_burden']:.3f})" for t in top[:3]))

    # Correlation: pathway-UMAP-1 vs process-space PC1 (from Inv.16b)
    print("\nCorrelation: pathway-UMAP-1 vs individual process-space PC1...", flush=True)
    ind_traj_path = RESULTS / "individual_trajectory.tsv"
    if ind_traj_path.exists():
        ind_rows = {r["sampleId"]: float(r["stat_pc1"])
                    for r in csv.DictReader(open(ind_traj_path), delimiter="\t")}
        proc_pc1 = np.array([ind_rows.get(s, np.nan) for s in sample_ids])
        valid = ~np.isnan(proc_pc1)
        rho, p = spearmanr(proc_pc1[valid], umap_coords[valid, 0])
        print(f"  Pathway-UMAP-1 vs Process-PC1: ρ={rho:.3f}, p={p:.2e}")
        rho2, p2 = spearmanr(proc_pc1[valid], svd_coords[valid, 0])
        print(f"  Pathway-SVD-1  vs Process-PC1: ρ={rho2:.3f}, p={p2:.2e}")

    # Step 6: Figure
    print("\nGenerating figure...", flush=True)
    make_figure(sample_ids, pop_ids, svd_coords, umap_coords, svd,
                mat_filt, pids_kept, pathway_meta, grp_top, ind_traj_path)

    # Step 7: Write outputs
    print("\nWriting outputs...", flush=True)
    rows = []
    for i, sid in enumerate(sample_ids):
        pop = pop_ids[i]
        rows.append({
            "sampleId":      sid,
            "population":    pop,
            "group":         POP_GROUP.get(pop, ""),
            "total_pw_burden": round(float(mat_filt[i].sum()), 2),
            "svd_pc1":       round(float(svd_coords[i, 0]), 4),
            "svd_pc2":       round(float(svd_coords[i, 1]), 4),
            "pw_umap1":      round(float(umap_coords[i, 0]), 4),
            "pw_umap2":      round(float(umap_coords[i, 1]), 4),
        })
    fields = ["sampleId","population","group","total_pw_burden",
              "svd_pc1","svd_pc2","pw_umap1","pw_umap2"]
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)

    out = {
        "n_samples":          n_samples,
        "n_pathways_total":   n_pathways,
        "n_pathways_kept":    n_kept,
        "min_carrier_frac":   MIN_CARRIER_FRAC,
        "svd_components":     n_svd,
        "svd_explained_var":  [round(float(v), 4) for v in svd.explained_variance_ratio_],
        "top_pathways_by_group": grp_top,
        "embeddings": {
            sid: {
                "population": pop_ids[i],
                "group":      POP_GROUP.get(pop_ids[i], ""),
                "pw_umap1":   round(float(umap_coords[i, 0]), 4),
                "pw_umap2":   round(float(umap_coords[i, 1]), 4),
            }
            for i, sid in enumerate(sample_ids)
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)

    elapsed = time.time() - t_start
    print(f"\nOutputs: {OUT_TSV}\n         {OUT_JSON}")
    print(f"Total time: {elapsed/60:.1f} min")

    # Print group SVD summary
    print("\n=== Group Summary (SVD-PC1, pathway space) ===")
    for grp in ["African","European","East_Asian","South_Asian","American"]:
        idx = [i for i, p in enumerate(pop_ids) if POP_GROUP.get(p) == grp]
        if idx:
            pc1_m = svd_coords[idx, 0].mean()
            tot_m = mat_filt[idx].sum(axis=1).mean()
            print(f"  {grp:12s}: SVD-PC1={pc1_m:+.3f}  "
                  f"total_pw_burden={tot_m:.1f}")


if __name__ == "__main__":
    main()
