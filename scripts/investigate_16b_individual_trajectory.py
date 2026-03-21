#!/usr/bin/env python3
"""Investigation 16b: Individual-Level Evolutionary Trajectory Map.

Companion to Inv.16 (population-level). Places all 3,202 individual samples in a
5-dimensional 'evolutionary process space' characterised by individual-level statistics.

  populations  ↔  cells (single-cell analogy)
  stat vectors ↔  expression vectors
  trajectory   ↔  differentiation trajectory

Features (per individual):
  1. froh_genome          — genome-wide FROH (sum ROH length / autosome length)
  2. n_roh_genome         — total ROH count across 22 autosomes
  3. het_rate_chr22       — heterozygosity rate at chr22 SNPs (hap1 ≠ hap2)
  4. rare_burden_chr22    — count of chr22 SNPs carried with global AF < 0.01
  5. rare_missense_chr22  — count of chr22 rare missense SNPs carried (has_missense + AF < 0.01)

Note: sweep_burden was not used because chr22 has no H12 GenomicWindow nodes
(Garud H scan did not cover chr22). rare_missense captures functional rare burden instead.

Data sources:
  - human_interpretation_results.json → roh_hmm (genome-wide ROH per individual)
  - data/hap_cache/chr22.npz + data/variants_cache/chr22.npz (chr22 haplotypes)
  - data/hap_cache/sample_meta.json (sample order / population assignment)
  - data/results/evo_trajectory.json (population-level PCs for panel F comparison)

Output:
  data/results/individual_trajectory.tsv / .json
  data/results/figures/fig20_individual_trajectory.png

Usage:
  conda run -n graphevo python -u scripts/investigate_16b_individual_trajectory.py
"""

import csv
import json
import time
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

# ── Config ────────────────────────────────────────────────────────────────────
RESULTS      = Path("data/results")
FIG_DIR      = RESULTS / "figures"
HAP_CACHE    = Path("data/hap_cache")
VAR_CACHE    = Path("data/variants_cache")
CHECKPOINT   = Path("human_interpretation_results.json")
OUT_TSV      = RESULTS / "individual_trajectory.tsv"
OUT_JSON     = RESULTS / "individual_trajectory.json"
CHUNK        = 50_000       # variants per hap-unpack chunk
RARE_AF_THRESH = 0.01       # global AF threshold for rare variant burden

# hg38 autosome lengths (bp) — from fix_froh.py
CHR_LENGTHS_HG38 = {
    "chr1":  248956422, "chr2":  242193529, "chr3":  198295559,
    "chr4":  190214555, "chr5":  181538259, "chr6":  170805979,
    "chr7":  159345973, "chr8":  145138636, "chr9":  138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16":  90338345, "chr17":  83257441, "chr18":  80373285,
    "chr19":  58617616, "chr20":  64444167, "chr21":  46709983,
    "chr22":  50818468,
}
AUTOSOME_LENGTH_BP = sum(CHR_LENGTHS_HG38.values())   # 2,875,001,001 bp

# Continental group colors (match Inv.16)
GROUP_COLORS = {
    "African":     "#E05C5C",
    "European":    "#5C8AE0",
    "South_Asian": "#E0A85C",
    "East_Asian":  "#5CC45C",
    "American":    "#A85CE0",
}

POP_GROUP = {
    "ACB": "African",  "ASW": "African",  "ESN": "African",  "GWD": "African",
    "LWK": "African",  "MSL": "African",  "YRI": "African",
    "CEU": "European", "FIN": "European", "GBR": "European", "IBS": "European",
    "TSI": "European",
    "BEB": "South_Asian", "GIH": "South_Asian", "ITU": "South_Asian",
    "PJL": "South_Asian", "STU": "South_Asian",
    "CDX": "East_Asian",  "CHB": "East_Asian",  "CHS": "East_Asian",
    "JPT": "East_Asian",  "KHV": "East_Asian",
    "CLM": "American",    "MXL": "American",    "PEL": "American",
    "PUR": "American",
}

FEATURE_NAMES = [
    "froh_genome", "n_roh_genome", "het_rate_chr22",
    "rare_burden_chr22", "rare_missense_chr22",
]

VELOCITY_COMPONENTS = {
    "drift_inbreeding": {
        "label": "Drift / Inbreeding",
        "color": "#8B4513",
        "features": {"froh_genome": +1.0, "n_roh_genome": +0.5, "het_rate_chr22": -1.0},
    },
    "functional_burden": {
        "label": "Functional Rare Burden",
        "color": "#2E8B57",
        "features": {"rare_missense_chr22": +1.0, "rare_burden_chr22": +0.5},
    },
    "population_diversity": {
        "label": "Population Diversity",
        "color": "#4169E1",
        "features": {"het_rate_chr22": +1.0, "froh_genome": -1.0},
    },
}


# ── Step 1: Genome-wide FROH from checkpoint ──────────────────────────────────

def load_froh_from_checkpoint(checkpoint_path: Path) -> dict:
    """Return {sampleId: {froh_genome, n_roh_genome}} from roh_hmm per-sample data."""
    print("Loading roh_hmm from checkpoint...", flush=True)
    t0 = time.time()
    with open(checkpoint_path) as f:
        data = json.load(f)

    roh_hmm = data.get("roh_hmm", {})
    if not roh_hmm:
        raise RuntimeError("roh_hmm key not found in checkpoint")

    # Accumulate per-sample across all pop|chr entries
    sample_total_bp = {}   # sampleId → total ROH length (bp)
    sample_n_roh    = {}   # sampleId → total ROH count

    for key, entry in roh_hmm.items():
        pop, chrom = key.split("|")
        if chrom not in CHR_LENGTHS_HG38:
            continue   # skip sex chromosomes if any
        for s in entry.get("per_sample", []):
            sid = s["sampleId"]
            sample_total_bp[sid] = sample_total_bp.get(sid, 0) + s["total_length"]
            sample_n_roh[sid]    = sample_n_roh.get(sid, 0)    + s["n_roh"]

    samples = {
        sid: {
            "froh_genome":  bp / AUTOSOME_LENGTH_BP,
            "n_roh_genome": sample_n_roh[sid],
        }
        for sid, bp in sample_total_bp.items()
    }

    print(f"  {len(samples)} samples with ROH data  ({time.time()-t0:.1f}s)", flush=True)
    return samples


# ── Step 2: Chunked haplotype computation on chr22 ────────────────────────────

def compute_individual_features_chr22(
    hap_cache_path: Path,
    var_cache_path: Path,
    n_samples: int = 3202,
) -> dict:
    """Compute het_rate, rare_burden, rare_missense per individual from chr22 SNPs."""
    print("Loading chr22 caches...", flush=True)
    hap_data = np.load(hap_cache_path)
    var_data  = np.load(var_cache_path)

    hap_packed = hap_data["hap"]           # (M, 801)
    vtype      = var_data["variant_type"]  # (M,) int8: 1=SNP
    ac_all     = var_data["ac"]            # (M, 26) int32
    an_all     = var_data["an"]            # (M, 26) int32
    mis_all    = var_data["has_missense"]  # (M,) bool

    snp_mask = vtype == 1
    n_snps   = snp_mask.sum()
    print(f"  chr22 SNPs: {n_snps:,} / {len(vtype):,} total variants", flush=True)

    hap_snp    = hap_packed[snp_mask]
    ac_snps    = ac_all[snp_mask].astype(np.int64)
    an_snps    = an_all[snp_mask].astype(np.int64)
    mis_snps   = mis_all[snp_mask]              # missense mask on SNPs

    # Global AF across all 26 populations
    global_ac = ac_snps.sum(axis=1)
    global_an = an_snps.sum(axis=1).clip(1)
    global_af = global_ac / global_an

    rare_mask          = global_af < RARE_AF_THRESH          # (n_snps,)
    rare_missense_mask = rare_mask & mis_snps                # (n_snps,)

    n_rare_sites     = rare_mask.sum()
    n_rare_mis_sites = rare_missense_mask.sum()
    print(f"  Rare SNPs (AF<{RARE_AF_THRESH}): {n_rare_sites:,}", flush=True)
    print(f"  Rare missense SNPs: {n_rare_mis_sites:,}", flush=True)

    het_count         = np.zeros(n_samples, dtype=np.int64)
    rare_count        = np.zeros(n_samples, dtype=np.int64)
    rare_missense_cnt = np.zeros(n_samples, dtype=np.int64)
    n_haplotypes = 2 * n_samples

    n_chunks = (n_snps + CHUNK - 1) // CHUNK
    print(f"  Processing {n_snps:,} SNPs in {n_chunks} chunks of {CHUNK:,}...", flush=True)
    t0 = time.time()

    for ci, chunk_start in enumerate(range(0, n_snps, CHUNK)):
        chunk_end = min(chunk_start + CHUNK, n_snps)

        # Unpack haplotypes: (chunk_size, 2*n_samples)
        chunk_hap = np.unpackbits(
            hap_snp[chunk_start:chunk_end], axis=1, bitorder="big"
        )[:, :n_haplotypes]

        h1 = chunk_hap[:, 0::2]   # (chunk, n_samples)
        h2 = chunk_hap[:, 1::2]

        het_count += (h1 != h2).sum(axis=0)

        alt_dose = h1.astype(np.int16) + h2    # 0/1/2

        cm_rare = rare_mask[chunk_start:chunk_end]
        if cm_rare.any():
            rare_count += (alt_dose[cm_rare] > 0).sum(axis=0)

        cm_mis = rare_missense_mask[chunk_start:chunk_end]
        if cm_mis.any():
            rare_missense_cnt += (alt_dose[cm_mis] > 0).sum(axis=0)

        if (ci + 1) % 5 == 0 or ci == n_chunks - 1:
            print(f"    chunk {ci+1}/{n_chunks}  {time.time()-t0:.1f}s", flush=True)

    het_rate_chr22      = het_count  / n_snps
    rare_burden_chr22   = rare_count.astype(float)
    rare_missense_chr22 = rare_missense_cnt.astype(float)

    print(f"  Het rate:       [{het_rate_chr22.min():.4f}, {het_rate_chr22.max():.4f}]", flush=True)
    print(f"  Rare burden:    [{rare_burden_chr22.min():.0f}, {rare_burden_chr22.max():.0f}]", flush=True)
    print(f"  Rare missense:  [{rare_missense_chr22.min():.0f}, {rare_missense_chr22.max():.0f}]", flush=True)

    return {
        "het_rate_chr22":      het_rate_chr22,
        "rare_burden_chr22":   rare_burden_chr22,
        "rare_missense_chr22": rare_missense_chr22,
    }


# ── Step 3: Build feature matrix ──────────────────────────────────────────────

def build_individual_matrix(
    sample_ids: list,
    froh_data: dict,
    chr22_features: dict,
) -> tuple:
    n = len(sample_ids)
    feat_idx = {f: i for i, f in enumerate(FEATURE_NAMES)}
    mat = np.full((n, len(FEATURE_NAMES)), np.nan)

    for fi, fname in enumerate(FEATURE_NAMES):
        if fname in chr22_features:
            mat[:, fi] = chr22_features[fname]

    for i, sid in enumerate(sample_ids):
        if sid in froh_data:
            mat[i, feat_idx["froh_genome"]]   = froh_data[sid]["froh_genome"]
            mat[i, feat_idx["n_roh_genome"]]  = froh_data[sid]["n_roh_genome"]

    # NaN → column median
    for j in range(mat.shape[1]):
        nan_mask = np.isnan(mat[:, j])
        if nan_mask.any():
            mat[nan_mask, j] = np.nanmedian(mat[:, j])

    scaler = StandardScaler()
    mat_z  = scaler.fit_transform(mat)
    return mat, mat_z, scaler


# ── Step 4: PCA + UMAP ────────────────────────────────────────────────────────

def compute_pca(mat_z: np.ndarray, n_components: int = 5):
    pca = PCA(n_components=min(n_components, mat_z.shape[1]))
    return pca.fit_transform(mat_z), pca


def compute_umap(mat_z: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.2):
    import umap as umap_lib
    reducer = umap_lib.UMAP(
        n_components=2, n_neighbors=n_neighbors, min_dist=min_dist,
        random_state=42, metric="euclidean",
    )
    return reducer.fit_transform(mat_z)


# ── Step 5: Velocity per individual ──────────────────────────────────────────

def compute_velocities(mat: np.ndarray) -> dict:
    feat_idx = {f: i for i, f in enumerate(FEATURE_NAMES)}
    velocities = {}
    for comp_name, comp_def in VELOCITY_COMPONENTS.items():
        scores = np.zeros(mat.shape[0])
        total_weight = 0.0
        for feat, weight in comp_def["features"].items():
            if feat in feat_idx:
                vals = mat[:, feat_idx[feat]]
                mu, sd = np.nanmean(vals), np.nanstd(vals)
                scores += weight * (vals - mu) / (sd + 1e-9)
                total_weight += abs(weight)
        if total_weight > 0:
            scores /= total_weight
        velocities[comp_name] = scores
    return velocities


# ── Visualization ─────────────────────────────────────────────────────────────

def draw_pop_ellipses(ax, coords2d, pop_ids, alpha=0.12):
    """Draw 1-sigma covariance ellipses per population."""
    for pop in sorted(set(pop_ids)):
        idx = [i for i, p in enumerate(pop_ids) if p == pop]
        if len(idx) < 5:
            continue
        xy   = coords2d[idx]
        mean = xy.mean(axis=0)
        cov  = np.cov(xy.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        order = eigvals.argsort()[::-1]
        eigvals, eigvecs = eigvals[order], eigvecs[:, order]
        angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))
        width, height = 2 * np.sqrt(np.maximum(eigvals, 0))
        color = GROUP_COLORS.get(POP_GROUP.get(pop, ""), "#888")
        ell = Ellipse(xy=mean, width=width, height=height, angle=angle,
                      facecolor=color, alpha=alpha, edgecolor=color,
                      linewidth=0.8, zorder=1)
        ax.add_patch(ell)


def make_figure(
    sample_ids, pop_ids, mat, mat_z, umap_coords, pca_coords, pca,
    velocities, pop_level_pca,
):
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(22, 16))
    fig.patch.set_facecolor("#0D0D0D")
    gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.35,
                          left=0.07, right=0.97, top=0.93, bottom=0.07)
    axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]
    for ax in axes:
        ax.set_facecolor("#161616")
        for spine in ax.spines.values():
            spine.set_color("#444")

    groups = [POP_GROUP.get(p, "Other") for p in pop_ids]
    colors = [GROUP_COLORS.get(g, "#888") for g in groups]
    feat_idx = {f: i for i, f in enumerate(FEATURE_NAMES)}
    ev = pca.explained_variance_ratio_

    # ── A: UMAP ───────────────────────────────────────────────────────────
    ax = axes[0]
    ax.set_title("A  Individual Process Space (UMAP, n=3,202)", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    draw_pop_ellipses(ax, umap_coords, pop_ids)
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
               c=colors, s=6, alpha=0.55, linewidths=0, zorder=2)
    for grp, col in GROUP_COLORS.items():
        ax.scatter([], [], c=col, s=30, label=grp.replace("_", " "))
    ax.legend(fontsize=7, loc="lower right", facecolor="#222",
              labelcolor="#CCC", framealpha=0.8)
    ax.set_xlabel("UMAP 1", color="#AAA", fontsize=8)
    ax.set_ylabel("UMAP 2", color="#AAA", fontsize=8)
    ax.tick_params(colors="#888", labelsize=7)

    # ── B: Statistical PCA ────────────────────────────────────────────────
    ax = axes[1]
    ax.set_title(f"B  Statistical PCA  (PC1={ev[0]*100:.1f}%, PC2={ev[1]*100:.1f}%)",
                 color="#EEE", fontsize=11, pad=6, loc="left")
    draw_pop_ellipses(ax, pca_coords[:, :2], pop_ids)
    ax.scatter(pca_coords[:, 0], pca_coords[:, 1],
               c=colors, s=6, alpha=0.55, linewidths=0, zorder=2)
    ax.set_xlabel(f"PC1 ({ev[0]*100:.1f}%)", color="#AAA", fontsize=8)
    ax.set_ylabel(f"PC2 ({ev[1]*100:.1f}%)", color="#AAA", fontsize=8)
    ax.axhline(0, color="#333", lw=0.5); ax.axvline(0, color="#333", lw=0.5)
    ax.tick_params(colors="#888", labelsize=7)

    # ── C: PC1 loadings ───────────────────────────────────────────────────
    ax = axes[2]
    ax.set_title("C  Stat-PC1 feature loadings", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    loadings = pca.components_[0]
    order = np.argsort(np.abs(loadings))[::-1]
    feat_labels = [FEATURE_NAMES[j].replace("_chr22", "").replace("_", "\n") for j in order]
    load_vals = loadings[order]
    bar_colors = ["#E05C5C" if v > 0 else "#5C8AE0" for v in load_vals]
    ax.barh(range(len(order)), load_vals[::-1], color=bar_colors[::-1], edgecolor="#333")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(feat_labels[::-1], fontsize=7, color="#CCC")
    ax.set_xlabel("PC1 loading", color="#AAA", fontsize=8)
    ax.axvline(0, color="#555", lw=0.8)
    ax.tick_params(axis="x", colors="#888", labelsize=7)

    # ── D: Violin plots by continental group ──────────────────────────────
    ax = axes[3]
    ax.set_title("D  FROH & Het Rate by Continental Group", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    group_order = ["African", "European", "South_Asian", "East_Asian", "American"]
    froh_by_grp = {g: [] for g in group_order}
    het_by_grp  = {g: [] for g in group_order}
    for i, pop in enumerate(pop_ids):
        grp = POP_GROUP.get(pop, "Other")
        if grp in froh_by_grp:
            froh_by_grp[grp].append(mat[i, feat_idx["froh_genome"]])
            het_by_grp[grp].append(mat[i, feat_idx["het_rate_chr22"]])

    for xi, grp in enumerate(group_order):
        col = GROUP_COLORS.get(grp, "#888")
        for offset, vals, color, lbl in [
            (-0.18, froh_by_grp[grp], col,     "FROH"),
            (+0.18, het_by_grp[grp],  "#CCCCCC", "Het rate"),
        ]:
            if len(vals) > 3:
                vp = ax.violinplot([vals], positions=[xi + offset], widths=0.3,
                                   showmedians=True, showextrema=False)
                for body in vp["bodies"]:
                    body.set_facecolor(color); body.set_alpha(0.5)
                    body.set_edgecolor("none")
                vp["cmedians"].set_color(color); vp["cmedians"].set_linewidth(2)

    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels([g.replace("_", "\n") for g in group_order], color="#CCC", fontsize=7)
    ax.tick_params(axis="y", colors="#888", labelsize=7)
    ax.set_ylabel("Value", color="#AAA", fontsize=8)
    ax.legend(
        handles=[mpatches.Patch(color="#888", alpha=0.6, label="FROH (genome-wide)"),
                 mpatches.Patch(color="#CCC", alpha=0.4, label="Het rate (chr22)")],
        fontsize=7, facecolor="#222", labelcolor="#CCC", framealpha=0.8,
    )

    # ── E: FROH vs Het rate scatter ───────────────────────────────────────
    ax = axes[4]
    froh_vals = mat[:, feat_idx["froh_genome"]]
    het_vals  = mat[:, feat_idx["het_rate_chr22"]]
    ax.scatter(froh_vals, het_vals, c=colors, s=6, alpha=0.5, linewidths=0, zorder=2)
    draw_pop_ellipses(ax, np.column_stack([froh_vals, het_vals]), pop_ids)
    from scipy.stats import spearmanr
    r, p = spearmanr(froh_vals, het_vals)
    ax.set_xlabel("FROH (genome-wide)", color="#AAA", fontsize=8)
    ax.set_ylabel("Het rate (chr22 SNPs)", color="#AAA", fontsize=8)
    ax.set_title(f"E  FROH vs Het Rate  (ρ={r:.3f}, p={p:.2e})",
                 color="#EEE", fontsize=10, pad=6, loc="left")
    ax.tick_params(colors="#888", labelsize=7)

    # ── F: Individual PC1 vs Population PC1 ──────────────────────────────
    ax = axes[5]
    if pop_level_pca:
        ind_pc1 = pca_coords[:, 0]
        pop_pc1 = np.array([pop_level_pca.get(p, {}).get("stat_pc1", np.nan)
                            for p in pop_ids])
        valid = ~np.isnan(pop_pc1)
        ax.scatter(pop_pc1[valid], ind_pc1[valid],
                   c=[colors[i] for i in np.where(valid)[0]],
                   s=6, alpha=0.5, linewidths=0, zorder=2)
        if valid.sum() > 10:
            r2, p2 = spearmanr(pop_pc1[valid], ind_pc1[valid])
            ax.set_title(f"F  Individual PC1 vs Pop Stat-PC1  (ρ={r2:.3f}, p={p2:.2e})",
                         color="#EEE", fontsize=10, pad=6, loc="left")
        ax.set_xlabel("Population Stat-PC1 (Inv.16)", color="#AAA", fontsize=8)
        ax.set_ylabel("Individual Stat-PC1", color="#AAA", fontsize=8)
        ax.axhline(0, color="#333", lw=0.5); ax.axvline(0, color="#333", lw=0.5)
    else:
        ax.text(0.5, 0.5, "Population PCA data not found",
                transform=ax.transAxes, ha="center", color="#888")
    ax.tick_params(colors="#888", labelsize=7)

    fig.suptitle(
        "Individual-Level Evolutionary Trajectory Map — 3,202 Human Samples in Process Space\n"
        "Dimensions: drift (FROH) · ROH count · diversity (het rate) · "
        "rare burden · rare missense burden (chr22)",
        color="#EEE", fontsize=12, y=0.97,
    )

    out_path = FIG_DIR / "fig20_individual_trajectory.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()
    print("=== Investigation 16b: Individual-Level Evolutionary Trajectory Map ===\n")

    # 1. Load genome-wide FROH from checkpoint
    froh_data = load_froh_from_checkpoint(CHECKPOINT)

    # 2. Sample metadata
    print("Loading sample metadata...", flush=True)
    meta = json.load(open(HAP_CACHE / "sample_meta.json"))
    sample_ids = meta["sample_ids"]
    pop_ids    = meta["pop_ids"]
    n_samples  = len(sample_ids)
    print(f"  {n_samples} samples, {len(set(pop_ids))} populations", flush=True)

    # 3. Chr22 individual features (chunked)
    print("\nComputing chr22 individual features...", flush=True)
    chr22_feats = compute_individual_features_chr22(
        HAP_CACHE / "chr22.npz",
        VAR_CACHE  / "chr22.npz",
        n_samples=n_samples,
    )

    # 4. Feature matrix
    print("\nBuilding feature matrix...", flush=True)
    mat, mat_z, _ = build_individual_matrix(sample_ids, froh_data, chr22_feats)
    n, d = mat.shape
    print(f"  Matrix: {n} samples × {d} features")
    for j, fname in enumerate(FEATURE_NAMES):
        v = mat[:, j]
        print(f"  {fname:25s}  mean={v.mean():.4f}  std={v.std():.4f}  "
              f"range=[{v.min():.4f}, {v.max():.4f}]")

    # 5. PCA
    print("\nComputing PCA...", flush=True)
    pca_coords, pca = compute_pca(mat_z)
    ev = pca.explained_variance_ratio_
    for i in range(len(ev)):
        top3 = sorted(range(d), key=lambda j: abs(pca.components_[i][j]), reverse=True)[:3]
        top3_str = ", ".join(f"{FEATURE_NAMES[j]}({pca.components_[i][j]:+.2f})" for j in top3)
        print(f"  PC{i+1}: {ev[i]*100:.1f}%  |  {top3_str}")

    # 6. UMAP
    print("\nComputing UMAP (n_neighbors=15, min_dist=0.2)...", flush=True)
    umap_coords = compute_umap(mat_z)
    print(f"  UMAP range: x=[{umap_coords[:,0].min():.2f},{umap_coords[:,0].max():.2f}]  "
          f"y=[{umap_coords[:,1].min():.2f},{umap_coords[:,1].max():.2f}]")

    # 7. Velocities
    print("\nComputing evolutionary velocities...", flush=True)
    velocities = compute_velocities(mat)
    for comp_name in VELOCITY_COMPONENTS:
        scores = velocities[comp_name]
        grp_means = {}
        for i, pop in enumerate(pop_ids):
            grp = POP_GROUP.get(pop, "Other")
            grp_means.setdefault(grp, []).append(scores[i])
        top3 = sorted(grp_means, key=lambda g: np.mean(grp_means[g]), reverse=True)[:3]
        bot3 = sorted(grp_means, key=lambda g: np.mean(grp_means[g]))[:3]
        print(f"  {comp_name}: high→ {', '.join(f'{g}({np.mean(grp_means[g]):+.3f})' for g in top3)}; "
              f"low→ {', '.join(f'{g}({np.mean(grp_means[g]):+.3f})' for g in bot3)}")

    # 8. Load population-level PCA from Inv.16 (panel F)
    pop_level_pca = {}
    inv16_path = RESULTS / "evo_trajectory.json"
    if inv16_path.exists():
        pop_level_pca = json.load(open(inv16_path)).get("embeddings", {})
        print(f"\nLoaded {len(pop_level_pca)} population PC coords from Inv.16")

    # 9. Figure
    print("\nGenerating figure...", flush=True)
    make_figure(sample_ids, pop_ids, mat, mat_z, umap_coords, pca_coords, pca,
                velocities, pop_level_pca)

    # 10. Write outputs
    print("\nWriting outputs...", flush=True)
    feat_idx = {f: i for i, f in enumerate(FEATURE_NAMES)}
    rows = []
    for i, sid in enumerate(sample_ids):
        pop = pop_ids[i]
        rows.append({
            "sampleId":           sid,
            "population":         pop,
            "group":              POP_GROUP.get(pop, ""),
            "froh_genome":        round(float(mat[i, feat_idx["froh_genome"]]),         6),
            "n_roh_genome":       int(mat[i, feat_idx["n_roh_genome"]]),
            "het_rate_chr22":     round(float(mat[i, feat_idx["het_rate_chr22"]]),       6),
            "rare_burden_chr22":  int(mat[i, feat_idx["rare_burden_chr22"]]),
            "rare_missense_chr22":int(mat[i, feat_idx["rare_missense_chr22"]]),
            "stat_pc1":           round(float(pca_coords[i, 0]), 4),
            "stat_pc2":           round(float(pca_coords[i, 1]), 4),
            "stat_pc3":           round(float(pca_coords[i, 2]), 4) if pca_coords.shape[1] > 2 else 0.0,
            "umap1":              round(float(umap_coords[i, 0]), 4),
            "umap2":              round(float(umap_coords[i, 1]), 4),
            "drift_inbreeding_score":   round(float(velocities["drift_inbreeding"][i]),   4),
            "functional_burden_score":  round(float(velocities["functional_burden"][i]),  4),
            "diversity_score":          round(float(velocities["population_diversity"][i]), 4),
        })

    fields = [
        "sampleId", "population", "group",
        "froh_genome", "n_roh_genome", "het_rate_chr22",
        "rare_burden_chr22", "rare_missense_chr22",
        "stat_pc1", "stat_pc2", "stat_pc3", "umap1", "umap2",
        "drift_inbreeding_score", "functional_burden_score", "diversity_score",
    ]
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Group summary stats
    group_stats = {}
    for grp in GROUP_COLORS:
        idx = [i for i, p in enumerate(pop_ids) if POP_GROUP.get(p) == grp]
        if not idx:
            continue
        group_stats[grp] = {
            feat: {"mean": round(float(mat[idx, fi].mean()), 6),
                   "std":  round(float(mat[idx, fi].std()),  6)}
            for feat, fi in feat_idx.items()
        }

    out_data = {
        "n_samples":              n,
        "n_features":             d,
        "features":               FEATURE_NAMES,
        "pca_explained_variance": [round(float(v), 4) for v in ev],
        "rare_af_threshold":      RARE_AF_THRESH,
        "group_statistics":       group_stats,
        "embeddings": {
            sid: {
                "population": pop_ids[i],
                "group":      POP_GROUP.get(pop_ids[i], ""),
                "stat_pc1":   round(float(pca_coords[i, 0]), 4),
                "stat_pc2":   round(float(pca_coords[i, 1]), 4),
                "umap1":      round(float(umap_coords[i, 0]), 4),
                "umap2":      round(float(umap_coords[i, 1]), 4),
            }
            for i, sid in enumerate(sample_ids)
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2)

    elapsed = time.time() - t_start
    print(f"\nOutputs: {OUT_TSV}")
    print(f"         {OUT_JSON}")
    print(f"Total time: {elapsed:.1f}s")
    print("\nVerification:")
    print(f"  TSV rows: {len(rows)}")
    for grp, stats in group_stats.items():
        print(f"  {grp:12s}  froh={stats['froh_genome']['mean']:.4f}  "
              f"het={stats['het_rate_chr22']['mean']:.4f}  "
              f"rare_mis={stats['rare_missense_chr22']['mean']:.1f}")


if __name__ == "__main__":
    main()
