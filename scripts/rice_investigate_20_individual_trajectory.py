#!/usr/bin/env python3
"""Rice Investigation 20: Individual-Level Evolutionary Trajectory Map.

Analogous to human Inv.16b. Places all 3,024 rice accessions in a
4-dimensional 'evolutionary process space' characterised by individual-level
statistics computed from the haplotype cache.

Rice is a predominantly self-pollinating species, so:
  - Heterozygosity rates will be very low (unlike human ~0.07)
  - Homozygosity rate will dominate (expected > 0.95 for most accessions)
  - Population structure should separate indica from japonica subgroups

Features (per individual, from Chr1 — the largest chromosome):
  1. het_rate       — fraction of heterozygous sites (hap1 != hap2)
  2. hom_alt_rate   — fraction of homozygous-alt sites (both haplotypes = 1)
  3. rare_burden    — count of rare alleles carried (global AF < 0.05)
  4. private_burden — count of private/ultra-rare alleles (global AF < 0.005)

Data sources:
  - data/rice_hap_cache/Chr1.npz (Chr1 haplotypes, 3M+ variants)
  - data/rice_hap_cache/sample_meta.json (sample order / population assignment)
  - rice_interpretation_results.json → roh per-sample data (sparse; supplemented by hom_alt_rate)

Output:
  data/results/rice_inv20_individual_trajectory.tsv / .json
  data/results/figures/rice_fig20_individual_trajectory.png

Usage:
  conda run -n graphevo python -u scripts/rice_investigate_20_individual_trajectory.py
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
HAP_CACHE    = Path("data/rice_hap_cache")
ROH_FILE     = Path("rice_interpretation_results.json")
OUT_TSV      = RESULTS / "rice_inv20_individual_trajectory.tsv"
OUT_JSON     = RESULTS / "rice_inv20_individual_trajectory.json"
CHUNK        = 50_000       # variants per hap-unpack chunk
RARE_AF_THRESH    = 0.05    # AF threshold for rare variant burden
PRIVATE_AF_THRESH = 0.005   # AF threshold for ultra-rare/private variants
CHR_TO_USE   = "Chr1"       # largest chromosome

# Rice chromosome lengths (IRGSP-1.0 / MSU7 reference)
CHR_LENGTHS_RICE = {
    "Chr1":  43270923, "Chr2":  35937250, "Chr3":  36413819,
    "Chr4":  35502694, "Chr5":  29958434, "Chr6":  31248787,
    "Chr7":  29697621, "Chr8":  28443022, "Chr9":  23012720,
    "Chr10": 23207287, "Chr11": 29021106, "Chr12": 27531856,
}
GENOME_LENGTH_BP = sum(CHR_LENGTHS_RICE.values())

# Population group mapping for rice 3K
# Major groups: Japonica (GJ-*), Indica (XI-*), cA-Aus, cB-Bas, admix
GROUP_COLORS = {
    "Japonica":  "#5C8AE0",    # blue
    "Indica":    "#E05C5C",    # red
    "Aus":       "#5CC45C",    # green
    "Basmati":   "#E0A85C",    # orange
    "Admixed":   "#A85CE0",    # purple
    "Unknown":   "#888888",    # gray
}

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

FEATURE_NAMES = [
    "het_rate", "hom_alt_rate", "rare_burden", "private_burden",
]


# ── Step 1: Load genome-wide FROH from ROH data (sparse) ─────────────────────

def load_froh_from_roh(roh_path: Path, sample_ids: list) -> np.ndarray:
    """Accumulate per-sample ROH across chromosomes from rice_interpretation_results.json.

    Returns froh array aligned to sample_ids (0.0 for samples with no ROH data).
    """
    print("Loading per-sample ROH data (sparse)...", flush=True)
    t0 = time.time()

    with open(roh_path) as f:
        data = json.load(f)

    roh = data.get("roh", {})
    sample_total_bp = {}
    sample_n_roh = {}

    for key, entry in roh.items():
        pop, chrom = key.split("|")
        if chrom not in CHR_LENGTHS_RICE:
            continue
        for s in entry.get("per_sample", []):
            sid = s["sampleId"]
            sample_total_bp[sid] = sample_total_bp.get(sid, 0) + s["total_length"]
            sample_n_roh[sid] = sample_n_roh.get(sid, 0) + s["n_roh"]

    # Build aligned arrays
    froh = np.zeros(len(sample_ids), dtype=np.float64)
    n_roh = np.zeros(len(sample_ids), dtype=np.int64)
    sid_to_idx = {sid: i for i, sid in enumerate(sample_ids)}

    matched = 0
    for sid, bp in sample_total_bp.items():
        if sid in sid_to_idx:
            froh[sid_to_idx[sid]] = bp / GENOME_LENGTH_BP
            n_roh[sid_to_idx[sid]] = sample_n_roh[sid]
            matched += 1

    print(f"  {matched} samples with ROH data out of {len(sample_ids)} total  ({time.time()-t0:.1f}s)",
          flush=True)
    print(f"  FROH range: [{froh.min():.6f}, {froh.max():.6f}]", flush=True)
    return froh, n_roh


# ── Step 2: Compute individual features from Chr1 hap cache ──────────────────

def compute_individual_features(
    hap_cache_path: Path,
    n_samples: int,
) -> dict:
    """Compute het_rate, hom_alt_rate, rare_burden, private_burden per individual."""
    print(f"Loading {CHR_TO_USE} haplotype cache...", flush=True)
    hap_data = np.load(hap_cache_path)
    hap_packed = hap_data["hap"]           # (M, ceil(2N/8))
    n_variants = hap_packed.shape[0]
    n_haplotypes = 2 * n_samples
    print(f"  {CHR_TO_USE}: {n_variants:,} variants, {n_samples:,} samples", flush=True)

    # First pass: compute global AF for rare variant classification
    # (need to unpack to count alt alleles per variant)
    print("  Computing global allele frequencies...", flush=True)
    t0 = time.time()

    global_ac = np.zeros(n_variants, dtype=np.int64)
    n_chunks = (n_variants + CHUNK - 1) // CHUNK

    for ci, cs in enumerate(range(0, n_variants, CHUNK)):
        ce = min(cs + CHUNK, n_variants)
        chunk_hap = np.unpackbits(
            hap_packed[cs:ce], axis=1, bitorder="big"
        )[:, :n_haplotypes]
        global_ac[cs:ce] = chunk_hap.sum(axis=1)

    global_af = global_ac / n_haplotypes
    print(f"  Global AF computed in {time.time()-t0:.1f}s", flush=True)
    print(f"  AF range: [{global_af.min():.6f}, {global_af.max():.6f}]", flush=True)
    print(f"  Monomorphic (AF=0): {(global_af == 0).sum():,}", flush=True)
    print(f"  Fixed alt (AF=1):   {(global_af == 1.0).sum():,}", flush=True)

    rare_mask    = (global_af > 0) & (global_af < RARE_AF_THRESH)
    private_mask = (global_af > 0) & (global_af < PRIVATE_AF_THRESH)
    print(f"  Rare variants (0 < AF < {RARE_AF_THRESH}): {rare_mask.sum():,}", flush=True)
    print(f"  Private variants (0 < AF < {PRIVATE_AF_THRESH}): {private_mask.sum():,}", flush=True)

    # Second pass: compute per-individual features
    print(f"  Computing per-individual features ({n_chunks} chunks)...", flush=True)
    het_count     = np.zeros(n_samples, dtype=np.int64)
    hom_alt_count = np.zeros(n_samples, dtype=np.int64)
    rare_count    = np.zeros(n_samples, dtype=np.int64)
    private_count = np.zeros(n_samples, dtype=np.int64)
    t0 = time.time()

    # Only count polymorphic sites for het/hom_alt rates
    poly_mask = (global_af > 0) & (global_af < 1.0)
    n_poly = poly_mask.sum()
    print(f"  Polymorphic variants: {n_poly:,}", flush=True)

    for ci, cs in enumerate(range(0, n_variants, CHUNK)):
        ce = min(cs + CHUNK, n_variants)
        chunk_hap = np.unpackbits(
            hap_packed[cs:ce], axis=1, bitorder="big"
        )[:, :n_haplotypes]

        h1 = chunk_hap[:, 0::2]   # (chunk, n_samples)
        h2 = chunk_hap[:, 1::2]

        # Only count at polymorphic sites
        cm_poly = poly_mask[cs:ce]
        if cm_poly.any():
            het_count     += ((h1[cm_poly] != h2[cm_poly])).sum(axis=0)
            hom_alt_count += ((h1[cm_poly] == 1) & (h2[cm_poly] == 1)).sum(axis=0)

        # Rare burden: count sites where individual carries at least one alt allele
        alt_dose = h1.astype(np.int16) + h2
        cm_rare = rare_mask[cs:ce]
        if cm_rare.any():
            rare_count += (alt_dose[cm_rare] > 0).sum(axis=0)

        cm_priv = private_mask[cs:ce]
        if cm_priv.any():
            private_count += (alt_dose[cm_priv] > 0).sum(axis=0)

        if (ci + 1) % 10 == 0 or ci == n_chunks - 1:
            print(f"    chunk {ci+1}/{n_chunks}  {time.time()-t0:.1f}s", flush=True)

    het_rate     = het_count / max(n_poly, 1)
    hom_alt_rate = hom_alt_count / max(n_poly, 1)

    print(f"  Het rate:       [{het_rate.min():.6f}, {het_rate.max():.6f}]  "
          f"mean={het_rate.mean():.6f}", flush=True)
    print(f"  Hom-alt rate:   [{hom_alt_rate.min():.6f}, {hom_alt_rate.max():.6f}]  "
          f"mean={hom_alt_rate.mean():.6f}", flush=True)
    print(f"  Rare burden:    [{rare_count.min()}, {rare_count.max()}]  "
          f"mean={rare_count.mean():.1f}", flush=True)
    print(f"  Private burden: [{private_count.min()}, {private_count.max()}]  "
          f"mean={private_count.mean():.1f}", flush=True)

    return {
        "het_rate":       het_rate,
        "hom_alt_rate":   hom_alt_rate,
        "rare_burden":    rare_count.astype(float),
        "private_burden": private_count.astype(float),
    }


# ── Step 3: Build feature matrix ─────────────────────────────────────────────

def build_individual_matrix(features: dict) -> tuple:
    n = len(features["het_rate"])
    mat = np.column_stack([features[f] for f in FEATURE_NAMES])
    assert mat.shape == (n, len(FEATURE_NAMES)), f"Shape mismatch: {mat.shape}"

    # NaN → column median (shouldn't happen but safety)
    for j in range(mat.shape[1]):
        nan_mask = np.isnan(mat[:, j])
        if nan_mask.any():
            mat[nan_mask, j] = np.nanmedian(mat[:, j])

    scaler = StandardScaler()
    mat_z = scaler.fit_transform(mat)
    return mat, mat_z, scaler


# ── Step 4: PCA + UMAP ───────────────────────────────────────────────────────

def compute_pca(mat_z: np.ndarray, n_components: int = 4):
    pca = PCA(n_components=min(n_components, mat_z.shape[1]))
    return pca.fit_transform(mat_z), pca


def compute_umap(mat_z: np.ndarray, n_neighbors: int = 15, min_dist: float = 0.2):
    import umap as umap_lib
    reducer = umap_lib.UMAP(
        n_components=2, n_neighbors=n_neighbors, min_dist=min_dist,
        random_state=42, metric="euclidean",
    )
    return reducer.fit_transform(mat_z)


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
        grp = POP_GROUP.get(pop, "Unknown")
        color = GROUP_COLORS.get(grp, "#888")
        ell = Ellipse(xy=mean, width=width, height=height, angle=angle,
                      facecolor=color, alpha=alpha, edgecolor=color,
                      linewidth=0.8, zorder=1)
        ax.add_patch(ell)


def make_figure(
    sample_ids, pop_ids, mat, mat_z, umap_coords, pca_coords, pca, froh,
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

    groups = [POP_GROUP.get(p, "Unknown") for p in pop_ids]
    colors = [GROUP_COLORS.get(g, "#888") for g in groups]
    feat_idx = {f: i for i, f in enumerate(FEATURE_NAMES)}
    ev = pca.explained_variance_ratio_

    # ── A: UMAP ───────────────────────────────────────────────────────────
    ax = axes[0]
    ax.set_title("A  Individual Process Space (UMAP, n=3,024)", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    draw_pop_ellipses(ax, umap_coords, pop_ids)
    ax.scatter(umap_coords[:, 0], umap_coords[:, 1],
               c=colors, s=6, alpha=0.55, linewidths=0, zorder=2)
    for grp, col in GROUP_COLORS.items():
        ax.scatter([], [], c=col, s=30, label=grp)
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
    feat_labels = [FEATURE_NAMES[j].replace("_", "\n") for j in order]
    load_vals = loadings[order]
    bar_colors = ["#E05C5C" if v > 0 else "#5C8AE0" for v in load_vals]
    ax.barh(range(len(order)), load_vals[::-1], color=bar_colors[::-1], edgecolor="#333")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(feat_labels[::-1], fontsize=7, color="#CCC")
    ax.set_xlabel("PC1 loading", color="#AAA", fontsize=8)
    ax.axvline(0, color="#555", lw=0.8)
    ax.tick_params(axis="x", colors="#888", labelsize=7)

    # ── D: Violin plots by major group ────────────────────────────────────
    ax = axes[3]
    ax.set_title("D  Het Rate & Hom-Alt Rate by Group", color="#EEE",
                 fontsize=11, pad=6, loc="left")
    group_order = ["Japonica", "Indica", "Aus", "Basmati", "Admixed"]
    het_by_grp   = {g: [] for g in group_order}
    hom_by_grp   = {g: [] for g in group_order}
    for i, pop in enumerate(pop_ids):
        grp = POP_GROUP.get(pop, "Unknown")
        if grp in het_by_grp:
            het_by_grp[grp].append(mat[i, feat_idx["het_rate"]])
            hom_by_grp[grp].append(mat[i, feat_idx["hom_alt_rate"]])

    for xi, grp in enumerate(group_order):
        col = GROUP_COLORS.get(grp, "#888")
        for offset, vals, color, lbl in [
            (-0.18, het_by_grp[grp], col,      "Het rate"),
            (+0.18, hom_by_grp[grp], "#CCCCCC", "Hom-alt rate"),
        ]:
            if len(vals) > 3:
                vp = ax.violinplot([vals], positions=[xi + offset], widths=0.3,
                                   showmedians=True, showextrema=False)
                for body in vp["bodies"]:
                    body.set_facecolor(color); body.set_alpha(0.5)
                    body.set_edgecolor("none")
                vp["cmedians"].set_color(color); vp["cmedians"].set_linewidth(2)

    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, color="#CCC", fontsize=7)
    ax.tick_params(axis="y", colors="#888", labelsize=7)
    ax.set_ylabel("Rate", color="#AAA", fontsize=8)
    ax.legend(
        handles=[mpatches.Patch(color="#888", alpha=0.6, label="Het rate (Chr1)"),
                 mpatches.Patch(color="#CCC", alpha=0.4, label="Hom-alt rate (Chr1)")],
        fontsize=7, facecolor="#222", labelcolor="#CCC", framealpha=0.8,
    )

    # ── E: Het rate vs Hom-alt rate scatter ───────────────────────────────
    ax = axes[4]
    het_vals = mat[:, feat_idx["het_rate"]]
    hom_vals = mat[:, feat_idx["hom_alt_rate"]]
    ax.scatter(hom_vals, het_vals, c=colors, s=6, alpha=0.5, linewidths=0, zorder=2)
    draw_pop_ellipses(ax, np.column_stack([hom_vals, het_vals]), pop_ids)
    from scipy.stats import spearmanr
    r, p = spearmanr(hom_vals, het_vals)
    ax.set_xlabel("Hom-alt rate (Chr1)", color="#AAA", fontsize=8)
    ax.set_ylabel("Het rate (Chr1)", color="#AAA", fontsize=8)
    ax.set_title(f"E  Hom-Alt vs Het Rate  (rho={r:.3f}, p={p:.2e})",
                 color="#EEE", fontsize=10, pad=6, loc="left")
    ax.tick_params(colors="#888", labelsize=7)

    # ── F: Rare burden vs Private burden scatter ──────────────────────────
    ax = axes[5]
    rare_vals = mat[:, feat_idx["rare_burden"]]
    priv_vals = mat[:, feat_idx["private_burden"]]
    ax.scatter(rare_vals, priv_vals, c=colors, s=6, alpha=0.5, linewidths=0, zorder=2)
    r2, p2 = spearmanr(rare_vals, priv_vals)
    ax.set_xlabel("Rare burden (AF<0.05, Chr1)", color="#AAA", fontsize=8)
    ax.set_ylabel("Private burden (AF<0.005, Chr1)", color="#AAA", fontsize=8)
    ax.set_title(f"F  Rare vs Private Burden  (rho={r2:.3f}, p={p2:.2e})",
                 color="#EEE", fontsize=10, pad=6, loc="left")
    ax.tick_params(colors="#888", labelsize=7)

    fig.suptitle(
        "Individual-Level Evolutionary Trajectory Map -- 3,024 Rice Accessions in Process Space\n"
        f"Dimensions: heterozygosity . homozygosity . "
        f"rare burden (AF<{RARE_AF_THRESH}) . private burden (AF<{PRIVATE_AF_THRESH})  [{CHR_TO_USE}]",
        color="#EEE", fontsize=12, y=0.97,
    )

    out_path = FIG_DIR / "rice_fig20_individual_trajectory.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    t_start = time.time()
    print("=== Rice Investigation 20: Individual-Level Evolutionary Trajectory Map ===\n")

    # 1. Sample metadata
    print("Loading sample metadata...", flush=True)
    meta = json.load(open(HAP_CACHE / "sample_meta.json"))
    sample_ids = meta["sample_ids"]
    pop_ids    = meta["pop_ids"]
    n_samples  = len(sample_ids)
    print(f"  {n_samples} samples, {len(set(pop_ids))} populations", flush=True)
    for grp in sorted(set(POP_GROUP.get(p, "Unknown") for p in pop_ids)):
        cnt = sum(1 for p in pop_ids if POP_GROUP.get(p, "Unknown") == grp)
        print(f"    {grp}: {cnt} samples")

    # 2. Load genome-wide FROH (sparse — most rice accessions have no ROH detected)
    froh, n_roh = load_froh_from_roh(ROH_FILE, sample_ids)

    # 3. Individual features from Chr1
    print(f"\nComputing {CHR_TO_USE} individual features...", flush=True)
    features = compute_individual_features(
        HAP_CACHE / f"{CHR_TO_USE}.npz",
        n_samples=n_samples,
    )

    # 4. Feature matrix
    print("\nBuilding feature matrix...", flush=True)
    mat, mat_z, _ = build_individual_matrix(features)
    n, d = mat.shape
    print(f"  Matrix: {n} samples x {d} features")
    for j, fname in enumerate(FEATURE_NAMES):
        v = mat[:, j]
        print(f"  {fname:20s}  mean={v.mean():.6f}  std={v.std():.6f}  "
              f"range=[{v.min():.6f}, {v.max():.6f}]")

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

    # 7. Figure
    print("\nGenerating figure...", flush=True)
    make_figure(sample_ids, pop_ids, mat, mat_z, umap_coords, pca_coords, pca, froh)

    # 8. Write outputs
    print("\nWriting outputs...", flush=True)
    feat_idx = {f: i for i, f in enumerate(FEATURE_NAMES)}
    rows = []
    for i, sid in enumerate(sample_ids):
        pop = pop_ids[i]
        grp = POP_GROUP.get(pop, "Unknown")
        rows.append({
            "sampleId":        sid,
            "population":      pop,
            "group":           grp,
            "het_rate":        round(float(mat[i, feat_idx["het_rate"]]),     6),
            "hom_alt_rate":    round(float(mat[i, feat_idx["hom_alt_rate"]]), 6),
            "rare_burden":     int(mat[i, feat_idx["rare_burden"]]),
            "private_burden":  int(mat[i, feat_idx["private_burden"]]),
            "froh_genome":     round(float(froh[i]), 6),
            "n_roh_genome":    int(n_roh[i]),
            "stat_pc1":        round(float(pca_coords[i, 0]), 4),
            "stat_pc2":        round(float(pca_coords[i, 1]), 4),
            "stat_pc3":        round(float(pca_coords[i, 2]), 4) if pca_coords.shape[1] > 2 else 0.0,
            "umap1":           round(float(umap_coords[i, 0]), 4),
            "umap2":           round(float(umap_coords[i, 1]), 4),
        })

    fields = [
        "sampleId", "population", "group",
        "het_rate", "hom_alt_rate", "rare_burden", "private_burden",
        "froh_genome", "n_roh_genome",
        "stat_pc1", "stat_pc2", "stat_pc3", "umap1", "umap2",
    ]
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Group summary stats
    group_stats = {}
    for grp in GROUP_COLORS:
        idx = [i for i, p in enumerate(pop_ids) if POP_GROUP.get(p, "Unknown") == grp]
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
        "chromosome":             CHR_TO_USE,
        "rare_af_threshold":      RARE_AF_THRESH,
        "private_af_threshold":   PRIVATE_AF_THRESH,
        "pca_explained_variance": [round(float(v), 4) for v in ev],
        "group_statistics":       group_stats,
        "embeddings": {
            sid: {
                "population": pop_ids[i],
                "group":      POP_GROUP.get(pop_ids[i], "Unknown"),
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
    for grp, stats in sorted(group_stats.items()):
        print(f"  {grp:12s}  het={stats['het_rate']['mean']:.6f}  "
              f"hom_alt={stats['hom_alt_rate']['mean']:.6f}  "
              f"rare={stats['rare_burden']['mean']:.1f}  "
              f"private={stats['private_burden']['mean']:.1f}")


if __name__ == "__main__":
    main()
