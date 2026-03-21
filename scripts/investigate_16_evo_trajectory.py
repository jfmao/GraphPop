#!/usr/bin/env python3
"""Investigation 16: Evolutionary Trajectory Map.

Question: Can we reveal the evolutionary trajectory of 26 human populations
by embedding them in a high-dimensional space of population genetic statistics?

Concept:
  Each population is characterized by ~17 statistics spanning diversity,
  divergence, selection, drift, inbreeding, and rare-variant burden.
  Unlike standard genotype PCA (which captures drift/distance), this
  'process space' captures the EVOLUTIONARY FORCES acting on each population.

  Inspired by single-cell trajectory inference:
    populations ↔ cells
    statistic vectors ↔ gene expression vectors
    evolutionary trajectory ↔ differentiation trajectory
    evolutionary velocity ↔ RNA velocity

Method:
  1. Build feature matrix: 26 pops × 17 statistics (z-scored)
  2. PCA for linear structure (interpretable axes)
  3. UMAP for nonlinear manifold (n=26 is small → careful tuning)
  4. Minimum spanning tree in full-dimensional space → trajectory backbone
  5. Evolutionary velocity vectors: deviation from neutral equilibrium in
     3 force components: drift/bottleneck, recent selection, expansion
  6. Compare statistical trajectory with genotype PCA (OoA gradient)

Output:
  data/results/evo_trajectory.tsv / .json
  data/results/figures/fig19_evo_trajectory.png

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_16_evo_trajectory.py
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
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.sparse import csr_matrix
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

# ── Config ────────────────────────────────────────────────────────────────────
RESULTS  = Path("data/results")
FIG_DIR  = RESULTS / "figures"
OUT_TSV  = RESULTS / "evo_trajectory.tsv"
OUT_JSON = RESULTS / "evo_trajectory.json"

# Continental group → color
GROUP_COLORS = {
    "African":     "#E05C5C",
    "European":    "#5C8AE0",
    "South_Asian": "#E0A85C",
    "East_Asian":  "#5CC45C",
    "American":    "#A85CE0",
}

# ── Feature definitions ────────────────────────────────────────────────────────
# (feature_name, data_source, interpretation_direction)
# direction: +1 means higher = more of this process, -1 means inverted
FEATURE_DEFS = [
    # Diversity
    ("mean_pi",         "roh_sweep",   +1, "diversity: long-term Ne"),
    ("mean_theta_w",    "roh_sweep",   +1, "diversity: segregating sites"),
    ("mean_tajima_d",   "roh_sweep",   -1, "expansion: more negative → recent expansion"),  # inverted
    ("mean_fis",        "roh_sweep",   -1, "inbreeding: FIS (inverted)"),
    # Divergence
    ("mean_fst",        "roh_sweep",   +1, "divergence: mean FST from other pops"),
    # Selection (sweep)
    ("n_hard_sweeps",   "roh_sweep",   +1, "selection: hard sweep count (genome-wide)"),
    ("mean_h12",        "roh_sweep",   +1, "selection: mean Garud H12"),
    ("h12_fraction",    "temporal",    +1, "selection: fraction windows H12>threshold"),
    ("h12_z",           "temporal",    +1, "selection: H12 enrichment z-score"),
    ("ihs_fraction",    "temporal",    +1, "selection: fraction variants |iHS|>2  (0 for AFR)"),
    ("ihs_z",           "temporal",    +1, "selection: iHS enrichment z-score (0 for AFR)"),
    # Inbreeding / bottleneck
    ("froh",            "temporal",    +1, "inbreeding: genome-wide ROH fraction"),
    ("froh_z",          "temporal",    +1, "bottleneck: FROH z-score vs expectation"),
    # Rare variants
    ("mean_rare",       "rare_burden", +1, "rare burden: private rare variants per gene"),
    # Geography (proxy for OoA drift)
    ("dist_ooa_km",     "ooa",         +1, "drift: OoA distance (serial founder effect)"),
    # OoA genetic PCA — for comparison
    ("pc1",             "ooa",         +1, "genetic PC1 (drift structure)"),
    ("pc2",             "ooa",         +1, "genetic PC2 (African substructure)"),
]

FEATURE_NAMES = [f[0] for f in FEATURE_DEFS]
FEATURE_LABELS = {f[0]: f[3] for f in FEATURE_DEFS}


# ── Velocity component definitions ────────────────────────────────────────────
# Each component is a weighted sum of features; direction = sign of weight
VELOCITY_COMPONENTS = {
    "drift_bottleneck": {
        "label": "Drift / Bottleneck",
        "color": "#8B4513",
        "features": {"froh_z": +1.0, "mean_rare": +0.5,
                     "dist_ooa_km": +0.5, "mean_pi": -1.0},
    },
    "recent_selection": {
        "label": "Recent Selection",
        "color": "#2E8B57",
        "features": {"ihs_fraction": +1.0, "ihs_z": +0.8,
                     "h12_fraction": +0.5, "n_hard_sweeps": +0.3},
    },
    "expansion": {
        "label": "Demographic Expansion",
        "color": "#4169E1",
        "features": {"mean_tajima_d": -1.0, "mean_pi": +0.5,
                     "froh": -0.5},
    },
}


# ── Load data ─────────────────────────────────────────────────────────────────

def load_features() -> dict:
    """Load all features into a pop → {feature: value} dict."""
    pops = {}

    # roh_sweep_correlation
    with open(RESULTS / "roh_sweep_correlation.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pop = row["pop"]
            pops.setdefault(pop, {"pop": pop, "group": row["group"]})
            for feat in ["mean_pi", "mean_theta_w", "mean_tajima_d",
                         "mean_fis", "mean_fst", "n_hard_sweeps", "mean_h12"]:
                try:
                    pops[pop][feat] = float(row[feat])
                except (ValueError, KeyError):
                    pops[pop][feat] = np.nan

    # temporal_selection
    with open(RESULTS / "temporal_selection.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pop = row["pop"]
            pops.setdefault(pop, {"pop": pop, "group": row.get("group", "")})
            for feat in ["froh", "froh_z", "h12_fraction", "h12_z",
                         "ihs_fraction", "ihs_z"]:
                v = row.get(feat, "")
                try:
                    pops[pop][feat] = float(v) if v not in ("", "nan", "None") else np.nan
                except (ValueError, TypeError):
                    pops[pop][feat] = np.nan
            pops[pop]["dominant_timescale"] = row.get("dominant_timescale", "")

    # ooa_gradient (geographic + genetic PCA)
    with open(RESULTS / "ooa_gradient.tsv") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            pop = row["pop"]
            pops.setdefault(pop, {"pop": pop})
            for feat in ["dist_ooa_km", "pc1", "pc2", "lat", "lon"]:
                try:
                    pops[pop][feat] = float(row[feat])
                except (ValueError, KeyError):
                    pops[pop][feat] = np.nan

    # rare_burden (mean_rare per population)
    rb = json.load(open(RESULTS / "rare_burden.json"))
    for pop, stats in rb.get("pop_summary", {}).items():
        pops.setdefault(pop, {"pop": pop})
        pops[pop]["mean_rare"] = stats.get("mean_rare", np.nan)

    return pops


def build_matrix(pops: dict):
    """Build (n_pops × n_features) matrix, z-scored, NaN-imputed."""
    pop_list = sorted(pops.keys())
    mat = np.zeros((len(pop_list), len(FEATURE_NAMES)), dtype=float)

    for i, pop in enumerate(pop_list):
        for j, feat in enumerate(FEATURE_NAMES):
            mat[i, j] = pops[pop].get(feat, np.nan)

    # Impute NaN with column median (conservative; African iHS → 0 via median)
    col_medians = np.nanmedian(mat, axis=0)
    for j in range(mat.shape[1]):
        nan_mask = np.isnan(mat[:, j])
        mat[nan_mask, j] = col_medians[j]

    # Z-score
    scaler = StandardScaler()
    mat_z = scaler.fit_transform(mat)

    return pop_list, mat, mat_z, scaler


# ── Trajectory analysis ───────────────────────────────────────────────────────

def compute_mst(mat_z: np.ndarray, pop_list: list) -> list:
    """Compute minimum spanning tree on z-scored feature matrix."""
    dist_mat = squareform(pdist(mat_z, metric="euclidean"))
    sparse_dist = csr_matrix(dist_mat)
    mst = minimum_spanning_tree(sparse_dist)
    mst_dense = mst.toarray()

    edges = []
    rows, cols = np.nonzero(mst_dense)
    for r, c in zip(rows, cols):
        edges.append({
            "pop1": pop_list[r],
            "pop2": pop_list[c],
            "dist": float(mst_dense[r, c]),
        })
    return edges


def compute_pca(mat_z: np.ndarray, n_components: int = 6):
    """PCA on z-scored feature matrix."""
    pca = PCA(n_components=min(n_components, mat_z.shape[1]))
    coords = pca.fit_transform(mat_z)
    return coords, pca


def compute_umap(mat_z: np.ndarray, n_neighbors: int = 8, min_dist: float = 0.4):
    """UMAP embedding (tuned for small n=26)."""
    import umap as umap_lib
    reducer = umap_lib.UMAP(
        n_components=2,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        random_state=42,
        metric="euclidean",
    )
    coords = reducer.fit_transform(mat_z)
    return coords, reducer


def compute_velocities(pops: dict, pop_list: list, mat: np.ndarray) -> dict:
    """Compute evolutionary velocity scores per population per component."""
    feat_idx = {f: i for i, f in enumerate(FEATURE_NAMES)}
    velocities = {comp: {} for comp in VELOCITY_COMPONENTS}

    for comp_name, comp_def in VELOCITY_COMPONENTS.items():
        # Collect feature values (raw, not z-scored, for interpretability)
        scores = np.zeros(len(pop_list))
        total_weight = 0.0
        for feat, weight in comp_def["features"].items():
            if feat in feat_idx:
                vals = mat[:, feat_idx[feat]]
                # z-score this feature across pops before weighting
                mu, sd = np.nanmean(vals), np.nanstd(vals)
                z = (vals - mu) / (sd + 1e-9)
                scores += weight * z
                total_weight += abs(weight)
        scores /= total_weight

        for i, pop in enumerate(pop_list):
            velocities[comp_name][pop] = float(scores[i])

    return velocities


def compare_with_genetic_pca(mat_z, pop_list, pops):
    """Correlate each statistical PC with genetic PC1/PC2."""
    from scipy.stats import spearmanr
    pca_coords, pca = compute_pca(mat_z)
    results = []
    for i in range(min(4, pca_coords.shape[1])):
        stat_pc = pca_coords[:, i]
        for genetic_feat in ["pc1", "pc2", "dist_ooa_km"]:
            gen_vals = np.array([pops[p].get(genetic_feat, np.nan) for p in pop_list])
            mask = ~np.isnan(gen_vals)
            r, p = spearmanr(stat_pc[mask], gen_vals[mask])
            results.append({
                "stat_pc": f"PC{i+1} ({pca.explained_variance_ratio_[i]*100:.1f}%)",
                "genetic_feature": genetic_feat,
                "spearman_r": round(float(r), 3),
                "p_value": round(float(p), 4),
            })
    return results, pca


# ── Visualization ─────────────────────────────────────────────────────────────

def make_figure(pop_list, pops, mat_z, mat, mst_edges, umap_coords,
                pca_coords, pca, velocities):

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(20, 16))
    fig.patch.set_facecolor("#0D0D0D")

    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.35,
                          left=0.07, right=0.97, top=0.93, bottom=0.07)
    axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]
    for ax in axes:
        ax.set_facecolor("#161616")
        for spine in ax.spines.values():
            spine.set_color("#444")

    groups = [pops[p]["group"] for p in pop_list]
    colors = [GROUP_COLORS.get(g, "#888") for g in groups]

    # ── Panel A: UMAP + MST + velocity arrows ─────────────────────────────
    ax = axes[0]
    ax.set_title("A  Evolutionary Process Space (UMAP + MST)", color="#EEE",
                 fontsize=11, pad=6, loc="left")

    # Draw MST edges
    pop_to_umap = {p: umap_coords[i] for i, p in enumerate(pop_list)}
    for edge in mst_edges:
        x1, y1 = pop_to_umap[edge["pop1"]]
        x2, y2 = pop_to_umap[edge["pop2"]]
        ax.plot([x1, x2], [y1, y2], color="#444", lw=1.2, zorder=1, alpha=0.7)

    # Scatter
    for i, pop in enumerate(pop_list):
        x, y = umap_coords[i]
        ax.scatter(x, y, c=colors[i], s=80, zorder=3, edgecolors="white",
                   linewidths=0.5)
        ax.annotate(pop, (x, y), textcoords="offset points", xytext=(4, 4),
                    fontsize=6, color="#CCC", zorder=4)

    # Velocity arrows (dominant component per pop)
    # Project velocity components onto UMAP axes using PCA loading proxy
    comp_directions = {}
    for comp_name, comp_def in VELOCITY_COMPONENTS.items():
        # Direction in UMAP space: average position of top-quartile pops
        scores = np.array([velocities[comp_name][p] for p in pop_list])
        top_idx = np.argsort(scores)[-7:]  # top 7 pops
        bot_idx = np.argsort(scores)[:7]
        top_center = umap_coords[top_idx].mean(axis=0)
        bot_center = umap_coords[bot_idx].mean(axis=0)
        comp_directions[comp_name] = top_center - bot_center  # direction of increasing score

    for i, pop in enumerate(pop_list):
        x, y = umap_coords[i]
        arrow_x, arrow_y = 0.0, 0.0
        for comp_name, comp_def in VELOCITY_COMPONENTS.items():
            score = velocities[comp_name][pop]
            direction = comp_directions[comp_name]
            norm = np.linalg.norm(direction)
            if norm > 0:
                arrow_x += score * direction[0] / norm * 0.3
                arrow_y += score * direction[1] / norm * 0.3
        if abs(arrow_x) + abs(arrow_y) > 0.05:
            ax.annotate("", xy=(x + arrow_x, y + arrow_y), xytext=(x, y),
                        arrowprops=dict(arrowstyle="->", color="#FFD700",
                                        lw=1.2, alpha=0.7))

    # Legend
    for grp, col in GROUP_COLORS.items():
        ax.scatter([], [], c=col, s=40, label=grp.replace("_", " "))
    ax.legend(fontsize=7, loc="lower right", facecolor="#222",
              labelcolor="#CCC", framealpha=0.8)
    ax.set_xlabel("UMAP 1", color="#AAA", fontsize=8)
    ax.set_ylabel("UMAP 2", color="#AAA", fontsize=8)
    ax.tick_params(colors="#888", labelsize=7)

    # ── Panel B: Statistical PC1/PC2 (interpretable axes) ─────────────────
    ax = axes[1]
    ax.set_title("B  Statistical PCA (process axes)", color="#EEE",
                 fontsize=11, pad=6, loc="left")

    ev = pca.explained_variance_ratio_
    pop_to_pca = {p: pca_coords[i] for i, p in enumerate(pop_list)}

    # Draw MST in PCA space
    for edge in mst_edges:
        x1, y1 = pop_to_pca[edge["pop1"]][:2]
        x2, y2 = pop_to_pca[edge["pop2"]][:2]
        ax.plot([x1, x2], [y1, y2], color="#444", lw=1.0, zorder=1, alpha=0.6)

    for i, pop in enumerate(pop_list):
        x, y = pca_coords[i, :2]
        ax.scatter(x, y, c=colors[i], s=80, zorder=3, edgecolors="white",
                   linewidths=0.5)
        ax.annotate(pop, (x, y), textcoords="offset points", xytext=(4, 4),
                    fontsize=6, color="#CCC", zorder=4)

    ax.set_xlabel(f"Stat-PC1 ({ev[0]*100:.1f}% var)", color="#AAA", fontsize=8)
    ax.set_ylabel(f"Stat-PC2 ({ev[1]*100:.1f}% var)", color="#AAA", fontsize=8)
    ax.tick_params(colors="#888", labelsize=7)
    ax.axhline(0, color="#333", lw=0.5)
    ax.axvline(0, color="#333", lw=0.5)

    # ── Panel C: PC1 loadings (feature contributions) ─────────────────────
    ax = axes[2]
    ax.set_title("C  Stat-PC1 feature loadings", color="#EEE",
                 fontsize=11, pad=6, loc="left")

    loadings = pca.components_[0]
    sorted_idx = np.argsort(np.abs(loadings))[::-1][:12]
    feat_names_short = [FEATURE_NAMES[j].replace("mean_", "").replace("_", "\n")
                        for j in sorted_idx]
    load_vals = loadings[sorted_idx]
    bar_colors = ["#E05C5C" if v > 0 else "#5C8AE0" for v in load_vals]
    bars = ax.barh(range(len(sorted_idx)), load_vals[::-1], color=bar_colors[::-1],
                   edgecolor="#333")
    ax.set_yticks(range(len(sorted_idx)))
    ax.set_yticklabels(feat_names_short[::-1], fontsize=6, color="#CCC")
    ax.set_xlabel("PC1 loading", color="#AAA", fontsize=8)
    ax.axvline(0, color="#555", lw=0.8)
    ax.tick_params(axis="x", colors="#888", labelsize=7)

    # ── Panel D: Velocity heatmap (26 pops × 3 components) ────────────────
    ax = axes[3]
    ax.set_title("D  Evolutionary velocity (force decomposition)", color="#EEE",
                 fontsize=11, pad=6, loc="left")

    # Sort populations by drift/bottleneck score
    drift_order = sorted(pop_list, key=lambda p: velocities["drift_bottleneck"][p])
    comp_names = list(VELOCITY_COMPONENTS.keys())
    vel_mat = np.array([[velocities[c][p] for c in comp_names]
                         for p in drift_order])

    # Normalize each column to [-1, 1]
    for col in range(vel_mat.shape[1]):
        mx = max(abs(vel_mat[:, col].max()), abs(vel_mat[:, col].min()))
        if mx > 0:
            vel_mat[:, col] /= mx

    im = ax.imshow(vel_mat, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(comp_names)))
    ax.set_xticklabels([VELOCITY_COMPONENTS[c]["label"]
                         for c in comp_names], fontsize=7, color="#CCC", rotation=15)
    ax.set_yticks(range(len(drift_order)))
    ax.set_yticklabels(drift_order, fontsize=6, color="#CCC")
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.03,
                 label="Normalized score").ax.yaxis.set_tick_params(color="#888",
                                                                      labelcolor="#888")

    # Color pop labels by group
    pop_label_colors = [GROUP_COLORS.get(pops[p]["group"], "#888")
                         for p in drift_order]
    for label, col in zip(ax.get_yticklabels(), pop_label_colors):
        label.set_color(col)

    # ── Panel E: Statistical PC1 vs Genetic PC1 (process vs drift) ────────
    ax = axes[4]
    ax.set_title("E  Statistical PC1 vs Genetic PC1\n(process space vs drift)", color="#EEE",
                 fontsize=11, pad=6, loc="left")

    from scipy.stats import spearmanr
    stat_pc1 = pca_coords[:, 0]
    gen_pc1 = np.array([pops[p].get("pc1", np.nan) for p in pop_list])
    mask = ~np.isnan(gen_pc1)

    for i, pop in enumerate(pop_list):
        if mask[i]:
            ax.scatter(gen_pc1[i], stat_pc1[i], c=colors[i], s=80,
                       edgecolors="white", linewidths=0.5, zorder=3)
            ax.annotate(pop, (gen_pc1[i], stat_pc1[i]),
                        textcoords="offset points", xytext=(4, 3),
                        fontsize=6, color="#CCC", zorder=4)

    r, p = spearmanr(gen_pc1[mask], stat_pc1[mask])
    ax.set_xlabel("Genetic PC1 (drift / OoA)", color="#AAA", fontsize=8)
    ax.set_ylabel(f"Statistical PC1 ({ev[0]*100:.1f}% var)", color="#AAA", fontsize=8)
    ax.set_title(f"E  Genetic PC1 vs Statistical PC1  (ρ={r:.3f}, p={p:.3f})",
                 color="#EEE", fontsize=10, pad=6, loc="left")
    ax.axhline(0, color="#333", lw=0.5)
    ax.axvline(0, color="#333", lw=0.5)
    ax.tick_params(colors="#888", labelsize=7)

    # ── Panel F: Feature correlation heatmap (what drives co-variation?) ──
    ax = axes[5]
    ax.set_title("F  Feature correlation matrix", color="#EEE",
                 fontsize=11, pad=6, loc="left")

    # Spearman correlation between all features across 26 pops
    from scipy.stats import spearmanr as sr
    n_feat = len(FEATURE_NAMES)
    corr_mat = np.zeros((n_feat, n_feat))
    for i in range(n_feat):
        for j in range(i, n_feat):
            vi = mat[:, i]
            vj = mat[:, j]
            mask2 = ~(np.isnan(vi) | np.isnan(vj))
            if mask2.sum() >= 5:
                r2, _ = sr(vi[mask2], vj[mask2])
                corr_mat[i, j] = corr_mat[j, i] = r2
            else:
                corr_mat[i, j] = corr_mat[j, i] = 0.0

    feat_short = [f.replace("mean_", "").replace("_", " ") for f in FEATURE_NAMES]
    im2 = ax.imshow(corr_mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(n_feat))
    ax.set_yticks(range(n_feat))
    ax.set_xticklabels(feat_short, fontsize=5, color="#CCC", rotation=45, ha="right")
    ax.set_yticklabels(feat_short, fontsize=5, color="#CCC")
    plt.colorbar(im2, ax=ax, fraction=0.03, pad=0.03,
                 label="Spearman ρ").ax.yaxis.set_tick_params(color="#888",
                                                               labelcolor="#888")

    fig.suptitle(
        "Evolutionary Trajectory Map — 26 Human Populations in Process Space\n"
        "Dimensions: diversity · divergence · selection · inbreeding · rare burden",
        color="#EEE", fontsize=13, y=0.97
    )

    out_path = FIG_DIR / "fig19_evo_trajectory.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    t0 = time.time()
    print("=== Investigation 16: Evolutionary Trajectory Map ===\n")

    # ── Load and build feature matrix ────────────────────────────────────
    print("Loading features...", flush=True)
    pops = load_features()
    pop_list, mat, mat_z, scaler = build_matrix(pops)
    n, d = mat_z.shape
    print(f"  Matrix: {n} populations × {d} features")
    print(f"  Features: {', '.join(FEATURE_NAMES)}")

    # ── PCA ───────────────────────────────────────────────────────────────
    print("\nComputing PCA...", flush=True)
    pca_coords, pca = compute_pca(mat_z, n_components=min(6, d))
    ev = pca.explained_variance_ratio_
    for i in range(len(ev)):
        # Top 3 loading features
        top3 = sorted(range(d), key=lambda j: abs(pca.components_[i][j]), reverse=True)[:3]
        top3_names = [f"{FEATURE_NAMES[j]}({pca.components_[i][j]:+.2f})" for j in top3]
        print(f"  PC{i+1}: {ev[i]*100:.1f}% var  |  top: {', '.join(top3_names)}")

    # ── UMAP ──────────────────────────────────────────────────────────────
    print("\nComputing UMAP (n_neighbors=8, min_dist=0.4)...", flush=True)
    umap_coords, _ = compute_umap(mat_z)
    print(f"  UMAP range: x=[{umap_coords[:,0].min():.2f}, {umap_coords[:,0].max():.2f}]  "
          f"y=[{umap_coords[:,1].min():.2f}, {umap_coords[:,1].max():.2f}]")

    # ── Minimum spanning tree ─────────────────────────────────────────────
    print("\nComputing minimum spanning tree...", flush=True)
    mst_edges = compute_mst(mat_z, pop_list)
    print(f"  MST edges: {len(mst_edges)}")
    print("  MST connections:")
    for e in sorted(mst_edges, key=lambda x: x["dist"]):
        p1, p2 = e["pop1"], e["pop2"]
        g1, g2 = pops[p1]["group"], pops[p2]["group"]
        cross = " ← cross-group" if g1 != g2 else ""
        print(f"    {p1:6s}–{p2:6s}  dist={e['dist']:.3f}{cross}")

    # ── Velocity vectors ──────────────────────────────────────────────────
    print("\nComputing evolutionary velocities...", flush=True)
    velocities = compute_velocities(pops, pop_list, mat)
    for comp_name in VELOCITY_COMPONENTS:
        scores = velocities[comp_name]
        top3 = sorted(scores, key=scores.get, reverse=True)[:3]
        bot3 = sorted(scores, key=scores.get)[:3]
        print(f"  {comp_name}:")
        print(f"    High: {', '.join(f'{p}({scores[p]:+.2f})' for p in top3)}")
        print(f"    Low:  {', '.join(f'{p}({scores[p]:+.2f})' for p in bot3)}")

    # ── Correlation with genetic PCA ──────────────────────────────────────
    print("\nCorrelation: statistical PCs vs genetic PC1/OoA distance:")
    corr_results, _ = compare_with_genetic_pca(mat_z, pop_list, pops)
    for r in corr_results:
        sig = "**" if r["p_value"] < 0.05 else ""
        print(f"  {r['stat_pc']:20s} vs {r['genetic_feature']:12s}  "
              f"ρ={r['spearman_r']:+.3f}  p={r['p_value']:.4f} {sig}")

    # ── Figure ────────────────────────────────────────────────────────────
    print("\nGenerating figure...", flush=True)
    make_figure(pop_list, pops, mat_z, mat, mst_edges, umap_coords,
                pca_coords, pca, velocities)

    # ── Write outputs ─────────────────────────────────────────────────────
    rows = []
    for i, pop in enumerate(pop_list):
        row = {
            "population":   pop,
            "group":        pops[pop].get("group", ""),
            "stat_pc1":     round(float(pca_coords[i, 0]), 4),
            "stat_pc2":     round(float(pca_coords[i, 1]), 4),
            "stat_pc3":     round(float(pca_coords[i, 2]), 4) if pca_coords.shape[1] > 2 else "",
            "umap1":        round(float(umap_coords[i, 0]), 4),
            "umap2":        round(float(umap_coords[i, 1]), 4),
            "drift_score":  round(velocities["drift_bottleneck"][pop], 4),
            "selection_score": round(velocities["recent_selection"][pop], 4),
            "expansion_score": round(velocities["expansion"][pop], 4),
        }
        rows.append(row)

    fields = ["population", "group", "stat_pc1", "stat_pc2", "stat_pc3",
              "umap1", "umap2", "drift_score", "selection_score", "expansion_score"]
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    out_data = {
        "n_populations": n,
        "n_features": d,
        "features": FEATURE_NAMES,
        "pca_explained_variance": [round(float(v), 4) for v in ev],
        "mst_edges": mst_edges,
        "velocity_scores": {comp: {p: round(velocities[comp][p], 4)
                                    for p in pop_list}
                             for comp in VELOCITY_COMPONENTS},
        "pc1_correlations": corr_results,
        "embeddings": {pop: {
            "stat_pc1": round(float(pca_coords[i, 0]), 4),
            "stat_pc2": round(float(pca_coords[i, 1]), 4),
            "umap1": round(float(umap_coords[i, 0]), 4),
            "umap2": round(float(umap_coords[i, 1]), 4),
        } for i, pop in enumerate(pop_list)},
    }

    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2)

    print(f"\nOutputs: {OUT_TSV}  {OUT_JSON}")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
