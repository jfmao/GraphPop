#!/usr/bin/env python3
"""Rice Inv.R19: Pathway Co-Selection Network.

Build a pathway co-selection network based on correlated FST profiles across
population pairs.  Pathways with similar differentiation patterns may reflect
shared selection pressures or co-regulation.  Analogous to human Inv.12.

Data source: data/results/rice/rice_inv17_pathway_fst.tsv
Output:      data/results/rice/rice_inv19_pathway_coselection.{tsv,json}
Figure:      data/results/rice/figures/rice_fig19_pathway_coselection.png
"""

import json
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

PATHWAY_FST_TSV = RESULTS_DIR / "rice_inv17_pathway_fst.tsv"

RHO_THRESHOLD = 0.7   # correlation threshold for co-selection edges


def load_pathway_fst():
    """Load pathway FST data and pivot to pathway x pair matrix."""
    df = pd.read_csv(PATHWAY_FST_TSV, sep="\t")
    print(f"  Loaded {len(df)} rows: {df['pair'].nunique()} pairs, "
          f"{df['pathwayId'].nunique()} pathways")
    sys.stdout.flush()

    # Pivot: pathways as rows, pairs as columns, FST as values
    pivot = df.pivot_table(index="pathwayId", columns="pair",
                           values="fst", aggfunc="first")
    # Drop pathways missing any pair
    pivot = pivot.dropna()
    print(f"  Pivot matrix: {pivot.shape[0]} pathways x {pivot.shape[1]} pairs "
          f"(after dropping NaN)")

    # Add pathway names
    name_map = dict(zip(df["pathwayId"], df["name"]))

    return pivot, name_map


def compute_pairwise_correlations(pivot):
    """Compute Spearman correlation between all pathway pairs."""
    t0 = time.time()
    pathways = pivot.index.tolist()
    n = len(pathways)
    fst_mat = pivot.values  # (n_pathways, n_pairs)

    # Check if we have enough pairs for meaningful correlation
    n_pairs = fst_mat.shape[1]
    print(f"  Computing pairwise Spearman correlations: {n} pathways, {n_pairs} pairs")
    sys.stdout.flush()

    if n_pairs < 3:
        print(f"  WARNING: Only {n_pairs} population pairs. "
              "Spearman correlation with <3 points has limited meaning.")

    # Compute full correlation matrix
    rho_mat = np.zeros((n, n))
    p_mat = np.ones((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            if n_pairs >= 3:
                rho, p = stats.spearmanr(fst_mat[i], fst_mat[j])
            else:
                # Pearson for 2 points (Spearman undefined for n<3)
                rho, p = stats.pearsonr(fst_mat[i], fst_mat[j])
            if np.isnan(rho):
                rho, p = 0.0, 1.0
            rho_mat[i, j] = rho
            rho_mat[j, i] = rho
            p_mat[i, j] = p
            p_mat[j, i] = p

    np.fill_diagonal(rho_mat, 1.0)
    np.fill_diagonal(p_mat, 0.0)

    elapsed = time.time() - t0
    print(f"  Correlation matrix computed in {elapsed:.1f}s")
    return rho_mat, p_mat, pathways


def build_network(rho_mat, p_mat, pathways, name_map, threshold=RHO_THRESHOLD):
    """Build co-selection edges from correlated pathway pairs."""
    n = len(pathways)
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if abs(rho_mat[i, j]) >= threshold:
                edges.append({
                    "pathway1": pathways[i],
                    "name1": name_map.get(pathways[i], pathways[i]),
                    "pathway2": pathways[j],
                    "name2": name_map.get(pathways[j], pathways[j]),
                    "rho": round(float(rho_mat[i, j]), 4),
                    "p_value": float(p_mat[i, j]),
                })

    print(f"  Network: {len(edges)} edges (|rho| >= {threshold}) "
          f"among {n} pathways")
    return edges


def detect_communities(rho_mat, pathways, name_map):
    """Detect communities using hierarchical clustering on distance = 1 - |rho|."""
    dist = 1.0 - np.abs(rho_mat)
    np.fill_diagonal(dist, 0)
    # Ensure symmetry and non-negative
    dist = np.clip((dist + dist.T) / 2, 0, 2)

    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")

    # Cut at a distance that produces reasonable clusters
    labels = fcluster(Z, t=0.3, criterion="distance")
    n_communities = len(set(labels))

    communities = {}
    for i, lab in enumerate(labels):
        lab = int(lab)
        if lab not in communities:
            communities[lab] = []
        communities[lab].append({
            "pathwayId": pathways[i],
            "name": name_map.get(pathways[i], pathways[i]),
        })

    print(f"  Detected {n_communities} communities "
          f"(sizes: {sorted([len(v) for v in communities.values()], reverse=True)[:10]})")
    return communities, labels, Z


def main():
    t0 = time.time()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Rice Inv.R19: Pathway Co-Selection Network ===")
    sys.stdout.flush()

    # Load data
    pivot, name_map = load_pathway_fst()

    # Compute correlations
    rho_mat, p_mat, pathways = compute_pairwise_correlations(pivot)

    # Build network
    edges = build_network(rho_mat, p_mat, pathways, name_map)

    # Detect communities
    communities, labels, Z = detect_communities(rho_mat, pathways, name_map)

    # ── Save TSV (edges) ──
    tsv_path = RESULTS_DIR / "rice_inv19_pathway_coselection.tsv"
    with open(tsv_path, "w") as f:
        f.write("pathway1\tname1\tpathway2\tname2\trho\tp_value\n")
        for e in sorted(edges, key=lambda x: -abs(x["rho"])):
            f.write(f"{e['pathway1']}\t{e['name1']}\t{e['pathway2']}\t"
                    f"{e['name2']}\t{e['rho']}\t{e['p_value']:.6e}\n")

    # ── Save JSON ──
    json_path = RESULTS_DIR / "rice_inv19_pathway_coselection.json"
    n_pos = sum(1 for e in edges if e["rho"] > 0)
    n_neg = sum(1 for e in edges if e["rho"] < 0)

    # Degree distribution
    degree = {}
    for e in edges:
        degree[e["pathway1"]] = degree.get(e["pathway1"], 0) + 1
        degree[e["pathway2"]] = degree.get(e["pathway2"], 0) + 1

    top_hubs = sorted(degree.items(), key=lambda x: -x[1])[:10]

    summary = {
        "n_pathways": len(pathways),
        "n_pairs": pivot.shape[1],
        "n_edges": len(edges),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "rho_threshold": RHO_THRESHOLD,
        "n_communities": len(communities),
        "communities": {
            str(k): {"n_pathways": len(v), "pathways": v[:5]}
            for k, v in sorted(communities.items())
        },
        "top_hubs": [
            {"pathwayId": pid, "name": name_map.get(pid, pid), "degree": d}
            for pid, d in top_hubs
        ],
        "top_edges": edges[:10],
        "total_time_s": time.time() - t0,
    }
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # ── Figure ──
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel A: Correlation heatmap (subsample if too large)
    ax = axes[0]
    n_show = min(50, len(pathways))
    if n_show < len(pathways):
        # Show top-degree pathways
        top_idx = sorted(range(len(pathways)),
                         key=lambda i: degree.get(pathways[i], 0),
                         reverse=True)[:n_show]
        sub_mat = rho_mat[np.ix_(top_idx, top_idx)]
        sub_names = [name_map.get(pathways[i], pathways[i])[:25] for i in top_idx]
    else:
        sub_mat = rho_mat
        sub_names = [name_map.get(p, p)[:25] for p in pathways]

    im = ax.imshow(sub_mat, cmap="RdBu_r", vmin=-1, vmax=1,
                   aspect="auto", interpolation="nearest")
    if n_show <= 30:
        ax.set_xticks(range(len(sub_names)))
        ax.set_xticklabels(sub_names, fontsize=5, rotation=90)
        ax.set_yticks(range(len(sub_names)))
        ax.set_yticklabels(sub_names, fontsize=5)
    ax.set_title(f"Pathway Correlation Matrix\n(top {n_show} by degree)",
                 fontweight="bold", fontsize=10)
    plt.colorbar(im, ax=ax, label="Spearman rho", shrink=0.6)

    # Panel B: Degree distribution
    ax = axes[1]
    if degree:
        deg_vals = list(degree.values())
        ax.hist(deg_vals, bins=min(30, max(deg_vals) + 1),
                color="#4393c3", edgecolor="black", linewidth=0.5)
        ax.set_xlabel("Degree (number of co-selected partners)")
        ax.set_ylabel("Number of pathways")
        ax.set_title(f"Degree Distribution\n({len(edges)} edges, |rho|>{RHO_THRESHOLD})",
                     fontweight="bold", fontsize=10)
    else:
        ax.text(0.5, 0.5, "No edges above threshold",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_title("Degree Distribution", fontweight="bold", fontsize=10)

    # Panel C: Dendrogram
    ax = axes[2]
    from scipy.cluster.hierarchy import dendrogram
    if len(pathways) <= 80:
        leaf_labels = [name_map.get(p, p)[:20] for p in pathways]
    else:
        leaf_labels = None  # too many to show
    dendrogram(Z, ax=ax, labels=leaf_labels, leaf_rotation=90,
               leaf_font_size=5, color_threshold=0.3)
    ax.set_ylabel("Distance (1 - |rho|)")
    ax.set_title(f"Hierarchical Clustering\n({len(communities)} communities)",
                 fontweight="bold", fontsize=10)
    ax.axhline(0.3, color="red", linestyle="--", alpha=0.5, label="cut=0.3")
    ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig19_pathway_coselection.png",
                dpi=200, bbox_inches="tight")
    plt.close()

    # ── Summary ──
    total_time = time.time() - t0
    print(f"\n=== Rice Inv.R19: Pathway Co-Selection Summary ===")
    print(f"  Total time: {total_time:.1f}s")
    print(f"  Pathways: {len(pathways)}")
    print(f"  Edges (|rho|>{RHO_THRESHOLD}): {len(edges)} "
          f"({n_pos} positive, {n_neg} negative)")
    print(f"  Communities: {len(communities)}")
    if top_hubs:
        hub = top_hubs[0]
        print(f"  Top hub: {hub[1]} ({name_map.get(hub[0], hub[0])}, degree={hub[1]})")
    if edges:
        top_e = edges[0] if edges[0]["rho"] >= 0 else max(edges, key=lambda x: x["rho"])
        print(f"  Strongest edge: {top_e['name1'][:30]} <-> "
              f"{top_e['name2'][:30]} (rho={top_e['rho']:.3f})")
    print(f"  Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
