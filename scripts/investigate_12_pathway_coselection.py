#!/usr/bin/env python3
"""Investigation 12: Pathway Co-selection Network.

Question: Which pairs of biological pathways are co-regulated by natural selection
across populations? Do co-selected pathway modules correspond to known physiological
systems (immune, metabolic, cardiac)?

Method:
  1. Load pathway_fst.json → 5-dimensional FST vectors per pathway
     (YRI/CEU, YRI/CHB, CEU/CHB, YRI/JPT, CEU/JPT)
  2. Pairwise Spearman correlation between all pathway FST vectors
  3. Build co-selection network: nodes=pathways, edges=strong co-selection (ρ > 0.7)
  4. Find network communities (greedy modularity via networkx)
  5. Identify most co-selected biological module

Why graph-native:
  Built from Inv.8 gene_fst + Neo4j Gene→Pathway traversal (which computed the
  pathway FST vectors). Higher-order network assembled from the primary graph.
  Impossible without integrated gene-to-pathway membership.
  Neo4j is used for enriching community members with pathway names and gene counts.

Output: data/results/pathway_coselection.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_12_pathway_coselection.py
"""

import csv
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import networkx as nx
from neo4j import GraphDatabase
from scipy.stats import spearmanr, rankdata
import os

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"

OUT_DIR  = Path("data/results")
OUT_TSV  = OUT_DIR / "pathway_coselection.tsv"
OUT_JSON = OUT_DIR / "pathway_coselection.json"

RHO_THRESHOLD   = 0.7    # minimum |Spearman ρ| for an edge
MIN_VARIANTS    = 10     # minimum variants for a pathway to be included
MIN_COMPLETE_PAIRS = 3   # minimum non-null FST pairs per pathway

FST_PAIR_KEYS = [
    ("fst_YRI_vs_CEU", "YRI/CEU"),
    ("fst_YRI_vs_CHB", "YRI/CHB"),
    ("fst_CEU_vs_CHB", "CEU/CHB"),
    ("fst_YRI_vs_JPT", "YRI/JPT"),
    ("fst_CEU_vs_JPT", "CEU/JPT"),
]


# ── Fast pairwise Spearman via rank matrices ──────────────────────────────────
def pairwise_spearman(mat):
    """Compute pairwise Spearman correlation for all rows of mat (n_pathways × n_dims).

    Uses scipy.stats.spearmanr with axis=1 for efficiency.
    Returns (n × n) correlation matrix.
    """
    n = mat.shape[0]
    # Rank each row
    ranked = np.apply_along_axis(rankdata, 1, mat)
    # Pearson on ranked = Spearman
    # Center rows
    ranked_c = ranked - ranked.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranked_c, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    ranked_n = ranked_c / norms
    corr = ranked_n @ ranked_n.T
    return np.clip(corr, -1.0, 1.0)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print("=== Investigation 12: Pathway Co-selection Network ===\n", flush=True)

    # ── Load pathway FST data ─────────────────────────────────────────────────
    pathway_fst_path = OUT_DIR / "pathway_fst.json"
    if not pathway_fst_path.exists():
        raise FileNotFoundError(f"Missing {pathway_fst_path} — run investigate_01 first")

    with open(pathway_fst_path) as f:
        pathways_raw = json.load(f)

    print(f"Loaded {len(pathways_raw):,} pathways from pathway_fst.json", flush=True)

    # Filter to pathways with enough variants and enough non-null FST pairs
    pathways = []
    for p in pathways_raw:
        if p.get("n_variants", 0) < MIN_VARIANTS:
            continue
        fst_vals = []
        for key, _ in FST_PAIR_KEYS:
            v = p.get(key)
            fst_vals.append(v if v is not None else np.nan)
        n_valid = sum(1 for v in fst_vals if not np.isnan(v))
        if n_valid >= MIN_COMPLETE_PAIRS:
            pathways.append({
                "pathway_id":   p["pathway_id"],
                "pathway_name": p.get("pathway_name", p["pathway_id"]),
                "n_variants":   p.get("n_variants", 0),
                "mean_fst":     p.get("mean_fst", np.nan),
                "fst_vec":      np.array(fst_vals),
            })

    print(f"  → {len(pathways):,} pathways pass filters "
          f"(≥{MIN_VARIANTS} variants, ≥{MIN_COMPLETE_PAIRS} non-null FST pairs)",
          flush=True)

    # Build FST matrix (impute missing with column median)
    n = len(pathways)
    ndim = len(FST_PAIR_KEYS)
    mat = np.array([p["fst_vec"] for p in pathways])  # n × ndim

    for d in range(ndim):
        col = mat[:, d]
        col_median = np.nanmedian(col)
        mat[np.isnan(mat[:, d]), d] = col_median

    print(f"FST matrix shape: {mat.shape}  (after median imputation)", flush=True)

    # ── Pairwise Spearman correlation ─────────────────────────────────────────
    print("Computing pairwise Spearman correlations...", flush=True)
    corr = pairwise_spearman(mat)   # n × n

    # ── Build co-selection network ────────────────────────────────────────────
    print(f"Building co-selection network (|ρ| > {RHO_THRESHOLD})...", flush=True)
    G = nx.Graph()
    path_id_to_idx = {p["pathway_id"]: i for i, p in enumerate(pathways)}

    for i, p in enumerate(pathways):
        G.add_node(p["pathway_id"],
                   name=p["pathway_name"],
                   n_variants=p["n_variants"],
                   mean_fst=p["mean_fst"])

    edge_rows = []
    for i in range(n):
        for j in range(i + 1, n):
            rho = corr[i, j]
            if abs(rho) >= RHO_THRESHOLD:
                p1_id = pathways[i]["pathway_id"]
                p2_id = pathways[j]["pathway_id"]
                G.add_edge(p1_id, p2_id, rho=float(rho))
                edge_rows.append({
                    "pathway1_id":   p1_id,
                    "pathway1_name": pathways[i]["pathway_name"],
                    "pathway2_id":   p2_id,
                    "pathway2_name": pathways[j]["pathway_name"],
                    "rho":           round(float(rho), 4),
                    "sign":          "positive" if rho > 0 else "negative",
                })

    print(f"  Network: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges",
          flush=True)

    # ── Community detection ───────────────────────────────────────────────────
    print("Finding communities (greedy modularity)...", flush=True)
    if G.number_of_edges() > 0:
        communities_gen = nx.community.greedy_modularity_communities(G)
        communities = [sorted(c) for c in communities_gen]
    else:
        communities = []

    print(f"  Found {len(communities)} communities", flush=True)

    # Annotate each community with mean FST and top pathways
    community_stats = []
    for idx, comm in enumerate(communities):
        fst_vals = [G.nodes[pid]["mean_fst"] for pid in comm
                    if G.nodes[pid]["mean_fst"] is not None
                    and not np.isnan(G.nodes[pid]["mean_fst"])]
        names = [G.nodes[pid]["name"] for pid in comm]
        mean_fst = float(np.mean(fst_vals)) if fst_vals else 0.0
        community_stats.append({
            "community_id": idx,
            "size": len(comm),
            "mean_fst": round(mean_fst, 4),
            "pathways": [{"id": pid, "name": G.nodes[pid]["name"]} for pid in comm],
            "top_pathway": names[0] if names else "",
        })

    community_stats.sort(key=lambda x: x["mean_fst"], reverse=True)
    print("\n  Top 5 communities by mean FST:")
    for c in community_stats[:5]:
        print(f"    community_{c['community_id']:2d}: size={c['size']:3d}  "
              f"mean_fst={c['mean_fst']:.4f}  top: {c['top_pathway'][:60]}")

    # ── Neo4j: enrich top community members with pathway gene counts ──────────
    top_comm = community_stats[0] if community_stats else None
    top_comm_gene_counts = {}
    if top_comm:
        top_ids = [p["id"] for p in top_comm["pathways"][:20]]
        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
        with driver.session(database=NEO4J_DB) as session:
            r = session.run(
                """
                MATCH (p:Pathway)<-[:IN_PATHWAY]-(g:Gene)
                WHERE p.pathwayId IN $ids
                RETURN p.pathwayId AS pid, count(DISTINCT g) AS n_genes
                """,
                ids=top_ids
            )
            for rec in r:
                top_comm_gene_counts[rec["pid"]] = rec["n_genes"]
        driver.close()

    # ── Write outputs ─────────────────────────────────────────────────────────
    edge_rows.sort(key=lambda r: abs(r["rho"]), reverse=True)

    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "pathway1_id", "pathway1_name",
            "pathway2_id", "pathway2_name",
            "rho", "sign"
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(edge_rows)

    # Summary stats
    pos_edges = [e for e in edge_rows if e["sign"] == "positive"]
    neg_edges = [e for e in edge_rows if e["sign"] == "negative"]

    out_data = {
        "n_pathways": len(pathways),
        "rho_threshold": RHO_THRESHOLD,
        "n_edges": len(edge_rows),
        "n_positive_edges": len(pos_edges),
        "n_negative_edges": len(neg_edges),
        "n_communities": len(communities),
        "top_communities": community_stats[:10],
        "top_edges": edge_rows[:50],
        "top_community_gene_counts": top_comm_gene_counts,
        "network_stats": {
            "n_nodes": G.number_of_nodes(),
            "n_edges": G.number_of_edges(),
            "density": nx.density(G) if G.number_of_nodes() > 1 else 0,
            "largest_component_size": (
                max(len(c) for c in nx.connected_components(G))
                if G.number_of_edges() > 0 else 0
            ),
        },
    }

    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2, default=str)

    print(f"\n=== Summary ===")
    print(f"  Pathways analyzed: {len(pathways):,}")
    print(f"  Co-selection edges (|ρ| > {RHO_THRESHOLD}): {len(edge_rows):,}")
    print(f"  Communities detected: {len(communities)}")
    print(f"\nOutputs: {OUT_TSV}  {OUT_JSON}")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
