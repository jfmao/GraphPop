#!/usr/bin/env python3
"""Rice Inv.R18: GO Enrichment in Sweep Genes.

Identify GO terms enriched among genes at Garud H sweep peaks, compared to the
background set of all genes with GO annotations.  Analogous to human Inv.05.

Data source: results/rice/rice_interpretation_results.json → sweep_annotations
             Neo4j (Gene → HAS_GO_TERM → GOTerm) for GO annotations
Output:      data/results/rice/rice_inv18_go_enrichment.{tsv,json}
Figure:      data/results/rice/figures/rice_fig18_go_enrichment.png
"""

import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from neo4j import GraphDatabase
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "data" / "results" / "rice"
FIG_DIR = RESULTS_DIR / "figures"

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

# Sweep annotation source
INTERP_JSON = ROOT / "results" / "rice" / "rice_interpretation_results.json"


def load_sweep_genes():
    """Extract unique gene IDs from sweep annotations, per population and overall."""
    with open(INTERP_JSON) as f:
        data = json.load(f)

    sweep_annot = data["annotate"]["sweep_annotations"]

    per_pop = {}   # pop -> set of geneIds
    all_genes = set()

    for pop, pdata in sweep_annot.items():
        genes = set()
        for sweep in pdata.get("sweeps", []):
            for g in sweep.get("genes", []):
                gid = g["geneId"]
                genes.add(gid)
                all_genes.add(gid)
        per_pop[pop] = genes

    return per_pop, all_genes


def fetch_go_annotations(driver):
    """Fetch Gene→GOTerm mapping from Neo4j for all rice genes."""
    t0 = time.time()
    gene_to_go = defaultdict(set)   # geneId -> set of goIds
    go_names = {}                    # goId -> name

    with driver.session() as session:
        result = session.run(
            "MATCH (g:Gene)-[:HAS_GO_TERM]->(go:GOTerm) "
            "WHERE g.geneId STARTS WITH 'LOC_Os' "
            "RETURN g.geneId AS geneId, go.goId AS goId, go.name AS goName"
        )
        for rec in result:
            gene_to_go[rec["geneId"]].add(rec["goId"])
            go_names[rec["goId"]] = rec["goName"]

    elapsed = time.time() - t0
    print(f"  Fetched GO annotations: {len(gene_to_go)} genes, "
          f"{len(go_names)} GO terms in {elapsed:.1f}s")
    sys.stdout.flush()
    return dict(gene_to_go), go_names


def fisher_enrichment(sweep_genes, gene_to_go, go_names, min_genes=2):
    """Run Fisher's exact test for each GO term: sweep genes vs background."""
    # Background = all genes with at least one GO annotation
    bg_genes = set(gene_to_go.keys())
    sweep_with_go = sweep_genes & bg_genes

    if len(sweep_with_go) < 2:
        return []

    # Count GO terms in sweep genes
    sweep_go_counts = Counter()
    for g in sweep_with_go:
        for go_id in gene_to_go[g]:
            sweep_go_counts[go_id] += 1

    # Count GO terms in background
    bg_go_counts = Counter()
    for g in bg_genes:
        for go_id in gene_to_go[g]:
            bg_go_counts[go_id] += 1

    n_sweep = len(sweep_with_go)
    n_bg = len(bg_genes)

    results = []
    for go_id, n_sweep_go in sweep_go_counts.items():
        if n_sweep_go < min_genes:
            continue
        n_bg_go = bg_go_counts[go_id]

        # 2x2 contingency table
        #                 In GO term    Not in GO term
        # Sweep genes     a             b
        # Background      c             d
        a = n_sweep_go
        b = n_sweep - n_sweep_go
        c = n_bg_go - n_sweep_go
        d = (n_bg - n_sweep) - c

        if d < 0:
            d = 0  # guard
        _, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

        fold = (a / n_sweep) / (n_bg_go / n_bg) if n_bg_go > 0 else float("inf")
        results.append({
            "goId": go_id,
            "goName": go_names.get(go_id, go_id),
            "n_sweep_genes": a,
            "n_bg_genes": n_bg_go,
            "n_sweep_total": n_sweep,
            "n_bg_total": n_bg,
            "fold_enrichment": round(fold, 3),
            "p_value": p_value,
        })

    # Sort by p-value
    results.sort(key=lambda x: x["p_value"])

    # BH correction
    n_tests = len(results)
    for i, r in enumerate(results):
        r["p_adjusted"] = min(1.0, r["p_value"] * n_tests / (i + 1))

    return results


def main():
    t0 = time.time()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Rice Inv.R18: GO Enrichment in Sweep Genes ===")
    sys.stdout.flush()

    # Load sweep genes
    per_pop, all_sweep_genes = load_sweep_genes()
    print(f"  Sweep genes: {len(all_sweep_genes)} unique across {len(per_pop)} populations")
    for pop, genes in sorted(per_pop.items(), key=lambda x: -len(x[1])):
        print(f"    {pop}: {len(genes)} genes")
    sys.stdout.flush()

    # Fetch GO annotations from Neo4j
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    gene_to_go, go_names = fetch_go_annotations(driver)
    driver.close()

    # ── Global enrichment (all sweep genes pooled) ──
    print("\n  Global enrichment (all sweep genes pooled)...")
    sys.stdout.flush()
    global_results = fisher_enrichment(all_sweep_genes, gene_to_go, go_names, min_genes=2)
    n_sig = sum(1 for r in global_results if r["p_adjusted"] < 0.05)
    print(f"    {len(global_results)} GO terms tested, {n_sig} significant (padj < 0.05)")

    # ── Per-population enrichment ──
    pop_results = {}
    for pop in sorted(per_pop.keys()):
        genes = per_pop[pop]
        if len(genes) < 3:
            continue
        res = fisher_enrichment(genes, gene_to_go, go_names, min_genes=2)
        n_sig_pop = sum(1 for r in res if r["p_adjusted"] < 0.05)
        if res:
            pop_results[pop] = res
            print(f"    {pop}: {len(genes)} sweep genes, {len(res)} GO terms, "
                  f"{n_sig_pop} sig (padj<0.05)")
    sys.stdout.flush()

    # ── Save TSV ──
    tsv_path = RESULTS_DIR / "rice_inv18_go_enrichment.tsv"
    with open(tsv_path, "w") as f:
        f.write("population\tgoId\tgoName\tn_sweep_genes\tn_bg_genes\t"
                "fold_enrichment\tp_value\tp_adjusted\n")
        # Global first
        for r in global_results:
            f.write(f"ALL\t{r['goId']}\t{r['goName']}\t{r['n_sweep_genes']}\t"
                    f"{r['n_bg_genes']}\t{r['fold_enrichment']:.3f}\t"
                    f"{r['p_value']:.6e}\t{r['p_adjusted']:.6e}\n")
        # Per-population
        for pop, res_list in sorted(pop_results.items()):
            for r in res_list:
                f.write(f"{pop}\t{r['goId']}\t{r['goName']}\t{r['n_sweep_genes']}\t"
                        f"{r['n_bg_genes']}\t{r['fold_enrichment']:.3f}\t"
                        f"{r['p_value']:.6e}\t{r['p_adjusted']:.6e}\n")

    # ── Save JSON ──
    json_path = RESULTS_DIR / "rice_inv18_go_enrichment.json"
    summary = {
        "n_sweep_genes": len(all_sweep_genes),
        "n_genes_with_go": len(gene_to_go),
        "n_go_terms": len(go_names),
        "global_enrichment": {
            "n_tested": len(global_results),
            "n_significant": n_sig,
            "top_10": global_results[:10],
        },
        "per_population": {},
        "total_time_s": time.time() - t0,
    }
    for pop, res_list in pop_results.items():
        n_s = sum(1 for r in res_list if r["p_adjusted"] < 0.05)
        summary["per_population"][pop] = {
            "n_sweep_genes": len(per_pop[pop]),
            "n_tested": len(res_list),
            "n_significant": n_s,
            "top_5": res_list[:5],
        }
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)

    # ── Figure ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    # Panel A: Global top GO terms
    ax = axes[0]
    top_global = global_results[:15]
    if top_global:
        names = [r["goName"][:40] for r in top_global]
        neg_log_p = [-np.log10(max(r["p_value"], 1e-300)) for r in top_global]
        colors = ["#d73027" if r["p_adjusted"] < 0.05 else "#4393c3" for r in top_global]
        ax.barh(range(len(names)), neg_log_p, color=colors, edgecolor="black", linewidth=0.5)
        ax.set_yticks(range(len(names)))
        ax.set_yticklabels(names, fontsize=7)
        ax.set_xlabel("-log10(p-value)")
        ax.axvline(-np.log10(0.05), color="gray", linestyle="--", alpha=0.5, label="p=0.05")
        ax.invert_yaxis()
    ax.set_title(f"Global GO Enrichment (all sweep genes)\n"
                 f"({len(all_sweep_genes)} genes, {n_sig} significant)",
                 fontweight="bold", fontsize=10)
    ax.legend(fontsize=7)

    # Panel B: Per-pop heatmap of top GO terms
    ax = axes[1]
    # Collect top GO terms across all populations
    top_gos = []
    seen = set()
    for pop, res_list in sorted(pop_results.items()):
        for r in res_list[:3]:
            if r["goId"] not in seen:
                top_gos.append(r["goId"])
                seen.add(r["goId"])
            if len(top_gos) >= 15:
                break
        if len(top_gos) >= 15:
            break

    # Also add top global
    for r in global_results[:5]:
        if r["goId"] not in seen:
            top_gos.append(r["goId"])
            seen.add(r["goId"])

    if top_gos and pop_results:
        pops = sorted(pop_results.keys())
        # Build -log10(p) matrix
        mat = np.zeros((len(top_gos), len(pops)))
        for j, pop in enumerate(pops):
            p_map = {r["goId"]: r["p_value"] for r in pop_results[pop]}
            for i, go_id in enumerate(top_gos):
                p = p_map.get(go_id, 1.0)
                mat[i, j] = -np.log10(max(p, 1e-300))

        im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", interpolation="nearest")
        ax.set_xticks(range(len(pops)))
        ax.set_xticklabels(pops, fontsize=6, rotation=45, ha="right")
        ax.set_yticks(range(len(top_gos)))
        ax.set_yticklabels([go_names.get(g, g)[:35] for g in top_gos], fontsize=6)
        plt.colorbar(im, ax=ax, label="-log10(p)", shrink=0.6)
        ax.set_title("Per-Population GO Enrichment", fontweight="bold", fontsize=10)
    else:
        ax.text(0.5, 0.5, "Insufficient data\nfor per-pop heatmap",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_title("Per-Population GO Enrichment", fontweight="bold", fontsize=10)

    plt.tight_layout()
    fig.savefig(FIG_DIR / "rice_fig18_go_enrichment.png", dpi=200, bbox_inches="tight")
    plt.close()

    # ── Summary ──
    total_time = time.time() - t0
    print(f"\n=== Rice Inv.R18: GO Enrichment Summary ===")
    print(f"  Total time: {total_time:.1f}s")
    print(f"  Sweep genes: {len(all_sweep_genes)}")
    print(f"  Genes with GO annotations: {len(gene_to_go)}")
    print(f"  Global: {len(global_results)} GO terms tested, {n_sig} significant")
    if global_results:
        top = global_results[0]
        print(f"    Top: {top['goName']} (fold={top['fold_enrichment']:.1f}, "
              f"p={top['p_value']:.2e}, padj={top['p_adjusted']:.2e})")
    for pop in sorted(pop_results.keys()):
        res = pop_results[pop]
        n_s = sum(1 for r in res if r["p_adjusted"] < 0.05)
        print(f"  {pop}: {len(res)} tested, {n_s} significant")
    print(f"  Saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
