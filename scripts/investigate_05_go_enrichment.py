#!/usr/bin/env python3
"""Investigation 5: GO Term Enrichment in Hard Sweep Genes (4-hop graph traversal).

Question: Which GO biological functions are enriched in genes under hard positive
selection, and does the enrichment pattern differ across human populations?

Method (4-hop traversal):
  GenomicWindow(sweep_type='hard', h12 > threshold)
    → positional overlap → Variant
    → [HAS_CONSEQUENCE] → Gene
    → [HAS_GO_TERM] → GOTerm

  Fisher's exact test for GO term enrichment:
    Foreground = genes in hard sweep windows for a population (or all pops)
    Background = all genes with any HAS_CONSEQUENCE variant

Why novel: 4-hop traversal across the full schema in one analysis.
  VCF-centric workflows need: sweep detection tool → VCF annotation → gene list
  → GO enrichment tool — four separate steps with manual joins.
  GraphPop integrates all layers: sweep windows, variants, genes, GO terms.

Output: data/results/go_enrichment_{population}.tsv  (one file per focal population)
        data/results/go_enrichment_summary.json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_05_go_enrichment.py
"""

import json
import time
from collections import defaultdict
from pathlib import Path

from neo4j import GraphDatabase
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB   = "neo4j"

OUT_DIR        = Path("data/results")
H12_THRESHOLD  = 0.3     # hard sweep window minimum h12
MIN_GENES_GO   = 5       # minimum genes per GO term to test
FDR_ALPHA      = 0.05    # Benjamini-Hochberg FDR threshold
MIN_SWEEP_GENES = 3      # minimum genes in foreground to report a population

# Focal populations (one representative per continental group + all combined)
FOCAL_POPS = [
    "ALL",         # combined across all populations
    "YRI",         # African
    "CEU",         # European
    "CHB",         # East Asian
    "GIH",         # South Asian
    "MXL",         # American
]


# ── Queries ───────────────────────────────────────────────────────────────────

SWEEP_GENE_CYPHER = """
MATCH (w:GenomicWindow)
WHERE w.h12 >= $threshold AND w.sweep_type = 'hard'
  AND ($pop = 'ALL' OR w.population = $pop)
WITH w
MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
WHERE v.chr = w.chr AND v.pos >= w.start AND v.pos <= w.end
RETURN DISTINCT g.geneId AS gene_id
"""

# All gene→GO pairs (background)
GENE_GO_CYPHER = """
MATCH (g:Gene)-[:HAS_GO_TERM]->(t:GOTerm)
RETURN g.geneId AS gene_id, t.goId AS go_id, t.name AS go_name, t.aspect AS aspect
"""

# Background: all genes with any HAS_CONSEQUENCE edge
BACKGROUND_CYPHER = """
MATCH (g:Gene)<-[:HAS_CONSEQUENCE]-(:Variant)
RETURN DISTINCT g.geneId AS gene_id
"""


# ── GO enrichment test ────────────────────────────────────────────────────────

def go_enrichment(foreground_genes: set, background_genes: set,
                  gene_go: dict, go_meta: dict,
                  min_genes: int = MIN_GENES_GO) -> list:
    """Fisher's exact test for GO term enrichment.

    foreground_genes: genes in sweep regions
    background_genes: all annotated genes
    gene_go: gene_id → set of go_ids
    go_meta: go_id → {name, aspect}

    Returns list of dicts sorted by adjusted p-value.
    """
    # Build GO → gene sets
    go_foreground = defaultdict(set)
    go_background = defaultdict(set)
    for g in background_genes:
        for go_id in gene_go.get(g, []):
            go_background[go_id].add(g)
    for g in foreground_genes:
        for go_id in gene_go.get(g, []):
            go_foreground[go_id].add(g)

    n_fg = len(foreground_genes)
    n_bg = len(background_genes)

    results = []
    for go_id, fg_genes in go_foreground.items():
        bg_genes = go_background[go_id]
        if len(bg_genes) < min_genes:
            continue
        a = len(fg_genes)                   # in foreground AND in GO term
        b = n_fg - a                        # in foreground, NOT in GO term
        c = len(bg_genes) - a              # NOT in foreground, in GO term
        d = n_bg - n_fg - c               # NOT in foreground, NOT in GO term
        if c < 0 or d < 0:
            continue
        _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        or_val = (a * d) / (b * c) if b > 0 and c > 0 else float("inf")
        meta = go_meta.get(go_id, {})
        results.append({
            "go_id":        go_id,
            "go_name":      meta.get("name", go_id),
            "aspect":       meta.get("aspect", ""),
            "n_fg":         a,
            "n_bg":         len(bg_genes),
            "n_fg_total":   n_fg,
            "pct_fg":       round(100 * a / n_fg, 1) if n_fg > 0 else 0,
            "pval":         pval,
            "odds_ratio":   round(or_val, 3),
            "genes":        sorted(fg_genes),
        })

    if not results:
        return results

    # Benjamini-Hochberg FDR correction
    pvals = [r["pval"] for r in results]
    _, padj, _, _ = multipletests(pvals, method="fdr_bh")
    for r, p in zip(results, padj):
        r["padj"] = round(float(p), 6)

    results.sort(key=lambda r: r["padj"])
    return results


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 5: GO Term Enrichment in Hard Sweep Genes ===\n", flush=True)

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()

    # Load background data (gene→GO and all annotated genes) once
    print("[1/3] Loading Gene→GO mapping and background gene set ...", flush=True)
    t0 = time.time()
    gene_go   = defaultdict(set)   # gene_id → set of go_ids
    go_meta   = {}                 # go_id → {name, aspect}
    with driver.session(database=NEO4J_DB) as session:
        result = session.run(GENE_GO_CYPHER)
        for rec in result:
            gene_go[rec["gene_id"]].add(rec["go_id"])
            go_meta[rec["go_id"]] = {"name": rec["go_name"], "aspect": rec["aspect"]}

        result = session.run(BACKGROUND_CYPHER)
        background_genes = {rec["gene_id"] for rec in result}

    print(f"  {len(gene_go):,} genes with GO annotations", flush=True)
    print(f"  {len(go_meta):,} unique GO terms", flush=True)
    print(f"  {len(background_genes):,} background genes (have HAS_CONSEQUENCE)", flush=True)
    print(f"  Done in {time.time()-t0:.1f}s", flush=True)

    # Run enrichment for each focal population
    print(f"\n[2/3] Running GO enrichment for {len(FOCAL_POPS)} populations ...", flush=True)
    summary = {}

    with driver.session(database=NEO4J_DB) as session:
        for pop in FOCAL_POPS:
            print(f"\n  Population: {pop} (h12 ≥ {H12_THRESHOLD}) ...", flush=True)
            t1 = time.time()

            # Get genes in hard sweep windows for this population
            result = session.run(SWEEP_GENE_CYPHER, threshold=H12_THRESHOLD, pop=pop)
            sweep_genes = {rec["gene_id"] for rec in result}
            print(f"    {len(sweep_genes):,} genes in hard sweep windows", flush=True)

            if len(sweep_genes) < MIN_SWEEP_GENES:
                print(f"    < {MIN_SWEEP_GENES} genes — skipping", flush=True)
                continue

            # Enrichment test
            enriched = go_enrichment(
                sweep_genes, background_genes, dict(gene_go), go_meta
            )
            sig = [r for r in enriched if r["padj"] < FDR_ALPHA]
            print(f"    {len(sig):,} significant GO terms (FDR < {FDR_ALPHA})", flush=True)
            print(f"    Done in {time.time()-t1:.1f}s", flush=True)

            # Write per-population TSV
            out_tsv = OUT_DIR / f"go_enrichment_{pop}.tsv"
            with open(out_tsv, "w") as f:
                f.write("go_id\tgo_name\taspect\tn_fg\tn_bg\tpct_fg\t"
                        "pval\tpadj\todds_ratio\tgenes\n")
                for r in enriched:
                    f.write(
                        f"{r['go_id']}\t{r['go_name']}\t{r['aspect']}\t"
                        f"{r['n_fg']}\t{r['n_bg']}\t{r['pct_fg']}\t"
                        f"{r['pval']:.2e}\t{r['padj']:.2e}\t{r['odds_ratio']}\t"
                        f"{','.join(r['genes'])}\n"
                    )

            summary[pop] = {
                "n_sweep_genes":   len(sweep_genes),
                "n_go_tested":     len(enriched),
                "n_significant":   len(sig),
                "top5": [{"go_id": r["go_id"], "go_name": r["go_name"],
                           "padj": r["padj"], "n_fg": r["n_fg"]}
                          for r in sig[:5]],
            }

            # Print top 10
            print(f"    Top 10 significant GO terms:")
            for r in sig[:10]:
                aspect_short = r["aspect"][:3].upper() if r["aspect"] else "?"
                print(f"      [{aspect_short}] {r['go_name'][:50]}  "
                      f"(n={r['n_fg']}, FDR={r['padj']:.2e})")

    driver.close()

    # Write summary JSON
    print(f"\n[3/3] Writing summary to {OUT_DIR}/go_enrichment_summary.json ...", flush=True)
    with open(OUT_DIR / "go_enrichment_summary.json", "w") as f:
        json.dump({
            "h12_threshold": H12_THRESHOLD,
            "fdr_alpha":     FDR_ALPHA,
            "min_genes_go":  MIN_GENES_GO,
            "populations":   summary,
        }, f, indent=2)

    print("\n=== Summary ===")
    for pop, s in summary.items():
        print(f"  {pop}: {s['n_sweep_genes']} sweep genes → {s['n_significant']} sig GO terms")
    print(f"\nOutput files: {OUT_DIR}/go_enrichment_*.tsv")


if __name__ == "__main__":
    main()
