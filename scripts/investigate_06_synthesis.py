#!/usr/bin/env python3
"""Investigation 6: Multi-Evidence Synthesis.

Question: Which genes and pathways have the strongest convergent evidence across
ALL layers of analysis — population differentiation (FST), convergent sweeps,
GO enrichment, and ROH/diversity extremes?

Method:
  Cross-links results from all 5 investigations:
    Inv.1: pathway_fst.json         — pathways ranked by mean FST
    Inv.2: convergent_sweeps.json   — genes with sweeps in 2+ continental groups
    Inv.3: roh_sweep_correlation.json — population-level ROH/sweep metrics
    Inv.5: go_enrichment_summary.json — GO terms enriched in sweep genes

  Builds a composite evidence score per pathway:
    evidence_score = (fst_rank_norm) + (n_convergent_genes_norm) +
                     (go_enrichment_norm) + (n_groups_norm)

  Assigns confidence tier: HIGH (≥3 lines), MEDIUM (2 lines), LOW (1 line)

Output: data/results/synthesis.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_06_synthesis.py
"""

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

# ── Config ────────────────────────────────────────────────────────────────────
OUT_DIR  = Path("data/results")
OUT_TSV  = OUT_DIR / "synthesis.tsv"
OUT_JSON = OUT_DIR / "synthesis.json"

# Input files
PATHWAY_FST_JSON       = OUT_DIR / "pathway_fst.json"
CONV_SWEEPS_JSON       = OUT_DIR / "convergent_sweeps.json"
GO_SUMMARY_JSON        = OUT_DIR / "go_enrichment_summary.json"
ROH_SWEEP_JSON         = OUT_DIR / "roh_sweep_correlation.json"
KCNE1_JSON             = OUT_DIR / "kcne1_deepdive.json"
OOA_JSON               = OUT_DIR / "ooa_gradient.json"

H12_HARD_THRESH = 0.3
MIN_EVIDENCE_LINES = 2   # minimum distinct evidence lines to report


# ── Helpers ───────────────────────────────────────────────────────────────────

def rank_normalize(values: list, reverse=False) -> np.ndarray:
    """Map values to [0, 1] via rank normalization."""
    arr = np.array(values, dtype=float)
    mask = np.isfinite(arr)
    if mask.sum() == 0:
        return np.zeros(len(arr))
    ranks = np.zeros(len(arr))
    finite_vals = arr[mask]
    order = np.argsort(finite_vals)
    if reverse:
        order = order[::-1]
    n = len(finite_vals)
    norm_ranks = np.arange(n) / max(n - 1, 1)
    rank_map = np.zeros(n)
    rank_map[order] = norm_ranks
    ranks[mask] = rank_map
    return ranks


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=== Investigation 6: Multi-Evidence Synthesis ===\n", flush=True)

    # ── Load Inv.1: Pathway FST ───────────────────────────────────────────────
    print("[1/6] Loading pathway FST results (Inv.1) ...", flush=True)
    with open(PATHWAY_FST_JSON) as f:
        pathway_fst_list = json.load(f)

    pathway_fst = {r["pathway_id"]: r for r in pathway_fst_list}
    pathway_names = {r["pathway_id"]: r["pathway_name"] for r in pathway_fst_list}
    print(f"  {len(pathway_fst):,} pathways with FST", flush=True)

    # ── Load Inv.2: Convergent sweeps ─────────────────────────────────────────
    print("[2/6] Loading convergent sweep results (Inv.2) ...", flush=True)
    with open(CONV_SWEEPS_JSON) as f:
        conv = json.load(f)

    conv_genes    = conv["convergent_genes"]         # gene_id → {groups, pops, max_h12, n_groups}
    conv_pathways = {r["pathway_id"]: r for r in conv["pathways"]}  # pathway_id → row
    print(f"  {len(conv_genes):,} convergent genes, {len(conv_pathways):,} pathways", flush=True)

    # ── Load Inv.5: GO enrichment ─────────────────────────────────────────────
    print("[3/6] Loading GO enrichment results (Inv.5) ...", flush=True)
    go_enrichment_available = GO_SUMMARY_JSON.exists()
    go_pop_summary = {}
    if go_enrichment_available:
        with open(GO_SUMMARY_JSON) as f:
            go_summary = json.load(f)
        go_pop_summary = go_summary.get("populations", {})
        print(f"  {len(go_pop_summary)} populations with GO enrichment", flush=True)
    else:
        print(f"  WARNING: {GO_SUMMARY_JSON} not found — skipping GO layer", flush=True)

    # Build gene → enriched GO terms (from per-population TSV files)
    gene_go_enriched = defaultdict(set)  # gene_id → set of enriched GO term names
    for pop in go_pop_summary:
        tsv_file = OUT_DIR / f"go_enrichment_{pop}.tsv"
        if not tsv_file.exists():
            continue
        with open(tsv_file) as f:
            for line in f:
                if line.startswith("go_id"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                go_name = parts[1]
                padj = float(parts[7]) if parts[7] not in ("nan", "") else 1.0
                if padj < 0.05 and len(parts) > 9:
                    genes = parts[9].split(",")
                    for g in genes:
                        gene_go_enriched[g.strip()].add(go_name)

    print(f"  {len(gene_go_enriched):,} genes in at least one enriched GO term", flush=True)

    # ── Load Inv.3: ROH/sweep correlation ─────────────────────────────────────
    print("[4/6] Loading ROH-sweep correlation results (Inv.3) ...", flush=True)
    roh_sweep = {}
    if ROH_SWEEP_JSON.exists():
        with open(ROH_SWEEP_JSON) as f:
            rs = json.load(f)
        roh_sweep = rs.get("per_population", {})
        correlations = rs.get("correlations", [])
        # Find key correlation: froh vs sweep intensity
        froh_h12 = next((c for c in correlations if c["x"] == "froh" and c["y"] == "mean_h12"), None)
        if froh_h12:
            print(f"  FROH vs mean_h12: ρ={froh_h12['rho']:.4f}  p={froh_h12['pval']:.4f}", flush=True)

    # ── Load KCNE1 deep-dive ──────────────────────────────────────────────────
    print("[5/6] Loading KCNE1 deep-dive (Inv.4) ...", flush=True)
    kcne1_data = {}
    if KCNE1_JSON.exists():
        with open(KCNE1_JSON) as f:
            kcne1_data = json.load(f)
        n_groups = kcne1_data.get("n_groups", 0)
        print(f"  KCNE1 sweep in {n_groups} continental groups", flush=True)
        locus_fst = kcne1_data.get("locus_fst", {})
        if locus_fst:
            ratios = [v["fst_ratio"] for v in locus_fst.values() if v.get("fst_ratio") is not None]
            if ratios:
                print(f"  KCNE1 FST ratio (locus/genome): {np.mean(ratios):.3f} "
                      f"(mean, {len(ratios)} pairs)", flush=True)

    # ── Build composite evidence per pathway ──────────────────────────────────
    print("\n[6/6] Building composite evidence scores ...", flush=True)

    all_pathway_ids = sorted(set(pathway_fst.keys()) | set(conv_pathways.keys()))
    print(f"  {len(all_pathway_ids):,} pathways to score", flush=True)

    # Collect raw scores for normalization
    fst_scores      = []
    conv_gene_scores = []
    conv_group_scores = []
    go_scores       = []

    for pid in all_pathway_ids:
        fst_data  = pathway_fst.get(pid, {})
        conv_data = conv_pathways.get(pid, {})

        mean_fst      = fst_data.get("mean_fst", None) or 0.0
        n_conv_genes  = conv_data.get("n_convergent_genes", 0)
        n_conv_groups = conv_data.get("n_groups", 0)

        # GO evidence: count genes in this pathway that appear in enriched GO terms
        pathway_conv_genes = conv_data.get("convergent_genes", [])
        n_go_genes = sum(1 for g in pathway_conv_genes if g in gene_go_enriched)

        fst_scores.append(mean_fst)
        conv_gene_scores.append(n_conv_genes)
        conv_group_scores.append(n_conv_groups)
        go_scores.append(n_go_genes)

    # Normalize all scores
    fst_norm       = rank_normalize(fst_scores,       reverse=False)
    conv_gene_norm = rank_normalize(conv_gene_scores,  reverse=False)
    conv_grp_norm  = rank_normalize(conv_group_scores, reverse=False)
    go_norm        = rank_normalize(go_scores,         reverse=False)

    # Build results
    pathway_results = []
    for i, pid in enumerate(all_pathway_ids):
        fst_data  = pathway_fst.get(pid, {})
        conv_data = conv_pathways.get(pid, {})
        pname = pathway_names.get(pid) or conv_data.get("pathway_name", pid)

        # Count evidence lines (non-zero contributions)
        evidence_lines = []
        if fst_scores[i] > 0:           evidence_lines.append("FST")
        if conv_gene_scores[i] > 0:     evidence_lines.append("Convergent_Sweep")
        if go_scores[i] > 0:            evidence_lines.append("GO_Enrichment")
        if conv_group_scores[i] >= 4:   evidence_lines.append("Pan_continental")

        # Composite score (equal weights)
        composite = (fst_norm[i] + conv_gene_norm[i] + conv_grp_norm[i] + go_norm[i]) / 4.0

        tier = "HIGH" if len(evidence_lines) >= 3 else \
               "MEDIUM" if len(evidence_lines) >= 2 else "LOW"

        pathway_results.append({
            "pathway_id":          pid,
            "pathway_name":        pname,
            "mean_fst":            round(fst_scores[i], 6),
            "n_convergent_genes":  conv_gene_scores[i],
            "n_conv_groups":       conv_group_scores[i],
            "n_go_enriched_genes": go_scores[i],
            "composite_score":     round(float(composite), 4),
            "evidence_lines":      evidence_lines,
            "n_evidence_lines":    len(evidence_lines),
            "tier":                tier,
            "continental_groups":  conv_data.get("continental_groups", []),
            "fisher_pval":         conv_data.get("fisher_pval"),
        })

    # Sort: tier → composite score
    tier_order = {"HIGH": 0, "MEDIUM": 1, "LOW": 2}
    pathway_results.sort(key=lambda r: (tier_order[r["tier"]], -r["composite_score"]))

    # ── Print top findings ─────────────────────────────────────────────────────
    high_med = [r for r in pathway_results if r["tier"] in ("HIGH", "MEDIUM")]
    print(f"\n  {len(high_med)} pathways with MEDIUM or HIGH evidence ({len([r for r in pathway_results if r['tier']=='HIGH'])} HIGH)")

    print("\n=== Top 30 Pathways by Composite Evidence ===")
    print(f"{'#':<4} {'Tier':<7} {'Score':>6} {'FST':>7} {'Conv':>5} {'Grps':>5} {'GO':>4}  Pathway")
    print("-" * 100)
    for i, r in enumerate(pathway_results[:30], 1):
        fst_s = f"{r['mean_fst']:.4f}" if r["mean_fst"] else "    --"
        grps  = ",".join(g[:3] for g in r["continental_groups"]) if r["continental_groups"] else "--"
        name  = r["pathway_name"][:55]
        print(f"  {i:<3} {r['tier']:<7} {r['composite_score']:>6.4f} {fst_s}  "
              f"{r['n_convergent_genes']:>5} {r['n_conv_groups']:>5} {r['n_go_enriched_genes']:>4}  {name}")

    # ── Notable gene-level synthesis ──────────────────────────────────────────
    print("\n=== Top Convergent Genes (2+ continental groups + GO enrichment) ===")
    print(f"{'Gene':<12} {'Groups':>6} {'max_h12':>8} {'GO_terms':>8}  Continental groups")
    print("-" * 80)
    gene_evidence = []
    for gid, info in conv_genes.items():
        n_go = len(gene_go_enriched.get(gid, set()))
        gene_evidence.append({
            "gene":    gid,
            "groups":  info["n_groups"],
            "max_h12": info["max_h12"],
            "n_go":    n_go,
            "group_names": info["groups"],
        })
    gene_evidence.sort(key=lambda x: (-x["groups"], -x["max_h12"], -x["n_go"]))
    for ge in gene_evidence[:20]:
        grp_str = ",".join(g[:3] for g in ge["group_names"])
        print(f"  {ge['gene']:<11} {ge['groups']:>6} {ge['max_h12']:>8.4f} {ge['n_go']:>8}  {grp_str}")

    # ── KCNE1 synthesis entry ─────────────────────────────────────────────────
    kcne1_entry = next((g for g in gene_evidence if g["gene"] == "KCNE1"), None)
    if kcne1_entry:
        print(f"\n=== KCNE1 Multi-Evidence Summary ===")
        print(f"  Convergent groups: {kcne1_entry['groups']} ({', '.join(kcne1_entry['group_names'])})")
        print(f"  Max h12:          {kcne1_entry['max_h12']:.4f}")
        print(f"  GO enriched terms: {kcne1_entry['n_go']}")
        if kcne1_data:
            locus_fst_vals = kcne1_data.get("locus_fst", {})
            ratios = [v["fst_ratio"] for v in locus_fst_vals.values()
                      if v.get("fst_ratio") is not None]
            if ratios:
                mean_ratio = np.mean(ratios)
                print(f"  FST locus/genome ratio: {mean_ratio:.3f}")
                if mean_ratio < 0.5:
                    print(f"  → ANCIENT sweep: haplotype fixed before OoA dispersal")
                elif mean_ratio > 1.5:
                    print(f"  → CONVERGENT sweep: independent selection events")
                else:
                    print(f"  → MIXED/NEUTRAL FST pattern")

    # ── OoA gradient synthesis ────────────────────────────────────────────────
    if OOA_JSON.exists():
        with open(OOA_JSON) as f:
            ooa = json.load(f)
        ooa_corr = ooa.get("correlations_with_ooa_dist", [])
        print(f"\n=== OoA Gradient Key Correlations (Inv.7) ===")
        for c in sorted(ooa_corr, key=lambda x: abs(x.get("rho", 0)), reverse=True)[:5]:
            star = "***" if c["pval"] < 0.001 else ("**" if c["pval"] < 0.01 else
                   ("*" if c["pval"] < 0.05 else ""))
            print(f"  {c['metric']:<18} ρ={c['rho']:+.4f}  p={c['pval']:.4f} {star}")

    # ── Write outputs ─────────────────────────────────────────────────────────
    with open(OUT_TSV, "w") as f:
        f.write("pathway_id\tpathway_name\ttier\tcomposite_score\tmean_fst\t"
                "n_convergent_genes\tn_conv_groups\tn_go_enriched_genes\tevidence_lines\n")
        for r in pathway_results:
            f.write(f"{r['pathway_id']}\t{r['pathway_name']}\t{r['tier']}\t"
                    f"{r['composite_score']}\t{r['mean_fst']}\t{r['n_convergent_genes']}\t"
                    f"{r['n_conv_groups']}\t{r['n_go_enriched_genes']}\t"
                    f"{','.join(r['evidence_lines'])}\n")

    with open(OUT_JSON, "w") as f:
        json.dump({
            "summary": {
                "n_pathways_total":  len(pathway_results),
                "n_high_evidence":   sum(1 for r in pathway_results if r["tier"] == "HIGH"),
                "n_medium_evidence": sum(1 for r in pathway_results if r["tier"] == "MEDIUM"),
                "n_low_evidence":    sum(1 for r in pathway_results if r["tier"] == "LOW"),
            },
            "top_pathways":     pathway_results[:50],
            "top_genes":        gene_evidence[:50],
        }, f, indent=2)

    print(f"\nOutput: {OUT_TSV}")
    print(f"        {OUT_JSON}")


if __name__ == "__main__":
    main()
