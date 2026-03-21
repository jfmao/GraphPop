#!/usr/bin/env python3
"""
Generate Source Data and Supplementary Table files for Nature journal submission.

Source Data files (paper/source_data/):
  - One TSV per figure containing the exact data plotted
  - Required by Nature for every figure panel showing quantitative data

Supplementary Tables (paper/supplementary_tables/):
  - ST01: 1000G population summary (26 populations)
  - ST02: Top 50 FST pathways (Inv.01)
  - ST03: Convergent sweep genes (Inv.02)
  - ST04: Top 100 genes by FST (Inv.08)
  - ST05: GO enrichment significant terms (Inv.05)
  - ST06: PBS top genes per population (Inv.09)
  - ST07: Temporal selection stratification (Inv.10)
  - ST08: Functional consequence bias by VEP impact (Inv.11)
  - ST09: Pathway co-selection top communities (Inv.12)
  - ST10: Cross-species functional convergence (Inv.13)
  - ST11: Top 50 genes by population-specific rare variant burden (Inv.14)
  - ST12: Convergent sweep genomic extents (Inv.15)
  - ST13: Benchmark performance summary

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/generate_tables.py
"""

import csv
import json
from pathlib import Path
from collections import defaultdict

RESULTS       = Path("data/results")
SOURCE_DIR    = Path("paper/source_data")
SUPP_DIR      = Path("paper/supplementary_tables")
BENCHMARKS    = Path("benchmarks/results")


# ── Helpers ───────────────────────────────────────────────────────────────────

def write_tsv(path: Path, rows: list, fieldnames: list, note: str = ""):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tag = f"  ({note})" if note else ""
    print(f"  {path.name:<60s}  {len(rows):>6,} rows{tag}")


def load_tsv(path: Path) -> list:
    if not path.exists():
        print(f"  WARNING: missing {path}")
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load_json(path: Path):
    if not path.exists():
        print(f"  WARNING: missing {path}")
        return {}
    with open(path) as f:
        return json.load(f)


def _fmt(v, decimals=4):
    """Format a numeric value for output."""
    if v is None or v == "" or v != v:  # None / empty / NaN
        return ""
    try:
        return str(round(float(v), decimals))
    except (TypeError, ValueError):
        return str(v)


# ── Source Data files (one per figure) ────────────────────────────────────────

def source_supp_fig01_pathway_fst():
    """Top 50 pathways by mean FST (backing fig1_pathway_fst.png)."""
    rows = load_tsv(RESULTS / "pathway_fst.tsv")
    # sort by mean_fst descending; keep top 50 with data
    rows = [r for r in rows if r.get("mean_fst", "")]
    rows.sort(key=lambda r: float(r["mean_fst"]) if r["mean_fst"] else 0, reverse=True)
    rows = rows[:50]
    fields = ["pathway_id", "pathway_name", "n_variants", "mean_fst",
              "fst_YRI_vs_CEU", "fst_YRI_vs_CHB", "fst_CEU_vs_CHB",
              "fst_YRI_vs_JPT", "fst_CEU_vs_JPT"]
    write_tsv(SOURCE_DIR / "SourceData_SuppFig01_pathway_fst.tsv", rows, fields,
              "top 50 of 2241 pathways")


def source_supp_fig02_gene_fst_landscape():
    """Gene FST distribution (backing fig2_fst_landscape.png)."""
    rows = load_tsv(RESULTS / "gene_fst.tsv")
    rows = [r for r in rows if r.get("mean_fst", "")]
    rows.sort(key=lambda r: float(r["mean_fst"]) if r["mean_fst"] else 0, reverse=True)
    fields = ["gene_id", "n_variants", "mean_fst", "max_fst",
              "fst_YRI_vs_CEU", "fst_YRI_vs_CHB", "fst_CEU_vs_CHB",
              "in_sweep", "n_conv_groups", "max_h12"]
    write_tsv(SOURCE_DIR / "SourceData_SuppFig02_gene_fst_landscape.tsv", rows, fields,
              "all 21194 genes")


def source_supp_fig03_convergent_sweeps():
    """Convergent sweep genes (backing fig3_convergent_kcne1.png)."""
    rows = load_tsv(RESULTS / "convergent_sweeps.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig03_convergent_sweeps.tsv", rows, fields)


def source_supp_fig04_roh_sweep():
    """ROH × sweep correlation (backing fig4_roh_sweep_correlation.png)."""
    rows = load_tsv(RESULTS / "roh_sweep_correlation.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig04_roh_sweep_correlation.tsv", rows, fields)


def source_supp_fig05_go_enrichment():
    """GO enrichment results (backing fig5_go_enrichment.png)."""
    rows = load_tsv(RESULTS / "go_enrichment_ALL.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig05_go_enrichment.tsv", rows, fields)


def source_supp_fig06_pca():
    """Population PCA coordinates (backing fig6_pca_populations.png)."""
    rows = load_tsv(RESULTS / "ooa_gradient.tsv")
    fields = ["pop", "group", "pc1", "pc2", "lat", "lon",
              "mean_pi", "mean_fst", "froh"]
    write_tsv(SOURCE_DIR / "SourceData_SuppFig06_pca_populations.tsv", rows, fields)


def source_supp_fig07_ooa_gradient():
    """OoA distance gradient (backing fig7_ooa_gradient.png)."""
    rows = load_tsv(RESULTS / "ooa_gradient.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig07_ooa_gradient.tsv", rows, fields)


def source_supp_fig08_synthesis():
    """Multi-evidence synthesis heatmap (backing fig8_synthesis_heatmap.png)."""
    rows = load_tsv(RESULTS / "synthesis.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig08_synthesis_heatmap.tsv", rows, fields)


def source_supp_fig09_summary_overview():
    """Summary overview — per-pop temporal stats (backing fig9_summary_overview.png)."""
    rows = load_tsv(RESULTS / "temporal_selection.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig09_summary_overview.tsv", rows, fields)


def source_supp_fig10_gene_fst():
    """Top 100 genes by FST (backing fig10_gene_fst.png)."""
    rows = load_tsv(RESULTS / "gene_fst.tsv")
    rows = [r for r in rows if r.get("mean_fst", "")]
    rows.sort(key=lambda r: float(r["mean_fst"]) if r["mean_fst"] else 0, reverse=True)
    rows = rows[:100]
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig10_gene_fst_top100.tsv", rows, fields,
              "top 100 of 21194")


def source_supp_fig11_pbs():
    """PBS scan results (backing fig11_pbs_scan.png)."""
    rows = load_tsv(RESULTS / "pbs_scan.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig11_pbs_scan.tsv", rows, fields)


def source_supp_fig12_temporal():
    """Temporal stratification (backing fig12_temporal_selection.png)."""
    rows = load_tsv(RESULTS / "temporal_selection.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig12_temporal_selection.tsv", rows, fields)


def source_supp_fig13_kcne1_overview():
    """KCNE1 deep-dive (backing fig13_new_overview.png)."""
    rows = load_tsv(RESULTS / "kcne1_deepdive.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig13_kcne1_deepdive.tsv", rows, fields)


def source_supp_fig14_consequence_bias():
    """Consequence bias (backing fig14_consequence_bias.png)."""
    rows = load_tsv(RESULTS / "consequence_bias.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig14_consequence_bias.tsv", rows, fields)


def source_supp_fig15_pathway_coselection():
    """Top 500 pathway co-selection edges (backing fig15_pathway_coselection.png)."""
    rows = load_tsv(RESULTS / "pathway_coselection.tsv")
    rows = rows[:500]  # top 500 by |rho|, already sorted
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig15_pathway_coselection_top500.tsv",
              rows, fields, "top 500 of 1.18M edges")


def source_supp_fig16_cross_species():
    """Cross-species convergence (backing fig16_cross_species.png)."""
    rows = load_tsv(RESULTS / "cross_species.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig16_cross_species.tsv", rows, fields)


def source_supp_fig17_rare_burden():
    """Rare burden top genes × populations (backing fig17_rare_burden.png)."""
    rows = load_tsv(RESULTS / "rare_burden.tsv")
    # top 50 genes by n_private_rare
    gene_totals = defaultdict(int)
    for r in rows:
        gene_totals[r["gene_id"]] += int(r["n_private_rare"])
    top_genes = {g for g, _ in sorted(gene_totals.items(),
                                       key=lambda x: x[1], reverse=True)[:50]}
    rows = [r for r in rows if r["gene_id"] in top_genes]
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig17_rare_burden_top50genes.tsv",
              rows, fields, "top 50 genes of 27899")


def source_supp_fig18_sweep_extent():
    """Sweep extent FST decay (backing fig18_sweep_extent.png)."""
    # Full sweep_extent.tsv + decay data from JSON
    rows = load_tsv(RESULTS / "sweep_extent.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig18_sweep_extent.tsv", rows, fields)

    # Also write the FST decay curves as a separate file
    # decay_curves is a dict: {gene_symbol: [{dist, fst, n_var}, ...]}
    d = load_json(RESULTS / "sweep_extent.json")
    decay_rows = []
    gene_meta = {g["gene_id"]: g for g in d.get("genes", [])}
    for symbol, bins in d.get("decay_curves", {}).items():
        for bin_d in bins:
            decay_rows.append({
                "gene_symbol":  symbol,
                "dist_kb":      _fmt(bin_d.get("dist", 0) / 1000, 0),
                "fst":          _fmt(bin_d.get("fst")),
                "n_variants":   bin_d.get("n_var", ""),
            })
    if decay_rows:
        write_tsv(SOURCE_DIR / "SourceData_SuppFig18_sweep_extent_decay_curves.tsv",
                  decay_rows,
                  ["gene_symbol", "dist_kb", "fst", "n_variants"],
                  "per-gene FST decay in 10 kb bins")


def source_supp_fig19_evo_trajectory():
    """Evolutionary trajectory embedding (backing fig19_evo_trajectory.png)."""
    rows = load_tsv(RESULTS / "evo_trajectory.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SOURCE_DIR / "SourceData_SuppFig19_evo_trajectory.tsv", rows, fields)


def source_benchmarks():
    """Benchmark performance table (backing main benchmarks figure)."""
    bench_file = BENCHMARKS / "benchmark_v4.jsonl"
    if not bench_file.exists():
        print(f"  WARNING: missing {bench_file}")
        return

    # Collect full-chromosome results (region == "full_chr22")
    rows = []
    with open(bench_file) as f:
        for line in f:
            rec = json.loads(line)
            if rec.get("region") != "full":
                continue
            base = {"region": "full_chr22", "chr": rec["chr"],
                    "start": rec.get("start", ""), "end": rec.get("end", "")}

            for stat_key, stat_label in [
                ("diversity",        "diversity_pi_thetaW_TajimaD"),
                ("differentiation",  "differentiation_fst_dxy"),
                ("sfs",              "sfs"),
                ("ihs",              "ihs"),
                ("xpehh",            "xpehh"),
            ]:
                stat = rec.get(stat_key, {})
                if not stat:
                    continue
                for tool, tdata in stat.items():
                    if not isinstance(tdata, dict) or "time_s" not in tdata:
                        continue
                    rows.append({
                        "statistic":  stat_label,
                        "tool":       tool,
                        "time_s":     _fmt(tdata["time_s"], 2),
                        "peak_mem_mb": _fmt(tdata.get("peak_mem_mb"), 1),
                        "n_variants": tdata.get("n_variants", ""),
                        "region":     rec["region"],
                    })

    fields = ["statistic", "tool", "time_s", "peak_mem_mb", "n_variants", "region"]
    write_tsv(SOURCE_DIR / "SourceData_Fig_benchmarks_full_chr22.tsv", rows, fields)


# ── Supplementary Tables ──────────────────────────────────────────────────────

def st01_population_summary():
    """ST01: 26-population summary table."""
    temporal = load_tsv(RESULTS / "temporal_selection.tsv")
    roh_sweep = {r["pop"]: r for r in load_tsv(RESULTS / "roh_sweep_correlation.tsv")}
    ooa = {r["pop"]: r for r in load_tsv(RESULTS / "ooa_gradient.tsv")}

    rows = []
    for r in temporal:
        pop = r["pop"]
        rs  = roh_sweep.get(pop, {})
        oo  = ooa.get(pop, {})
        rows.append({
            "population":          pop,
            "continental_group":   r.get("group", ""),
            "latitude":            _fmt(oo.get("lat"), 2),
            "longitude":           _fmt(oo.get("lon"), 2),
            "dist_ooa_km":         _fmt(oo.get("dist_ooa_km"), 0),
            "n_hard_sweeps":       rs.get("n_hard_sweeps", ""),
            "mean_pi":             _fmt(rs.get("mean_pi")),
            "mean_tajima_d":       _fmt(rs.get("mean_tajima_d"), 4),
            "mean_fst":            _fmt(rs.get("mean_fst")),
            "froh":                _fmt(r.get("froh")),
            "froh_z":              _fmt(r.get("froh_z"), 3),
            "h12_fraction":        _fmt(r.get("h12_fraction")),
            "h12_z":               _fmt(r.get("h12_z"), 3),
            "ihs_fraction":        _fmt(r.get("ihs_fraction")),
            "ihs_z":               _fmt(r.get("ihs_z"), 3),
            "dominant_timescale":  r.get("dominant_timescale", ""),
        })

    fields = ["population", "continental_group", "latitude", "longitude",
              "dist_ooa_km", "n_hard_sweeps", "mean_pi", "mean_tajima_d",
              "mean_fst", "froh", "froh_z", "h12_fraction", "h12_z",
              "ihs_fraction", "ihs_z", "dominant_timescale"]
    write_tsv(SUPP_DIR / "SupplementaryTable01_population_summary.tsv", rows, fields)


def st02_top_pathways():
    """ST02: Top 50 pathways ranked by mean Hudson FST."""
    rows = load_tsv(RESULTS / "pathway_fst.tsv")
    rows = [r for r in rows if r.get("mean_fst", "")]
    rows.sort(key=lambda r: float(r["mean_fst"]) if r["mean_fst"] else 0, reverse=True)
    rows = rows[:50]
    # Add rank
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    fields = ["rank", "pathway_id", "pathway_name", "n_variants", "mean_fst",
              "fst_YRI_vs_CEU", "fst_YRI_vs_CHB", "fst_CEU_vs_CHB",
              "fst_YRI_vs_JPT", "fst_CEU_vs_JPT"]
    write_tsv(SUPP_DIR / "SupplementaryTable02_top50_fst_pathways.tsv", rows, fields)


def st03_convergent_sweeps():
    """ST03: Convergent sweep genes (9 genes)."""
    rows = load_tsv(RESULTS / "convergent_sweeps.tsv")
    # Enrich with gene FST
    gene_fst = {r["gene_id"]: r for r in load_tsv(RESULTS / "gene_fst.tsv")}
    enriched = []
    for r in rows:
        # convergent_sweeps.tsv uses "pathway_id" as key, but gene data
        # may be in gene_fst; let's add gene-level FST if available
        gid = r.get("convergent_genes", "").split(",")[0].strip() if r.get("convergent_genes") else ""
        gf  = gene_fst.get(gid, {})
        enriched.append({
            "pathway_id":       r.get("pathway_id", ""),
            "pathway_name":     r.get("pathway_name", ""),
            "n_convergent_genes": r.get("n_convergent_genes", ""),
            "n_groups":         r.get("n_groups", ""),
            "continental_groups": r.get("continental_groups", ""),
            "convergent_genes": r.get("convergent_genes", ""),
            "fisher_pval":      _fmt(r.get("fisher_pval"), 5),
            "odds_ratio":       _fmt(r.get("odds_ratio"), 3),
            "gene_mean_fst":    _fmt(gf.get("mean_fst")),
            "gene_max_fst":     _fmt(gf.get("max_fst")),
            "max_h12":          _fmt(gf.get("max_h12"), 3),
        })
    fields = ["pathway_id", "pathway_name", "n_convergent_genes",
              "n_groups", "continental_groups", "convergent_genes",
              "fisher_pval", "odds_ratio", "gene_mean_fst", "gene_max_fst", "max_h12"]
    write_tsv(SUPP_DIR / "SupplementaryTable03_convergent_sweep_genes.tsv",
              enriched, fields)


def st04_top_genes_fst():
    """ST04: Top 100 genes by mean FST."""
    rows = load_tsv(RESULTS / "gene_fst.tsv")
    rows = [r for r in rows if r.get("mean_fst", "")]
    rows.sort(key=lambda r: float(r["mean_fst"]) if r["mean_fst"] else 0, reverse=True)
    rows = rows[:100]
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    fields = ["rank", "gene_id", "n_variants", "mean_fst", "max_fst",
              "fst_YRI_vs_CEU", "fst_YRI_vs_CHB", "fst_CEU_vs_CHB",
              "in_sweep", "n_conv_groups", "max_h12"]
    write_tsv(SUPP_DIR / "SupplementaryTable04_top100_genes_fst.tsv", rows, fields)


def st05_go_enrichment():
    """ST05: Significant GO terms (padj < 0.05) across all populations."""
    rows = load_tsv(RESULTS / "go_enrichment_ALL.tsv")
    sig = [r for r in rows if r.get("padj", "1") and float(r.get("padj", 1)) < 0.05]
    sig.sort(key=lambda r: float(r.get("padj", 1)))
    fields = list(sig[0].keys()) if sig else list(rows[0].keys()) if rows else []
    write_tsv(SUPP_DIR / "SupplementaryTable05_go_enrichment_significant.tsv",
              sig, fields, "padj < 0.05")


def st06_pbs_top_genes():
    """ST06: PBS top genes per population."""
    rows = load_tsv(RESULTS / "pbs_scan.tsv")
    # Unpack the wide-format top_gene_1/2/3 into long format
    long_rows = []
    for r in rows:
        pop   = r["pop"]
        group = r["group"]
        for k in range(1, 6):
            gene = r.get(f"top_gene_{k}", "")
            pbs  = r.get(f"pbs_{k}", "")
            if gene:
                long_rows.append({
                    "population":   pop,
                    "group":        group,
                    "rank":         k,
                    "gene_symbol":  gene,
                    "pbs_score":    _fmt(pbs, 5),
                    "pbs_threshold": _fmt(r.get("pbs_threshold"), 5),
                    "n_top_windows": r.get("n_top_windows", ""),
                })
    fields = ["population", "group", "rank", "gene_symbol", "pbs_score",
              "pbs_threshold", "n_top_windows"]
    write_tsv(SUPP_DIR / "SupplementaryTable06_pbs_top_genes_per_pop.tsv",
              long_rows, fields)


def st07_temporal_stratification():
    """ST07: Temporal selection stratification — all 26 populations."""
    rows = load_tsv(RESULTS / "temporal_selection.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SUPP_DIR / "SupplementaryTable07_temporal_stratification.tsv",
              rows, fields)


def st08_consequence_bias():
    """ST08: Functional consequence bias by VEP impact category."""
    rows = load_tsv(RESULTS / "consequence_bias.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SUPP_DIR / "SupplementaryTable08_consequence_bias.tsv", rows, fields)


def st09_pathway_communities():
    """ST09: Pathway co-selection top communities."""
    d = load_json(RESULTS / "pathway_coselection.json")
    rows = []
    for comm in d.get("top_communities", []):
        for p in comm.get("pathways", []):
            rows.append({
                "community_id":   comm["community_id"],
                "community_size": comm["size"],
                "community_mean_fst": _fmt(comm.get("mean_fst")),
                "pathway_id":     p["id"],
                "pathway_name":   p["name"],
            })
    fields = ["community_id", "community_size", "community_mean_fst",
              "pathway_id", "pathway_name"]
    write_tsv(SUPP_DIR / "SupplementaryTable09_pathway_coselection_communities.tsv",
              rows, fields)


def st10_cross_species():
    """ST10: Cross-species functional convergence by module."""
    rows = load_tsv(RESULTS / "cross_species.tsv")
    fields = list(rows[0].keys()) if rows else []
    write_tsv(SUPP_DIR / "SupplementaryTable10_cross_species_convergence.tsv",
              rows, fields)


def st11_rare_burden():
    """ST11: Top 50 genes by total population-specific rare variant burden."""
    d = load_json(RESULTS / "rare_burden.json")
    rows = []
    pop_summary = d.get("pop_summary", {})
    for entry in d.get("top_genes", [])[:50]:
        rows.append({
            "rank":               len(rows) + 1,
            "gene_id":            entry["gene_id"],
            "gene_symbol":        entry["symbol"],
            "total_private_rare": entry["total_private_rare"],
        })
    fields = ["rank", "gene_id", "gene_symbol", "total_private_rare"]
    write_tsv(SUPP_DIR / "SupplementaryTable11_rare_burden_top50_genes.tsv",
              rows, fields)

    # Also write population-level burden summary
    pop_rows = []
    for pop, stats in pop_summary.items():
        pop_rows.append({
            "population":   pop,
            "mean_rare_per_gene": _fmt(stats.get("mean_rare"), 2),
            "froh":         _fmt(stats.get("froh")),
            "froh_z":       _fmt(stats.get("froh_z"), 3),
        })
    pop_rows.sort(key=lambda r: float(r["mean_rare_per_gene"]) if r["mean_rare_per_gene"] else 0,
                  reverse=True)
    write_tsv(SUPP_DIR / "SupplementaryTable11b_rare_burden_by_population.tsv",
              pop_rows, ["population", "mean_rare_per_gene", "froh", "froh_z"])


def st12_sweep_extent():
    """ST12: Convergent sweep genomic extents."""
    rows = load_tsv(RESULTS / "sweep_extent.tsv")
    d    = load_json(RESULTS / "sweep_extent.json")
    gene_notes = {g["gene_id"]: g for g in d.get("genes", [])}

    enriched = []
    for r in rows:
        note = gene_notes.get(r.get("gene_id", ""), {})
        enriched.append({
            "gene_id":        r.get("gene_id", ""),
            "gene_symbol":    r.get("gene_symbol", ""),
            "chr":            r.get("chr", ""),
            "gene_start":     r.get("gene_start", ""),
            "gene_end":       r.get("gene_end", ""),
            "sweep_pop":      r.get("sweep_pop", ""),
            "ref_pop":        r.get("ref_pop", ""),
            "peak_fst":       _fmt(r.get("peak_fst")),
            "extent_bp":      r.get("extent_bp", ""),
            "extent_kb":      _fmt(int(r["extent_bp"]) / 1000, 0) if r.get("extent_bp") else "",
            "max_h12":        _fmt(r.get("max_h12"), 3),
            "n_sweep_pops":   note.get("n_sweep_pops", ""),
        })
    fields = ["gene_id", "gene_symbol", "chr", "gene_start", "gene_end",
              "sweep_pop", "ref_pop", "peak_fst", "extent_bp", "extent_kb",
              "max_h12", "n_sweep_pops"]
    write_tsv(SUPP_DIR / "SupplementaryTable12_sweep_extents.tsv", enriched, fields)


def st13_benchmark_summary():
    """ST13: Benchmark performance summary (speedup vs scikit-allel)."""
    bench_file = BENCHMARKS / "benchmark_v4.jsonl"
    if not bench_file.exists():
        print(f"  WARNING: missing {bench_file}")
        return

    # Collect full-chr22 results
    stat_times = defaultdict(dict)  # stat_key → tool → time_s
    stat_mem   = defaultdict(dict)

    with open(bench_file) as f:
        for line in f:
            rec = json.loads(line)
            if rec.get("region") != "full":
                continue
            for stat_key in ["diversity", "differentiation", "sfs", "ihs", "xpehh"]:
                stat = rec.get(stat_key, {})
                if not stat:
                    continue
                for tool, tdata in stat.items():
                    if isinstance(tdata, dict) and "time_s" in tdata:
                        stat_times[stat_key][tool] = tdata["time_s"]
                        stat_mem[stat_key][tool]   = tdata.get("peak_mem_mb", None)

    rows = []
    stat_labels = {
        "diversity":       "π / θ_W / Tajima's D",
        "differentiation": "FST / Dxy",
        "sfs":             "SFS",
        "ihs":             "iHS",
        "xpehh":           "XP-EHH",
    }

    for stat_key, label in stat_labels.items():
        times = stat_times.get(stat_key, {})
        mems  = stat_mem.get(stat_key, {})
        sa_time = times.get("scikit-allel")
        for tool in ["GraphPop", "GraphPop-numpy", "scikit-allel", "vcftools",
                     "PLINK2", "PLINK2-pgen", "bcftools", "selscan"]:
            t = times.get(tool)
            if t is None:
                continue
            speedup = round(sa_time / t, 1) if sa_time and t > 0 else ""
            rows.append({
                "statistic":    label,
                "tool":         tool,
                "time_s":       _fmt(t, 2),
                "peak_mem_mb":  _fmt(mems.get(tool), 1),
                "speedup_vs_scikit_allel": str(speedup) if speedup != "" else "",
                "chromosome":   "chr22 (full)",
            })

    fields = ["statistic", "tool", "time_s", "peak_mem_mb",
              "speedup_vs_scikit_allel", "chromosome"]
    write_tsv(SUPP_DIR / "SupplementaryTable13_benchmark_performance.tsv", rows, fields)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=== Generating Source Data and Supplementary Tables ===\n")

    print("── Source Data files ──────────────────────────────────")
    source_supp_fig01_pathway_fst()
    source_supp_fig02_gene_fst_landscape()
    source_supp_fig03_convergent_sweeps()
    source_supp_fig04_roh_sweep()
    source_supp_fig05_go_enrichment()
    source_supp_fig06_pca()
    source_supp_fig07_ooa_gradient()
    source_supp_fig08_synthesis()
    source_supp_fig09_summary_overview()
    source_supp_fig10_gene_fst()
    source_supp_fig11_pbs()
    source_supp_fig12_temporal()
    source_supp_fig13_kcne1_overview()
    source_supp_fig14_consequence_bias()
    source_supp_fig15_pathway_coselection()
    source_supp_fig16_cross_species()
    source_supp_fig17_rare_burden()
    source_supp_fig18_sweep_extent()
    source_supp_fig19_evo_trajectory()
    source_benchmarks()

    print("\n── Supplementary Tables ───────────────────────────────")
    st01_population_summary()
    st02_top_pathways()
    st03_convergent_sweeps()
    st04_top_genes_fst()
    st05_go_enrichment()
    st06_pbs_top_genes()
    st07_temporal_stratification()
    st08_consequence_bias()
    st09_pathway_communities()
    st10_cross_species()
    st11_rare_burden()
    st12_sweep_extent()
    st13_benchmark_summary()

    print("\n── Summary ─────────────────────────────────────────────")
    src_files  = list(SOURCE_DIR.glob("*.tsv"))
    supp_files = list(SUPP_DIR.glob("*.tsv"))
    print(f"  Source Data files:     {len(src_files):>3}  →  {SOURCE_DIR}")
    print(f"  Supplementary Tables:  {len(supp_files):>3}  →  {SUPP_DIR}")


if __name__ == "__main__":
    main()
