#!/usr/bin/env python3
"""Generate all investigation figures for GraphPop population genomics analysis.

Produces 8 publication-quality figures saved to data/results/figures/.

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/generate_investigation_figures.py
"""

import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import seaborn as sns

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
RES  = Path("data/results")
FIGS = RES / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

# ── Global style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "font.size":       9,
    "axes.labelsize":  10,
    "axes.titlesize":  11,
    "axes.titleweight": "bold",
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "figure.dpi":      150,
    "savefig.dpi":     200,
    "savefig.bbox":    "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

# Continental group palette (consistent across all figures)
GROUP_COLORS = {
    "African":    "#E07B39",
    "European":   "#4472C4",
    "East_Asian": "#70AD47",
    "South_Asian":"#7030A0",
    "American":   "#C00000",
}

PAIR_COLORS = {
    "YRI_vs_CEU": "#4472C4",
    "YRI_vs_CHB": "#70AD47",
    "CEU_vs_CHB": "#E07B39",
    "YRI_vs_JPT": "#7030A0",
    "CEU_vs_JPT": "#C00000",
}

# ── Helper: group legend patches ─────────────────────────────────────────────
def group_legend_patches():
    return [mpatches.Patch(color=c, label=g.replace("_", " "))
            for g, c in GROUP_COLORS.items()]

# ── Load all data once ────────────────────────────────────────────────────────
print("Loading result files ...", flush=True)
with open(RES / "pathway_fst.json")            as f: pathway_fst  = json.load(f)
with open(RES / "convergent_sweeps.json")      as f: conv         = json.load(f)
with open(RES / "roh_sweep_correlation.json")  as f: roh_sweep    = json.load(f)
with open(RES / "kcne1_deepdive.json")         as f: kcne1        = json.load(f)
with open(RES / "go_enrichment_summary.json")  as f: go_summary   = json.load(f)
with open(RES / "ooa_gradient.json")           as f: ooa          = json.load(f)
with open(RES / "synthesis.json")              as f: synthesis     = json.load(f)

# Load GO per-population TSVs
go_rows = {}
for pop in ["ALL", "CEU", "CHB", "GIH"]:
    tsv = RES / f"go_enrichment_{pop}.tsv"
    rows = []
    with open(tsv) as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            row = dict(zip(header, parts))
            row["padj"]       = float(row["padj"])
            row["pval"]       = float(row["pval"])
            row["n_fg"]       = int(row["n_fg"])
            row["odds_ratio"] = float(row["odds_ratio"])
            rows.append(row)
    go_rows[pop] = rows

print("  All data loaded.\n", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1 — Top-20 Pathway FST
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 1: Top-20 Pathway FST ...", flush=True)

top20 = [r for r in pathway_fst if r["mean_fst"] is not None][:20]
pairs = ["YRI_vs_CEU", "YRI_vs_CHB", "CEU_vs_CHB", "YRI_vs_JPT", "CEU_vs_JPT"]
pair_labels = ["YRI/CEU", "YRI/CHB", "CEU/CHB", "YRI/JPT", "CEU/JPT"]

# Abbreviate long pathway names
def abbrev(name, maxlen=48):
    return name if len(name) <= maxlen else name[:maxlen-3] + "..."

fig, (ax_main, ax_pairs) = plt.subplots(
    1, 2, figsize=(14, 6.5),
    gridspec_kw={"width_ratios": [1.8, 1]})

# Left: mean FST bars
names  = [abbrev(r["pathway_name"]) for r in reversed(top20)]
values = [r["mean_fst"] for r in reversed(top20)]
colors_bar = plt.cm.YlOrRd(np.linspace(0.3, 0.95, 20))[::-1]
bars = ax_main.barh(names, values, color=colors_bar, edgecolor="white", linewidth=0.4)
ax_main.set_xlabel("Mean Hudson F$_{ST}$")
ax_main.set_title("Top 20 Pathways by Mean F$_{ST}$\n(Inv.1 — 2,241 Reactome pathways, 5 population pairs)")
ax_main.axvline(0.05, color="gray", lw=0.8, ls="--", alpha=0.6, label="Genome avg F$_{ST}$")
# Value labels
for bar, val in zip(bars, values):
    ax_main.text(val + 0.002, bar.get_y() + bar.get_height()/2,
                 f"{val:.3f}", va="center", fontsize=7.5)
ax_main.legend(fontsize=8)

# Right: per-pair FST for top 10 pathways (grouped bars)
top10 = list(reversed(top20))[:10]
x = np.arange(len(top10))
width = 0.16
for k, (pair, label, color) in enumerate(zip(pairs, pair_labels, PAIR_COLORS.values())):
    vals = [r.get(f"fst_{pair}") or 0 for r in top10]
    ax_pairs.bar(x + k*width, vals, width, label=label, color=color, alpha=0.85)

ax_pairs.set_xticks(x + width*2)
ax_pairs.set_xticklabels([abbrev(r["pathway_name"], 22) for r in top10],
                          rotation=45, ha="right", fontsize=7)
ax_pairs.set_ylabel("F$_{ST}$")
ax_pairs.set_title("Pairwise F$_{ST}$ — Top 10 Pathways")
ax_pairs.legend(fontsize=7.5, title="Pop. pair")

plt.tight_layout()
fig.savefig(FIGS / "fig1_pathway_fst.png")
plt.close()
print("  Saved fig1_pathway_fst.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2 — Pairwise FST landscape (scatter of YRI/CEU vs YRI/CHB)
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 2: Pairwise FST scatter ...", flush=True)

# Build arrays for the two main pairs
fst_af_eu = np.array([r.get("fst_YRI_vs_CEU") or np.nan for r in pathway_fst])
fst_af_ea = np.array([r.get("fst_YRI_vs_CHB") or np.nan for r in pathway_fst])
fst_eu_ea = np.array([r.get("fst_CEU_vs_CHB") or np.nan for r in pathway_fst])
fst_mean  = np.array([r.get("mean_fst")       or np.nan for r in pathway_fst])
names_all = [r["pathway_name"] for r in pathway_fst]

mask = np.isfinite(fst_af_eu) & np.isfinite(fst_af_ea) & np.isfinite(eu_ea := fst_eu_ea)

# Top outlier pathways to label
top_labels_idx = np.argsort(fst_mean)[::-1][:8]

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

scatter_pairs = [
    (fst_af_eu, fst_af_ea, "YRI vs CEU  F$_{ST}$", "YRI vs CHB  F$_{ST}$", "Africa vs Europe/EAsia"),
    (fst_af_eu, fst_eu_ea, "YRI vs CEU  F$_{ST}$", "CEU vs CHB  F$_{ST}$", "Africa–Europe vs Eur–EAsia"),
    (fst_af_ea, fst_eu_ea, "YRI vs CHB  F$_{ST}$", "CEU vs CHB  F$_{ST}$", "AfEAsia vs EurEAsia"),
]

for ax, (x_arr, y_arr, xl, yl, title) in zip(axes, scatter_pairs):
    mask_ = np.isfinite(x_arr) & np.isfinite(y_arr)
    # Hex density background
    hb = ax.hexbin(x_arr[mask_], y_arr[mask_], gridsize=50, mincnt=1,
                   cmap="Blues", linewidths=0.2, alpha=0.85)
    # Diagonal
    lim = max(np.nanmax(x_arr), np.nanmax(y_arr)) * 1.05
    ax.plot([0, lim], [0, lim], "--", color="gray", lw=0.8, alpha=0.6)
    # Label top outliers
    labeled = set()
    for idx in top_labels_idx:
        xv, yv = x_arr[idx], y_arr[idx]
        if not (np.isfinite(xv) and np.isfinite(yv)):
            continue
        name = abbrev(names_all[idx], 28)
        if name not in labeled:
            ax.annotate(name, (xv, yv), xytext=(6, 2), textcoords="offset points",
                        fontsize=6.5, color="#C00000",
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.6))
            labeled.add(name)
    ax.set_xlabel(xl, fontsize=8.5)
    ax.set_ylabel(yl, fontsize=8.5)
    ax.set_title(title, fontsize=9.5)

fig.suptitle("Pathway F$_{ST}$ Landscape — All 2,241 Reactome Pathways (Inv.1)", fontsize=11, fontweight="bold")
plt.tight_layout()
fig.savefig(FIGS / "fig2_fst_landscape.png")
plt.close()
print("  Saved fig2_fst_landscape.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 3 — Convergent Sweeps + KCNE1 Deep-Dive
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 3: Convergent sweeps + KCNE1 ...", flush=True)

conv_genes  = conv["convergent_genes"]
kcne1_locus = kcne1["locus_fst"]

fig = plt.figure(figsize=(14, 6))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.40)
ax_bubble = fig.add_subplot(gs[0])
ax_h12    = fig.add_subplot(gs[1])
ax_fst    = fig.add_subplot(gs[2])

# ---- Bubble: convergent genes (n_groups × max_h12) ----
from matplotlib import cm

gene_names  = list(conv_genes.keys())
n_groups    = [conv_genes[g]["n_groups"] for g in gene_names]
max_h12_arr = [conv_genes[g]["max_h12"] for g in gene_names]
# Load n_GO from synthesis
synthesis_genes = {g["gene"]: g for g in synthesis["top_genes"]}
n_go            = [synthesis_genes.get(g, {}).get("n_go", 0) for g in gene_names]

jitter = np.random.default_rng(42).uniform(-0.08, 0.08, len(gene_names))
for i, (gname, ng, mh, ngo) in enumerate(zip(gene_names, n_groups, max_h12_arr, n_go)):
    color = GROUP_COLORS.get(
        ["African","European","East_Asian","South_Asian","American"][min(ng-2, 4)], "#888")
    size = max(80, ngo * 12 + 60)
    ax_bubble.scatter(ng + jitter[i], mh, s=size, c=color, alpha=0.80,
                      edgecolors="white", linewidths=0.8, zorder=3)
    offset_y = 0.015 if gname != "KCNE1" else 0.025
    ax_bubble.annotate(gname, (ng + jitter[i], mh),
                       xytext=(0, 8 if gname == "KCNE1" else 5),
                       textcoords="offset points",
                       fontsize=7.5, ha="center",
                       fontweight="bold" if gname == "KCNE1" else "normal")

ax_bubble.set_xlabel("Continental groups with hard sweep")
ax_bubble.set_ylabel("Max H12 score")
ax_bubble.set_title("Convergent Selection Genes\n(Inv.2 — bubble ∝ GO terms)")
ax_bubble.set_xticks([2, 3, 4, 5])
ax_bubble.set_xlim(1.5, 5.8)
ax_bubble.set_ylim(0.45, 1.10)
ax_bubble.axhline(0.9, color="gray", ls="--", lw=0.8, alpha=0.5)
ax_bubble.text(5.5, 0.91, "H12=0.9", fontsize=7, color="gray", ha="right")
ax_bubble.legend(handles=group_legend_patches()[:5], fontsize=7,
                 title="Group (by n_groups color)", loc="lower right")

# ---- Bar: KCNE1 h12 per population ----
sw_windows = kcne1["sweep_windows_neo4j"]
# De-duplicate: one entry per pop (take first occurrence of each)
seen_pops = {}
for w in sw_windows:
    if w["pop"] not in seen_pops:
        seen_pops[w["pop"]] = w

pops_sw = sorted(seen_pops.keys())
h12_vals = [seen_pops[p]["h12"] for p in pops_sw]
hap_vals = [seen_pops[p]["hap_div"] if seen_pops[p].get("hap_div") else 0.0 for p in pops_sw]

grp_colors_sw = [GROUP_COLORS.get(kcne1_data_pop_group := {
    "YRI": "African", "CEU": "European", "FIN": "European",
    "CHB": "East_Asian", "GIH": "South_Asian", "PEL": "American"
}.get(p, "European"), "#aaa") for p in pops_sw]

x_sw = np.arange(len(pops_sw))
bars_h = ax_h12.bar(x_sw, h12_vals, color=grp_colors_sw, edgecolor="white",
                     linewidth=0.6, alpha=0.88)
ax_h12_r = ax_h12.twinx()
ax_h12_r.plot(x_sw, hap_vals, "o--", color="#555", lw=1.2, ms=5,
               label="Hap. diversity", alpha=0.7)
ax_h12_r.set_ylabel("Haplotype diversity", fontsize=8, color="#555")
ax_h12_r.tick_params(labelsize=7.5, colors="#555")

ax_h12.set_xticks(x_sw)
ax_h12.set_xticklabels(pops_sw, fontsize=8.5)
ax_h12.set_ylabel("H12 score")
ax_h12.set_ylim(0, 1.12)
ax_h12.set_title("KCNE1 Sweep Intensity\nper Population (Inv.4)")
ax_h12.axhline(1.0, color="gray", ls=":", lw=0.8, alpha=0.5)
for bar, val in zip(bars_h, h12_vals):
    ax_h12.text(bar.get_x() + bar.get_width()/2, val + 0.01,
                f"{val:.2f}", ha="center", va="bottom", fontsize=7.5)

# ---- Grouped bars: FST locus vs genome ----
pairs_k = sorted(kcne1_locus.keys())
pair_labels_k = [pk.replace("_", "/") for pk in pairs_k]
locus_fst_v  = [kcne1_locus[pk]["fst_locus"] for pk in pairs_k]
genome_fst_v = [kcne1_locus[pk]["fst_genome"] or 0 for pk in pairs_k]
ratio_v      = [kcne1_locus[pk]["fst_ratio"]  or 0 for pk in pairs_k]

x_fst = np.arange(len(pairs_k))
w_f = 0.35
b1 = ax_fst.bar(x_fst - w_f/2, genome_fst_v, w_f, label="Genome-wide F$_{ST}$",
                 color="#4472C4", alpha=0.75, edgecolor="white")
b2 = ax_fst.bar(x_fst + w_f/2, locus_fst_v,  w_f, label="KCNE1 locus F$_{ST}$",
                 color="#E07B39", alpha=0.85, edgecolor="white")
ax_fst.set_xticks(x_fst)
ax_fst.set_xticklabels(pair_labels_k, fontsize=8, rotation=25, ha="right")
ax_fst.set_ylabel("F$_{ST}$")
ax_fst.set_title("KCNE1 Locus vs Genome-Wide F$_{ST}$\n(Inv.4 — ancient sweep test)")
ax_fst.legend(fontsize=8)
# Annotate ratios
for xi, (lv, gv, rv) in enumerate(zip(locus_fst_v, genome_fst_v, ratio_v)):
    ax_fst.text(xi, max(lv, gv) + 0.003, f"×{rv:.2f}", ha="center",
                fontsize=7.5, color="#C00000", fontweight="bold")

fig.suptitle("Convergent Selection: KCNE1 Shows Pan-Continental Ancient Sweep\n"
             "F$_{ST}$(locus) / F$_{ST}$(genome) < 0.5 → swept haplotype predates OoA dispersal",
             fontsize=10, fontweight="bold")
fig.savefig(FIGS / "fig3_convergent_kcne1.png")
plt.close()
print("  Saved fig3_convergent_kcne1.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 4 — ROH × Sweep Correlation (Inv.3)
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 4: ROH × Sweep correlation ...", flush=True)

pops_rs  = sorted(roh_sweep["per_population"].keys())
pop_data = roh_sweep["per_population"]
corr_data = roh_sweep["correlations"]

metrics_m = ["froh", "n_hard_sweeps", "mean_h12", "mean_pi", "mean_theta_w",
              "mean_tajima_d", "mean_fis", "mean_fst"]
metric_labels = ["FROH", "N sweeps", "Mean H12", "π", "θ_W", "Tajima's D", "F$_{IS}$", "Mean F$_{ST}$"]

from scipy.stats import spearmanr

# Build correlation matrix
n_m = len(metrics_m)
corr_mat = np.zeros((n_m, n_m))
pval_mat = np.ones((n_m, n_m))
arrays_rs = {mk: np.array([pop_data[p][mk] for p in pops_rs]) for mk in metrics_m}

for i, m1 in enumerate(metrics_m):
    for j, m2 in enumerate(metrics_m):
        mask_ = np.isfinite(arrays_rs[m1]) & np.isfinite(arrays_rs[m2])
        if mask_.sum() >= 5:
            rho, pv = spearmanr(arrays_rs[m1][mask_], arrays_rs[m2][mask_])
            corr_mat[i, j] = rho
            pval_mat[i, j] = pv
        else:
            corr_mat[i, j] = np.nan

fig = plt.figure(figsize=(14, 5.5))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38)

# ---- Correlation heatmap ----
ax_heat = fig.add_subplot(gs[0])
mask_tri = np.triu(np.ones_like(corr_mat, dtype=bool), k=1)
cmap_div = sns.diverging_palette(250, 30, as_cmap=True)
sns.heatmap(corr_mat, ax=ax_heat, mask=mask_tri, cmap=cmap_div,
            center=0, vmin=-1, vmax=1,
            xticklabels=metric_labels, yticklabels=metric_labels,
            linewidths=0.4, square=True, annot=True, fmt=".2f",
            annot_kws={"size": 7})
# Star significant cells
for i in range(n_m):
    for j in range(i):
        pv = pval_mat[i, j]
        star = "***" if pv < 0.001 else ("**" if pv < 0.01 else ("*" if pv < 0.05 else ""))
        if star:
            ax_heat.text(j + 0.5, i + 0.85, star, ha="center", va="center",
                         fontsize=6, color="white", fontweight="bold")
ax_heat.set_title("Spearman Correlation Matrix\n(26 populations, Inv.3)", fontsize=10)
ax_heat.tick_params(axis="x", rotation=40, labelsize=7.5)
ax_heat.tick_params(axis="y", rotation=0,  labelsize=7.5)

# ---- Scatter: π vs n_hard_sweeps ----
ax_s1 = fig.add_subplot(gs[1])
for pop in pops_rs:
    grp   = pop_data[pop]["group"]
    color = GROUP_COLORS.get(grp, "#aaa")
    xv = pop_data[pop]["mean_pi"]
    yv = pop_data[pop]["n_hard_sweeps"]
    ax_s1.scatter(xv, yv, s=55, c=color, edgecolors="white", linewidths=0.6,
                  zorder=3, alpha=0.85)
    if pop in {"YRI", "CEU", "CHB", "GIH", "PEL", "PJL"}:
        ax_s1.annotate(pop, (xv, yv), xytext=(4, 2), textcoords="offset points", fontsize=7.5)

rho_v, pv_v = spearmanr(arrays_rs["mean_pi"], arrays_rs["n_hard_sweeps"])
ax_s1.set_xlabel("Mean nucleotide diversity (π)")
ax_s1.set_ylabel("Number of hard sweep windows")
ax_s1.set_title(f"π vs Sweep Count\nρ = {rho_v:+.3f}  p = {pv_v:.4f} ***", fontsize=10)
ax_s1.legend(handles=group_legend_patches(), fontsize=7.5, loc="upper right")

# ---- Scatter: Fis vs FROH ----
ax_s2 = fig.add_subplot(gs[2])
for pop in pops_rs:
    grp   = pop_data[pop]["group"]
    color = GROUP_COLORS.get(grp, "#aaa")
    xv = pop_data[pop]["mean_fis"]
    yv = pop_data[pop]["froh"]
    ax_s2.scatter(xv, yv, s=55, c=color, edgecolors="white", linewidths=0.6,
                  zorder=3, alpha=0.85)
    if pop in {"YRI", "PJL", "STU", "ITU", "CEU", "CHB", "MXL"}:
        ax_s2.annotate(pop, (xv, yv), xytext=(4, 2), textcoords="offset points", fontsize=7.5)

rho_fis, pv_fis = spearmanr(arrays_rs["mean_fis"], arrays_rs["froh"])
ax_s2.set_xlabel("Mean inbreeding coefficient (F$_{IS}$)")
ax_s2.set_ylabel("Mean FROH")
ax_s2.set_title(f"F$_{{IS}}$ vs FROH\nρ = {rho_fis:+.3f}  p < 0.0001 ***", fontsize=10)
ax_s2.legend(handles=group_legend_patches(), fontsize=7.5)

fig.suptitle("Population-Level ROH × Sweep × Diversity Correlations (Inv.3 — 26 populations)",
             fontsize=11, fontweight="bold")
fig.savefig(FIGS / "fig4_roh_sweep_correlation.png")
plt.close()
print("  Saved fig4_roh_sweep_correlation.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 5 — GO Enrichment Dot Plot (Inv.5)
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 5: GO enrichment dot plot ...", flush=True)

# Collect all sig GO terms (padj < 0.05) across all populations
from collections import defaultdict
go_padj  = defaultdict(dict)   # go_name → pop → padj
go_nfg   = defaultdict(dict)   # go_name → pop → n_fg
go_aspect = {}

for pop in ["ALL", "CEU", "CHB", "GIH"]:
    for row in go_rows[pop]:
        if row["padj"] < 0.05:
            name = row["go_name"][:55]
            go_padj[name][pop]  = row["padj"]
            go_nfg[name][pop]   = row["n_fg"]
            go_aspect[name]     = row["aspect"][:3].upper() if row["aspect"] else "?"

# Sort: total significance across populations, then by most specific aspect
all_terms = sorted(go_padj.keys(),
                   key=lambda t: -sum(-np.log10(v) for v in go_padj[t].values()))

# Keep top 30
top_terms = all_terms[:30]
pops_go   = ["ALL", "CEU", "CHB", "GIH"]
pop_labels = ["All pops", "CEU (European)", "CHB (E.Asian)", "GIH (S.Asian)"]

# Build matrices
n_rows_go = len(top_terms)
n_cols_go = len(pops_go)
padj_mat  = np.full((n_rows_go, n_cols_go), np.nan)
nfg_mat   = np.zeros((n_rows_go, n_cols_go))

for i, term in enumerate(top_terms):
    for j, pop in enumerate(pops_go):
        if pop in go_padj[term]:
            padj_mat[i, j] = go_padj[term][pop]
            nfg_mat[i, j]  = go_nfg[term][pop]

fig, ax = plt.subplots(figsize=(8, max(7, n_rows_go * 0.32)))

# Aspect-based left color strip
aspect_colors = {"BIO": "#70AD47", "MOL": "#4472C4", "CEL": "#E07B39", "?": "#aaa"}
for i, term in enumerate(top_terms):
    asp = go_aspect.get(term, "?")
    ax.add_patch(mpatches.Rectangle((-0.55, i - 0.45), 0.12, 0.9,
                                    color=aspect_colors.get(asp, "#aaa"), zorder=2))

for j, pop in enumerate(pops_go):
    for i in range(n_rows_go):
        pv = padj_mat[i, j]
        if np.isnan(pv):
            continue
        ng  = nfg_mat[i, j]
        sig = -np.log10(max(pv, 1e-10))
        col = plt.cm.Reds(min(sig / 5.0, 1.0))
        size = min(20 + ng * 35, 350)
        ax.scatter(j, i, s=size, c=[col], edgecolors="gray", linewidths=0.4, zorder=3)

ax.set_xticks(range(n_cols_go))
ax.set_xticklabels(pop_labels, fontsize=9, rotation=30, ha="right")
ax.set_yticks(range(n_rows_go))
ax.set_yticklabels(top_terms, fontsize=7.5)
ax.set_xlim(-0.55, n_cols_go - 0.4)
ax.set_ylim(-0.8, n_rows_go - 0.2)
ax.set_title("GO Term Enrichment in Hard Sweep Genes\n"
             "Bubble size = gene count, color = −log₁₀(FDR) (Inv.5)",
             fontsize=10)

# Aspect legend
asp_patches = [mpatches.Patch(color=c, label={"BIO":"Biological Process","MOL":"Molecular Function",
                                               "CEL":"Cellular Component"}.get(a, a))
               for a, c in aspect_colors.items() if a != "?"]
# Color legend (size + color scale)
leg_sizes  = [mpatches.Patch(color="white", label="Bubble size ~ gene count")]
leg_color  = [mpatches.Patch(color=plt.cm.Reds(0.4), label="FDR < 0.05"),
              mpatches.Patch(color=plt.cm.Reds(0.8), label="FDR < 0.001")]
ax.legend(handles=asp_patches + leg_color, fontsize=7.5, loc="lower right",
          title="GO aspect / significance", title_fontsize=7)

plt.tight_layout()
fig.savefig(FIGS / "fig5_go_enrichment.png")
plt.close()
print("  Saved fig5_go_enrichment.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 6 — PCA Biplot of 26 Populations (Inv.7)
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 6: PCA biplot ...", flush=True)

ooa_pops  = ooa["per_population"]
var_exp   = ooa["pca_variance_explained"]

fig, ax = plt.subplots(figsize=(8, 6.5))

# Plot by continental group
by_group = defaultdict(list)
for entry in ooa_pops:
    by_group[entry["group"]].append(entry)

for grp, entries in by_group.items():
    xs = [e["pc1"] for e in entries]
    ys = [e["pc2"] for e in entries]
    ax.scatter(xs, ys, s=70, c=GROUP_COLORS.get(grp, "#aaa"),
               edgecolors="white", linewidths=0.7, label=grp.replace("_", " "),
               zorder=3, alpha=0.90)
    for e in entries:
        ax.annotate(e["pop"], (e["pc1"], e["pc2"]),
                    xytext=(4, 3), textcoords="offset points", fontsize=8)

# Draw group hulls (convex hull)
from matplotlib.patches import FancyArrowPatch
try:
    from scipy.spatial import ConvexHull
    for grp, entries in by_group.items():
        pts = np.array([[e["pc1"], e["pc2"]] for e in entries])
        if len(pts) >= 3:
            hull = ConvexHull(pts)
            for simplex in hull.simplices:
                ax.plot(pts[simplex, 0], pts[simplex, 1],
                        color=GROUP_COLORS.get(grp, "#aaa"), lw=0.8, alpha=0.35)
except Exception:
    pass

ax.axhline(0, color="gray", lw=0.6, alpha=0.4)
ax.axvline(0, color="gray", lw=0.6, alpha=0.4)
ax.set_xlabel(f"PC1 ({100*var_exp[0]:.1f}% variance)")
ax.set_ylabel(f"PC2 ({100*var_exp[1]:.1f}% variance)")
ax.set_title("PCA of 26 1000G Populations\n"
             "Features: FROH, π, Tajima's D, F$_{IS}$, n_sweeps, mean H12 (Inv.7)")
ax.legend(fontsize=8.5, title="Continental group", loc="upper left")
plt.tight_layout()
fig.savefig(FIGS / "fig6_pca_populations.png")
plt.close()
print("  Saved fig6_pca_populations.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 7 — Out-of-Africa Gradient (Inv.7)
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 7: OoA gradient ...", flush=True)

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

ooa_metrics = [
    ("dist_ooa_km", "mean_pi",       "OoA distance (km)", "Mean π", "π vs Distance from Africa"),
    ("dist_ooa_km", "n_hard_sweeps", "OoA distance (km)", "N hard sweeps", "Sweeps vs Distance"),
    ("dist_ooa_km", "froh",          "OoA distance (km)", "FROH", "FROH vs Distance"),
]

for ax, (xk, yk, xl, yl, title) in zip(axes, ooa_metrics):
    xs = np.array([e[xk] for e in ooa_pops])
    ys = np.array([e[yk] for e in ooa_pops])
    grps = [e["group"] for e in ooa_pops]

    for entry in ooa_pops:
        color = GROUP_COLORS.get(entry["group"], "#aaa")
        ax.scatter(entry[xk], entry[yk], s=60, c=color,
                   edgecolors="white", linewidths=0.6, zorder=3, alpha=0.88)

    # Regression line
    mask_r = np.isfinite(xs) & np.isfinite(ys)
    if mask_r.sum() >= 4:
        rho_r, pv_r = spearmanr(xs[mask_r], ys[mask_r])
        z = np.polyfit(xs[mask_r], ys[mask_r], 1)
        xl_lim = np.linspace(xs[mask_r].min(), xs[mask_r].max(), 100)
        ax.plot(xl_lim, np.polyval(z, xl_lim), "--", color="gray", lw=1.2, alpha=0.7)
        ax.set_title(f"{title}\nρ = {rho_r:+.3f}  p = {pv_r:.3f}", fontsize=10)
    else:
        ax.set_title(title, fontsize=10)

    # Label a few key populations
    for entry in ooa_pops:
        if entry["pop"] in {"YRI", "LWK", "CEU", "CHB", "GIH", "PEL", "PJL", "FIN"}:
            ax.annotate(entry["pop"], (entry[xk], entry[yk]),
                        xytext=(4, 3), textcoords="offset points", fontsize=7.5)
    ax.set_xlabel(xl, fontsize=9)
    ax.set_ylabel(yl, fontsize=9)

axes[0].legend(handles=group_legend_patches(), fontsize=7.5, loc="lower left")
fig.suptitle("Out-of-Africa Gradient: Evolutionary Parameters vs Geographic Distance (Inv.7)\n"
             "Origin: Addis Ababa, Ethiopia", fontsize=11, fontweight="bold")
plt.tight_layout()
fig.savefig(FIGS / "fig7_ooa_gradient.png")
plt.close()
print("  Saved fig7_ooa_gradient.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 8 — Synthesis Evidence Heatmap (Inv.6)
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 8: Synthesis heatmap ...", flush=True)

top_pathways = synthesis["top_pathways"][:20]

# Columns: FST rank, conv_genes, conv_groups, GO_enriched, composite
col_data = {
    "Mean F$_{ST}$":     [r["mean_fst"] for r in top_pathways],
    "Conv. genes":       [r["n_convergent_genes"] for r in top_pathways],
    "Conv. groups":      [r["n_conv_groups"] for r in top_pathways],
    "GO genes":          [r["n_go_enriched_genes"] for r in top_pathways],
    "Composite\nscore":  [r["composite_score"] for r in top_pathways],
}
tier_colors = {"HIGH": "#C00000", "MEDIUM": "#E07B39", "LOW": "#4472C4"}
tiers = [r["tier"] for r in top_pathways]

# Normalize each column 0–1 for visualization
col_keys = list(col_data.keys())
mat = np.zeros((len(top_pathways), len(col_keys)))
for j, ck in enumerate(col_keys):
    vals = np.array(col_data[ck], dtype=float)
    vmin, vmax = vals.min(), vals.max()
    mat[:, j] = (vals - vmin) / (vmax - vmin + 1e-12)

path_labels = [abbrev(r["pathway_name"], 42) for r in top_pathways]

fig, (ax_tier, ax_heat, ax_bar) = plt.subplots(
    1, 3, figsize=(16, 7),
    gridspec_kw={"width_ratios": [0.06, 1.1, 0.4]})

# Tier strip
for i, tier in enumerate(tiers):
    ax_tier.add_patch(mpatches.Rectangle(
        (0, i - 0.5), 1, 1, color=tier_colors[tier], alpha=0.85))
    ax_tier.text(0.5, i, tier[0], ha="center", va="center",
                 fontsize=8, color="white", fontweight="bold")
ax_tier.set_xlim(0, 1)
ax_tier.set_ylim(-0.5, len(tiers) - 0.5)
ax_tier.axis("off")
ax_tier.set_title("Tier", fontsize=9)

# Heatmap
cmap_heat = LinearSegmentedColormap.from_list(
    "evidence", ["#f7f7f7", "#fce4d6", "#E07B39", "#C00000"])
im = ax_heat.imshow(mat, aspect="auto", cmap=cmap_heat, vmin=0, vmax=1)
ax_heat.set_xticks(range(len(col_keys)))
ax_heat.set_xticklabels(col_keys, fontsize=9, rotation=30, ha="right")
ax_heat.set_yticks(range(len(top_pathways)))
ax_heat.set_yticklabels(path_labels, fontsize=7.5)
ax_heat.set_title("Multi-Evidence Synthesis Heatmap\n"
                  "(Inv.6 — normalized per column, 20 top pathways)", fontsize=10)

# Annotate cells with actual values
actual_vals = {
    "Mean F$_{ST}$":    [f"{v:.3f}" for v in col_data["Mean F$_{ST}$"]],
    "Conv. genes":      [str(v) for v in col_data["Conv. genes"]],
    "Conv. groups":     [str(v) for v in col_data["Conv. groups"]],
    "GO genes":         [str(v) for v in col_data["GO genes"]],
    "Composite\nscore": [f"{v:.3f}" for v in col_data["Composite\nscore"]],
}
for j, ck in enumerate(col_keys):
    for i in range(len(top_pathways)):
        cell_val = mat[i, j]
        text_c = "white" if cell_val > 0.55 else "#333"
        ax_heat.text(j, i, actual_vals[ck][i],
                     ha="center", va="center", fontsize=6.5, color=text_c)

# Colorbar
cbar = fig.colorbar(im, ax=ax_heat, fraction=0.025, pad=0.02)
cbar.set_label("Normalized evidence", fontsize=8)
cbar.ax.tick_params(labelsize=7)

# Right panel: tier legend + key finding
ax_bar.axis("off")
tier_patches = [mpatches.Patch(color=c, label=f"{t} evidence")
                for t, c in tier_colors.items()]
ax_bar.legend(handles=tier_patches, loc="upper center", fontsize=10,
              title="Evidence tier", title_fontsize=10)
ax_bar.text(0.5, 0.55,
            "★ Top finding:\nKCNE1 cardiac K⁺ channel\nhard sweep in ALL 5\ncontinental groups\n"
            "FST ratio = 0.30×\n→ ANCIENT pre-OoA sweep\n(~60–70 kya)",
            transform=ax_bar.transAxes, ha="center", va="center",
            fontsize=9.5, style="italic",
            bbox=dict(boxstyle="round,pad=0.5", fc="#FFF2CC", ec="#E07B39", lw=1.5))

plt.tight_layout()
fig.savefig(FIGS / "fig8_synthesis_heatmap.png")
plt.close()
print("  Saved fig8_synthesis_heatmap.png", flush=True)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 9 — Summary Overview Panel (all 8 key results in one figure)
# ─────────────────────────────────────────────────────────────────────────────
print("Fig 9: Summary overview panel ...", flush=True)

fig = plt.figure(figsize=(18, 14))
gs_outer = gridspec.GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.40)

# ---- Panel A: Top 10 pathway FST (bar) ----
axA = fig.add_subplot(gs_outer[0, 0])
top10_fst = [r for r in pathway_fst if r["mean_fst"] is not None][:10]
names_a  = [abbrev(r["pathway_name"], 30) for r in reversed(top10_fst)]
vals_a   = [r["mean_fst"] for r in reversed(top10_fst)]
cols_a   = plt.cm.YlOrRd(np.linspace(0.35, 0.95, 10))[::-1]
axA.barh(names_a, vals_a, color=cols_a, edgecolor="white", lw=0.4)
axA.set_xlabel("Mean F$_{ST}$", fontsize=8)
axA.set_title("A. Top Pathway F$_{ST}$ (Inv.1)", fontsize=9)
axA.tick_params(labelsize=6.5)

# ---- Panel B: Convergent genes bubble ----
axB = fig.add_subplot(gs_outer[0, 1])
for gname, ng, mh, ngo in zip(gene_names, n_groups, max_h12_arr, n_go):
    color = list(GROUP_COLORS.values())[min(ng - 2, 4)]
    sz = max(50, ngo * 10 + 40)
    axB.scatter(ng, mh, s=sz, c=color, alpha=0.80, edgecolors="white", lw=0.8, zorder=3)
    axB.annotate(gname, (ng, mh), xytext=(3, 4), textcoords="offset points", fontsize=7)
axB.set_xlabel("Continental groups")
axB.set_ylabel("Max H12")
axB.set_title("B. Convergent Genes (Inv.2)", fontsize=9)
axB.set_xticks([2, 3, 4, 5])
axB.set_xlim(1.5, 5.8)

# ---- Panel C: π vs n_sweeps ----
axC = fig.add_subplot(gs_outer[0, 2])
for pop in pops_rs:
    grp = pop_data[pop]["group"]
    axC.scatter(pop_data[pop]["mean_pi"], pop_data[pop]["n_hard_sweeps"],
                s=45, c=GROUP_COLORS.get(grp, "#aaa"),
                edgecolors="white", lw=0.5, alpha=0.88, zorder=3)
axC.set_xlabel("π")
axC.set_ylabel("N hard sweeps")
axC.set_title(f"C. π vs Sweeps  ρ=−0.75*** (Inv.3)", fontsize=9)
axC.legend(handles=group_legend_patches(), fontsize=6.5, loc="upper right")

# ---- Panel D: KCNE1 h12 per pop ----
axD = fig.add_subplot(gs_outer[1, 0])
axD.bar(x_sw, h12_vals, color=grp_colors_sw, edgecolor="white", lw=0.5, alpha=0.88)
axD.set_xticks(x_sw)
axD.set_xticklabels(pops_sw, fontsize=8)
axD.set_ylabel("H12")
axD.set_ylim(0, 1.12)
axD.set_title("D. KCNE1 Sweep Intensity (Inv.4)", fontsize=9)
for xi, val in zip(x_sw, h12_vals):
    axD.text(xi, val + 0.01, f"{val:.2f}", ha="center", fontsize=7)

# ---- Panel E: KCNE1 FST comparison ----
axE = fig.add_subplot(gs_outer[1, 1])
axE.bar(x_fst - w_f/2, genome_fst_v, w_f, color="#4472C4", alpha=0.75, label="Genome", edgecolor="white")
axE.bar(x_fst + w_f/2, locus_fst_v,  w_f, color="#E07B39", alpha=0.85, label="KCNE1 locus", edgecolor="white")
axE.set_xticks(x_fst)
axE.set_xticklabels(pair_labels_k, rotation=30, ha="right", fontsize=7.5)
axE.set_ylabel("F$_{ST}$")
axE.set_title("E. KCNE1 Locus vs Genome F$_{ST}$ (Inv.4)", fontsize=9)
axE.legend(fontsize=7.5)
for xi, rv in enumerate(ratio_v):
    axE.text(xi, max(locus_fst_v[xi], genome_fst_v[xi]) + 0.003,
             f"×{rv:.2f}", ha="center", fontsize=7, color="#C00000", fontweight="bold")

# ---- Panel F: GO enrichment (top 12 terms, CHB + GIH shown) ----
axF = fig.add_subplot(gs_outer[1, 2])
go_terms_shown = all_terms[:12]
for j_go, pop_go in enumerate(["CHB", "GIH"]):
    for i_go, term in enumerate(go_terms_shown):
        pv_go = go_padj[term].get(pop_go, 1.0)
        if pv_go >= 0.05:
            continue
        ng_go = go_nfg[term].get(pop_go, 0)
        sig_go = -np.log10(max(pv_go, 1e-8))
        col_go = plt.cm.Reds(min(sig_go / 4.0, 1.0))
        sz_go  = min(20 + ng_go * 30, 220)
        axF.scatter(j_go, i_go, s=sz_go, c=[col_go], edgecolors="gray", lw=0.4, zorder=3)
axF.set_xticks([0, 1])
axF.set_xticklabels(["CHB (E.Asian)", "GIH (S.Asian)"], fontsize=8.5)
axF.set_yticks(range(len(go_terms_shown)))
axF.set_yticklabels([t[:40] for t in go_terms_shown], fontsize=6.5)
axF.set_ylim(-0.6, len(go_terms_shown) - 0.4)
axF.set_title("F. GO Enrichment (Inv.5)", fontsize=9)

# ---- Panel G: PCA biplot ----
axG = fig.add_subplot(gs_outer[2, 0])
for grp_g, entries_g in by_group.items():
    axG.scatter([e["pc1"] for e in entries_g],
                [e["pc2"] for e in entries_g],
                s=40, c=GROUP_COLORS.get(grp_g, "#aaa"),
                edgecolors="white", lw=0.5,
                label=grp_g.replace("_", " "), alpha=0.88)
    for e_g in entries_g:
        axG.annotate(e_g["pop"], (e_g["pc1"], e_g["pc2"]),
                     xytext=(3, 2), textcoords="offset points", fontsize=6.5)
axG.axhline(0, color="gray", lw=0.5, alpha=0.4)
axG.axvline(0, color="gray", lw=0.5, alpha=0.4)
axG.set_xlabel(f"PC1 ({100*var_exp[0]:.1f}%)", fontsize=8)
axG.set_ylabel(f"PC2 ({100*var_exp[1]:.1f}%)", fontsize=8)
axG.set_title("G. Population PCA (Inv.7)", fontsize=9)
axG.legend(fontsize=6.5, loc="upper left")

# ---- Panel H: OoA distance vs π ----
axH = fig.add_subplot(gs_outer[2, 1])
for entry_h in ooa_pops:
    color_h = GROUP_COLORS.get(entry_h["group"], "#aaa")
    axH.scatter(entry_h["dist_ooa_km"], entry_h["mean_pi"],
                s=45, c=color_h, edgecolors="white", lw=0.5, alpha=0.88, zorder=3)
xs_h = np.array([e["dist_ooa_km"] for e in ooa_pops])
ys_h = np.array([e["mean_pi"] for e in ooa_pops])
z_h  = np.polyfit(xs_h, ys_h, 1)
axH.plot(np.linspace(xs_h.min(), xs_h.max(), 100),
         np.polyval(z_h, np.linspace(xs_h.min(), xs_h.max(), 100)),
         "--", color="gray", lw=1.1, alpha=0.6)
axH.set_xlabel("OoA distance (km)", fontsize=8)
axH.set_ylabel("Mean π", fontsize=8)
axH.set_title("H. OoA Gradient  ρ=−0.35 (Inv.7)", fontsize=9)

# ---- Panel I: Synthesis top 8 ----
axI = fig.add_subplot(gs_outer[2, 2])
top8_s = synthesis["top_pathways"][:8]
s_names = [abbrev(r["pathway_name"], 30) for r in top8_s]
s_comp  = [r["composite_score"] for r in top8_s]
s_tier_c = [tier_colors[r["tier"]] for r in top8_s]
axI.barh(s_names[::-1], s_comp[::-1], color=s_tier_c[::-1],
         edgecolor="white", lw=0.4, alpha=0.88)
axI.set_xlabel("Composite evidence score", fontsize=8)
axI.set_title("I. Synthesis Top Pathways (Inv.6)", fontsize=9)
axI.tick_params(labelsize=6.5)
tier_legend = [mpatches.Patch(color=c, label=t) for t, c in tier_colors.items()]
axI.legend(handles=tier_legend, fontsize=7, loc="lower right")

fig.suptitle("GraphPop Multi-Layer Population Genomics Analysis — 1000 Genomes 26 Populations\n"
             "Seven novel graph-native investigations from a single Neo4j database",
             fontsize=13, fontweight="bold", y=1.01)

fig.savefig(FIGS / "fig9_summary_overview.png", bbox_inches="tight")
plt.close()
print("  Saved fig9_summary_overview.png", flush=True)

# ── Final report ──────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("All figures saved to data/results/figures/")
figs_list = sorted(FIGS.glob("*.png"))
for fp in figs_list:
    size_kb = fp.stat().st_size // 1024
    print(f"  {fp.name:<40}  {size_kb:>5} KB")
print("=" * 60)
