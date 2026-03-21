#!/usr/bin/env python3
"""Generate supplementary figures for Investigations 8, 9, and 10.

Produces:
  fig10_gene_fst.png           - Gene-level FST heatmap (top 30 genes)
  fig11_pbs_scan.png           - Population-specific PBS sweep genes
  fig12_temporal_selection.png - Temporal stratification fingerprint (if data available)

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/generate_supp_figures.py
"""

import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
RES  = Path("data/results")
FIGS = RES / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

# ── Global style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        9,
    "axes.labelsize":   10,
    "axes.titlesize":   11,
    "axes.titleweight": "bold",
    "legend.fontsize":  8,
    "xtick.labelsize":  8,
    "ytick.labelsize":  8,
    "figure.dpi":       150,
    "savefig.dpi":      200,
    "savefig.bbox":     "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

GROUP_COLORS = {
    "African":    "#E07B39",
    "European":   "#4472C4",
    "East_Asian": "#70AD47",
    "South_Asian":"#7030A0",
    "American":   "#C00000",
}
POP_TO_GROUP = {
    "ACB":"African","ASW":"African","ESN":"African","GWD":"African",
    "LWK":"African","MSL":"African","YRI":"African",
    "CEU":"European","FIN":"European","GBR":"European","IBS":"European","TSI":"European",
    "CDX":"East_Asian","CHB":"East_Asian","CHS":"East_Asian","JPT":"East_Asian","KHV":"East_Asian",
    "BEB":"South_Asian","GIH":"South_Asian","ITU":"South_Asian","PJL":"South_Asian","STU":"South_Asian",
    "CLM":"American","MXL":"American","PEL":"American","PUR":"American",
}


# ── Figure 10: Gene-Level FST Heatmap ─────────────────────────────────────────
def fig10_gene_fst():
    p = RES / "gene_fst.json"
    if not p.exists():
        print("  [skip] gene_fst.json not found")
        return

    with open(p) as f:
        d = json.load(f)

    genes = d["top_genes"][:30]
    pairs = d["focal_pairs"]  # e.g. ["YRI_vs_CEU", ...]
    pair_labels = [p.replace("_vs_", "/") for p in pairs]

    gene_names = [g["gene_id"] for g in genes]
    mean_fst   = [g["mean_fst"] for g in genes]
    in_sweep   = [g.get("in_convergent_sweep", False) for g in genes]

    # Build FST matrix (rows = genes, cols = pairs)
    mat = np.zeros((len(genes), len(pairs)))
    for i, g in enumerate(genes):
        for j, tag in enumerate(pairs):
            val = g.get(f"fst_{tag}")
            mat[i, j] = val if val is not None else 0.0

    fig, (ax_heat, ax_bar) = plt.subplots(
        1, 2, figsize=(13, 9),
        gridspec_kw={"width_ratios": [3, 1], "wspace": 0.35}
    )
    fig.suptitle("Investigation 8: Gene-Level Population Differentiation (Hudson FST)\n"
                 "Top 30 genes by mean FST across 5 population pairs", fontsize=12, y=1.01)

    # Heatmap
    im = ax_heat.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.9)
    ax_heat.set_xticks(range(len(pairs)))
    ax_heat.set_xticklabels(pair_labels, rotation=30, ha="right", fontsize=8)
    ax_heat.set_yticks(range(len(genes)))

    yticklabels = []
    for name, sw in zip(gene_names, in_sweep):
        label = f"★ {name}" if sw else name
        yticklabels.append(label)
    ax_heat.set_yticklabels(yticklabels, fontsize=7.5)
    ax_heat.set_title("Per-pair Hudson FST")

    # Add value annotations
    for i in range(len(genes)):
        for j in range(len(pairs)):
            v = mat[i, j]
            if v > 0.01:
                ax_heat.text(j, i, f"{v:.2f}", ha="center", va="center",
                             fontsize=6, color="white" if v > 0.5 else "black")

    plt.colorbar(im, ax=ax_heat, label="Hudson FST", shrink=0.6, pad=0.02)

    # Highlight known genes
    known = {"SLC24A5": "#4472C4", "MKRN2": "#70AD47", "ID1": "#7030A0",
             "RNF135": "#E07B39", "INTS11": "#C00000", "ETFDH": "#888888"}
    for i, name in enumerate(gene_names):
        if name in known:
            ax_heat.axhline(i - 0.5, color="gray", linewidth=0.4, linestyle="--", alpha=0.5)
            ax_heat.axhline(i + 0.5, color="gray", linewidth=0.4, linestyle="--", alpha=0.5)

    # Bar chart — mean FST with highlighted genes
    colors = ["#E07B39" if name == "SLC24A5" else
              "#4472C4" if name in ("BTF3", "FMC1") else
              "#888888" for name in gene_names]
    ax_bar.barh(range(len(genes)), mean_fst, color=colors, edgecolor="white", linewidth=0.3)
    ax_bar.set_yticks(range(len(genes)))
    ax_bar.set_yticklabels(gene_names, fontsize=7.5)
    ax_bar.set_xlabel("Mean Hudson FST")
    ax_bar.set_xlim(0, 0.85)
    ax_bar.set_title("Mean FST\n(all pairs)")
    ax_bar.invert_yaxis()

    # Mark SLC24A5 explicitly
    if "SLC24A5" in gene_names:
        idx = gene_names.index("SLC24A5")
        ax_bar.text(mean_fst[idx] + 0.01, idx, "SLC24A5\n(skin pigmentation)",
                    va="center", fontsize=7, color="#E07B39", fontweight="bold")

    patches = [
        mpatches.Patch(color="#E07B39", label="SLC24A5 (skin pigmentation)"),
        mpatches.Patch(color="#4472C4", label="High YRI/CEU differentiation"),
        mpatches.Patch(color="#888888", label="Other top genes"),
    ]
    fig.legend(handles=patches, loc="lower center", ncol=3, fontsize=8,
               bbox_to_anchor=(0.5, -0.03))

    out = FIGS / "fig10_gene_fst.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Figure 11: PBS Population-Specific Sweep Genes ────────────────────────────
def fig11_pbs_scan():
    p = RES / "pbs_scan.json"
    if not p.exists():
        print("  [skip] pbs_scan.json not found")
        return

    with open(p) as f:
        d = json.load(f)

    pops  = list(d["populations"].keys())
    pdata = d["populations"]

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Investigation 9: Population-Specific Selection (PBS Scan)\n"
                 "Top 15 genes above 99th-percentile PBS per population", fontsize=12)

    for ax, pop in zip(axes.flat, pops):
        r    = pdata[pop]
        top  = r["top_genes"][:15]
        gids = [g["gene_id"] for g in top]
        pbs  = [g["max_pbs"] for g in top]
        n_w  = [g["n_windows"] for g in top]

        grp   = r["group"]
        color = GROUP_COLORS.get(grp, "#888888")

        # Color by n_windows (depth of signal)
        cmap  = plt.cm.get_cmap("Blues")
        norm  = plt.Normalize(1, max(n_w) + 1)
        bar_colors = [cmap(norm(n)) for n in n_w]

        bars = ax.barh(range(len(gids)), pbs, color=bar_colors, edgecolor="white", lw=0.3)
        ax.set_yticks(range(len(gids)))
        ax.set_yticklabels(gids, fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("Max PBS score")
        ax.set_title(f"{pop} ({grp})\nthreshold={r['pbs_threshold']:.4f} | "
                     f"{r['n_selected_genes']} genes selected",
                     color=color, fontweight="bold")

        # Annotate known biology
        bio_labels = {
            "SLC24A5": "skin pigmentation",
            "HLA-DPA1": "MHC/immune",
            "HLA-DPB1": "MHC/immune",
            "SH2B1": "insulin/obesity",
            "ATXN2L": "neurodegeneration",
            "BCAR3": "breast cancer",
            "BICD2": "dynein",
        }
        for i, gid in enumerate(gids):
            if gid in bio_labels:
                ax.text(pbs[i] * 1.02, i, bio_labels[gid],
                        va="center", fontsize=6.5, color="#C00000", style="italic")

    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    out = FIGS / "fig11_pbs_scan.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Figure 12: Temporal Stratification Fingerprint ────────────────────────────
def fig12_temporal_selection():
    p = RES / "temporal_selection.json"
    if not p.exists():
        print("  [skip] temporal_selection.json not found (Inv.10 still running)")
        return

    with open(p) as f:
        d = json.load(f)

    rows = d["populations"]
    pops_with_all = [r for r in rows if r.get("froh_z") is not None
                     and r.get("ihs_z") is not None and r.get("h12_z") is not None]

    if not pops_with_all:
        print("  [skip] No populations with all 3 temporal layers")
        return

    fig = plt.figure(figsize=(15, 10))
    fig.suptitle("Investigation 10: Temporal Stratification of Selection Signals\n"
                 "Three complementary timescales: FROH (very recent), iHS (recent), H12 (intermediate)",
                 fontsize=12)

    gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.35)
    ax_scatter = fig.add_subplot(gs[0, :])   # full-width top
    ax_bar     = fig.add_subplot(gs[1, 0])
    ax_ts      = fig.add_subplot(gs[1, 1])

    # ── Panel A: FROH_z vs iHS_z scatter (size = H12_z) ──
    for r in rows:
        fz  = r.get("froh_z")
        iz  = r.get("ihs_z")
        hz  = r.get("h12_z")
        if fz is None or iz is None:
            continue
        grp   = r["group"]
        color = GROUP_COLORS.get(grp, "#888888")
        sz    = max(20, 80 + (hz or 0) * 30)
        ax_scatter.scatter(fz, iz, c=color, s=sz, alpha=0.8, zorder=3,
                           edgecolors="white", linewidth=0.5)
        ax_scatter.annotate(r["pop"], (fz, iz), fontsize=6.5,
                            xytext=(3, 2), textcoords="offset points", color=color)

    ax_scatter.axhline(0, color="gray", lw=0.8, ls="--", alpha=0.5)
    ax_scatter.axvline(0, color="gray", lw=0.8, ls="--", alpha=0.5)
    ax_scatter.set_xlabel("FROH z-score (very recent <5 kya)")
    ax_scatter.set_ylabel("iHS fraction z-score (recent ~5-30 kya)")
    ax_scatter.set_title("Temporal Selection Fingerprint\n(bubble size = H12 z-score; intermediate ~10-100 kya)",
                         loc="left")

    # Group legend
    for grp, col in GROUP_COLORS.items():
        ax_scatter.scatter([], [], c=col, s=50, label=grp)
    ax_scatter.legend(loc="upper right", frameon=False, fontsize=8)

    # ── Panel B: Dominant timescale bar chart ──
    from collections import Counter
    by_grp = {}
    ts_set = set()
    for r in rows:
        g  = r["group"]
        ts = r.get("dominant_timescale", "unknown")
        ts_set.add(ts)
        if g not in by_grp:
            by_grp[g] = []
        by_grp[g].append(ts)

    ts_list = sorted(ts_set)
    ts_colors = {
        "very_recent(FROH)":    "#C00000",
        "recent(iHS)":          "#4472C4",
        "intermediate(H12)":    "#70AD47",
        "unknown":              "#AAAAAA",
    }

    grps  = sorted(by_grp.keys())
    x     = np.arange(len(grps))
    bot   = np.zeros(len(grps))
    for ts in ts_list:
        counts = [sum(1 for t in by_grp[g] if t == ts) for g in grps]
        ax_bar.bar(x, counts, bottom=bot, color=ts_colors.get(ts, "#AAAAAA"),
                   label=ts, edgecolor="white", linewidth=0.5)
        bot += np.array(counts, dtype=float)

    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels([g.replace("_", "\n") for g in grps], fontsize=8)
    ax_bar.set_ylabel("Number of populations")
    ax_bar.set_title("Dominant Timescale\nby Continental Group")
    ax_bar.legend(loc="upper right", fontsize=7, frameon=False)

    # ── Panel C: Heatmap of z-scores across all populations with iHS ──
    pop_labels = [r["pop"] for r in pops_with_all]
    mat = np.array([[r["froh_z"], r["ihs_z"], r["h12_z"]] for r in pops_with_all])
    col_grp = [GROUP_COLORS.get(r["group"], "#888888") for r in pops_with_all]

    im = ax_ts.imshow(mat.T, aspect="auto", cmap="RdBu_r",
                      vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax_ts.set_yticks([0, 1, 2])
    ax_ts.set_yticklabels(["FROH_z\n(very recent)", "iHS_z\n(recent)", "H12_z\n(intermediate)"], fontsize=8)
    ax_ts.set_xticks(range(len(pop_labels)))
    ax_ts.set_xticklabels(pop_labels, rotation=45, ha="right", fontsize=7)
    ax_ts.set_title("Temporal Layer Z-scores\nper Population")
    plt.colorbar(im, ax=ax_ts, label="Z-score", shrink=0.8)

    # Color-code x-tick labels by group
    for tick, grp_color in zip(ax_ts.get_xticklabels(), col_grp):
        tick.set_color(grp_color)

    out = FIGS / "fig12_temporal_selection.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Overview panel for new investigations ─────────────────────────────────────
def fig13_new_overview():
    """4-panel summary for Inv.8–10 combined."""
    gene_p = RES / "gene_fst.json"
    pbs_p  = RES / "pbs_scan.json"
    if not gene_p.exists() or not pbs_p.exists():
        print("  [skip] gene_fst or pbs_scan JSON not found")
        return

    with open(gene_p) as f:
        gd = json.load(f)
    with open(pbs_p) as f:
        pd = json.load(f)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.suptitle("Summary: Investigations 8–10 — Gene-Level Selection Landscape",
                 fontsize=13, fontweight="bold")

    ax1, ax2, ax3, ax4 = axes.flat

    # ── A: Top genes by mean FST (bar chart) ──
    top20 = gd["top_genes"][:20]
    gnames = [g["gene_id"] for g in top20]
    gmean  = [g["mean_fst"]  for g in top20]
    gmax   = [g["max_fst"]   for g in top20]
    highlight = {"SLC24A5": "#E07B39", "MKRN2": "#70AD47", "ID1": "#7030A0",
                 "BTF3": "#4472C4", "RNF135": "#C00000"}
    bar_col = [highlight.get(g, "#AAAAAA") for g in gnames]

    ax1.barh(range(len(gnames)), gmean, color=bar_col, alpha=0.85, label="Mean FST")
    ax1.barh(range(len(gnames)), gmax, color=bar_col, alpha=0.3, label="Max FST")
    ax1.set_yticks(range(len(gnames)))
    ax1.set_yticklabels(gnames, fontsize=7)
    ax1.invert_yaxis()
    ax1.set_xlabel("Hudson FST")
    ax1.set_title("A. Top 20 Genes by Mean FST\n(opaque=mean, transparent=max)")
    ax1.legend(fontsize=7, frameon=False)

    # Highlight SLC24A5
    if "SLC24A5" in gnames:
        i = gnames.index("SLC24A5")
        ax1.axhline(i, color="#E07B39", lw=1.5, ls="--", alpha=0.6)

    # ── B: N_selected_genes per pop for PBS ──
    pop_order = sorted(pd["populations"].keys(),
                       key=lambda p: -pd["populations"][p]["n_selected_genes"])
    n_genes  = [pd["populations"][p]["n_selected_genes"] for p in pop_order]
    grp_cols = [GROUP_COLORS.get(pd["populations"][p]["group"], "#888888") for p in pop_order]

    ax2.bar(range(len(pop_order)), n_genes, color=grp_cols, edgecolor="white", lw=0.5)
    ax2.set_xticks(range(len(pop_order)))
    ax2.set_xticklabels(pop_order, rotation=30, ha="right", fontsize=9)
    ax2.set_ylabel("Genes above P99 PBS")
    ax2.set_title("B. Population-Specific Selection (PBS)\nGenes above 99th percentile")

    # Group legend
    for grp, col in GROUP_COLORS.items():
        ax2.bar([], [], color=col, label=grp)
    ax2.legend(fontsize=7, frameon=False)

    # ── C: SLC24A5 FST across pairs (bar) ──
    slc_data = next((g for g in gd["top_genes"] if g["gene_id"] == "SLC24A5"), None)
    if slc_data:
        pairs = gd["focal_pairs"]
        pair_labels = [p.replace("_vs_", "/") for p in pairs]
        fst_vals = [slc_data.get(f"fst_{p}") or 0.0 for p in pairs]
        pair_colors = ["#E07B39" if "YRI" in p else "#4472C4" if "CEU" in p
                       else "#70AD47" for p in pairs]
        ax3.bar(range(len(pairs)), fst_vals, color=pair_colors, edgecolor="white", lw=0.5)
        ax3.set_xticks(range(len(pairs)))
        ax3.set_xticklabels(pair_labels, rotation=30, ha="right")
        ax3.set_ylabel("Hudson FST")
        ax3.set_ylim(0, 1.0)
        ax3.axhline(0.15, color="gray", ls="--", lw=0.8, label="Background FST")
        ax3.set_title("C. SLC24A5 (Skin Pigmentation)\nFST by Population Pair")
        ax3.legend(fontsize=7, frameon=False)
        for i, v in enumerate(fst_vals):
            if v > 0.01:
                ax3.text(i, v + 0.02, f"{v:.3f}", ha="center", fontsize=8, fontweight="bold")
    else:
        ax3.text(0.5, 0.5, "SLC24A5 not in top 100\ngene FST results",
                 ha="center", va="center", transform=ax3.transAxes, fontsize=10)
        ax3.set_visible(False)

    # ── D: PBS top-gene biology summary ──
    bio_summary = {
        "CEU": [("TMCC1", "cytoskeleton"), ("BCAR3", "breast cancer"), ("DNAH7", "dynein")],
        "CHB": [("IPPK", "inositol kinase"), ("BICD2", "dynein/Golgi"), ("CENPP", "centromere")],
        "FIN": [("PSD3", "PH domain"), ("DNAH7", "dynein"), ("SMAD3", "TGF-β")],
        "GIH": [("SLC24A5", "skin pigment"), ("MYEF2", "myelin"), ("CTXN2", "cortexin")],
        "PEL": [("SH2B1", "insulin/leptin"), ("ATXN2L", "spinocerebellar"), ("SLC7A5", "amino acid transport")],
        "YRI": [("HLA-DPA1", "MHC class II"), ("HLA-DPB1", "MHC class II"), ("HLA-DPB2", "MHC class II")],
    }
    ax4.axis("off")
    y = 0.97
    ax4.text(0.0, y, "D. PBS Top Genes — Biological Annotation", fontsize=10, fontweight="bold",
             transform=ax4.transAxes, va="top")
    y -= 0.07
    for pop, entries in bio_summary.items():
        grp = pd["populations"].get(pop, {}).get("group", "")
        color = GROUP_COLORS.get(grp, "#888888")
        ax4.text(0.0, y, f"{pop} ({grp}):", fontsize=9, fontweight="bold",
                 color=color, transform=ax4.transAxes, va="top")
        for gene, func in entries:
            y -= 0.055
            ax4.text(0.05, y, f"• {gene}: {func}", fontsize=8,
                     transform=ax4.transAxes, va="top", color="#333333")
        y -= 0.065

    out = FIGS / "fig13_new_overview.png"
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Figure 14: Consequence Bias — FST by VEP Impact Level ─────────────────────
def fig14_consequence_bias():
    p = RES / "consequence_bias.json"
    if not p.exists():
        print("  [skip] consequence_bias.json not found")
        return

    with open(p) as f:
        d = json.load(f)

    impact_stats = d["impact_stats"]
    mw = d.get("mann_whitney", {})
    focal_pairs = d.get("focal_pairs", ["YRI/CEU", "YRI/CHB", "CEU/CHB"])
    ihs_pops    = d.get("ihs_pops", ["CEU", "CHB"])

    IMPACTS = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
    IMPACT_COLORS = {"HIGH": "#C00000", "MODERATE": "#E07B39",
                     "LOW": "#4472C4", "MODIFIER": "#AAAAAA"}

    # Build dictionaries for quick lookup
    stats = {r["impact"]: r for r in impact_stats}

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Investigation 11: Functional Consequence Bias in Selection\n"
                 "Population differentiation (Hudson FST) and positive selection "
                 "signal (|iHS| > 2) by VEP impact category",
                 fontsize=12)

    # Panel A: FST per impact per population pair
    ax = axes[0]
    pair_labels = focal_pairs
    x = np.arange(len(pair_labels))
    w = 0.2
    offsets = np.linspace(-0.3, 0.3, len(IMPACTS))
    for i, impact in enumerate(IMPACTS):
        if impact not in stats:
            continue
        s = stats[impact]
        fst_vals = []
        for pair in focal_pairs:
            p1, p2 = pair.split("/")
            key = f"fst_{p1}_{p2}"
            fst_vals.append(s.get(key) or 0.0)
        ax.bar(x + offsets[i], fst_vals, w, color=IMPACT_COLORS[impact],
               label=impact, edgecolor="white", lw=0.4)

    ax.set_xticks(x)
    ax.set_xticklabels(pair_labels)
    ax.set_ylabel("Hudson FST (ratio-of-averages)")
    ax.set_title("A. Mean FST by Impact Category\nand Population Pair")
    ax.legend(fontsize=7, frameon=False)

    # Panel B: iHS enrichment (|iHS| > 2) per impact × population
    ax = axes[1]
    ihs_pops_plot = ihs_pops[:3]  # show top 3
    x2 = np.arange(len(ihs_pops_plot))
    offsets2 = np.linspace(-0.3, 0.3, len(IMPACTS))
    for i, impact in enumerate(IMPACTS):
        if impact not in stats:
            continue
        s = stats[impact]
        enrich_vals = [s.get(f"ihs_enrich_{pop}") or 0.0 for pop in ihs_pops_plot]
        ax.bar(x2 + offsets2[i], enrich_vals, 0.2,
               color=IMPACT_COLORS[impact], label=impact, edgecolor="white", lw=0.4)

    ax.set_xticks(x2)
    ax.set_xticklabels(ihs_pops_plot)
    ax.set_ylabel("Fraction of variants with |iHS| > 2")
    ax.set_title("B. Positive Selection Enrichment\n(|iHS| > 2) by Impact and Population")
    ax.legend(fontsize=7, frameon=False)

    # Panel C: N variants per category + interpretation text
    ax = axes[2]
    n_vars = [stats.get(imp, {}).get("n_variants", 0) for imp in IMPACTS]
    colors = [IMPACT_COLORS[imp] for imp in IMPACTS]
    bars = ax.bar(range(len(IMPACTS)), n_vars, color=colors, edgecolor="white", lw=0.4)
    ax.set_xticks(range(len(IMPACTS)))
    ax.set_xticklabels(IMPACTS)
    ax.set_ylabel("Number of variant-consequence pairs")
    ax.set_title("C. Variant Count per Impact Category")
    for bar, n in zip(bars, n_vars):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 100,
                f"{n:,}", ha="center", va="bottom", fontsize=7)

    # Add Mann-Whitney annotation
    mw_yri_ceu = mw.get("YRI/CEU", {}).get("HIGH_vs_LOW", {})
    if mw_yri_ceu:
        pval = mw_yri_ceu.get("pval", "")
        h_mean = mw_yri_ceu.get("mean1", "")
        l_mean = mw_yri_ceu.get("mean2", "")
        ax.text(0.05, 0.95,
                f"Mann-Whitney HIGH vs LOW\n(YRI/CEU FST)\n"
                f"HIGH mean={h_mean:.4f}\nLOW mean={l_mean:.4f}\np={pval}",
                transform=ax.transAxes, fontsize=7, va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8))

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGS / "fig14_consequence_bias.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Figure 15: Pathway Co-selection Network ────────────────────────────────────
def fig15_pathway_coselection():
    p = RES / "pathway_coselection.json"
    if not p.exists():
        print("  [skip] pathway_coselection.json not found")
        return

    with open(p) as f:
        d = json.load(f)

    import networkx as nx

    top_communities = d.get("top_communities", [])[:3]
    top_edges       = d.get("top_edges", [])[:200]

    fig, axes = plt.subplots(1, 2, figsize=(15, 7))
    fig.suptitle("Investigation 12: Pathway Co-selection Network\n"
                 "Pathways with correlated selection signals (Spearman |ρ| > 0.7)",
                 fontsize=12)

    # Panel A: Spring-layout network of top communities (subsample)
    ax = axes[0]
    ax.set_title(f"A. Co-selection Network\n"
                 f"({d['n_edges']:,} edges, {d['n_communities']} communities detected)")

    # Build a small subgraph from top edges for visualization
    G_sub = nx.Graph()
    for edge in top_edges[:100]:
        p1 = edge["pathway1_id"]
        p2 = edge["pathway2_id"]
        # Truncate names
        n1 = edge["pathway1_name"][:30]
        n2 = edge["pathway2_name"][:30]
        G_sub.add_node(p1, name=n1)
        G_sub.add_node(p2, name=n2)
        G_sub.add_edge(p1, p2, rho=abs(edge["rho"]))

    if len(G_sub.nodes) > 0:
        pos = nx.spring_layout(G_sub, seed=42, k=1.5)
        # Color by community (assign by degree rank as proxy)
        degrees = dict(G_sub.degree())
        max_deg = max(degrees.values()) if degrees else 1
        node_colors = [plt.cm.tab10(degrees[n] / max_deg) for n in G_sub.nodes()]
        edge_weights = [G_sub[u][v]["rho"] for u, v in G_sub.edges()]

        nx.draw_networkx_edges(G_sub, pos, ax=ax, alpha=0.3,
                               width=[w * 2 for w in edge_weights],
                               edge_color="#888888")
        nx.draw_networkx_nodes(G_sub, pos, ax=ax, node_color=node_colors,
                               node_size=60, alpha=0.8)

        # Label top-degree nodes only
        top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:8]
        labels = {n: G_sub.nodes[n]["name"] for n, _ in top_nodes}
        nx.draw_networkx_labels(G_sub, pos, labels=labels, ax=ax, font_size=5.5)

    ax.axis("off")

    # Panel B: Community bar chart
    ax2 = axes[1]
    comm_sizes  = [c["size"] for c in top_communities]
    comm_fsts   = [c["mean_fst"] for c in top_communities]
    comm_labels = [f"Community {c['community_id']}\n{c['top_pathway'][:25]}"
                   for c in top_communities]

    x = np.arange(len(comm_labels))
    ax2_twin = ax2.twinx()
    bars = ax2.bar(x, comm_sizes, color=["#4472C4", "#E07B39", "#70AD47"],
                   alpha=0.7, label="Community size", edgecolor="white")
    ax2_twin.plot(x, comm_fsts, "D--", color="#C00000", markersize=8,
                  label="Mean FST", zorder=3)
    ax2.set_xticks(x)
    ax2.set_xticklabels(comm_labels, fontsize=8)
    ax2.set_ylabel("Community size (# pathways)")
    ax2_twin.set_ylabel("Mean FST", color="#C00000")
    ax2.set_title("B. Top 3 Communities\nSize and Mean FST")

    for bar, n in zip(bars, comm_sizes):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
                 str(n), ha="center", va="bottom", fontsize=9, fontweight="bold")

    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, fontsize=8, frameon=False)

    # Network stats annotation
    stats_text = (f"Total: {d['n_pathways']:,} pathways\n"
                  f"Edges (|ρ|>0.7): {d['n_edges']:,}\n"
                  f"Communities: {d['n_communities']}\n"
                  f"Largest community: {d['network_stats']['largest_component_size']}")
    ax2.text(0.98, 0.02, stats_text, transform=ax2.transAxes,
             fontsize=8, va="bottom", ha="right",
             bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8))

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGS / "fig15_pathway_coselection.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Figure 16: Cross-Species Convergence ──────────────────────────────────────
def fig16_cross_species():
    p = RES / "cross_species.json"
    if not p.exists():
        print("  [skip] cross_species.json not found")
        return

    with open(p) as f:
        d = json.load(f)

    rows = d.get("module_convergence", [])
    if not rows:
        print("  [skip] No module convergence data")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Investigation 13: Cross-Species Parallel Adaptation\n"
                 "Human population genomics × Rice domestication: convergent functional modules",
                 fontsize=12)

    # Panel A: Convergence score bar chart
    ax = axes[0]
    modules    = [r["module"] for r in rows]
    h_signal   = [r["human_signal_score"] for r in rows]
    r_npops    = [r["rice_n_pops_total"] for r in rows]
    conv_score = [r["convergence_score"] for r in rows]

    MODULE_COLORS = {
        "development_growth":  "#4472C4",
        "metabolism":          "#70AD47",
        "stress_response":     "#E07B39",
        "immunity_defense":    "#C00000",
        "signal_transduction": "#7030A0",
        "ion_transport":       "#00B0F0",
    }
    colors = [MODULE_COLORS.get(m, "#888888") for m in modules]

    y = range(len(modules))
    ax.barh(y, conv_score, color=colors, edgecolor="white", lw=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels([m.replace("_", " ") for m in modules], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Convergence Score (human × rice)")
    ax.set_title("A. Module-Level Convergence Score\n(Human signal × Rice populations sweeping)")

    # Annotate rice genes
    for i, r in enumerate(rows):
        if r["rice_genes"]:
            rice_label = r["rice_genes"][:30] + ("..." if len(r["rice_genes"]) > 30 else "")
            ax.text(conv_score[i] * 0.02, i, rice_label,
                    va="center", fontsize=6.5, color="#555555")

    # Panel B: 2D scatter — human signal vs rice populations
    ax2 = axes[1]
    ax2.scatter(h_signal, r_npops, c=colors, s=150, zorder=3,
                edgecolors="white", linewidths=0.8)
    for i, m in enumerate(modules):
        ax2.annotate(m.replace("_", "\n"), (h_signal[i], r_npops[i]),
                     fontsize=7.5, xytext=(4, 4), textcoords="offset points",
                     color=MODULE_COLORS.get(m, "#888888"))

    ax2.set_xlabel("Human selection signal score\n(# genes × enrichment)")
    ax2.set_ylabel("Rice populations sweeping (total)")
    ax2.set_title("B. Human vs Rice Selection Intensity\nper Functional Module")

    # Add diagonal reference
    xlim = ax2.get_xlim()
    ylim = ax2.get_ylim()
    ax2.axhline(0, color="gray", lw=0.7, ls="--", alpha=0.4)
    ax2.axvline(0, color="gray", lw=0.7, ls="--", alpha=0.4)

    # Summary annotation
    n_human = d["n_human_sweep_genes"]
    n_rice  = d["n_rice_domestication_genes"]
    ax2.text(0.98, 0.02,
             f"Human sweep genes: {n_human}\n"
             f"Rice domestication genes: {n_rice}\n"
             f"Functional modules: {len(rows)}",
             transform=ax2.transAxes, fontsize=8, va="bottom", ha="right",
             bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", alpha=0.8))

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGS / "fig16_cross_species.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Figure 17: Rare Variant Burden Heatmap ────────────────────────────────────
def fig17_rare_burden():
    p = RES / "rare_burden.json"
    if not p.exists():
        print("  [skip] rare_burden.json not found")
        return

    with open(p) as f:
        d = json.load(f)

    top_genes = [g["gene_id"] for g in d.get("top_genes", [])[:20]
                 if g.get("total_private_rare", 0) > 0]
    if not top_genes:
        print("  [skip] No top genes with rare burden (Inv.14 may still be running)")
        return

    pop_summary = d.get("pop_summary", {})
    pops = sorted(pop_summary.keys())

    # Build matrix: gene × pop → n_private_rare
    # Load TSV for full data
    tsv_path = RES / "rare_burden.tsv"
    if not tsv_path.exists():
        print("  [skip] rare_burden.tsv not found")
        return

    import csv as csv_mod
    gene_pop_mat = {g: {p: 0 for p in pops} for g in top_genes}
    with open(tsv_path) as f:
        for row in csv_mod.DictReader(f, delimiter="\t"):
            gid = row["gene_id"]
            pop = row["population"]
            if gid in gene_pop_mat and pop in gene_pop_mat[gid]:
                gene_pop_mat[gid][pop] = int(row["n_private_rare"])

    mat = np.array([[gene_pop_mat[g][p] for p in pops] for g in top_genes])

    # Sort populations by FROH z-score
    froh_z = np.array([pop_summary.get(p, {}).get("froh_z") or 0.0 for p in pops])
    pop_order = np.argsort(froh_z)[::-1]
    pops_sorted  = [pops[i] for i in pop_order]
    mat_sorted   = mat[:, pop_order]
    froh_z_sorted = froh_z[pop_order]

    fig, (ax_heat, ax_froh) = plt.subplots(
        2, 1, figsize=(16, 10),
        gridspec_kw={"height_ratios": [8, 1], "hspace": 0.05}
    )
    fig.suptitle("Investigation 14: Population-Specific Rare Variant Burden\n"
                 "Top 20 genes × 26 populations (private AF<0.5% variants)",
                 fontsize=12)

    # Heatmap
    im = ax_heat.imshow(mat_sorted, aspect="auto", cmap="YlOrRd",
                        vmin=0, vmax=max(1, mat_sorted.max()))
    ax_heat.set_xticks(range(len(pops_sorted)))
    ax_heat.set_xticklabels(pops_sorted, rotation=45, ha="right", fontsize=8)
    ax_heat.set_yticks(range(len(top_genes)))
    ax_heat.set_yticklabels(top_genes, fontsize=8)
    ax_heat.set_title("Private Rare Variant Count (per gene × population)")
    plt.colorbar(im, ax=ax_heat, label="n private rare variants", shrink=0.6,
                 orientation="horizontal", pad=0.01)

    # FROH z-score bar below
    froh_colors = ["#C00000" if z > 1 else "#4472C4" if z < -1 else "#AAAAAA"
                   for z in froh_z_sorted]
    ax_froh.bar(range(len(pops_sorted)), froh_z_sorted, color=froh_colors,
                edgecolor="white", lw=0.3)
    ax_froh.set_xticks([])
    ax_froh.set_ylabel("FROH z")
    ax_froh.axhline(0, color="gray", lw=0.7)
    ax_froh.set_xlim(-0.5, len(pops_sorted) - 0.5)

    # Correlation annotation
    froh_corr = d.get("froh_correlation", {})
    if froh_corr.get("spearman_r") is not None:
        ax_heat.text(
            0.99, 0.01,
            f"FROH × burden: ρ={froh_corr['spearman_r']:.3f}  "
            f"p={froh_corr['spearman_p']}",
            transform=ax_heat.transAxes, fontsize=8, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", alpha=0.9)
        )

    out = FIGS / "fig17_rare_burden.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Figure 18: Sweep Spatial Extent (FST Decay Curves) ────────────────────────
def fig18_sweep_extent():
    p = RES / "sweep_extent.json"
    if not p.exists():
        print("  [skip] sweep_extent.json not found")
        return

    with open(p) as f:
        d = json.load(f)

    genes_data = d.get("genes", [])
    decay_curves = d.get("decay_curves", {})

    if not genes_data:
        print("  [skip] No sweep extent data")
        return

    n_genes = len(genes_data)
    ncols = min(3, n_genes)
    nrows = (n_genes + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows))
    if n_genes == 1:
        axes = [axes]
    else:
        axes = axes.flat

    fig.suptitle("Investigation 15: Sweep Spatial Extent via FST Decay\n"
                 "Hudson FST in 10 kb bins radiating from convergent sweep gene centers",
                 fontsize=12)

    for ax, gene_info in zip(axes, genes_data):
        gid  = gene_info["gene_id"]
        bins = decay_curves.get(gid, [])

        if not bins:
            ax.text(0.5, 0.5, f"{gid}\nNo data", ha="center", va="center",
                    transform=ax.transAxes)
            ax.axis("off")
            continue

        dists = [b["dist"] / 1000 for b in bins]  # convert to kb
        fsts  = [b["fst"] if b["fst"] is not None else np.nan for b in bins]

        ax.plot(dists, fsts, "o-", markersize=2, lw=1.2, color="#4472C4", alpha=0.8)
        ax.fill_between(dists, 0, [f if not np.isnan(f) else 0 for f in fsts],
                        alpha=0.15, color="#4472C4")

        # Mark thresholds
        peak_fst = gene_info.get("peak_fst", 0)
        threshold = gene_info.get("fst_threshold", 0)
        if threshold:
            ax.axhline(threshold, color="#C00000", ls="--", lw=0.9,
                       label=f"Threshold={threshold:.3f}")

        # Mark extent boundaries
        left_kb  = -gene_info.get("extent_left_kb", 0)
        right_kb =  gene_info.get("extent_right_kb", 0)
        if left_kb:
            ax.axvline(left_kb, color="#E07B39", ls=":", lw=1.0, alpha=0.7)
        if right_kb:
            ax.axvline(right_kb, color="#E07B39", ls=":", lw=1.0, alpha=0.7)

        ax.set_xlabel("Distance from gene center (kb)")
        ax.set_ylabel("Hudson FST")
        sweep_pop  = gene_info.get("sweep_pop", "")
        ref_pop    = gene_info.get("ref_pop", "")
        extent_kb  = gene_info.get("total_extent_kb", 0)
        max_h12    = gene_info.get("max_h12", 0)
        ax.set_title(f"{gid}\n{ref_pop} vs {sweep_pop}  |  "
                     f"peak={peak_fst:.3f}  extent={extent_kb:.0f} kb  "
                     f"h12={max_h12:.2f}",
                     fontsize=8)
        ax.legend(fontsize=6.5, frameon=False)
        ax.set_xlim(-1000, 1000)

    # Hide unused subplots
    for ax in list(axes)[n_genes:]:
        ax.set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIGS / "fig18_sweep_extent.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}  ({out.stat().st_size//1024} KB)")


# ── Main ──────────────────────────────────────────────────────────────────────
def fig19_evo_trajectory():
    """fig19: Evolutionary trajectory map — re-plots from evo_trajectory.json/.tsv.

    The full multi-panel version is generated by investigate_16_evo_trajectory.py.
    This function produces a focused 2-panel figure (UMAP + velocity heatmap)
    suitable for inclusion in the supplementary PDF alongside fig14–fig18.
    """
    import csv

    tsv_path  = RES / "evo_trajectory.tsv"
    json_path = RES / "evo_trajectory.json"
    if not tsv_path.exists() or not json_path.exists():
        print("    Skipping: evo_trajectory.tsv/.json not found — run investigate_16 first")
        return

    with open(tsv_path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    d = json.load(open(json_path))

    GROUP_COLORS = {
        "African":     "#E05C5C",
        "European":    "#5C8AE0",
        "South_Asian": "#E0A85C",
        "East_Asian":  "#5CC45C",
        "American":    "#A85CE0",
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: UMAP + MST
    ax = axes[0]
    ax.set_title("A  Evolutionary Process Space (UMAP)", fontsize=11, pad=6, loc="left")

    # MST edges
    pop_xy = {r["population"]: (float(r["umap1"]), float(r["umap2"])) for r in rows}
    for edge in d.get("mst_edges", []):
        x1, y1 = pop_xy[edge["pop1"]]
        x2, y2 = pop_xy[edge["pop2"]]
        ax.plot([x1, x2], [y1, y2], color="#444", lw=1.2, alpha=0.7, zorder=1)

    for r in rows:
        x, y = float(r["umap1"]), float(r["umap2"])
        col = GROUP_COLORS.get(r["group"], "#888")
        ax.scatter(x, y, c=col, s=70, zorder=3, edgecolors="white", linewidths=0.4)
        ax.annotate(r["population"], (x, y), textcoords="offset points",
                    xytext=(3, 3), fontsize=6, zorder=4)

    for grp, col in GROUP_COLORS.items():
        ax.scatter([], [], c=col, s=35, label=grp.replace("_", " "))
    ax.legend(fontsize=7, loc="lower right")
    ax.set_xlabel("UMAP 1", fontsize=8)
    ax.set_ylabel("UMAP 2", fontsize=8)

    # Panel B: velocity heatmap
    ax = axes[1]
    ax.set_title("B  Evolutionary velocity (force decomposition)",
                 fontsize=11, pad=6, loc="left")

    vel_scores = d.get("velocity_scores", {})
    comp_labels = {"drift_bottleneck": "Drift/\nBottleneck",
                   "recent_selection": "Recent\nSelection",
                   "expansion":        "Expansion"}
    comp_keys = list(comp_labels.keys())

    drift_scores = vel_scores.get("drift_bottleneck", {})
    pop_order = sorted([r["population"] for r in rows],
                       key=lambda p: drift_scores.get(p, 0))
    pop_groups = {r["population"]: r["group"] for r in rows}

    vel_mat = np.array([[vel_scores.get(c, {}).get(p, 0) for c in comp_keys]
                         for p in pop_order])
    for col in range(vel_mat.shape[1]):
        mx = max(abs(vel_mat[:, col].max()), abs(vel_mat[:, col].min()))
        if mx > 0:
            vel_mat[:, col] /= mx

    im = ax.imshow(vel_mat, aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(comp_keys)))
    ax.set_xticklabels([comp_labels[c] for c in comp_keys], fontsize=8)
    ax.set_yticks(range(len(pop_order)))
    ax.set_yticklabels(pop_order, fontsize=6)
    for lbl, pop in zip(ax.get_yticklabels(), pop_order):
        lbl.set_color(GROUP_COLORS.get(pop_groups.get(pop, ""), "#888"))
    plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04, label="Norm. score")

    fig.suptitle("Evolutionary Trajectory Map — 26 Human Populations in Process Space",
                 fontsize=12, y=1.01)
    fig.tight_layout()

    out = FIGS / "fig19_evo_trajectory_summary.png"
    fig.savefig(out, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"    Saved: {out}")


def main():
    print("=== Generating Supplementary Figures (Inv.8–16) ===\n")
    print("[fig10] Gene-level FST heatmap ...")
    fig10_gene_fst()
    print("[fig11] PBS population-specific sweeps ...")
    fig11_pbs_scan()
    print("[fig12] Temporal selection fingerprint ...")
    fig12_temporal_selection()
    print("[fig13] Summary overview (Inv.8–10) ...")
    fig13_new_overview()
    print("[fig14] Consequence bias (Inv.11) ...")
    fig14_consequence_bias()
    print("[fig15] Pathway co-selection network (Inv.12) ...")
    fig15_pathway_coselection()
    print("[fig16] Cross-species convergence (Inv.13) ...")
    fig16_cross_species()
    print("[fig17] Rare variant burden heatmap (Inv.14) ...")
    fig17_rare_burden()
    print("[fig18] Sweep spatial extent FST decay (Inv.15) ...")
    fig18_sweep_extent()
    print("[fig19] Evolutionary trajectory map (Inv.16) ...")
    fig19_evo_trajectory()
    print("\nDone. Figures saved to data/results/figures/")


if __name__ == "__main__":
    main()
