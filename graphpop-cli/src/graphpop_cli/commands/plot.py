"""graphpop plot — generate standard population genomics figures from TSV results."""
from __future__ import annotations

import csv
import re
from pathlib import Path

import click

from ..cli import pass_ctx

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

try:
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


# ---------------------------------------------------------------------------
# Nature-style settings
# ---------------------------------------------------------------------------
WONG_PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#D55E00",
    "#56B4E9", "#CC79A7", "#F0E442", "#000000",
]


def _apply_style():
    """Apply Nature Methods figure style."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "lines.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
    })


def _check_matplotlib():
    if not HAS_MPL:
        click.echo(
            "Error: matplotlib is required for graphpop plot.\n"
            "Install with: pip install matplotlib numpy",
            err=True,
        )
        raise SystemExit(1)


def _read_tsv(path: str) -> list[dict]:
    """Read a TSV file, skipping comment lines."""
    rows = []
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    return list(reader)


def _read_tsv_dir(directory: str, pattern: str = "*.tsv") -> list[dict]:
    """Read all TSV files in a directory."""
    rows = []
    for p in sorted(Path(directory).glob(pattern)):
        rows.extend(_read_tsv(str(p)))
    return rows


def _save_fig(fig, output: str, dpi: int = 300):
    """Save figure in the requested format."""
    fig.savefig(output, dpi=dpi, bbox_inches="tight", facecolor="white")
    click.echo(f"Saved: {output}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot group
# ---------------------------------------------------------------------------
@click.group()
def plot():
    """Generate standard population genomics figures from TSV results.

    \b
    Plot types:
      diversity-bar   Per-population diversity ranking
      fst-heatmap     Pairwise Fst matrix with clustering
      manhattan       Genome-wide statistic scan
      pinpis          piN/piS ratios across populations
      sfs-plot        Site frequency spectrum
      roh-landscape   Per-population FROH distribution
    """
    pass


# ---------------------------------------------------------------------------
# diversity-bar
# ---------------------------------------------------------------------------
@plot.command("diversity-bar")
@click.argument("input_dir", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file (PNG/PDF)")
@click.option("--stat", default="pi", help="Statistic to plot (pi, theta_w, tajima_d, fis)")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2, help="Figure width in inches")
@click.option("--height", type=float, default=3.5, help="Figure height in inches")
def diversity_bar(input_dir, output, stat, title, width, height):
    """Plot per-population diversity as a horizontal bar chart.

    INPUT_DIR should contain per-population TSV files from graphpop diversity
    or graphpop run-all (e.g., results/diversity/).

    \b
    Examples:
      graphpop plot diversity-bar results/diversity/ -o fig_diversity.png
      graphpop plot diversity-bar results/diversity/ --stat tajima_d -o fig_tajima.png
    """
    _check_matplotlib()
    _apply_style()

    rows = _read_tsv_dir(input_dir)
    if not rows:
        click.echo("No data found.", err=True)
        return

    # Aggregate by population (mean across chromosomes)
    pop_vals = {}
    for r in rows:
        pop = r.get("population", r.get("file_pop", "unknown"))
        val = float(r.get(stat, 0))
        pop_vals.setdefault(pop, []).append(val)

    pops = sorted(pop_vals.keys(), key=lambda p: np.mean(pop_vals[p]))
    means = [np.mean(pop_vals[p]) for p in pops]
    colors = [WONG_PALETTE[i % len(WONG_PALETTE)] for i in range(len(pops))]

    fig, ax = plt.subplots(figsize=(width, height))
    y = range(len(pops))
    ax.barh(y, means, color=colors, height=0.7, edgecolor="none")
    ax.set_yticks(y)
    ax.set_yticklabels(pops)
    ax.set_xlabel(stat.replace("_", " ").title() if stat != "pi" else "Nucleotide diversity (π)")
    ax.set_title(title or f"Per-population {stat}", fontweight="bold")

    for i, v in enumerate(means):
        ax.text(v + max(means) * 0.01, i, f"{v:.4f}", va="center", fontsize=5)

    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# fst-heatmap
# ---------------------------------------------------------------------------
@plot.command("fst-heatmap")
@click.argument("input_dir", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file")
@click.option("--stat", default="fst_wc", help="Fst statistic (fst_hudson, fst_wc)")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2)
@click.option("--height", type=float, default=6.0)
def fst_heatmap(input_dir, output, stat, title, width, height):
    """Plot pairwise Fst as a heatmap matrix.

    INPUT_DIR should contain pairwise TSV files from graphpop divergence
    or graphpop run-all (e.g., results/divergence/).

    \b
    Examples:
      graphpop plot fst-heatmap results/divergence/ -o fig_fst.png
      graphpop plot fst-heatmap results/divergence/ --stat fst_hudson -o fig_fst_hudson.pdf
    """
    _check_matplotlib()
    _apply_style()

    rows = _read_tsv_dir(input_dir)
    if not rows:
        click.echo("No data found.", err=True)
        return

    # Build pairwise Fst matrix
    pair_vals = {}
    for r in rows:
        p1 = r.get("pop1", "")
        p2 = r.get("pop2", "")
        if not p1 or not p2:
            # Try parsing from filename pattern
            continue
        val = float(r.get(stat, 0))
        pair_vals.setdefault((p1, p2), []).append(val)

    if not pair_vals:
        click.echo("No pairwise data found. Check TSV format.", err=True)
        return

    # Get unique populations
    all_pops = sorted(set(p for pair in pair_vals for p in pair))
    n = len(all_pops)
    pop_idx = {p: i for i, p in enumerate(all_pops)}

    matrix = np.zeros((n, n))
    for (p1, p2), vals in pair_vals.items():
        i, j = pop_idx[p1], pop_idx[p2]
        mean_val = np.mean(vals)
        matrix[i, j] = mean_val
        matrix[j, i] = mean_val

    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(matrix, cmap="YlOrRd", aspect="equal")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(all_pops, rotation=45, ha="right")
    ax.set_yticklabels(all_pops)

    # Annotate cells if small matrix
    if n <= 15:
        for i in range(n):
            for j in range(n):
                if i != j:
                    val = matrix[i, j]
                    color = "white" if val > np.max(matrix) * 0.6 else "black"
                    ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                            fontsize=4, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, label=stat.replace("_", " "))
    ax.set_title(title or f"Pairwise {stat}", fontweight="bold")
    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# manhattan
# ---------------------------------------------------------------------------
@plot.command("manhattan")
@click.argument("input_file", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file")
@click.option("--stat", default="ihs", help="Statistic column name")
@click.option("--threshold", type=float, help="Significance threshold line")
@click.option("--abs-value/--raw-value", default=True, help="Plot absolute values")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2)
@click.option("--height", type=float, default=3.0)
def manhattan(input_file, output, stat, threshold, abs_value, title, width, height):
    """Plot a Manhattan plot of per-variant or per-window statistics.

    INPUT_FILE should be a TSV with columns: pos (or start) and the statistic.
    For multi-chromosome input, include a chr column.

    \b
    Examples:
      graphpop plot manhattan ihs_results.tsv --stat ihs --threshold 2.5 -o fig_ihs.png
      graphpop plot manhattan windows.tsv --stat fst --threshold 0.5 -o fig_fst_scan.png
      graphpop plot manhattan xpehh.tsv --stat xpehh --raw-value -o fig_xpehh.pdf
    """
    _check_matplotlib()
    _apply_style()

    rows = _read_tsv(input_file)
    if not rows:
        click.echo("No data found.", err=True)
        return

    # Extract positions and values
    positions = []
    values = []
    chroms = []
    for r in rows:
        pos = int(r.get("pos", r.get("start", 0)))
        val = float(r.get(stat, r.get(f"{stat}_unstd", 0)))
        if abs_value:
            val = abs(val)
        chrom = r.get("chr", r.get("chromosome", ""))
        positions.append(pos)
        values.append(val)
        chroms.append(chrom)

    fig, ax = plt.subplots(figsize=(width, height))

    # Color by chromosome
    unique_chrs = sorted(set(chroms), key=lambda c: (len(c), c))
    if len(unique_chrs) > 1:
        chr_colors = {c: WONG_PALETTE[i % 2] for i, c in enumerate(unique_chrs)}
        # Make positions additive
        chr_offsets = {}
        offset = 0
        for c in unique_chrs:
            chr_offsets[c] = offset
            chr_positions = [p for p, ch in zip(positions, chroms) if ch == c]
            if chr_positions:
                offset += max(chr_positions) + max(chr_positions) * 0.05

        adj_pos = [p + chr_offsets.get(c, 0) for p, c in zip(positions, chroms)]
        colors = [chr_colors[c] for c in chroms]
        ax.scatter(adj_pos, values, c=colors, s=1, alpha=0.5, rasterized=True)

        # Chromosome labels
        for c in unique_chrs:
            c_positions = [p + chr_offsets[c] for p, ch in zip(positions, chroms) if ch == c]
            if c_positions:
                mid = (min(c_positions) + max(c_positions)) / 2
                label = c.replace("chr", "").replace("Chr", "")
                ax.text(mid, -max(values) * 0.05, label, ha="center", fontsize=4)
        ax.set_xlabel("Chromosome")
    else:
        ax.scatter(positions, values, c=WONG_PALETTE[0], s=1, alpha=0.5, rasterized=True)
        ax.set_xlabel(f"Position on {unique_chrs[0] if unique_chrs else 'chromosome'} (bp)")

    if threshold is not None:
        ax.axhline(threshold, color="#D55E00", linestyle="--", linewidth=0.8, alpha=0.7)

    ylabel = f"|{stat}|" if abs_value else stat
    ax.set_ylabel(ylabel)
    ax.set_title(title or f"Manhattan plot: {stat}", fontweight="bold")
    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# pinpis
# ---------------------------------------------------------------------------
@plot.command("pinpis")
@click.argument("input_file", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2)
@click.option("--height", type=float, default=4.0)
def pinpis(input_file, output, title, width, height):
    """Plot piN/piS ratios across populations.

    INPUT_FILE should be a TSV with columns: population, piN_piS (or piN and piS
    columns to compute the ratio).

    \b
    Generate input:
      # For each population, compute piN and piS:
      graphpop diversity chr1 1 43270923 POP --consequence missense_variant -o piN.tsv
      graphpop diversity chr1 1 43270923 POP --consequence synonymous_variant -o piS.tsv
      # Combine into a single TSV with columns: population, piN_piS

    Examples:
      graphpop plot pinpis pinpis_ratios.tsv -o fig_pinpis.png
    """
    _check_matplotlib()
    _apply_style()

    rows = _read_tsv(input_file)
    if not rows:
        click.echo("No data found.", err=True)
        return

    # Try different column name patterns
    pops = []
    ratios = []
    for r in rows:
        pop = r.get("population", r.get("pop", ""))
        ratio = r.get("piN_piS", r.get("pinpis", r.get("ratio", None)))
        if ratio is None:
            piN = float(r.get("piN", r.get("pi_N", r.get("pi_missense", 0))))
            piS = float(r.get("piS", r.get("pi_S", r.get("pi_synonymous", 1))))
            ratio = piN / piS if piS > 0 else 0
        else:
            ratio = float(ratio)
        pops.append(pop)
        ratios.append(ratio)

    # Sort by ratio
    sorted_pairs = sorted(zip(pops, ratios), key=lambda x: x[1])
    pops = [p for p, _ in sorted_pairs]
    ratios = [r for _, r in sorted_pairs]

    fig, ax = plt.subplots(figsize=(width, height))
    colors = [WONG_PALETTE[0] if r <= 1.0 else WONG_PALETTE[3] for r in ratios]
    bars = ax.barh(range(len(pops)), ratios, color=colors, height=0.7, edgecolor="none")

    ax.axvline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.5,
               label="Neutral expectation (πN/πS = 1)")
    ax.set_yticks(range(len(pops)))
    ax.set_yticklabels(pops)
    ax.set_xlabel("πN/πS ratio")
    ax.set_title(title or "Cost of domestication: πN/πS across populations", fontweight="bold")

    for i, v in enumerate(ratios):
        ax.text(v + max(ratios) * 0.01, i, f"{v:.3f}", va="center", fontsize=5)

    ax.legend(fontsize=5, loc="lower right")
    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# sfs-plot
# ---------------------------------------------------------------------------
@plot.command("sfs-plot")
@click.argument("input_file", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file")
@click.option("--title", help="Figure title")
@click.option("--log-scale/--linear", default=False, help="Use log scale for y-axis")
@click.option("--width", type=float, default=5.0)
@click.option("--height", type=float, default=3.5)
def sfs_plot(input_file, output, title, log_scale, width, height):
    """Plot a site frequency spectrum.

    INPUT_FILE should be a TSV from graphpop sfs with an 'sfs' column
    containing comma-separated counts.

    \b
    Examples:
      graphpop sfs chr22 1 51304566 EUR -o sfs.tsv
      graphpop plot sfs-plot sfs.tsv -o fig_sfs.png
      graphpop plot sfs-plot sfs.tsv --log-scale -o fig_sfs_log.png
    """
    _check_matplotlib()
    _apply_style()

    rows = _read_tsv(input_file)
    if not rows:
        click.echo("No data found.", err=True)
        return

    sfs_str = rows[0].get("sfs", "")
    counts = [int(x) for x in sfs_str.split(",") if x.strip()]

    fig, ax = plt.subplots(figsize=(width, height))
    x = range(len(counts))
    ax.bar(x, counts, color=WONG_PALETTE[0], edgecolor="none", width=0.8)

    if log_scale:
        ax.set_yscale("log")
        ax.set_ylabel("Count (log scale)")
    else:
        ax.set_ylabel("Count")

    ax.set_xlabel("Allele count")
    ax.set_title(title or "Site frequency spectrum", fontweight="bold")

    # Label first and last bins
    if len(counts) > 2:
        ax.set_xticks([0, len(counts) // 4, len(counts) // 2,
                        3 * len(counts) // 4, len(counts) - 1])

    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# roh-landscape
# ---------------------------------------------------------------------------
@plot.command("roh-landscape")
@click.argument("input_dir", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2)
@click.option("--height", type=float, default=4.0)
def roh_landscape(input_dir, output, title, width, height):
    """Plot per-population FROH distribution as violin/box plots.

    INPUT_DIR should contain per-population ROH TSV files from graphpop roh
    or graphpop run-all (e.g., results/roh/).

    \b
    Examples:
      graphpop plot roh-landscape results/roh/ -o fig_roh.png
    """
    _check_matplotlib()
    _apply_style()

    rows = _read_tsv_dir(input_dir)
    if not rows:
        click.echo("No data found.", err=True)
        return

    # Group FROH by population
    pop_froh = {}
    for r in rows:
        pop = r.get("population", r.get("file_pop", "unknown"))
        froh = float(r.get("froh", 0))
        pop_froh.setdefault(pop, []).append(froh)

    pops = sorted(pop_froh.keys(), key=lambda p: np.median(pop_froh[p]))
    data = [pop_froh[p] for p in pops]

    fig, ax = plt.subplots(figsize=(width, height))
    parts = ax.violinplot(data, positions=range(len(pops)), showmedians=True,
                           showextrema=False)

    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(WONG_PALETTE[i % len(WONG_PALETTE)])
        body.set_alpha(0.7)
        body.set_edgecolor("none")
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.0)

    # Add mean markers
    means = [np.mean(d) for d in data]
    ax.scatter(range(len(pops)), means, color="black", s=15, zorder=3, marker="D")

    ax.set_xticks(range(len(pops)))
    ax.set_xticklabels(pops, rotation=45, ha="right")
    ax.set_ylabel("FROH (fraction of genome in ROH)")
    ax.set_title(title or "Inbreeding landscape: FROH by population", fontweight="bold")

    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# gene-zoom
# ---------------------------------------------------------------------------
@plot.command("gene-zoom")
@click.argument("target")
@click.option("--pop", "population", required=True, help="Population name")
@click.option("--pop2", help="Second population (for Fst track)")
@click.option("-o", "--output", required=True, help="Output figure file (PNG/PDF)")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2, help="Figure width in inches")
@click.option("--height", type=float, default=6.0, help="Figure height in inches")
@pass_ctx
def gene_zoom(ctx, target, population, pop2, output, title, width, height):
    """Multi-track regional plot for a gene or genomic region.

    TARGET is either a gene name (e.g., KCNE1) or a region in chr:start-end
    format (e.g., chr6:9000000-9600000). The command resolves gene names to
    coordinates via the Gene node in the graph.

    \b
    Tracks (top to bottom):
      1. Fst (from GenomicWindow or per-variant)
      2. |iHS| (from Variant properties)
      3. Gene model (from HAS_CONSEQUENCE edges)

    \b
    Examples:
      graphpop plot gene-zoom KCNE1 --pop EUR -o fig_kcne1.png
      graphpop plot gene-zoom chr6:9000000-9600000 --pop GJ-tmp -o fig_hd1.png
      graphpop plot gene-zoom GW5 --pop GJ-tmp --pop2 GJ-trop -o fig_gw5.pdf
    """
    _check_matplotlib()
    _apply_style()

    # --- Resolve target to chr, start, end ---
    region_match = re.match(r'^(chr\w+|Chr\w+):(\d+)-(\d+)$', target)
    if region_match:
        chrom = region_match.group(1)
        reg_start = int(region_match.group(2))
        reg_end = int(region_match.group(3))
        region_label = f"{chrom}:{reg_start}-{reg_end}"
    else:
        # Resolve gene name
        recs = ctx.run(
            "MATCH (g:Gene) "
            "WHERE g.symbol = $target OR g.geneId = $target "
            "RETURN g.chr AS chr, g.start AS start, g.end AS end, "
            "g.symbol AS symbol LIMIT 1",
            {"target": target},
        )
        if not recs:
            click.echo(f"Gene '{target}' not found in the graph.", err=True)
            raise SystemExit(1)
        g = recs[0]
        chrom = g["chr"]
        # Pad 20% on each side for context
        gene_len = (g["end"] or 0) - (g["start"] or 0)
        pad = max(gene_len * 0.2, 10000)
        reg_start = max(0, int((g["start"] or 0) - pad))
        reg_end = int((g["end"] or 0) + pad)
        region_label = f"{g['symbol']} ({chrom}:{reg_start}-{reg_end})"

    click.echo(f"Region: {region_label}", err=True)

    # --- Query Fst from GenomicWindow ---
    fst_pos = []
    fst_vals = []
    fst_query = (
        "MATCH (w:GenomicWindow) "
        "WHERE w.chr = $chrom AND w.population = $population "
        "AND w.start >= $reg_start AND w.end <= $reg_end "
        "RETURN w.start AS start, w.end AS end, "
        "w.fst AS fst "
        "ORDER BY w.start"
    )
    region_params = {
        "chrom": chrom, "population": population,
        "reg_start": reg_start, "reg_end": reg_end,
    }
    try:
        fst_recs = ctx.run(fst_query, region_params)
        for r in fst_recs:
            mid = ((r["start"] or 0) + (r["end"] or 0)) / 2
            val = r.get("fst")
            if val is not None:
                fst_pos.append(mid)
                fst_vals.append(float(val))
    except SystemExit:
        click.echo("Warning: no GenomicWindow Fst data.", err=True)

    # --- Query |iHS| from Variant nodes ---
    ihs_prop = f"ihs_{population}"
    ihs_pos = []
    ihs_vals = []
    ihs_query = (
        f"MATCH (v:Variant) "
        f"WHERE v.chr = $chrom AND v.pos >= $reg_start AND v.pos <= $reg_end "
        f"AND v.{ihs_prop} IS NOT NULL "
        f"RETURN v.pos AS pos, v.{ihs_prop} AS ihs "
        f"ORDER BY v.pos"
    )
    try:
        ihs_recs = ctx.run(ihs_query, region_params)
        for r in ihs_recs:
            ihs_pos.append(r["pos"])
            ihs_vals.append(abs(float(r["ihs"])))
    except SystemExit:
        click.echo("Warning: no iHS data.", err=True)

    # --- Query gene models ---
    gene_query = (
        "MATCH (v:Variant)-[hc:HAS_CONSEQUENCE]->(g:Gene) "
        "WHERE v.chr = $chrom AND v.pos >= $reg_start AND v.pos <= $reg_end "
        "RETURN DISTINCT g.symbol AS gene, g.start AS start, g.end AS end "
        "ORDER BY g.start"
    )
    gene_recs = ctx.run(gene_query, region_params)

    # --- Build figure ---
    fig, axes = plt.subplots(3, 1, figsize=(width, height), sharex=True,
                              gridspec_kw={"height_ratios": [2, 2, 1]})

    # Track 1: Fst
    ax_fst = axes[0]
    if fst_pos:
        ax_fst.fill_between(fst_pos, fst_vals, alpha=0.3, color=WONG_PALETTE[0])
        ax_fst.plot(fst_pos, fst_vals, color=WONG_PALETTE[0], linewidth=0.8)
    ax_fst.set_ylabel("Fst")
    ax_fst.set_title(title or f"Gene zoom: {region_label}", fontweight="bold")

    # Track 2: |iHS|
    ax_ihs = axes[1]
    if ihs_pos:
        ax_ihs.scatter(ihs_pos, ihs_vals, s=2, color=WONG_PALETTE[3],
                        alpha=0.6, rasterized=True)
        # Threshold line at 2.0
        ax_ihs.axhline(2.0, color="grey", linestyle="--", linewidth=0.6, alpha=0.5)
    ax_ihs.set_ylabel("|iHS|")

    # Track 3: Gene models
    ax_gene = axes[2]
    if gene_recs:
        y_pos = 0.5
        for i, g in enumerate(gene_recs):
            g_start = g.get("start") or reg_start
            g_end = g.get("end") or reg_end
            color = WONG_PALETTE[i % len(WONG_PALETTE)]
            ax_gene.barh(y_pos, g_end - g_start, left=g_start, height=0.3,
                          color=color, edgecolor="black", linewidth=0.3)
            mid = (g_start + g_end) / 2
            ax_gene.text(mid, y_pos + 0.25, g.get("gene", ""),
                          ha="center", va="bottom", fontsize=5, style="italic")
            y_pos += 0.5
    ax_gene.set_ylabel("Genes")
    ax_gene.set_yticks([])
    ax_gene.set_xlabel(f"Position on {chrom} (bp)")
    ax_gene.set_xlim(reg_start, reg_end)

    # Add vertical lines at peak iHS positions
    if ihs_vals:
        peak_thresh = max(ihs_vals) * 0.9 if max(ihs_vals) > 0 else 999
        for pos, val in zip(ihs_pos, ihs_vals):
            if val >= peak_thresh:
                for ax in axes:
                    ax.axvline(pos, color=WONG_PALETTE[3], linestyle=":",
                               linewidth=0.5, alpha=0.4)

    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# pop-tree
# ---------------------------------------------------------------------------
@plot.command("pop-tree")
@click.argument("input_dir", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file (PNG/PDF)")
@click.option("--method", default="upgma", type=click.Choice(["upgma", "nj"]),
              help="Tree method: upgma (default) or nj")
@click.option("--stat", default="fst_wc", help="Fst statistic column (fst_wc or fst_hudson)")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2, help="Figure width in inches")
@click.option("--height", type=float, default=5.0, help="Figure height in inches")
def pop_tree(input_dir, output, method, stat, title, width, height):
    """Build a UPGMA or neighbor-joining tree from pairwise Fst data.

    INPUT_DIR should contain pairwise divergence TSV files (same format as
    fst-heatmap input, from graphpop divergence or graphpop run-all).

    Uses scipy.cluster.hierarchy for UPGMA clustering. Neighbor-joining is
    approximated via the 'weighted' linkage method.

    \b
    Examples:
      graphpop plot pop-tree results/divergence/ -o fig_tree.png
      graphpop plot pop-tree results/divergence/ --method nj --stat fst_hudson -o fig_nj.png
    """
    _check_matplotlib()
    _apply_style()

    if not HAS_SCIPY:
        click.echo(
            "Error: scipy is required for pop-tree.\n"
            "Install with: pip install scipy",
            err=True,
        )
        raise SystemExit(1)

    rows = _read_tsv_dir(input_dir)
    if not rows:
        click.echo("No data found.", err=True)
        return

    # Build pairwise Fst matrix
    pair_vals = {}
    for r in rows:
        p1 = r.get("pop1", "")
        p2 = r.get("pop2", "")
        if not p1 or not p2:
            continue
        val = float(r.get(stat, 0))
        pair_vals.setdefault((p1, p2), []).append(val)

    if not pair_vals:
        click.echo(f"No pairwise data found. Check TSV format and --stat={stat}.", err=True)
        return

    # Build distance matrix
    all_pops = sorted(set(p for pair in pair_vals for p in pair))
    n = len(all_pops)
    pop_idx = {p: i for i, p in enumerate(all_pops)}

    matrix = np.zeros((n, n))
    for (p1, p2), vals in pair_vals.items():
        i, j = pop_idx[p1], pop_idx[p2]
        mean_val = max(0, np.mean(vals))  # Clamp negative Fst to 0
        matrix[i, j] = mean_val
        matrix[j, i] = mean_val

    # Convert to condensed distance form for scipy
    dist_condensed = squareform(matrix, checks=False)

    # Linkage method
    if method == "upgma":
        linkage_method = "average"
    else:
        # NJ approximation via weighted linkage
        linkage_method = "weighted"

    Z = linkage(dist_condensed, method=linkage_method)

    # Plot dendrogram
    fig, ax = plt.subplots(figsize=(width, height))
    dendrogram(
        Z,
        labels=all_pops,
        ax=ax,
        leaf_rotation=0,
        orientation="left",
        leaf_font_size=7,
        color_threshold=0,
        above_threshold_color=WONG_PALETTE[0],
    )

    ax.set_xlabel(f"Genetic distance ({stat.replace('_', ' ')})")
    tree_label = "UPGMA" if method == "upgma" else "Neighbor-joining"
    ax.set_title(title or f"Population tree ({tree_label}, {stat})", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# chromosome — multi-track chromosome view
# ---------------------------------------------------------------------------
@plot.command("chromosome")
@click.option("--chr", "chrom", required=True, help="Chromosome name (e.g., chr22)")
@click.option("--pop", "population", required=True, help="Population name")
@click.option("--stats", default="fst,ihs",
              help="Comma-separated statistics to plot (e.g., fst,ihs,pi,tajima_d,xpehh)")
@click.option("-o", "--output", required=True, help="Output figure file (PNG/PDF)")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2, help="Figure width in inches")
@click.option("--height", type=float, default=2.0,
              help="Height per track in inches (total = n_tracks * height)")
@pass_ctx
def chromosome(ctx, chrom, population, stats, output, title, width, height):
    """Multi-track chromosome view of population statistics.

    Draws stacked tracks for each requested statistic along the chromosome.
    Window-level stats (fst, pi, tajima_d) are queried from GenomicWindow nodes.
    Variant-level stats (ihs, xpehh) are queried from Variant nodes.

    \b
    Examples:
      graphpop plot chromosome --chr chr22 --pop EUR --stats fst,ihs,pi -o fig_chr22.png
      graphpop plot chromosome --chr Chr01 --pop GJ-tmp --stats fst,pi -o fig_chr01.png
    """
    _check_matplotlib()
    _apply_style()

    stat_list = [s.strip() for s in stats.split(",") if s.strip()]
    if not stat_list:
        click.echo("No statistics specified.", err=True)
        raise SystemExit(1)

    window_stats = {"fst", "pi", "theta_w", "tajima_d"}
    variant_stats = {"ihs", "xpehh"}

    n_tracks = len(stat_list)
    fig_height = height * n_tracks
    fig, axes = plt.subplots(n_tracks, 1, figsize=(width, fig_height), sharex=True,
                              squeeze=False)
    axes = axes.flatten()

    track_colors = [WONG_PALETTE[i % len(WONG_PALETTE)] for i in range(n_tracks)]

    for idx, stat in enumerate(stat_list):
        ax = axes[idx]
        color = track_colors[idx]

        if stat.lower() in window_stats:
            query = (
                f"MATCH (w:GenomicWindow) "
                f"WHERE w.chr = $chrom AND w.population = $population "
                f"AND w.{stat} IS NOT NULL "
                f"RETURN w.start AS start, w.end AS end, w.{stat} AS value "
                f"ORDER BY w.start"
            )
            recs = ctx.run(query, {"chrom": chrom, "population": population})
            if recs:
                positions = [((r["start"] or 0) + (r["end"] or 0)) / 2 for r in recs]
                values = [float(r["value"]) for r in recs]
                ax.fill_between(positions, values, alpha=0.3, color=color)
                ax.plot(positions, values, color=color, linewidth=0.6)
            else:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=7, color="grey")

        elif stat.lower() in variant_stats:
            prop = f"{stat}_{population}"
            query = (
                f"MATCH (v:Variant) "
                f"WHERE v.chr = $chrom AND v.{prop} IS NOT NULL "
                f"RETURN v.pos AS pos, v.{prop} AS value "
                f"ORDER BY v.pos"
            )
            recs = ctx.run(query, {"chrom": chrom})
            if recs:
                positions = [r["pos"] for r in recs]
                values = [abs(float(r["value"])) for r in recs]
                ax.scatter(positions, values, s=1, alpha=0.4, color=color,
                           rasterized=True)
            else:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=7, color="grey")
        else:
            ax.text(0.5, 0.5, f"Unknown stat: {stat}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=7, color="red")

        label = f"|{stat}|" if stat.lower() in variant_stats else stat
        ax.set_ylabel(label, fontsize=7)

        # Alternating background
        if idx % 2 == 1:
            ax.set_facecolor("#f8f8f8")

    axes[-1].set_xlabel(f"Position on {chrom} (bp)")
    axes[0].set_title(
        title or f"Chromosome view: {chrom} ({population})", fontweight="bold"
    )

    fig.tight_layout()
    _save_fig(fig, output)


# ---------------------------------------------------------------------------
# pca-scatter
# ---------------------------------------------------------------------------
@plot.command("pca-scatter")
@click.argument("input_file", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file (PNG/PDF)")
@click.option("--color-by", "color_by", default="population",
              help="Column name to color points by (default: population)")
@click.option("--pc", "pc_axes", default="1,2",
              help="Which PCs to plot, comma-separated (default: 1,2)")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=5.0, help="Figure width in inches")
@click.option("--height", type=float, default=5.0, help="Figure height in inches")
def pca_scatter(input_file, output, color_by, pc_axes, title, width, height):
    """PCA scatter plot from a TSV with PC columns.

    INPUT_FILE should be a TSV with columns like pc1, pc2 (or PC1, PC2) and
    a grouping column (default: population) for coloring.

    \b
    Examples:
      graphpop plot pca-scatter pca_results.tsv -o fig_pca.png
      graphpop plot pca-scatter pca.tsv --color-by superpopulation --pc 1,3 -o pca13.png
    """
    _check_matplotlib()
    _apply_style()

    rows = _read_tsv(input_file)
    if not rows:
        click.echo("No data found.", err=True)
        return

    # Parse PC axes
    try:
        pc_a, pc_b = [int(x.strip()) for x in pc_axes.split(",")]
    except (ValueError, IndexError):
        click.echo("Invalid --pc format. Use e.g. '1,2'.", err=True)
        raise SystemExit(1)

    # Find PC column names (case-insensitive)
    sample_keys = list(rows[0].keys())
    pc_col_a = _find_pc_col(sample_keys, pc_a)
    pc_col_b = _find_pc_col(sample_keys, pc_b)
    if not pc_col_a or not pc_col_b:
        click.echo(
            f"Could not find PC{pc_a} and PC{pc_b} columns in: {sample_keys}",
            err=True,
        )
        raise SystemExit(1)

    # Group by color column
    groups = {}
    for r in rows:
        group = r.get(color_by, "unknown")
        x = float(r[pc_col_a])
        y = float(r[pc_col_b])
        groups.setdefault(group, ([], []))
        groups[group][0].append(x)
        groups[group][1].append(y)

    fig, ax = plt.subplots(figsize=(width, height))
    for i, (group_name, (xs, ys)) in enumerate(sorted(groups.items())):
        color = WONG_PALETTE[i % len(WONG_PALETTE)]
        ax.scatter(xs, ys, s=8, alpha=0.7, color=color, label=group_name,
                   edgecolors="none")

    ax.set_xlabel(f"PC{pc_a}")
    ax.set_ylabel(f"PC{pc_b}")
    ax.set_title(title or f"PCA: PC{pc_a} vs PC{pc_b}", fontweight="bold")

    # Legend outside if many groups
    if len(groups) > 8:
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=5,
                  markerscale=1.5, frameon=False)
    else:
        ax.legend(fontsize=6, markerscale=1.5, frameon=False)

    fig.tight_layout()
    _save_fig(fig, output)


def _find_pc_col(keys: list[str], pc_num: int) -> str | None:
    """Find the column name for a given PC number (case-insensitive)."""
    candidates = [f"pc{pc_num}", f"PC{pc_num}", f"Pc{pc_num}",
                  f"pc_{pc_num}", f"PC_{pc_num}"]
    for c in candidates:
        if c in keys:
            return c
    # Fallback: partial match
    for k in keys:
        if k.lower().replace("_", "") == f"pc{pc_num}":
            return k
    return None


# ---------------------------------------------------------------------------
# heatmap — general-purpose heatmap from a matrix TSV
# ---------------------------------------------------------------------------
@plot.command("heatmap")
@click.argument("input_file", type=click.Path(exists=True))
@click.option("-o", "--output", required=True, help="Output figure file (PNG/PDF)")
@click.option("--cmap", default="viridis", help="Matplotlib colormap (default: viridis)")
@click.option("--annotate", is_flag=True, help="Add numeric values to cells")
@click.option("--title", help="Figure title")
@click.option("--width", type=float, default=7.2, help="Figure width in inches")
@click.option("--height", type=float, default=6.0, help="Figure height in inches")
def heatmap(input_file, output, cmap, annotate, title, width, height):
    """General-purpose heatmap from a matrix TSV.

    INPUT_FILE should be a TSV where the first column contains row labels
    and the header row contains column labels. All other cells are numeric.

    \b
    Examples:
      graphpop plot heatmap matrix.tsv -o fig_heatmap.png --cmap YlOrRd --annotate
      graphpop plot heatmap fst_matrix.tsv -o fig_fst.png --title "Pairwise Fst"
    """
    _check_matplotlib()
    _apply_style()

    rows = _read_tsv(input_file)
    if not rows:
        click.echo("No data found.", err=True)
        return

    # First column = row labels, rest = numeric matrix
    col_keys = list(rows[0].keys())
    label_col = col_keys[0]
    value_cols = col_keys[1:]

    row_labels = [r[label_col] for r in rows]
    matrix = np.zeros((len(rows), len(value_cols)))
    for i, r in enumerate(rows):
        for j, c in enumerate(value_cols):
            try:
                matrix[i, j] = float(r[c])
            except (ValueError, TypeError):
                matrix[i, j] = np.nan

    fig, ax = plt.subplots(figsize=(width, height))
    im = ax.imshow(matrix, cmap=cmap, aspect="auto")

    ax.set_xticks(range(len(value_cols)))
    ax.set_yticks(range(len(row_labels)))
    ax.set_xticklabels(value_cols, rotation=45, ha="right")
    ax.set_yticklabels(row_labels)

    if annotate:
        for i in range(len(row_labels)):
            for j in range(len(value_cols)):
                val = matrix[i, j]
                if not np.isnan(val):
                    text_color = ("white"
                                  if val > np.nanmax(matrix) * 0.6
                                  else "black")
                    ax.text(j, i, f"{val:.3g}", ha="center", va="center",
                            fontsize=5, color=text_color)

    fig.colorbar(im, ax=ax, shrink=0.8)
    ax.set_title(title or "Heatmap", fontweight="bold")
    fig.tight_layout()
    _save_fig(fig, output)
