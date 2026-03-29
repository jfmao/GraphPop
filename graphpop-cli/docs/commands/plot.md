# graphpop plot

Generate standard population genomics figures from TSV results

## Description

`graphpop plot` provides built-in visualization for common population genomics
analyses, producing publication-ready figures directly from GraphPop TSV output.
All plots follow Nature Methods figure guidelines: white background, Arial font,
Wong/Okabe-Ito colorblind-safe palette, and clean axis styling.

Six plot types are available, each designed for a specific analysis output:

| Plot type | Input | Output |
|-----------|-------|--------|
| `diversity-bar` | Directory of diversity TSVs | Per-population diversity ranking |
| `fst-heatmap` | Directory of divergence TSVs | Pairwise Fst matrix with color scale |
| `manhattan` | Single TSV with pos + stat columns | Genome-wide scan plot |
| `pinpis` | TSV with population + piN/piS | Cost of domestication bar chart |
| `sfs-plot` | TSV from `graphpop sfs` | Site frequency spectrum histogram |
| `roh-landscape` | Directory of ROH TSVs | Per-population FROH violin plots |

## Usage

```
graphpop plot <plot-type> [OPTIONS] INPUT -o OUTPUT
```

## Subcommands

### graphpop plot diversity-bar

```
graphpop plot diversity-bar INPUT_DIR -o OUTPUT [--stat pi] [--title TEXT]
```

**Arguments:**

| Name | Type | Description |
|------|------|-------------|
| `INPUT_DIR` | path | Directory containing per-population diversity TSV files |

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o, --output` | path | required | Output figure file (PNG, PDF, SVG) |
| `--stat` | string | `pi` | Statistic to plot: `pi`, `theta_w`, `tajima_d`, `fis` |
| `--title` | string | auto | Figure title |
| `--width` | float | 7.2 | Figure width in inches |
| `--height` | float | 3.5 | Figure height in inches |

**Value:** A horizontal bar chart ranking populations by the chosen statistic.
Values are annotated on each bar. Populations sorted from lowest (top) to
highest (bottom).

### graphpop plot fst-heatmap

```
graphpop plot fst-heatmap INPUT_DIR -o OUTPUT [--stat fst_wc] [--title TEXT]
```

**Arguments:**

| Name | Type | Description |
|------|------|-------------|
| `INPUT_DIR` | path | Directory containing pairwise divergence TSV files |

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o, --output` | path | required | Output figure file |
| `--stat` | string | `fst_wc` | Fst statistic: `fst_hudson`, `fst_wc` |
| `--title` | string | auto | Figure title |

**Value:** A square heatmap matrix with YlOrRd color scale. For matrices with
≤15 populations, individual Fst values are annotated in each cell. Populations
are ordered alphabetically.

### graphpop plot manhattan

```
graphpop plot manhattan INPUT_FILE -o OUTPUT [--stat ihs] [--threshold FLOAT]
                        [--abs-value | --raw-value]
```

**Arguments:**

| Name | Type | Description |
|------|------|-------------|
| `INPUT_FILE` | path | TSV with `pos` (or `start`) and statistic columns |

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o, --output` | path | required | Output figure file |
| `--stat` | string | `ihs` | Column name for the statistic |
| `--threshold` | float | none | Horizontal dashed line at this value |
| `--abs-value` | flag | yes | Plot absolute values (default) |
| `--raw-value` | flag | no | Plot raw signed values |
| `--title` | string | auto | Figure title |

**Value:** A scatter plot with genomic position on x-axis and statistic on y-axis.
For multi-chromosome input (with `chr` column), positions are made additive and
chromosomes are colored in alternating two-tone. Points are rasterized for
efficient file sizes with dense data.

### graphpop plot pinpis

```
graphpop plot pinpis INPUT_FILE -o OUTPUT [--title TEXT]
```

**Arguments:**

| Name | Type | Description |
|------|------|-------------|
| `INPUT_FILE` | path | TSV with `population` and `piN_piS` columns |

**Value:** A horizontal bar chart with a dashed reference line at piN/piS = 1.0
(neutral expectation). Bars are colored blue (≤1.0) or vermillion (>1.0).
Populations sorted by ratio. Values annotated on bars.

### graphpop plot sfs-plot

```
graphpop plot sfs-plot INPUT_FILE -o OUTPUT [--log-scale | --linear]
```

**Arguments:**

| Name | Type | Description |
|------|------|-------------|
| `INPUT_FILE` | path | TSV from `graphpop sfs` with comma-separated `sfs` column |

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--log-scale` | flag | no | Use log scale for y-axis |

**Value:** A bar chart (histogram) of the site frequency spectrum. X-axis is
allele count, y-axis is number of sites.

### graphpop plot roh-landscape

```
graphpop plot roh-landscape INPUT_DIR -o OUTPUT [--title TEXT]
```

**Arguments:**

| Name | Type | Description |
|------|------|-------------|
| `INPUT_DIR` | path | Directory containing per-population ROH TSV files |

**Value:** Violin plots showing the distribution of FROH across individuals for
each population. Diamond markers show population means. Populations sorted by
median FROH. Wong palette colors.

## Details

All plots use the Nature Methods figure style:
- Font: Arial/Helvetica, 7pt default, 8pt titles
- Colors: Wong/Okabe-Ito colorblind-safe palette
- Background: white, no gridlines
- Axes: top and right spines removed, outward ticks
- Output: 300 DPI for PNG, vector for PDF/SVG

The `--width` and `--height` options (in inches) default to Nature Methods
double-column width (7.2 inches = 183mm). For single-column figures, use
`--width 3.5`.

Requires `matplotlib` and `numpy`. Install with:
```bash
pip install matplotlib numpy
```

## Examples

```bash
# Full workflow: compute → plot
graphpop run-all --phase 1 -d results/
graphpop plot diversity-bar results/diversity/ -o fig_diversity.png
graphpop plot fst-heatmap results/divergence/ -o fig_fst.png
graphpop plot roh-landscape results/roh/ -o fig_roh.png

# Selection scan visualization
graphpop ihs chr22 EUR --persist -o ihs_chr22.tsv
graphpop plot manhattan ihs_chr22.tsv --stat ihs --threshold 2.5 -o fig_ihs.png

# SFS comparison
graphpop sfs chr22 1 51304566 YRI -o sfs_yri.tsv
graphpop plot sfs-plot sfs_yri.tsv --title "SFS: YRI" -o fig_sfs_yri.png

# piN/piS visualization
graphpop plot pinpis pinpis_ratios.tsv -o fig_pinpis.pdf
```

## See Also

- `graphpop diversity` — compute diversity statistics
- `graphpop divergence` — compute pairwise Fst
- `graphpop ihs`, `graphpop xpehh`, `graphpop nsl` — selection scans
- `graphpop roh` — runs of homozygosity
- `graphpop sfs` — site frequency spectrum
- `graphpop run-all` — orchestrate full-genome analysis
- `graphpop aggregate` — generate summary tables
