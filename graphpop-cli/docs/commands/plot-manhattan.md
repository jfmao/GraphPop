# Plot Manhattan plot for genome-wide selection scans

## Description

`graphpop plot manhattan` produces a Manhattan plot: a scatter plot with genomic
position on the x-axis and a selection or summary statistic on the y-axis. This
is the standard visualization for identifying outlier peaks across a chromosome
or genome. Points exceeding a significance threshold are highlighted, and for
multi-chromosome input, chromosomes are displayed in alternating colors.

This is a standalone reference for the `manhattan` subcommand. For shared plot
style details (font, palette, DPI, figure guidelines), see
[graphpop plot](plot.md).

## Usage

```
graphpop plot manhattan INPUT_FILE -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `INPUT_FILE` | path | yes | -- | TSV or CSV file with at minimum a position column (`pos` or `start`) and a statistic column. Optionally includes a `chr` column for multi-chromosome input. |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--stat` | string | `ihs` | Column name for the statistic to plot. |
| `--threshold` | float | none | Draw a horizontal dashed line at this value. Points above are highlighted in vermillion. |
| `--abs-value` | flag | yes | Plot absolute values of the statistic (default). |
| `--raw-value` | flag | no | Plot raw signed values instead of absolute values. |
| `--title` | string | auto | Figure title. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | 3.0 | Figure height in inches. |
| `--point-size` | float | 2 | Point size in matplotlib units. |
| `--chr-col` | string | `chr` | Column name for chromosome identifier. |
| `--pos-col` | string | `pos` | Column name for position. Also accepts `start`. |

## Value

A scatter plot. For single-chromosome input, points are drawn in a single color
with position on the x-axis. For multi-chromosome input, positions are converted
to genome-wide additive coordinates and chromosomes alternate between two shades
(blue-grey / steel-blue). Points above `--threshold` are drawn in vermillion.
Points are rasterized for efficient file sizes when the input exceeds 100,000
rows.

## Details

### Multi-chromosome handling

When the input contains a `chr` column, the command sorts chromosomes
numerically (chr1, chr2, ..., chr22, chrX) and computes cumulative offsets.
Chromosome labels are placed at the midpoint of each chromosome's span on the
x-axis.

### Rasterization

To keep file sizes manageable for genome-wide plots with millions of points,
scatter points are rasterized at 300 DPI while axes and text remain vector.
This produces clean output in both PNG and PDF formats.

For full style details, see [graphpop plot](plot.md).

## Examples

```bash
# iHS Manhattan plot with significance threshold
graphpop plot manhattan ihs_chr22.tsv -o fig_ihs_manhattan.png \
    --stat ihs --threshold 2.5

# XP-EHH across the genome (multi-chromosome file)
graphpop plot manhattan xpehh_genome.tsv -o fig_xpehh.png \
    --stat xpehh --threshold 3.0 --title "XP-EHH: EUR vs AFR"

# Raw (signed) values for Tajima's D
graphpop plot manhattan tajima_windows.tsv -o fig_tajima.png \
    --stat tajima_d --raw-value --threshold -2.0

# Custom column names
graphpop plot manhattan scan_results.tsv -o fig_scan.png \
    --stat fst_wc --pos-col start --chr-col chromosome --threshold 0.4

# Compact figure for supplementary material
graphpop plot manhattan nsl_results.tsv -o sfig_nsl.pdf \
    --stat nsl --threshold 2.5 --width 7.2 --height 2.5
```

## See Also

- `graphpop plot` -- full plot documentation with all subcommands and style details
- `graphpop plot gene-zoom` -- zoomed multi-track view of a specific region
- `graphpop plot chromosome` -- multi-track single-chromosome view
- `graphpop ihs` -- compute iHS scores
- `graphpop genome-scan` -- compute windowed statistics
