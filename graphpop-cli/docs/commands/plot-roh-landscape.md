# Plot ROH landscape violin plots

## Description

`graphpop plot roh-landscape` visualizes the distribution of genomic inbreeding
coefficients (FROH) across individuals within each population, using violin
plots. FROH is the fraction of the genome contained in runs of homozygosity
(ROH), and its distribution reveals population history: bottlenecked or inbred
populations show higher and more dispersed FROH values, while large outbred
populations cluster near zero.

This is a standalone reference for the `roh-landscape` subcommand. For shared
plot style details (font, palette, DPI, figure guidelines), see
[graphpop plot](plot.md).

## Usage

```
graphpop plot roh-landscape INPUT_DIR -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `INPUT_DIR` | path | yes | -- | Directory containing per-population ROH TSV files (as produced by `graphpop roh` or `graphpop run-all`). Each file must contain per-individual FROH values. |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--title` | string | auto | Figure title. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | 4.0 | Figure height in inches. Auto-adjusts for many populations. |
| `--sort-by` | choice | `median` | Sort populations by: `median`, `mean`, `name`. |

## Value

A figure with one violin plot per population, arranged horizontally. The y-axis
shows FROH (0.0 to 1.0). Each violin shows the kernel density estimate of the
FROH distribution for that population. Diamond markers indicate population means.
Populations are sorted by median FROH (default) so that the most inbred
populations appear on the right.

## Details

### Input file discovery

The command scans `INPUT_DIR` for TSV files matching `roh_*.tsv` or
`*_roh.tsv`. Each file must contain columns `sampleId` and `froh`. The
population name is extracted from the filename. Files without a `froh` column
are skipped with a warning.

### Violin rendering

Violins are drawn with the Wong/Okabe-Ito palette, cycling colors across
populations. The violin width is normalized so that the widest violin spans the
full allocated width. Interior quartile lines (25th, 50th, 75th percentiles)
are drawn as thin horizontal lines within each violin.

### Sorting

Default sorting by median FROH places populations with the strongest
bottleneck/inbreeding signals on the right side of the plot. Use `--sort-by
name` for alphabetical order or `--sort-by mean` for mean-based ordering.

For full style details, see [graphpop plot](plot.md).

## Examples

```bash
# Standard ROH landscape violin plot
graphpop plot roh-landscape results/roh/ -o fig_roh.png

# Sort by population name instead of median FROH
graphpop plot roh-landscape results/roh/ -o fig_roh_alpha.png \
    --sort-by name

# Rice dataset with custom title
graphpop plot roh-landscape rice_results/roh/ -o rice_roh.pdf \
    --title "Rice 3K: Genomic Inbreeding (FROH) by Subpopulation"

# Compact figure for supplementary
graphpop plot roh-landscape results/roh/ -o sfig_roh.png \
    --width 5 --height 3

# Publication figure at single-column width
graphpop plot roh-landscape results/roh/ -o fig5_roh.pdf \
    --width 3.5 --height 3.0 --sort-by median
```

## See Also

- `graphpop plot` -- full plot documentation with all subcommands and style details
- `graphpop roh` -- compute runs of homozygosity
- `graphpop plot diversity-bar` -- diversity ranking bar chart
- `graphpop plot fst-heatmap` -- pairwise Fst visualization
