# Plot diversity bar chart

## Description

`graphpop plot diversity-bar` produces a horizontal bar chart ranking populations
by a chosen diversity statistic. Each bar represents one population, with the
numeric value annotated directly on the bar. Populations are sorted from lowest
(top) to highest (bottom), providing an immediate visual comparison of genetic
diversity across the dataset.

This is a standalone reference for the `diversity-bar` subcommand. For shared
plot style details (font, palette, DPI, figure guidelines), see
[graphpop plot](plot.md).

## Usage

```
graphpop plot diversity-bar INPUT_DIR -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `INPUT_DIR` | path | yes | -- | Directory containing per-population diversity TSV files (one file per population, as produced by `graphpop pop-summary` or `graphpop run-all`). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--stat` | string | `pi` | Statistic to plot. One of: `pi`, `theta_w`, `tajima_d`, `fis`. |
| `--title` | string | auto | Figure title. Defaults to the statistic name. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | 3.5 | Figure height in inches. Auto-adjusts for many populations. |

## Value

A horizontal bar chart. Each bar is colored from the Wong/Okabe-Ito palette.
Numeric values are annotated at the end of each bar using the format appropriate
for the statistic (4 decimal places for pi and theta_w, 2 decimal places for
Tajima's D and F_IS). An x-axis label identifies the statistic.

## Details

### Input file discovery

The command scans `INPUT_DIR` for TSV files matching the pattern
`diversity_*.tsv` or `*_diversity.tsv`. Each file must contain a header row with
the statistic column names (pi, theta_w, tajima_d, fis). The population name is
extracted from the filename.

### Sorting

Populations are sorted by the chosen statistic in ascending order (lowest at
top). For Tajima's D, which can be negative, this places the most negative
values at top.

For full style details, see [graphpop plot](plot.md).

## Examples

```bash
# Bar chart of nucleotide diversity
graphpop plot diversity-bar results/diversity/ -o fig_pi.png

# Compare Watterson's theta across populations
graphpop plot diversity-bar results/diversity/ -o fig_theta.png --stat theta_w

# Tajima's D comparison with custom title
graphpop plot diversity-bar results/diversity/ -o fig_tajima.pdf \
    --stat tajima_d --title "Tajima's D by Population"

# Inbreeding coefficient comparison, compact figure
graphpop plot diversity-bar results/diversity/ -o fig_fis.png \
    --stat fis --width 5 --height 4

# Rice dataset diversity ranking
graphpop plot diversity-bar rice_results/diversity/ -o rice_pi.png \
    --stat pi --title "Rice 3K: Nucleotide Diversity"
```

## See Also

- `graphpop plot` -- full plot documentation with all subcommands and style details
- `graphpop diversity` -- compute diversity statistics
- `graphpop pop-summary` -- whole-chromosome diversity with persistence
- `graphpop plot fst-heatmap` -- pairwise Fst visualization
