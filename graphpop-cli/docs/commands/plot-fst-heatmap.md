# Plot pairwise Fst heatmap

## Description

`graphpop plot fst-heatmap` produces a square heatmap of pairwise Fst values
between all populations, providing a visual summary of population
differentiation. Cells are colored on a YlOrRd scale, and numeric values are
annotated in each cell when the number of populations is 15 or fewer.

This is a standalone reference for the `fst-heatmap` subcommand. For shared
plot style details (font, palette, DPI, figure guidelines), see
[graphpop plot](plot.md).

## Usage

```
graphpop plot fst-heatmap INPUT_DIR -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `INPUT_DIR` | path | yes | -- | Directory containing pairwise divergence TSV files (as produced by `graphpop divergence` or `graphpop run-all`). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--stat` | string | `fst_wc` | Fst statistic to plot: `fst_wc` (Weir-Cockerham) or `fst_hudson` (Hudson). |
| `--title` | string | auto | Figure title. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | 7.2 | Figure height in inches. |
| `--cmap` | string | `YlOrRd` | Matplotlib colormap. |
| `--annotate` | flag | auto | Force cell annotation on/off. Default: on for <=15 populations. |
| `--fmt` | string | `.3f` | Format string for cell annotations. |

## Value

A square heatmap with populations on both axes, ordered alphabetically. The
diagonal is blank (Fst with self is undefined). The color bar on the right
indicates the Fst scale. For small matrices, numeric Fst values are displayed
in each cell (white text on dark cells for readability).

## Details

### Matrix assembly

The command reads all pairwise divergence TSV files from the input directory.
Each file is expected to contain a `pop1`, `pop2`, and `fst_wc` (or
`fst_hudson`) column. The command assembles the symmetric matrix from all
available pairs. Missing pairs produce an error.

### Negative Fst handling

Negative Fst estimates (which can occur for very closely related populations
due to sampling variance) are displayed as-is in annotations but clipped to
0.0 for color mapping.

For full style details, see [graphpop plot](plot.md).

## Examples

```bash
# Standard Fst heatmap
graphpop plot fst-heatmap results/divergence/ -o fig_fst.png

# Hudson's Fst with custom color map
graphpop plot fst-heatmap results/divergence/ -o fig_fst_hudson.png \
    --stat fst_hudson --cmap viridis

# PDF for publication
graphpop plot fst-heatmap results/divergence/ -o fig_fst.pdf \
    --width 3.5 --height 3.5 --fmt ".2f"

# Rice population Fst matrix with title
graphpop plot fst-heatmap rice_results/divergence/ -o rice_fst.png \
    --title "Rice 3K Pairwise Fst (Weir-Cockerham)"

# Large matrix without annotations
graphpop plot fst-heatmap results/divergence/ -o fst_large.png \
    --annotate --width 10 --height 10
```

## See Also

- `graphpop plot` -- full plot documentation with all subcommands and style details
- `graphpop divergence` -- compute pairwise Fst
- `graphpop plot pop-tree` -- population tree from Fst matrix
- `graphpop plot diversity-bar` -- diversity ranking bar chart
