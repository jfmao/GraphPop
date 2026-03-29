# Plot general-purpose heatmap

## Description

`graphpop plot heatmap` renders a matrix of numeric values as a color-scaled
heatmap. It is a general-purpose visualization tool for any rectangular data:
correlation matrices, gene-by-population statistic tables, window-by-statistic
comparison matrices, or any TSV/CSV with row labels, column headers, and numeric
cells.

The first column is used as row labels. All remaining columns are treated as
numeric data columns. The heatmap supports cell annotation, optional row/column
clustering, and configurable color maps.

## Usage

```
graphpop plot heatmap INPUT_FILE -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `INPUT_FILE` | path | yes | -- | TSV or CSV file with row labels in the first column and numeric values in remaining columns. |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--cmap` | string | `YlOrRd` | Matplotlib colormap name (e.g., `viridis`, `RdBu_r`, `coolwarm`, `YlOrRd`). |
| `--annotate` | flag | false | Display numeric values in each cell. |
| `--cluster-rows` | flag | false | Hierarchically cluster rows and add a dendrogram. |
| `--cluster-cols` | flag | false | Hierarchically cluster columns and add a dendrogram. |
| `--vmin` | float | auto | Minimum value for color scale. |
| `--vmax` | float | auto | Maximum value for color scale. |
| `--fmt` | string | `.2f` | Format string for cell annotations (Python format spec). |
| `--title` | string | none | Figure title. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | none | Figure height auto-calculated from matrix dimensions. |

## Value

A heatmap figure with row labels on the y-axis and column labels on the x-axis.
A color bar is shown on the right side. When `--annotate` is set, numeric values
are printed in each cell (white text on dark cells, black text on light cells
for readability). When clustering is enabled, dendrograms are drawn on the
corresponding axis.

## Details

### Color map selection

- **Diverging** (`RdBu_r`, `coolwarm`): best for data centered around zero
  (e.g., correlation, Tajima's D, delta Fst). Set `--vmin` and `--vmax`
  symmetrically.
- **Sequential** (`YlOrRd`, `viridis`): best for data with a natural zero
  (e.g., Fst, pi, allele frequency).

### Clustering

Row and column clustering use Ward's minimum variance method on Euclidean
distances. The resulting dendrograms are drawn to the left of (rows) or above
(columns) the heatmap. Clustering reorders the matrix to group similar
rows/columns together.

### Large matrices

For matrices with more than 50 rows or columns, cell annotations are
automatically suppressed (even if `--annotate` is set) to maintain readability.
Override with `--fmt` and a shorter format string if needed.

### Figure style

Follows Nature Methods guidelines. See `graphpop plot` for full style details.

## Examples

```bash
# Fst correlation matrix with annotations
graphpop plot heatmap fst_matrix.tsv -o fst_heatmap.png \
    --cmap YlOrRd --annotate --fmt ".3f"

# Gene-by-population statistic table with clustering
graphpop plot heatmap gene_stats.tsv -o gene_heatmap.png \
    --cluster-rows --cluster-cols --cmap viridis

# Diverging color map for correlation matrix
graphpop plot heatmap stat_correlations.tsv -o corr_heatmap.pdf \
    --cmap RdBu_r --vmin -1 --vmax 1 --annotate

# Custom dimensions for a wide matrix
graphpop plot heatmap pathway_scores.tsv -o pathway_heatmap.png \
    --width 10 --height 6 --title "Pathway Selection Scores"

# Minimal heatmap without annotations
graphpop plot heatmap window_matrix.tsv -o windows.png --cmap inferno
```

## See Also

- `graphpop plot fst-heatmap` -- specialized Fst matrix heatmap
- `graphpop plot diversity-bar` -- bar chart for diversity statistics
- `graphpop plot` -- overview of all plot subcommands and shared style options
