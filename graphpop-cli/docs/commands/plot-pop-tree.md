# Plot population tree from Fst matrix

## Description

`graphpop plot pop-tree` constructs and visualizes a population tree from a
pairwise Fst distance matrix. The tree is built using either UPGMA (default) or
neighbor-joining (NJ) clustering, and rendered as a dendrogram with branch
lengths proportional to genetic distance. This provides a quick overview of
population relationships based on genome-wide differentiation.

The input is a directory of pairwise divergence TSV files (as produced by
`graphpop divergence` or `graphpop run-all`). The command reads Fst values,
constructs a symmetric distance matrix, applies the chosen clustering method,
and renders the tree.

## Usage

```
graphpop plot pop-tree INPUT_DIR -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `INPUT_DIR` | path | yes | -- | Directory containing pairwise divergence TSV files with `fst_wc` or `fst_hudson` columns. |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--method` | choice | `upgma` | Clustering method: `upgma` (Unweighted Pair Group Method with Arithmetic Mean) or `nj` (Neighbor-Joining). |
| `--stat` | string | `fst_wc` | Fst statistic to use: `fst_wc` (Weir-Cockerham) or `fst_hudson`. |
| `--color-by` | string | none | Color leaf labels by a grouping variable. Accepts a TSV file with columns `population` and `group` (e.g., super-population assignments). |
| `--title` | string | auto | Figure title. |
| `--width` | float | 5.0 | Figure width in inches. |
| `--height` | float | none | Figure height in inches. Auto-calculated from number of populations. |
| `--orientation` | choice | `horizontal` | Tree orientation: `horizontal` (root left, leaves right) or `vertical` (root top, leaves bottom). |

## Value

A dendrogram figure with population names as leaf labels and branch lengths
proportional to Fst distance. When `--color-by` is specified, leaf labels are
colored by group assignment.

## Details

### Tree construction

- **UPGMA** assumes a molecular clock (equal rates across lineages) and produces
  an ultrametric tree. Suitable when populations diverged under relatively
  uniform drift.
- **NJ** (Saitou and Nei 1987) does not assume a clock and produces an unrooted
  tree rendered with midpoint rooting. More appropriate when drift rates differ
  substantially across populations.

### Distance matrix

The command reads all pairwise TSV files in the input directory, extracts the
specified Fst statistic, and assembles a symmetric matrix. Missing pairs produce
an error. Negative Fst values (which can occur for very closely related
populations) are set to 0 before tree construction.

### Figure style

Follows Nature Methods guidelines. Branch lines are drawn in grey. Leaf labels
use 7pt Arial. When `--color-by` is used, colors are drawn from the
Wong/Okabe-Ito palette. See `graphpop plot` for full style details.

## Examples

```bash
# UPGMA tree from Fst matrix
graphpop plot pop-tree results/divergence/ -o pop_tree.png

# Neighbor-joining tree colored by super-population
graphpop plot pop-tree results/divergence/ -o pop_tree_nj.pdf \
    --method nj --color-by superpop_groups.tsv

# Vertical orientation with Hudson's Fst
graphpop plot pop-tree results/divergence/ -o pop_tree_vert.png \
    --stat fst_hudson --orientation vertical

# Custom dimensions for a large number of populations
graphpop plot pop-tree rice_divergence/ -o rice_tree.png \
    --width 6 --height 8 --title "Rice 3K Population Tree"

# PDF for publication
graphpop plot pop-tree results/divergence/ -o fig_tree.pdf \
    --method nj --width 3.5 --title "Population relationships (NJ, Fst)"
```

## See Also

- `graphpop plot fst-heatmap` -- Fst matrix as a heatmap
- `graphpop divergence` -- compute pairwise Fst
- `graphpop plot pca-scatter` -- PCA-based population visualization
- `graphpop plot` -- overview of all plot subcommands and shared style options
