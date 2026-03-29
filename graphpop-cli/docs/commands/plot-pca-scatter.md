# Plot PCA scatter colored by population

## Description

`graphpop plot pca-scatter` produces a scatter plot of principal component
analysis (PCA) results, with points colored by population label. PCA is the
standard method for visualizing population structure: each point represents an
individual, and the axes represent the top principal components of genetic
variation. Clusters in the scatter plot correspond to genetically distinct
populations.

The input is a TSV or CSV file with sample identifiers, population labels, and
PC coordinates. This file is typically produced by external PCA tools (PLINK,
smartpca) or by `graphpop query` extracting PCA results from the graph.

## Usage

```
graphpop plot pca-scatter INPUT_FILE -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `INPUT_FILE` | path | yes | -- | TSV or CSV file with columns for sample ID, population, and PC coordinates. |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--pc` | string | `1,2` | Comma-separated pair of PC numbers to plot (e.g., `1,2` or `2,3`). |
| `--color-by` | string | `population` | Column name to use for point coloring. |
| `--sample-col` | string | `sampleId` | Column name for sample identifiers. |
| `--pc-prefix` | string | `PC` | Column name prefix for PC coordinates (columns are expected as `PC1`, `PC2`, etc.). |
| `--highlight` | string | none | Comma-separated list of population names to highlight (other populations shown in grey). |
| `--alpha` | float | 0.7 | Point transparency (0.0 = fully transparent, 1.0 = opaque). |
| `--point-size` | float | 15 | Point size in matplotlib units. |
| `--legend` | choice | `auto` | Legend placement: `auto`, `right`, `none`. |
| `--title` | string | auto | Figure title. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | 5.0 | Figure height in inches. |

## Value

A scatter plot with the chosen PCs on x and y axes. Points are colored by the
`--color-by` column using the Wong/Okabe-Ito palette (up to 8 colors) or a
categorical colormap for larger numbers of groups. Axis labels show the PC
number and variance explained (if a `variance_explained` column is present in
the input).

## Details

### Input format

The input file must contain at minimum:
- A sample identifier column (default name: `sampleId`).
- A grouping column for coloring (default name: `population`).
- Two or more PC coordinate columns (default names: `PC1`, `PC2`, ...).

Column names can be overridden with `--sample-col`, `--color-by`, and
`--pc-prefix`. The file may be tab-separated or comma-separated (auto-detected).

### Variance explained

If the input file contains a header comment line starting with `#` that includes
variance explained percentages, or if columns `variance_explained_1`,
`variance_explained_2`, etc. are present, the axis labels will display the
percentage (e.g., "PC1 (8.3%)").

### Highlighting

When `--highlight` is specified, only the named populations are drawn with full
color. All other populations are drawn in light grey with reduced alpha. This is
useful for emphasizing specific populations in a crowded plot.

### Figure style

Follows Nature Methods guidelines. See `graphpop plot` for full style details.

## Examples

```bash
# Basic PCA scatter of PC1 vs PC2
graphpop plot pca-scatter pca_results.tsv -o pca_scatter.png

# Plot PC2 vs PC3, highlight East Asian populations
graphpop plot pca-scatter pca_results.tsv -o pca_23.png \
    --pc 2,3 --highlight CHB,JPT,CHS,CDX,KHV

# Custom column names from PLINK output
graphpop plot pca-scatter plink_pca.eigenvec -o pca_plink.png \
    --sample-col IID --color-by FID --pc-prefix PC

# PDF for publication at single-column width
graphpop plot pca-scatter pca_results.tsv -o fig_pca.pdf \
    --width 3.5 --height 3.5 --point-size 8 --alpha 0.6

# Color by super-population instead of population
graphpop plot pca-scatter pca_results.tsv -o pca_superpop.png \
    --color-by super_pop --legend right
```

## See Also

- `graphpop plot pop-tree` -- population tree from Fst distances
- `graphpop plot fst-heatmap` -- Fst matrix heatmap
- `graphpop plot heatmap` -- general-purpose heatmap
- `graphpop plot` -- overview of all plot subcommands and shared style options
