# Plot site frequency spectrum

## Description

`graphpop plot sfs-plot` creates a bar chart of the site frequency spectrum
(SFS) from a TSV file produced by `graphpop sfs`. Each bar represents the
count of variants at a given allele count bin.

## Usage

```
graphpop plot sfs-plot INPUT_FILE -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Description |
|------|------|----------|-------------|
| `INPUT_FILE` | path | yes | TSV file from `graphpop sfs` containing an `sfs` column with comma-separated counts. |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | *required* | Output figure file (PNG or PDF). |
| `--title` | string | auto | Figure title. |
| `--log-scale` / `--linear` | flag | `--linear` | Use log scale for the y-axis. |
| `--width` | float | 5.0 | Figure width in inches. |
| `--height` | float | 3.5 | Figure height in inches. |

## Details

The input TSV should have an `sfs` column containing comma-separated integer
counts (one per allele count bin). This is the default output format of
`graphpop sfs`. Both folded and unfolded SFS are supported — the plot
automatically adjusts the x-axis labels.

For large sample sizes, the `--log-scale` option is recommended to make
rare-variant bins visible alongside the dominant low-frequency class.

## Examples

```bash
# Compute SFS then plot it
graphpop sfs chr22 1 51304566 EUR -o sfs.tsv
graphpop plot sfs-plot sfs.tsv -o fig_sfs.png

# Log-scale SFS plot
graphpop plot sfs-plot sfs.tsv --log-scale -o fig_sfs_log.png

# Unfolded SFS
graphpop sfs chr22 1 51304566 EUR --unfolded -o sfs_unfolded.tsv
graphpop plot sfs-plot sfs_unfolded.tsv --title "Unfolded SFS (EUR)" -o fig_usfs.png
```

## See Also

- `graphpop sfs` -- compute site frequency spectrum
- `graphpop joint-sfs` -- compute 2D joint SFS
- `graphpop plot manhattan` -- genome-wide Manhattan plot
- `graphpop plot diversity-bar` -- per-population diversity bar chart
