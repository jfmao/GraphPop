# Plot multi-track chromosome view

## Description

`graphpop plot chromosome` produces a multi-track figure spanning an entire
chromosome, showing one or more statistics as parallel tracks with a shared
genomic position axis. This is the standard view for visualizing how selection
and diversity signals vary along a chromosome, identifying outlier peaks, and
correlating patterns across statistics.

Each track displays windowed or variant-level statistics as a scatter or line
plot. Gene density can be shown as an additional annotation track. Use this
command after running `graphpop genome-scan` and selection scans with `--persist`
to visualize all results in a single figure.

## Usage

```
graphpop plot chromosome --chr CHR --pop POP -o OUTPUT [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--chr` | string | required | Chromosome to plot. |
| `--pop` | string | required | Population for within-population statistics. |
| `--pop2` | string | none | Second population for pairwise statistics (Fst, XP-EHH). |
| `--stats` | string | `fst,ihs` | Comma-separated list of statistics to plot as tracks. Available: `fst`, `pi`, `tajima_d`, `ihs`, `xpehh`, `nsl`, `h12`, `dxy`, `theta_w`. |
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--thresholds` | string | none | Comma-separated threshold values, one per statistic in `--stats` order. Draws horizontal dashed lines. |
| `--gene-track` | flag | false | Add a gene density track at the bottom. |
| `--smoothing` | integer | none | Smoothing window (number of points) for a running average overlay on each track. |
| `--title` | string | auto | Figure title. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | none | Figure height in inches. Auto-calculated from number of tracks (1.5 inches per track). |

## Value

A vertically stacked multi-panel figure. All tracks share the same x-axis
(genomic position in Mb). Each track shows the named statistic with points
colored by magnitude. If `--thresholds` are provided, horizontal dashed lines
mark the thresholds and points exceeding them are highlighted.

## Details

### Track layout

Tracks are stacked top to bottom in the order specified by `--stats`. Each track
has its own y-axis with an appropriate label and scale. The x-axis shows
positions in megabases with tick marks every 10 Mb.

### Data sources

- Window-level statistics (fst, pi, tajima_d, theta_w, dxy) are read from
  GenomicWindow nodes for the specified chromosome and population(s).
- Variant-level statistics (ihs, xpehh, nsl, h12) are read from Variant nodes
  with persisted score properties. Absolute values are plotted for iHS, XP-EHH,
  and nSL.

### Smoothing

When `--smoothing N` is specified, a running average with a window of N data
points is overlaid as a colored line on each track. The raw data points remain
visible as semi-transparent dots underneath.

### Figure style

Follows Nature Methods guidelines. See `graphpop plot` for full style details.
Chromosome plots default to double-column width (7.2 inches) with height
auto-scaled to the number of tracks.

## Examples

```bash
# Two-track chromosome view: Fst and iHS
graphpop plot chromosome --chr chr22 --pop EUR --pop2 AFR \
    --stats fst,ihs --thresholds 0.3,2.5 -o chr22_tracks.png

# Four-track view with gene density
graphpop plot chromosome --chr chr1 --pop GJ-tmp --pop2 XI-1A \
    --stats fst,pi,ihs,tajima_d --gene-track -o chr1_rice.png

# Smoothed diversity landscape
graphpop plot chromosome --chr chr22 --pop YRI \
    --stats pi,theta_w --smoothing 20 -o diversity_landscape.pdf

# Single-statistic view with custom threshold
graphpop plot chromosome --chr chr22 --pop EUR \
    --stats ihs --thresholds 3.0 --title "iHS scan: EUR chr22" \
    -o ihs_chr22.png

# PDF with custom dimensions
graphpop plot chromosome --chr chr1 --pop CEU --pop2 YRI \
    --stats fst,xpehh,h12 --width 10 --height 6 -o chr1_selection.pdf
```

## See Also

- `graphpop plot gene-zoom` -- zoomed multi-track view centered on a gene
- `graphpop plot manhattan` -- genome-wide Manhattan plot
- `graphpop genome-scan` -- compute windowed statistics
- `graphpop plot` -- overview of all plot subcommands and shared style options
