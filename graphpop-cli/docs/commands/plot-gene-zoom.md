# Plot multi-track gene zoom view

## Description

`graphpop plot gene-zoom` produces a multi-track figure centered on a gene or
genomic interval, showing Fst, |iHS|, and gene models in vertically stacked
panels. This is the primary visualization for examining selection signals in
their gene-level context. The top track shows Fst (per-window or per-variant),
the middle track shows |iHS| (per-variant), and the bottom track shows gene
model rectangles with exon/intron structure.

The target can be specified as a gene name (resolved to coordinates via Gene
nodes in the graph) or as an explicit `chr:start-end` interval. When a gene
name is used, the plot is padded by `--flank` base pairs on each side to show
the surrounding genomic context.

## Usage

```
graphpop plot gene-zoom TARGET --pop POP -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `TARGET` | string | yes | -- | Gene symbol (e.g., `LCT`, `Os01g0100100`) or genomic interval in `chr:start-end` format (e.g., `chr22:20000000-21000000`). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop` | string | required | Population for iHS and Fst statistics. |
| `--pop2` | string | none | Second population for Fst. If omitted, Fst track is hidden. |
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--flank` | integer | 100000 | Base pairs of flanking sequence on each side when targeting a gene. |
| `--tracks` | string | `fst,ihs,genes` | Comma-separated tracks to display. Available: `fst`, `ihs`, `xpehh`, `nsl`, `h12`, `genes`. |
| `--threshold` | float | 2.0 | Horizontal threshold line on selection statistic tracks. |
| `--title` | string | auto | Figure title. Defaults to gene name or interval. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | 5.0 | Figure height in inches. |
| `--highlight` | string | none | Highlight a sub-region within the plot (format: `start-end`). Draws a shaded rectangle. |

## Value

A multi-panel figure file. Each track shares the same x-axis (genomic position).
Selection statistic tracks use scatter plots with points colored by significance
(grey below threshold, colored above). The gene track shows rectangles for gene
bodies with arrows indicating strand direction.

## Details

### Data sources

- **Fst track**: reads `fst_wc_{pop1}_{pop2}` from GenomicWindow nodes
  overlapping the region, or per-variant Fst if window data is unavailable.
- **iHS track**: reads `ihs_{pop}` from Variant nodes. Values are plotted as
  |iHS|.
- **Gene track**: reads Gene nodes with start/end coordinates overlapping the
  region. Exon positions are read from Gene node properties when available.

All data must be persisted before plotting. If a requested track has no data,
it is omitted with a warning.

### Figure style

Follows Nature Methods guidelines: Arial font, white background, colorblind-safe
palette, removed top/right spines, outward ticks. Points exceeding the threshold
are highlighted in vermillion. The gene track uses blue rectangles for gene
bodies. See `graphpop plot` for full style details.

### Performance

Queries are limited to the plotted region plus flanks. Typical rendering time
is 1-3 seconds for a 1 Mb region.

## Examples

```bash
# Zoom into LCT gene with Fst and iHS tracks
graphpop plot gene-zoom LCT --pop EUR --pop2 AFR \
    -o fig_lct_zoom.png

# Rice gene with custom flank and threshold
graphpop plot gene-zoom Os01g0100100 --pop GJ-tmp --pop2 XI-1A \
    --flank 200000 --threshold 2.5 -o rice_gene_zoom.pdf

# Genomic interval with additional XP-EHH track
graphpop plot gene-zoom chr22:20000000-21000000 --pop EUR --pop2 EAS \
    --tracks fst,ihs,xpehh,genes -o region_zoom.png

# Highlight a specific sub-region
graphpop plot gene-zoom EDAR --pop EAS --pop2 AFR \
    --highlight 109510000-109530000 -o edar_zoom.png

# PDF for publication, single-column width
graphpop plot gene-zoom SLC24A5 --pop EUR --pop2 AFR \
    --width 3.5 --height 4.0 -o slc24a5_zoom.pdf
```

## See Also

- `graphpop plot chromosome` -- full chromosome multi-track view
- `graphpop plot manhattan` -- genome-wide Manhattan plot
- `graphpop lookup gene` -- text-based gene annotation lookup
- `graphpop plot` -- overview of all plot subcommands and shared style options
