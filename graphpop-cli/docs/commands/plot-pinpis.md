# Plot piN/piS ratio bar chart

## Description

`graphpop plot pinpis` produces a horizontal bar chart of the ratio of
nonsynonymous to synonymous nucleotide diversity (piN/piS) across populations. A
dashed reference line at piN/piS = 1.0 marks the neutral expectation: ratios
below 1.0 indicate purifying selection removing deleterious nonsynonymous
variation, while ratios above 1.0 suggest relaxed constraint or positive
selection elevating nonsynonymous diversity.

This visualization is central to the "cost of domestication" hypothesis in crop
genetics, where domesticated populations often show elevated piN/piS relative to
wild relatives due to reduced effective population size and relaxed purifying
selection.

This is a standalone reference for the `pinpis` subcommand. For shared plot
style details (font, palette, DPI, figure guidelines), see
[graphpop plot](plot.md).

## Usage

```
graphpop plot pinpis INPUT_FILE -o OUTPUT [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `INPUT_FILE` | path | yes | -- | TSV or CSV file with `population` and `piN_piS` columns. |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | path | required | Output figure file (PNG, PDF, SVG). |
| `--title` | string | auto | Figure title. |
| `--width` | float | 7.2 | Figure width in inches. |
| `--height` | float | 3.5 | Figure height in inches. |

## Value

A horizontal bar chart with population names on the y-axis and piN/piS ratio on
the x-axis. A vertical dashed line marks piN/piS = 1.0. Bars are colored blue
(sky blue) when the ratio is at or below 1.0 and vermillion when above 1.0.
Populations are sorted by ratio (lowest at top). Numeric values are annotated on
each bar.

## Details

### Input format

The input file must contain at minimum two columns:
- `population`: population identifier string.
- `piN_piS`: numeric piN/piS ratio.

Additional columns are ignored. The file is auto-detected as TSV or CSV based
on the delimiter.

### Color logic

The two-color scheme (blue for <= 1.0, vermillion for > 1.0) provides an
immediate visual signal for populations under strong purifying selection versus
those with relaxed constraint. Colors are drawn from the Wong/Okabe-Ito
colorblind-safe palette.

For full style details, see [graphpop plot](plot.md).

## Examples

```bash
# Basic piN/piS bar chart
graphpop plot pinpis pinpis_ratios.tsv -o fig_pinpis.png

# Custom title for rice dataset
graphpop plot pinpis rice_pinpis.tsv -o rice_pinpis.pdf \
    --title "Cost of Domestication: piN/piS by Subpopulation"

# Compact figure for supplementary material
graphpop plot pinpis pinpis_ratios.tsv -o sfig_pinpis.png \
    --width 5 --height 3

# Publication PDF at single-column width
graphpop plot pinpis pinpis_ratios.tsv -o fig3b_pinpis.pdf \
    --width 3.5 --height 3.0
```

## See Also

- `graphpop plot` -- full plot documentation with all subcommands and style details
- `graphpop diversity` -- compute nucleotide diversity (pi)
- `graphpop plot diversity-bar` -- diversity ranking bar chart
