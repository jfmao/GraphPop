# Export high-signal regions as BED format

## Description

`graphpop export-bed` exports genomic regions from the graph in BED format,
suitable for intersection with external tools (bedtools, IGV, UCSC Genome
Browser). Two modes are supported:

- **Window mode** (default): queries GenomicWindow nodes where a statistic
  exceeds a threshold, outputting one BED interval per window.
- **Variant mode** (`--level variant`): queries Variant nodes with persisted
  scores exceeding a threshold, then merges nearby variants into contiguous
  intervals using `--merge-distance`.

The BED output follows the standard 4-column format (chr, start, end, name),
where the name column contains the peak statistic value and gene annotation
when available. An optional 5-column format adds a score column for BED tools
that support it.

## Usage

```
graphpop export-bed --stat STATISTIC --threshold VALUE [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--stat` | string | required | Statistic to threshold. For variant-level: `ihs`, `xpehh`, `nsl`, `h12`. For window-level: `fst_wc`, `pi`, `tajima_d`, `dxy`. |
| `--threshold` | float | required | Minimum value for inclusion. For absolute-value statistics (iHS, nSL, XP-EHH), applies to the absolute value. |
| `--pop` | string | required | Population for within-population statistics. |
| `--pop2` | string | none | Second population for pairwise statistics. |
| `--chr` | string | all | Restrict to a single chromosome. |
| `--level` | choice | `window` | Query level: `window` (GenomicWindow nodes) or `variant` (Variant nodes). |
| `--merge-distance` | integer | 50000 | When `--level variant`, merge variants within this distance (bp) into a single interval. |
| `--bed-format` | choice | `bed4` | Output columns: `bed4` (chr, start, end, name) or `bed5` (adds integer score 0-1000). |
| `-o`, `--output` | path | stdout | Output file. |

## Value

Standard BED format output (0-based, half-open intervals):

```
chr22   17000000    17050000    fst=0.42;gene=CECR1
chr22   20150000    20200000    fst=0.38;gene=DGCR6
```

The name column encodes the peak statistic value and overlapping gene symbol(s).
When `--bed-format bed5`, a score column (0-1000) is added, computed as
min(1000, int(value / threshold * 500)).

## Details

### Window mode

Queries GenomicWindow nodes with the property `{stat}_{pop}` (or
`{stat}_{pop1}_{pop2}` for pairwise statistics) exceeding the threshold. Each
window becomes one BED interval. Coordinates are converted from 1-based
inclusive (GraphPop convention) to 0-based half-open (BED convention).

### Variant mode

Queries Variant nodes with persisted score properties. Passing variants are
sorted by position and merged into contiguous intervals using
`--merge-distance`. The name field reports the peak value within the merged
interval and any overlapping genes.

### Gene annotation

Gene symbols are added to the name field by positional overlap with Gene node
boundaries. If multiple genes overlap a single interval, they are
semicolon-separated.

### Performance

Window-mode export completes in 1-3 seconds for a full genome. Variant-mode
export depends on the number of scored variants and may take 5-15 seconds for
a genome-wide query with millions of scored variants.

## Examples

```bash
# Export windows with Fst > 0.3 as BED
graphpop export-bed --stat fst_wc --threshold 0.3 \
    --pop EUR --pop2 EAS -o high_fst_regions.bed

# Export iHS outlier regions (variant-level, merged)
graphpop export-bed --stat ihs --threshold 2.5 --pop EUR \
    --level variant --merge-distance 100000 -o ihs_peaks.bed

# BED5 format for UCSC track upload
graphpop export-bed --stat xpehh --threshold 3.0 \
    --pop EUR --pop2 AFR --bed-format bed5 -o xpehh_track.bed

# Single chromosome export
graphpop export-bed --stat tajima_d --threshold -2.0 \
    --pop YRI --chr chr22 -o tajima_d_sweeps.bed

# Rice: export differentiated regions between indica and japonica
graphpop export-bed --stat fst_wc --threshold 0.4 \
    --pop GJ-tmp --pop2 XI-1A --level window -o rice_diff_regions.bed
```

## See Also

- `graphpop converge` -- multi-statistic convergence query
- `graphpop filter` -- filter and display persisted statistics
- `graphpop export-windows` -- export all GenomicWindow data as TSV
- `graphpop genome-scan` -- compute and persist windowed statistics
