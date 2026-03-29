# Sliding-window genome scan

## Description

`graphpop genome-scan` performs a sliding-window scan across a chromosome,
computing diversity and divergence statistics for each window. It returns one
row per window, enabling visualization of how summary statistics vary across the
genome. This is the primary tool for identifying genomic regions with unusual
patterns of variation -- selective sweeps, balancing selection, introgression
tracts, or regions of low recombination.

The scan always **writes GenomicWindow nodes to the graph**, creating permanent,
queryable records of windowed statistics. These nodes can later be queried
directly in Cypher, exported with `graphpop export-windows`, or used by
downstream analyses. The `--persist` flag is accepted for explicitness but
persistence is the default behavior.

This command operates on the **FAST PATH**, computing all statistics from allele
count arrays on Variant nodes. When `--pop2` is provided, divergence statistics
(F_ST, D_xy, PBS) are computed alongside diversity statistics in the same pass.

Use `genome-scan` for genome-wide visualization of variation, for identifying
outlier windows under selection, or for generating windowed statistics that
feed into downstream analyses (e.g., correlation of F_ST with recombination
rate, or overlap with annotation features).

## Usage

```
graphpop genome-scan CHR POPULATION WINDOW_SIZE STEP_SIZE [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome name (e.g., `chr22`, `chr1`). The scan covers the entire chromosome. |
| `POPULATION` | string | yes | -- | Primary population identifier. Diversity statistics are computed for this population. |
| `WINDOW_SIZE` | integer | yes | -- | Window size in base pairs (e.g., `100000` for 100 kb windows) |
| `STEP_SIZE` | integer | yes | -- | Step size in base pairs (e.g., `50000` for 50% overlap). Set equal to WINDOW_SIZE for non-overlapping windows. |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop2` | string | none | Second population for computing divergence statistics (F_ST, D_xy, PBS) alongside diversity. When omitted, divergence columns are null. |
| `--persist` | flag | false | Explicitly enable persistence of GenomicWindow nodes (this is the default behavior; the flag is accepted for clarity) |
| `--consequence` | string | none | Filter variants by VEP consequence type |
| `--pathway` | string | none | Filter variants to those in genes belonging to the named pathway |
| `--gene` | string | none | Filter variants to those within the named gene |
| `--min-af` | float | none | Exclude variants with allele frequency below this threshold |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json` |

## Value

Returns one row per window with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `window_id` | string | Unique identifier for the window (format: `{chr}_{start}_{end}_{population}`) |
| `chr` | string | Chromosome name |
| `start` | integer | Window start position (1-based, inclusive) |
| `end` | integer | Window end position (1-based, inclusive) |
| `population` | string | Primary population identifier |
| `n_variants` | integer | Number of variant sites in the window (after any filtering) |
| `n_segregating` | integer | Number of segregating sites in the primary population |
| `pi` | float | Nucleotide diversity |
| `theta_w` | float | Watterson's theta |
| `tajima_d` | float | Tajima's D |
| `fay_wu_h` | float | Fay and Wu's H (requires ancestral allele annotation) |
| `fst` | float or null | Hudson's F_ST between POPULATION and POP2. Null if `--pop2` is not provided. |
| `fst_wc` | float or null | Weir and Cockerham's F_ST. Null if `--pop2` is not provided. |
| `dxy` | float or null | Absolute divergence (D_xy). Null if `--pop2` is not provided. |
| `pbs` | float or null | Population branch statistic. Null unless a three-population design is used via the underlying procedure. |

## Details

### Algorithm

The procedure iterates across the chromosome in steps of STEP_SIZE, collecting
all Variant nodes within each [start, start + WINDOW_SIZE) interval. For each
window:

1. Allele count arrays are read from all Variant nodes in the window.
2. Diversity statistics (pi, theta_w, Tajima's D, Fay and Wu's H) are computed
   for the primary population.
3. If `--pop2` is specified, divergence statistics (F_ST, D_xy) are computed
   for the population pair.
4. A GenomicWindow node is created in the graph with all computed properties.

Windows at the end of the chromosome that extend beyond the last variant are
included if they contain at least one variant.

### Window sizing guidance

Common configurations:

- **100 kb windows, 50 kb step**: Standard for genome-wide scans. Provides
  good resolution while averaging over enough variants for stable estimates.
- **50 kb windows, 25 kb step**: Higher resolution, useful for fine-mapping
  selection signals.
- **10 kb windows, 5 kb step**: Gene-scale resolution, but statistics become
  noisy when variant density is low.
- **1 Mb windows, 500 kb step**: Broad patterns, useful for visualizing
  chromosome-wide trends.

Set STEP_SIZE equal to WINDOW_SIZE for non-overlapping windows (faster, fewer
output rows). Use STEP_SIZE = WINDOW_SIZE / 2 for standard 50% overlap.

### Conditioning behavior

When `--consequence`, `--pathway`, or `--gene` is specified, only matching
variants are included in each window. This produces an annotation-conditioned
genome scan. For example, `--consequence missense_variant` creates a scan based
only on nonsynonymous variants, which can reveal windows where functional
variation is unusually depleted or enriched.

Note that conditioning reduces the number of variants per window, which
increases variance in the estimated statistics. Use larger windows when
conditioning on rare annotation categories.

### Persistence behavior

GenomicWindow nodes are **always written to the graph** as a side effect of
this command. Each node stores the chromosome, start, end, population, and all
computed statistics as properties. These nodes are permanent and can be queried
with Cypher:

```cypher
MATCH (w:GenomicWindow)
WHERE w.population = 'EUR' AND w.fst > 0.5
RETURN w.chr, w.start, w.end, w.fst, w.pi
ORDER BY w.fst DESC
```

Use `graphpop export-windows` to batch-export GenomicWindow nodes with filters.
Running `genome-scan` again with the same parameters will create new window
nodes; manage duplicates through database queries or by clearing previous runs.

### Performance

FAST PATH: O(W x V_w x K) where W is the number of windows and V_w is the
average number of variants per window. A genome scan of chromosome 22 (~50 Mb)
with 100 kb windows and 50 kb step produces ~1000 windows and typically
completes in 10-30 seconds. The persistence step (writing GenomicWindow nodes)
adds a small overhead per window.

Adding `--pop2` for divergence statistics increases computation by approximately
50% per window but is still done in a single pass over the data.

## Examples

```bash
# Basic: 100 kb non-overlapping windows for EUR on chromosome 22
graphpop genome-scan chr22 EUR 100000 100000 -o scan_eur.tsv

# Overlapping windows with divergence statistics
graphpop genome-scan chr22 EUR 100000 50000 \
    --pop2 EAS \
    -o scan_eur_eas.tsv

# Conditioned on missense variants with 200 kb windows
graphpop genome-scan chr22 CEU 200000 100000 \
    --consequence missense_variant \
    -o scan_missense.tsv

# Pathway-conditioned scan for rice subpopulation
graphpop genome-scan chr1 GJ-tmp 100000 50000 \
    --pathway "Starch and sucrose metabolism" \
    --pop2 XI-1A \
    --format json \
    -o scan_starch.json

# Fine-resolution scan for a specific region (use small windows)
graphpop genome-scan chr22 YRI 10000 5000 -o scan_yri_fine.tsv
```

## See Also

- `graphpop diversity` -- single-region diversity (when you need one region, not a scan)
- `graphpop divergence` -- single-region divergence statistics
- `graphpop export-windows` -- batch export of GenomicWindow nodes with filters
- `graphpop pop-summary` -- whole-chromosome summary (complements per-window scan)
- `graphpop filter` -- query persisted FULL PATH statistics by annotation
