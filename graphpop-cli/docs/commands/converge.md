# Find regions where multiple selection statistics converge

## Description

`graphpop converge` identifies genomic regions where two or more persisted
statistics simultaneously exceed user-defined thresholds. Selection leaves
signatures across multiple dimensions -- iHS captures ongoing sweeps, XP-EHH
captures completed sweeps, Fst captures population differentiation, H12 captures
haplotype homozygosity -- and the strongest candidates are those that show
extreme values in several statistics at once. This command performs a single
graph query that intersects threshold exceedances and returns annotated
convergence regions.

The query operates on persisted Variant-level or GenomicWindow-level statistics.
Statistics must have been previously computed and written to the graph (e.g., via
`graphpop ihs --persist`, `graphpop genome-scan`, `graphpop xpehh --persist`).
Results are joined to Gene annotations when available.

## Usage

```
graphpop converge --stats STAT1,STAT2[,...] --thresholds T1,T2[,...] [OPTIONS]
```

## Arguments

This command takes no positional arguments. All configuration is via options.

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--stats` | string | required | Comma-separated list of statistic names to intersect (e.g., `ihs,xpehh,h12,fst_wc`). |
| `--thresholds` | string | required | Comma-separated thresholds, one per statistic in the same order. For absolute-value statistics (iHS, nSL), the threshold applies to the absolute value. |
| `--pop` | string | required | Population for variant-level statistics. |
| `--pop2` | string | none | Second population for pairwise statistics (XP-EHH, Fst). Required if any pairwise statistic is included. |
| `--chr` | string | all | Restrict to a single chromosome. |
| `--level` | choice | `variant` | Query level: `variant` (per-variant scores) or `window` (GenomicWindow statistics). |
| `--merge-distance` | integer | 50000 | When `--level variant`, merge nearby hits within this distance (bp) into contiguous regions. |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file. |

## Value

Returns one row per convergence region (after merging nearby variants):

| Column | Type | Description |
|--------|------|-------------|
| `chr` | string | Chromosome. |
| `start` | integer | Region start position. |
| `end` | integer | Region end position. |
| `n_variants` | integer | Number of variants in the region exceeding all thresholds. |
| `max_ihs` | float | Maximum |iHS| in the region (if iHS was queried). |
| `max_xpehh` | float | Maximum |XP-EHH| in the region (if XP-EHH was queried). |
| `max_h12` | float | Maximum H12 in the region (if H12 was queried). |
| `max_fst` | float | Maximum Fst in the region (if Fst was queried). |
| `genes` | string | Comma-separated gene symbols overlapping the region. |
| `top_variant` | string | Variant with the highest composite rank across all queried statistics. |

## Details

### Threshold interpretation

Each statistic is compared against its threshold independently. A variant or
window passes the filter only if it exceeds the threshold for every listed
statistic simultaneously. For iHS and nSL, the comparison uses absolute values
(|iHS| >= threshold). For Fst and H12, raw values are used. For XP-EHH,
absolute values are used by default; use negative thresholds to capture
selection in pop2.

### Merge strategy

When `--level variant`, individual variants passing all thresholds are merged
into contiguous regions using `--merge-distance`. Two passing variants closer
than the merge distance are joined into a single region. The summary statistics
(max values) span the merged region.

### Gene annotation

If Gene nodes are present in the graph, the `genes` column lists all genes
whose boundaries overlap the convergence region. This join is performed via
positional overlap, not via IN_GENE edges, so it works even if variant-to-gene
edges are incomplete.

### Performance

The query issues one Cypher statement per chromosome with a WHERE clause
combining all thresholds. Performance scales with the number of variants that
have persisted scores. Typical runtime for a single chromosome with ~1M scored
variants is 2-5 seconds.

## Examples

```bash
# Find regions with strong iHS and Fst signals in Europeans
graphpop converge --stats ihs,fst_wc --thresholds 2.5,0.3 \
    --pop EUR --pop2 EAS --chr chr22 -o convergence_chr22.tsv

# Three-way convergence: iHS + XP-EHH + H12
graphpop converge --stats ihs,xpehh,h12 --thresholds 2.0,2.0,0.3 \
    --pop EUR --pop2 AFR -o triple_convergence.tsv

# Window-level convergence across the whole genome
graphpop converge --stats fst_wc,tajima_d --thresholds 0.4,-2.0 \
    --pop CEU --pop2 YRI --level window --format json

# Tight merge distance for fine-mapping
graphpop converge --stats ihs,nsl --thresholds 2.5,2.5 \
    --pop GJ-tmp --chr chr1 --merge-distance 10000

# Export convergence regions for a rice population
graphpop converge --stats ihs,xpehh,fst_wc --thresholds 2.0,2.0,0.2 \
    --pop GJ-tmp --pop2 XI-1A -o rice_convergence.tsv
```

## See Also

- `graphpop rank-genes` -- rank genes by composite selection score
- `graphpop filter` -- single-statistic filtering of persisted scores
- `graphpop export-bed` -- export high-signal regions as BED format
- `graphpop ihs` -- compute and persist iHS scores
- `graphpop xpehh` -- compute and persist XP-EHH scores
- `graphpop garud-h` -- compute and persist Garud's H statistics
