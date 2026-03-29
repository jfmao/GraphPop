# graphpop export-windows

## Title

Batch Export GenomicWindow Nodes with Filters

## Description

Queries persisted `GenomicWindow` nodes from the graph database and exports them
as TSV, CSV, or JSON. Supports filtering by chromosome, population, and
statistic thresholds (Fst, pi, Tajima's D). This is the primary way to extract
genome-scan results for downstream analysis or visualization.

## Usage

```
graphpop export-windows [CHR] [POPULATION] [OPTIONS]
```

## Arguments

| Argument | Type | Required | Description |
|---|---|---|---|
| `CHR` | TEXT | No | Chromosome filter (e.g., `chr1`, `chr22`). If omitted, all chromosomes are returned. |
| `POPULATION` | TEXT | No | Population filter (e.g., `EUR`, `GJ-tmp`). If omitted, all populations are returned. |

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `-o`, `--output` | PATH | stdout | Output file path. |
| `--format` | CHOICE | `tsv` | Output format: `tsv`, `csv`, or `json`. |
| `--min-pi` | FLOAT | *(none)* | Minimum pi (nucleotide diversity) filter. |
| `--max-pi` | FLOAT | *(none)* | Maximum pi filter. |
| `--min-fst` | FLOAT | *(none)* | Minimum Fst filter. |
| `--min-tajima-d` | FLOAT | *(none)* | Minimum Tajima's D filter. |
| `--max-tajima-d` | FLOAT | *(none)* | Maximum Tajima's D filter. |
| `--run-id` | TEXT | *(none)* | Filter by a specific run ID. |
| `--limit` | INT | *(none)* | Maximum number of windows to return. |

## Value

Returns a table with the following columns:

| Column | Description |
|---|---|
| `window_id` | Window identifier (unique) |
| `chr` | Chromosome |
| `start` | Window start position (bp) |
| `end` | Window end position (bp) |
| `population` | Population ID |
| `run_id` | Run identifier |
| `n_variants` | Total variants in window |
| `n_segregating` | Segregating sites in window |
| `pi` | Nucleotide diversity |
| `theta_w` | Watterson's theta |
| `tajima_d` | Tajima's D |
| `fst` | Hudson's Fst |
| `fst_wc` | Weir & Cockerham Fst |
| `dxy` | Absolute divergence |
| `pbs` | Population branch statistic |
| `fay_wu_h` | Fay & Wu's H |

Columns with no data for a given window appear as `null` (JSON) or empty (TSV).
Results are ordered by chromosome and start position.

## Details

GenomicWindow nodes are created by `graphpop genome-scan` with `--persist`
(the default). Each node stores the summary statistics computed for that
genomic window and population.

All filters are combined with AND logic. For example, `--min-fst 0.5
--max-tajima-d -2` returns only windows with Fst >= 0.5 AND Tajima's D <= -2.

When neither `CHR` nor `POPULATION` is specified and no filters are given, the
command exports all GenomicWindow nodes in the database. This can be a large
result set; use `--limit` to cap the output.

### Use Cases

- **Outlier detection** -- Export high-Fst or low-Tajima's D windows as
  candidates for selection.
- **Visualization** -- Export to TSV for plotting in R or Python (Manhattan
  plots, genome scans).
- **Cross-population comparison** -- Export windows for multiple populations
  and merge externally.
- **Quality control** -- Check that genome-scan results are present and
  reasonable.

## Examples

### Export all windows for one population

```bash
graphpop export-windows chr22 EUR -o windows_eur_chr22.tsv
```

### Export high-differentiation outlier windows

```bash
graphpop export-windows --min-fst 0.5 -o high_fst_outliers.tsv
```

### Export windows with negative Tajima's D (possible sweeps)

```bash
graphpop export-windows chr1 AFR --max-tajima-d -2.0 -o sweep_candidates.tsv
```

### Export all windows as JSON

```bash
graphpop export-windows --format json -o all_windows.json
```

### Export with multiple filters

```bash
graphpop export-windows \
    --min-fst 0.3 \
    --max-tajima-d -1.5 \
    --min-pi 0.001 \
    --limit 500 \
    -o selection_hotspots.tsv
```

### Pipe to analysis tools

```bash
# Count windows per population
graphpop export-windows --min-fst 0.5 | awk -F'\t' 'NR>1 {print $5}' | sort | uniq -c

# Quick look at top differentiated windows
graphpop export-windows --min-fst 0.5 --limit 20
```

## See Also

- `graphpop genome-scan` -- Create GenomicWindow nodes (must run first).
- `graphpop filter` -- Filter per-variant statistics by annotation.
- `graphpop aggregate` -- Generate summary tables from run-all output.
- `graphpop query` -- Run arbitrary Cypher for custom window queries.
