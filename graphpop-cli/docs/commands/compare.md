# Compare statistics between two populations

## Description

`graphpop compare` computes the per-window difference in a statistic between two
populations along a chromosome. It reads persisted GenomicWindow statistics and
returns, for each window, the value in each population and both the signed delta
(pop1 - pop2) and absolute delta. This is useful for identifying regions where
two populations diverge in diversity, haplotype structure, or other windowed
metrics.

Use `compare` after running `graphpop genome-scan` for both populations. The
command does not compute new statistics; it compares values already stored on
GenomicWindow nodes.

## Usage

```
graphpop compare POP1 POP2 CHR --stat STATISTIC [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `POP1` | string | yes | -- | First population identifier. |
| `POP2` | string | yes | -- | Second population identifier. |
| `CHR` | string | yes | -- | Chromosome to compare along. |
| `--stat` | string | yes | -- | Statistic to compare. Must be a property persisted on GenomicWindow nodes (e.g., `pi`, `theta_w`, `tajima_d`, `het_exp`, `fis`). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--sort-by` | choice | `position` | Sort output by: `position`, `delta`, `abs_delta`. |
| `--top` | integer | all | Return only the top N rows (by `--sort-by` column). |
| `--min-abs-delta` | float | none | Filter to windows where |delta| exceeds this value. |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Output file path. |

## Value

Returns one row per genomic window:

| Column | Type | Description |
|--------|------|-------------|
| `chr` | string | Chromosome. |
| `start` | integer | Window start position. |
| `end` | integer | Window end position. |
| `{stat}_{pop1}` | float | Statistic value in population 1. |
| `{stat}_{pop2}` | float | Statistic value in population 2. |
| `delta` | float | Signed difference: pop1 value minus pop2 value. |
| `abs_delta` | float | Absolute value of delta. |
| `genes` | string | Genes overlapping the window (if Gene nodes are present). |

## Details

### Window alignment

Both populations must have been scanned with the same window size and step size.
The command matches windows by chromosome and start position. Windows present in
one population but missing in the other are reported with `NA` for the missing
value and excluded from delta computation.

### Statistic resolution

The `--stat` option refers to the property name prefix on GenomicWindow nodes.
For example, `--stat pi` looks for `pi_{pop}` properties. Population-specific
suffixes are appended automatically based on POP1 and POP2 arguments.

### Performance

The query reads GenomicWindow nodes for a single chromosome and performs an
in-memory join. Typical runtime is under 2 seconds for chromosomes with up to
10,000 windows.

## Examples

```bash
# Compare nucleotide diversity between EUR and AFR on chr22
graphpop compare EUR AFR chr22 --stat pi -o pi_compare_chr22.tsv

# Find windows with largest diversity difference, top 20
graphpop compare CEU YRI chr22 --stat pi --sort-by abs_delta --top 20

# Compare Tajima's D, filtering for large absolute differences
graphpop compare GJ-tmp XI-1A chr1 --stat tajima_d \
    --min-abs-delta 1.0 -o tajima_diff_rice.tsv

# JSON output for downstream analysis
graphpop compare EUR EAS chr22 --stat theta_w --format json \
    -o theta_compare.json

# Compare inbreeding coefficients across a chromosome
graphpop compare GJ-tmp XI-1A chr1 --stat fis --sort-by abs_delta \
    -o fis_compare.tsv
```

## See Also

- `graphpop genome-scan` -- compute and persist windowed statistics
- `graphpop divergence` -- compute pairwise divergence (Fst, Dxy) between populations
- `graphpop export-windows` -- export raw GenomicWindow data
- `graphpop converge` -- find multi-statistic convergence regions
