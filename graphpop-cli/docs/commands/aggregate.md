# graphpop aggregate

## Title

Generate Summary Tables from run-all Results

## Description

Reads the per-procedure TSV output files from a `graphpop run-all` execution
and produces publication-ready summary tables. Aggregates per-population
diversity, pairwise Fst, ROH statistics, selection scan peaks, and sweep
windows into consolidated TSV files.

## Usage

```
graphpop aggregate [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `-d`, `--results-dir` | PATH | *(required)* | Directory containing per-procedure TSV results (from `run-all`). Must exist. |
| `-j`, `--json-results` | PATH | *(none)* | JSON results file from `run-all` (optional; provides additional metadata). |
| `-o`, `--output-dir` | PATH | `graphpop_tables` | Output directory for the generated summary tables. |

## Value

Creates the following summary tables in the output directory:

| File | Content | Source |
|---|---|---|
| `population_summary.tsv` | Per-population diversity, theta, Tajima's D, F_IS, per chromosome | `diversity/` |
| `fst_matrix.tsv` | Pairwise Fst (Hudson and W&C), Dxy, Da for each population pair per chromosome | `divergence/` |
| `roh_summary.tsv` | Per-population ROH summary: mean FROH, mean number of ROH segments, max FROH | `roh/` |
| `ihs_peaks.tsv` | Top 100 iHS peaks across all populations and chromosomes | `ihs/` |
| `nsl_peaks.tsv` | Top 100 nSL peaks across all populations and chromosomes | `nsl/` |
| `xpehh_peaks.tsv` | Top 100 XP-EHH peaks across all pairs and chromosomes | `xpehh/` |
| `sweep_windows.tsv` | Garud's H12 windows exceeding threshold (0.1) | `garud_h/` |

## Details

### Table Generation

**population_summary.tsv** -- Reads each single-row TSV from `diversity/` and
extracts the population name and chromosome from the filename (e.g.,
`EUR_chr22.tsv` yields `population=EUR`, `chr=chr22`). Columns: `population`,
`chr`, `pi`, `theta_w`, `tajima_d`, `het_exp`, `het_obs`, `fis`, `n_variants`,
`n_segregating`.

**fst_matrix.tsv** -- Reads from `divergence/`. Filenames follow the pattern
`POP1_vs_POP2_CHR.tsv`. Columns: `pop1`, `pop2`, `chr`, `fst_hudson`,
`fst_wc`, `dxy`, `da`.

**roh_summary.tsv** -- Reads multi-row TSV files from `roh/` (one row per
sample). Computes per-population aggregates: `mean_froh`, `mean_n_roh`,
`max_froh`, `n_samples`.

**ihs_peaks.tsv / nsl_peaks.tsv / xpehh_peaks.tsv** -- Reads all per-variant
score files, computes absolute scores, sorts descending, and takes the top 100.

**sweep_windows.tsv** -- Reads Garud's H windows and filters for H12 >= 0.1,
sorted by H12 descending.

### Expected Input Structure

The results directory should have the structure created by `graphpop run-all`:

```
results/
  diversity/
    EUR_chr22.tsv
    AFR_chr22.tsv
    ...
  divergence/
    EUR_vs_AFR_chr22.tsv
    ...
  roh/
    EUR_chr22.tsv
    ...
  ihs/
    EUR_chr22.tsv
    ...
```

Only directories that exist are processed. Missing directories are silently
skipped.

## Examples

### Standard aggregation after run-all

```bash
graphpop run-all -d results/
graphpop aggregate -d results/ -o tables/
```

### With JSON results for additional metadata

```bash
graphpop aggregate -d results/ -j results/results.json -o tables/
```

### Aggregate and inspect

```bash
graphpop aggregate -d results/ -o tables/

# View population summary
column -t -s $'\t' tables/population_summary.tsv | head -20

# Count selection peaks
wc -l tables/ihs_peaks.tsv tables/nsl_peaks.tsv tables/xpehh_peaks.tsv
```

### Generate tables for a subset (re-run on filtered input)

```bash
# Copy only chr22 results
mkdir -p results_chr22/diversity results_chr22/divergence
cp results/diversity/*_chr22.tsv results_chr22/diversity/
cp results/divergence/*_chr22.tsv results_chr22/divergence/

graphpop aggregate -d results_chr22/ -o tables_chr22/
```

## See Also

- `graphpop run-all` -- Generate the input results for aggregation.
- `graphpop export-windows` -- Alternative export of GenomicWindow nodes.
- `graphpop filter` -- Query persisted results with annotation filters.
