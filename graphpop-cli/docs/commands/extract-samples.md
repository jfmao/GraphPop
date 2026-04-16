# Export Sample nodes and population metadata for a population

## Description

`graphpop extract samples` returns all Sample nodes belonging to a population,
together with per-sample metadata and, when available, population-level summary
statistics joined from the Population node.

Output includes the sample identifier, population label, the sample's packed
index (its position in the bit-packed genotype arrays), and the FROH value if
runs-of-homozygosity analysis has been run. Population-level fields
(`pop_n_samples`, `pop_mean_froh`) are included when the corresponding
properties are present on the Population node; they are empty otherwise.

This command is useful for generating sample sheets for external tools
(PLINK, ADMIXTURE, R), for verifying that a population loaded correctly, and for
cross-referencing sample identifiers between GraphPop and upstream VCF files.

## Usage

```
graphpop extract samples --pop POPULATION [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop` | string | *required* | Population name as stored in the graph (e.g., `EUR`, `GJ-tmp`, `CEU`). |
| `-o`, `--output` | path | stdout | Output file path. |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json`. |

## Value

One row per Sample node in the population:

| Column | Type | Description |
|--------|------|-------------|
| `sampleId` | string | Sample identifier (matches VCF sample column headers). |
| `population` | string | Population label as stored in the graph. |
| `packed_index` | integer | Zero-based index of the sample in the packed genotype bit arrays. Used internally by CARRIES edge traversal and bit-packing operations. |
| `froh` | float | Fraction of the autosomal genome covered by runs of homozygosity (if `graphpop roh` has been run; empty otherwise). |
| `pop_n_samples` | integer | Total number of samples in the population (from the Population node; empty if not yet computed). |
| `pop_mean_froh` | float | Mean FROH across all samples in the population (from the Population node; empty if not yet computed). |

## Details

### Population node join

The `pop_n_samples` and `pop_mean_froh` columns are sourced from the Population
node, not computed on the fly. They reflect whatever was written to the
Population node when `graphpop pop-summary` was last run. If `pop-summary` has
not been run, or if FROH has not been computed for the population, these columns
will be empty in the output.

### packed_index

The packed index is assigned at import time and determines a sample's bit
position in all packed genotype arrays. It is stable for the lifetime of a
database and can be used to address genotype data in `extract genotypes --format-gt raw`
output without going through CARRIES edge traversal.

### Performance

This command issues a single Cypher query over Sample nodes filtered by
population label. For populations of up to 10,000 samples, the query typically
completes in under 1 second. For larger populations, expect 1-3 seconds.

## Examples

```bash
# Export all EUR samples to a TSV file
graphpop extract samples --pop EUR -o eur_samples.tsv

# Export GJ-tmp samples (rice) as JSON
graphpop extract samples --pop GJ-tmp --format json -o gj_tmp_samples.json

# Quick check: count samples in a population
graphpop extract samples --pop AFR | wc -l

# Generate a PLINK keep-file (sample IDs, one per line)
graphpop extract samples --pop EAS \
    | awk 'NR>1 {print $1, $1}' > eas_keep.txt
```

## See Also

- `graphpop extract variants` -- export Variant nodes with annotation and AF filters
- `graphpop extract genotypes` -- export a sample x variant genotype matrix for a region
- `graphpop extract` -- overview of all extract subcommands
- `graphpop export-bed` -- export high-signal regions in BED format
- `graphpop pop-summary` -- compute and persist whole-population summary statistics
- `graphpop roh` -- compute runs of homozygosity and populate FROH on Sample nodes
