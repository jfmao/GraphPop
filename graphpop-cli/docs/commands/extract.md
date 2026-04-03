# Extract variants, samples, or genotypes from the graph

## Description

`graphpop extract` exports data from the graph in tabular form. Three
subcommands cover the most common extraction needs:

| Subcommand | Purpose |
|------------|---------|
| `variants` | Export variant records with flexible filtering by chromosome, region, population, consequence, pathway, gene, and allele frequency. Select which fields to include. |
| `samples` | Export per-population sample lists with metadata. |
| `genotypes` | Export a dosage or genotype matrix for a genomic region. |

These commands are read-only data exports. They do not compute new statistics;
they extract what is stored in the graph. Use `extract` to generate input for
external tools (R, Python, PLINK) or to inspect the data underlying a
computation.

## Usage

```
graphpop extract variants [OPTIONS]
graphpop extract samples --pop POPULATION [OPTIONS]
graphpop extract genotypes --chr CHR --start START --end END --pop POPULATION [OPTIONS]
```

## Options

### Shared options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Output file path. |

### extract variants options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--chr` | string | none | Chromosome filter. |
| `--start` | integer | none | Start position filter. |
| `--end` | integer | none | End position filter. |
| `--pop` | string | none | Population name (for AF lookup and filtering). |
| `--min-af` | float | none | Minimum allele frequency in `--pop`. |
| `--max-af` | float | none | Maximum allele frequency in `--pop`. |
| `--consequence` | string | none | Filter by VEP consequence (e.g., `missense_variant`). |
| `--pathway` | string | none | Filter to variants in genes belonging to this pathway (substring match). |
| `--gene` | string | none | Filter to variants within a specific gene (symbol or ID). |
| `--fields` | string | `variantId,pos,ref,alt,af` | Comma-separated list of variant properties to return. |
| `--limit` | integer | 10000 | Maximum rows to return. |

### extract samples options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop` | string | *required* | Population name. |

### extract genotypes options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--chr` | string | *required* | Chromosome. |
| `--start` | integer | *required* | Start position. |
| `--end` | integer | *required* | End position. |
| `--pop` | string | *required* | Population name. |
| `--format-gt` | choice | `dosage` | Output mode: `dosage` (0/1/2), `gt` (0/0, 0/1, 1/1), `raw` (hex gt_packed). |
| `--limit` | integer | 1000 | Maximum number of variants. |

## Value

### extract variants

One row per variant with the requested `--fields`. When `--pop` is specified
and `af` is in the field list, the column is named `af_{population}`.

### extract samples

| Column | Type | Description |
|--------|------|-------------|
| `sampleId` | string | Sample identifier. |
| `population` | string | Population label. |
| `packed_index` | integer | Sample index in packed genotype arrays. |
| `froh` | float | Fraction of genome in runs of homozygosity (if computed). |
| `pop_n_samples` | integer | Total samples in the population. |
| `pop_mean_froh` | float | Population mean FROH (if computed). |

### extract genotypes

In `dosage` or `gt` mode: one row per sample-variant pair with columns
`sampleId`, `variantId`, `pos`, `genotype`.

In `raw` mode: one row per variant with columns `variantId`, `pos`, `ref`,
`alt`, `gt_packed_hex`, `af`.

## Details

### Variant field resolution

The `--fields` option accepts variant node properties (pos, ref, alt, etc.).
When `--pop` is specified, `af` is resolved to the population-specific allele
frequency from the `af[]` array. Unknown property names produce `null` values.

### Genotype reconstruction

The `dosage` and `gt` modes query individual CARRIES edges between Sample and
Variant nodes. Samples without a CARRIES edge to a variant are homozygous
reference (implicit in the graph schema). The `raw` mode returns the packed
genotype byte array as a hex string per variant.

## Examples

```bash
# Export all missense variants on chr22 for EUR
graphpop extract variants --chr chr22 --pop EUR \
    --consequence missense_variant -o missense_chr22.tsv

# Export variants in a gene
graphpop extract variants --gene KCNE1 --pop EUR -o kcne1_variants.tsv

# Export sample list for a population
graphpop extract samples --pop EUR -o eur_samples.tsv

# Extract dosage matrix for a region
graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 \
    --pop EUR -o geno.tsv

# Extract raw packed genotypes
graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 \
    --pop EUR --format-gt raw -o geno_raw.tsv
```

## See Also

- `graphpop filter` -- filter persisted statistics with threshold queries
- `graphpop export-bed` -- export regions as BED format
- `graphpop dump` -- full database export
- `graphpop query` -- arbitrary Cypher queries
- `graphpop export-windows` -- export GenomicWindow statistics
