# Extract variants, samples, or genotypes from the graph

## Description

`graphpop extract` exports data from the graph in tabular form. Three
subcommands cover the most common extraction needs:

| Subcommand | Purpose |
|------------|---------|
| `variants` | Export variant records with flexible filtering by chromosome, region, population, consequence, pathway, and allele frequency. Select which fields to include. |
| `samples` | Export per-population sample lists with metadata. |
| `genotypes` | Export a dosage or genotype matrix for a genomic region: rows are variants, columns are samples. |

These commands are read-only data exports. They do not compute new statistics;
they extract what is stored in the graph. Use `extract` to generate input for
external tools (R, Python, PLINK) or to inspect the data underlying a
computation.

## Usage

```
graphpop extract variants [FILTERS] [OPTIONS]
graphpop extract samples [OPTIONS]
graphpop extract genotypes CHR START END [OPTIONS]
```

## Arguments

### extract variants

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--chr` | string | no | all | Chromosome filter. |
| `--start` | integer | no | -- | Start position (requires `--chr`). |
| `--end` | integer | no | -- | End position (requires `--chr`). |

### extract genotypes

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome. |
| `START` | integer | yes | -- | Start position (1-based, inclusive). |
| `END` | integer | yes | -- | End position (1-based, inclusive). |

## Options

### Shared options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Output file path. |

### extract variants options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop` | string | none | Filter to variants segregating in this population. |
| `--consequence` | string | none | Filter by VEP consequence (e.g., `missense_variant`). |
| `--pathway` | string | none | Filter to variants in genes belonging to this pathway. |
| `--gene` | string | none | Filter to variants within a specific gene. |
| `--min-af` | float | none | Minimum allele frequency in `--pop`. |
| `--max-af` | float | none | Maximum allele frequency in `--pop`. |
| `--fields` | string | `variantId,pos,ref,alt,consequence,gene` | Comma-separated list of fields to include. Available: `variantId`, `pos`, `ref`, `alt`, `consequence`, `gene`, `af_{pop}`, `ac_{pop}`, `an_{pop}`, `ihs_{pop}`, `xpehh_{pop1}_{pop2}`. |

### extract samples options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop` | string | all | Export samples for a specific population only. |
| `--include-metadata` | flag | false | Include sample-level properties (sex, super_pop, etc.). |

### extract genotypes options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop` | string | all | Restrict to samples in this population. |
| `--mode` | choice | `dosage` | Output mode: `dosage` (0/1/2 matrix), `gt` (0/0, 0/1, 1/1 strings). |
| `--max-variants` | integer | 10000 | Safety limit on number of variants. Override with `0` for unlimited. |

## Value

### extract variants

One row per variant with the requested fields.

### extract samples

| Column | Type | Description |
|--------|------|-------------|
| `sampleId` | string | Sample identifier. |
| `population` | string | Population label. |
| `super_pop` | string | Super-population (if `--include-metadata`). |

### extract genotypes

A matrix with variant IDs as row labels, sample IDs as column headers, and
dosage values (0, 1, 2) or genotype strings as cell values. The first column
is `variantId`, followed by one column per sample.

## Details

### Variant field resolution

The `--fields` option accepts both static properties (pos, ref, alt) and
population-qualified statistics (af_EUR, ihs_CEU). Statistics must have been
persisted to appear in the output. Unknown fields produce an empty column with
a warning.

### Genotype reconstruction

The `extract genotypes` command reconstructs the genotype matrix from CARRIES
relationships. Samples without a CARRIES edge to a variant are assigned
dosage 0 (homozygous reference). The `--max-variants` safety limit prevents
accidental extraction of millions of variants; set to 0 for large exports.

### Performance

Variant extraction scales with the result set size: ~100k variants/second for
simple filters. Genotype extraction is slower due to CARRIES traversal:
expect 1-5 seconds per 1,000 variants depending on sample count.

## Examples

```bash
# Export all missense variants on chr22 for EUR
graphpop extract variants --chr chr22 --pop EUR \
    --consequence missense_variant --fields variantId,pos,af_EUR,gene \
    -o missense_chr22.tsv

# Export sample list for CEU population
graphpop extract samples --pop CEU -o ceu_samples.tsv

# Extract dosage matrix for a 100 kb region
graphpop extract genotypes chr22 20000000 20100000 --pop EUR \
    --mode dosage -o genotypes_region.tsv

# Extract genotype strings for a gene region
graphpop extract genotypes chr22 23522552 23838780 --pop YRI \
    --mode gt --format csv -o ewsr1_gt.csv

# Export rare variants in a pathway
graphpop extract variants --pathway "DNA repair" --max-af 0.01 \
    --fields variantId,pos,consequence,gene,af_EUR,af_AFR \
    -o rare_dna_repair.tsv
```

## See Also

- `graphpop filter` -- filter persisted statistics with threshold queries
- `graphpop export-bed` -- export regions as BED format
- `graphpop dump` -- full database export
- `graphpop query` -- arbitrary Cypher queries
- `graphpop export-windows` -- export GenomicWindow statistics
