# Extract a sample x variant genotype matrix for a genomic region

## Description

`graphpop extract genotypes` exports genotype data for every sample in a
population across a specified genomic region. Three output modes are supported:

- **dosage** (default): integer allele dosage per sample per variant (0, 1, 2).
  Missing genotypes are reported as `-1`.
- **gt**: explicit genotype string per sample per variant (`0/0`, `0/1`, `1/1`).
  Missing genotypes are reported as `./.`.
- **raw**: one row per variant containing the hex-encoded packed genotype byte
  array (`gt_packed_hex`). This is the compact binary representation stored
  in the graph, suitable for downstream HPC pipelines that can decode the
  bit-packed format directly.

The `dosage` and `gt` modes reconstruct genotypes by traversing CARRIES edges
between Sample and Variant nodes. Samples without a CARRIES edge to a given
variant are homozygous reference (the implicit representation used throughout
GraphPop to avoid storing ~90 % of genotype calls explicitly).

The `raw` mode bypasses per-sample CARRIES traversal entirely and returns the
packed bit array from each Variant node. This is the fastest extraction path
and produces the most compact output. Use it when you need to feed data to a
tool that understands the GraphPop packed format or when you want to avoid the
overhead of reconstructing individual genotype strings.

## Usage

```
graphpop extract genotypes --chr CHR --start START --end END --pop POPULATION [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--chr` | string | *required* | Chromosome (e.g., `chr22`, `chr1`). |
| `--start` | integer | *required* | Start position of the region (1-based, inclusive). |
| `--end` | integer | *required* | End position of the region (1-based, inclusive). |
| `--pop` | string | *required* | Population name as stored in the graph. |
| `--format-gt` | choice | `dosage` | Genotype output mode: `dosage` (0/1/2), `gt` (0/0, 0/1, 1/1), or `raw` (hex gt_packed per variant). |
| `--limit` | integer | `1000` | Maximum number of variants to include. |
| `-o`, `--output` | path | stdout | Output file path. |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json`. |

## Value

### dosage and gt modes

One row per (sample, variant) pair:

| Column | Type | Description |
|--------|------|-------------|
| `sampleId` | string | Sample identifier. |
| `variantId` | string | Variant identifier (e.g., `chr22:16051249:A:G`). |
| `pos` | integer | Variant position (1-based). |
| `genotype` | integer or string | Dosage (0/1/2/-1) in `dosage` mode; genotype string (`0/0`, `0/1`, `1/1`, `./.`) in `gt` mode. |

### raw mode

One row per variant:

| Column | Type | Description |
|--------|------|-------------|
| `variantId` | string | Variant identifier. |
| `pos` | integer | Variant position (1-based). |
| `ref` | string | Reference allele. |
| `alt` | string | Alternate allele. |
| `gt_packed_hex` | string | Hex-encoded packed genotype byte array for all samples in the population. Bit pairs encode genotype: `00` = hom-ref, `01` = het, `10` = hom-alt, `11` = missing. |
| `af` | float | Population allele frequency. |

## Details

### Implicit homozygous reference

CARRIES edges exist only for non-reference genotypes (heterozygous or
homozygous alternate). A sample absent from the CARRIES neighborhood of a
Variant is therefore homozygous reference. The `dosage` and `gt` modes
reconstruct the full N×V matrix, filling in reference calls for all
sample-variant combinations where no CARRIES edge exists.

### Packed bit array format

In `raw` mode, the `gt_packed_hex` string encodes all samples in the population
in the order given by their `packed_index` values. Each sample occupies two
bits: `00` = homozygous reference, `01` = heterozygous, `10` = homozygous
alternate, `11` = missing. The array length is `ceil(N * 2 / 8)` bytes, where N
is the number of samples. Trailing bits in the last byte are zero-padded. The
`packed_index` values for a population can be retrieved with
`graphpop extract samples`.

### Limit behavior

The default `--limit 1000` restricts the number of variants returned (not the
number of rows). In `dosage` or `gt` mode, the actual row count is
`n_samples × n_variants`. For a population of 2,504 samples and 1,000 variants,
this produces up to 2,504,000 rows. Use `--limit` to stay within memory
constraints, or use `raw` mode for compact bulk extraction.

### Performance

`dosage` and `gt` modes traverse CARRIES edges for up to `limit` variants and
all samples in the population. Performance depends on the number of CARRIES
edges in the region (proportional to the mean allele frequency and the number of
samples). For a 1 Mb region at mean AF ~0.10 with 2,504 samples, expect 5-15
seconds. `raw` mode reads packed arrays directly from Variant node properties
and is typically 5-10x faster.

## Examples

```bash
# Dosage matrix for a 1 Mb region on chr22, EUR population
graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 \
    --pop EUR -o geno_dosage.tsv

# Explicit genotype strings for the same region
graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 \
    --pop EUR --format-gt gt -o geno_gt.tsv

# Raw packed genotypes for HPC downstream processing
graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 \
    --pop EUR --format-gt raw --limit 5000 -o geno_raw.tsv

# Rice: dosage matrix for a chr1 region in GJ-tmp, JSON output
graphpop extract genotypes --chr chr1 --start 1000000 --end 2000000 \
    --pop GJ-tmp --format-gt dosage --format json -o rice_geno.json

# Pipe dosage output directly into Python for PCA
graphpop extract genotypes --chr chr22 --start 16000000 --end 17000000 \
    --pop EUR --format-gt dosage --format csv \
    | python scripts/run_pca.py
```

## See Also

- `graphpop extract variants` -- export Variant nodes with annotation and AF filters
- `graphpop extract samples` -- export Sample nodes and population metadata
- `graphpop extract` -- overview of all extract subcommands
- `graphpop export-bed` -- export high-signal regions in BED format
- `graphpop ld` -- compute pairwise linkage disequilibrium in a region
