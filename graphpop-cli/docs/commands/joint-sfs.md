# Compute the joint site frequency spectrum between two populations

## Description

`graphpop joint-sfs` computes the two-dimensional joint site frequency spectrum
(2D-SFS) between two populations over a genomic region. The joint SFS is a
matrix where entry (i, j) counts the number of variant sites with allele count i
in population 1 and allele count j in population 2. It captures the full joint
distribution of allele frequencies and is the sufficient statistic for many
demographic inference methods.

By default, the command returns the **folded** joint SFS (using minor allele
counts). When `--unfolded` is specified, it returns the **unfolded** joint SFS
based on derived allele counts, which requires ancestral allele annotations.

This command operates on the **FAST PATH** with O(V x K) complexity, reading
allele count arrays directly from Variant nodes.

Use `joint-sfs` when you need the full two-population frequency distribution
for demographic inference (input to dadi, moments, fastsimcoal2), when computing
custom divergence statistics, or when you need more information than the scalar
summaries (F_ST, D_xy) provided by `graphpop divergence`.

## Usage

```
graphpop joint-sfs CHR START END POP1 POP2 [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome name (e.g., `chr22`, `chr1`) |
| `START` | integer | yes | -- | Start position of the genomic region (1-based, inclusive) |
| `END` | integer | yes | -- | End position of the genomic region (1-based, inclusive) |
| `POP1` | string | yes | -- | First population identifier |
| `POP2` | string | yes | -- | Second population identifier |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--unfolded` | flag | false | Compute the unfolded (derived allele) joint SFS. Requires ancestral allele annotations. Variants without ancestral state are excluded. |
| `--consequence` | string | none | Filter variants by VEP consequence type |
| `--pathway` | string | none | Filter variants to those in genes belonging to the named pathway |
| `--gene` | string | none | Filter variants to those within the named gene |
| `--min-af` | float | none | Exclude variants with allele frequency below this threshold in either population |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json` |

## Value

Returns a single row with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `joint_sfs` | list of integers | The joint SFS as a flattened 1D array in row-major order. Element at logical position (i, j) is stored at index `i * dim2 + j`. Reshape to a `dim1 x dim2` matrix for analysis. |
| `dim1` | integer | Number of rows in the joint SFS matrix (max_ac1 + 1 for unfolded, floor(max_ac1 / 2) + 1 for folded) |
| `dim2` | integer | Number of columns in the joint SFS matrix (max_ac2 + 1 for unfolded, floor(max_ac2 / 2) + 1 for folded) |
| `n_variants` | integer | Total number of variant sites used in the computation |
| `max_ac1` | integer | Maximum allele count for POP1 (2 x sample_size_1, i.e., total chromosomes) |
| `max_ac2` | integer | Maximum allele count for POP2 (2 x sample_size_2, i.e., total chromosomes) |

## Details

### Algorithm

For each variant in the region, the allele counts in POP1 and POP2 are read
from the allele count arrays on the Variant node. The (ac1, ac2) pair
determines which cell of the joint SFS matrix is incremented:

- **Folded**: Minor allele counts are used for both populations. The matrix
  dimensions are (floor(n1/2) + 1) x (floor(n2/2) + 1).

- **Unfolded**: Derived allele counts are used, requiring ancestral allele
  polarization. The matrix dimensions are (n1 + 1) x (n2 + 1). Variants
  without ancestral annotation are excluded.

The output is a flattened array in row-major order. To reconstruct the matrix
in Python:

```python
import numpy as np
matrix = np.array(joint_sfs).reshape(dim1, dim2)
```

### Relationship to divergence statistics

The joint SFS contains all the information needed to compute F_ST, D_xy, and
other pairwise statistics. The `graphpop divergence` command computes scalar
summaries; `joint-sfs` provides the full distribution. If you need both, use
`joint-sfs` and derive the summaries, or run both commands (both are FAST PATH
and complete in seconds).

### Conditioning behavior

Annotation conditioning works identically to other FAST PATH commands. When
`--consequence` or `--pathway` is specified, only matching variants contribute
to the joint SFS. This is useful for comparing the joint frequency distribution
at functional versus neutral sites, or for generating annotation-stratified
inputs to demographic inference software.

### Performance

FAST PATH: O(V x K). The computation itself is fast, but the output size scales
with sample size: for two populations of 1000 diploid individuals each, the
unfolded joint SFS has 2001 x 2001 = ~4M entries. The folded spectrum is
roughly 4x smaller. For very large sample sizes, consider using `--format json`
for more compact output, or pipe to downstream tools directly.

### Persistence

This command does not write results to the graph.

## Examples

```bash
# Basic: folded joint SFS between CEU and YRI on chromosome 22
graphpop joint-sfs chr22 1 50818468 CEU YRI

# Unfolded joint SFS for demographic inference, output to file
graphpop joint-sfs chr22 1 50818468 CEU YRI \
    --unfolded \
    -o joint_sfs_ceu_yri.tsv

# Joint SFS conditioned on missense variants
graphpop joint-sfs chr22 1 50818468 EUR EAS \
    --consequence missense_variant \
    -o joint_sfs_missense.tsv

# Pathway-conditioned joint SFS in JSON format (for dadi/moments input)
graphpop joint-sfs chr1 1 43270923 GJ-tmp XI-1A \
    --pathway "Starch and sucrose metabolism" \
    --format json \
    -o joint_sfs_starch.json

# Common-variant joint SFS with minimum allele frequency filter
graphpop joint-sfs chr22 1 50818468 CEU CHB \
    --min-af 0.05 \
    -o joint_sfs_common.tsv
```

## See Also

- `graphpop sfs` -- single-population site frequency spectrum
- `graphpop divergence` -- scalar divergence statistics (F_ST, D_xy, PBS) derived from the joint SFS
- `graphpop diversity` -- within-population summary statistics
