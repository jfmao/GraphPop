# Compute nucleotide diversity and related summary statistics

## Description

`graphpop diversity` computes a suite of within-population summary statistics
for a genomic region: nucleotide diversity (pi), Watterson's theta, Tajima's D,
Fay and Wu's H (raw and normalized), expected and observed heterozygosity, and
the inbreeding coefficient F_IS. These are the core descriptors of genetic
variation within a single population.

This command operates on the **FAST PATH**, reading pre-computed allele count
arrays (ac[], an[], af[]) stored directly on Variant nodes. Computation scales
as O(V x K) where V is the number of variants in the region and K is the number
of populations, making it independent of sample size. A region with 100,000
variants returns in seconds regardless of whether the dataset contains 100 or
100,000 individuals.

Use `diversity` when you need a quick characterization of genetic variation in a
region, when comparing diversity across populations, or as input to downstream
analyses such as identifying regions of reduced diversity (selective sweeps) or
elevated diversity (balancing selection).

## Usage

```
graphpop diversity CHR START END POPULATION [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome name (e.g., `chr22`, `chr1`) |
| `START` | integer | yes | -- | Start position of the genomic region (1-based, inclusive) |
| `END` | integer | yes | -- | End position of the genomic region (1-based, inclusive) |
| `POPULATION` | string | yes | -- | Population identifier as stored in the graph (e.g., `EUR`, `GJ-tmp`, `CEU`) |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--consequence` | string | none | Filter variants by VEP consequence type (e.g., `missense_variant`, `synonymous_variant`, `stop_gained`) |
| `--pathway` | string | none | Filter variants to those in genes belonging to the named biological pathway |
| `--gene` | string | none | Filter variants to those within the named gene |
| `--min-af` | float | none | Exclude variants with allele frequency below this threshold (in the target population) |
| `--max-af` | float | none | Exclude variants with allele frequency above this threshold |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json` |

## Value

Returns a single row with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `pi` | float | Nucleotide diversity (average pairwise differences per site) |
| `theta_w` | float | Watterson's theta (expected diversity under neutral model, based on number of segregating sites) |
| `tajima_d` | float | Tajima's D statistic (deviation of pi from theta_w; negative values suggest purifying selection or population expansion, positive values suggest balancing selection or population contraction) |
| `fay_wu_h` | float | Fay and Wu's H (sensitive to high-frequency derived alleles; requires ancestral allele polarization; strongly negative values indicate positive selection) |
| `fay_wu_h_norm` | float | Normalized Fay and Wu's H (variance-adjusted, comparable across regions of different sizes) |
| `het_exp` | float | Expected heterozygosity (mean 2pq across sites) |
| `het_obs` | float | Observed heterozygosity (fraction of heterozygous genotypes, averaged across sites) |
| `fis` | float | Inbreeding coefficient (1 - H_obs / H_exp); positive values indicate excess homozygosity |
| `n_variants` | integer | Total number of variant sites in the region (after any filtering) |
| `n_segregating` | integer | Number of segregating (polymorphic) sites in the target population |
| `n_polarized` | integer | Number of variants with known ancestral allele (used for Fay and Wu's H) |

## Details

### Algorithm

All statistics are derived from allele count arrays stored on Variant nodes.
For a region with V variants and K populations tracked in the database:

- **pi** is computed as the sum of `2 * p * (1-p) * n/(n-1)` across all sites,
  divided by the number of callable sites, where p is the allele frequency and n
  is the allele number (number of observed chromosomes).
- **theta_w** uses the number of segregating sites S divided by the harmonic
  number a1 = sum(1/i for i = 1..n-1).
- **tajima_d** is the standardized difference (pi - theta_w) / sd, using
  Tajima's (1989) variance formula.
- **fay_wu_h** requires ancestral allele annotation. Variants without
  polarization information are excluded from H computation. The `n_polarized`
  column reports how many variants contributed.
- **het_obs** is computed from genotype counts on Variant nodes when available,
  or estimated from allele frequencies when individual genotype summaries are
  stored.

### Conditioning behavior

When `--consequence`, `--pathway`, or `--gene` is specified, the underlying
`VariantQuery` restricts which Variant nodes are traversed before computing
statistics. This is a true conditioned analysis: all statistics reflect only the
filtered variant set. For example, `--consequence missense_variant` computes
diversity using only nonsynonymous sites, which is useful for comparing
functional versus neutral diversity.

Multiple conditioning options can be combined. When both `--consequence` and
`--pathway` are specified, only variants matching both criteria are included.

The `--min-af` and `--max-af` options filter on population-specific allele
frequency and are applied before statistic computation.

### Performance

FAST PATH complexity: O(V x K) where V is the number of variants in the region
and K is the number of populations in the allele count arrays. Typical
performance for a 50 Mb chromosome region with ~1M variants is 1-3 seconds.
Conditioning by annotation may add index lookup time but typically remains under
5 seconds. No individual genotype traversal is required.

### Persistence

This command does not write results to the graph. Output is returned to stdout
or to the file specified with `-o`. To persist diversity statistics to the graph,
use `graphpop pop-summary` (whole chromosome) or `graphpop genome-scan` (per
window).

## Examples

```bash
# Basic: diversity statistics for EUR on chromosome 22
graphpop diversity chr22 1 50818468 EUR

# Conditioned on missense variants only, output to file
graphpop diversity chr22 1 50818468 EUR \
    --consequence missense_variant \
    -o diversity_missense.tsv

# Diversity for a specific gene in JSON format
graphpop diversity chr22 23522552 23838780 CEU \
    --gene EWSR1 \
    --format json \
    -o ewsr1_diversity.json

# Pathway-conditioned diversity with allele frequency filter
graphpop diversity chr1 1 43270923 GJ-tmp \
    --pathway "Starch and sucrose metabolism" \
    --min-af 0.01 \
    -o starch_diversity.tsv

# Compare common vs rare variant diversity (pipe to downstream tools)
graphpop diversity chr22 1 50818468 YRI --max-af 0.05 --format csv
graphpop diversity chr22 1 50818468 YRI --min-af 0.05 --format csv
```

## See Also

- `graphpop pop-summary` -- whole-chromosome diversity with persistence to Population nodes
- `graphpop genome-scan` -- sliding-window diversity scan with persistence to GenomicWindow nodes
- `graphpop divergence` -- between-population divergence (Fst, Dxy)
- `graphpop sfs` -- site frequency spectrum (the distribution underlying these summary statistics)
