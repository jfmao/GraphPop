# Compute population divergence statistics

## Description

`graphpop divergence` computes between-population divergence statistics for a
genomic region: Hudson's F_ST, Weir and Cockerham's F_ST, absolute divergence
(D_xy), net divergence (D_a), and optionally the population branch statistic
(PBS) when a third outgroup population is provided. These statistics quantify
how genetically differentiated two populations are and are fundamental to
detecting local adaptation, population structure, and barriers to gene flow.

This command operates on the **FAST PATH**, reading pre-computed allele count
arrays (ac[], an[], af[]) stored on Variant nodes. Computation scales as
O(V x K) and is independent of sample size. Divergence between two populations
across an entire chromosome typically completes in seconds.

Use `divergence` for pairwise population comparisons, to identify genomic
regions of elevated differentiation, or to compute PBS for detecting
population-specific selection with a three-population design.

## Usage

```
graphpop divergence CHR START END POP1 POP2 [OPTIONS]
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
| `--pop3` | string | none | Third (outgroup) population for PBS computation. When provided, PBS values are returned for POP1 relative to POP2 with POP3 as outgroup. |
| `--consequence` | string | none | Filter variants by VEP consequence type (e.g., `missense_variant`, `synonymous_variant`) |
| `--pathway` | string | none | Filter variants to those in genes belonging to the named biological pathway |
| `--gene` | string | none | Filter variants to those within the named gene |
| `--min-af` | float | none | Exclude variants with allele frequency below this threshold (in either population) |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json` |

## Value

Returns a single row with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `fst_hudson` | float | Hudson's F_ST estimator (ratio of averages). Ranges from 0 (no differentiation) to 1 (fixed differences). Robust for unequal sample sizes. |
| `fst_wc` | float | Weir and Cockerham's (1984) F_ST estimator. Accounts for sample size differences through variance component decomposition. The standard estimator for most population genetics applications. |
| `dxy` | float | Absolute divergence (average number of pairwise differences between populations per site). Unlike F_ST, D_xy is not affected by within-population diversity. |
| `da` | float | Net divergence (D_xy minus mean within-population diversity). Approximates divergence time under a simple isolation model. |
| `pbs` | float or null | Population branch statistic for POP1. Only computed when `--pop3` is provided. Positive values indicate POP1-specific differentiation; extreme values suggest population-specific selection. Null when `--pop3` is not given. |
| `n_variants` | integer | Number of variant sites used in computation (after any filtering) |

## Details

### Algorithm

All statistics are computed from allele frequency arrays on Variant nodes:

- **Hudson's F_ST** is computed as the ratio of between-population variance to
  total variance, averaged across sites (ratio of averages, not average of
  ratios). This is the estimator recommended by Bhatia et al. (2013) for its
  robustness to unequal sample sizes.

- **Weir and Cockerham's F_ST** decomposes variance into between-population (a),
  between-individual-within-population (b), and within-individual (c) components,
  returning a / (a + b + c) summed across sites.

- **D_xy** is the mean number of pairwise differences between one chromosome
  sampled from POP1 and one from POP2, per site. Computed as
  sum(p1 * (1-p2) + p2 * (1-p1)) / L where L is the region length.

- **D_a** = D_xy - (pi_1 + pi_2) / 2, removing within-population diversity.

- **PBS** transforms pairwise F_ST values between three populations using
  T = -log(1 - F_ST), then computes PBS_1 = (T_12 + T_13 - T_23) / 2. This
  isolates the branch length specific to POP1 and is the statistic introduced by
  Yi et al. (2010) for detecting high-altitude adaptation in Tibetans.

### PBS requirements

PBS requires three populations. If `--pop3` is not provided, the `pbs` column
is null. The choice of outgroup population affects PBS values: the outgroup
should be a population not subject to the selective pressure being tested. For
example, when testing for selection in East Asians, a common design uses
East Asian (POP1), European (POP2), and African (POP3) populations.

### Conditioning behavior

When `--consequence`, `--pathway`, or `--gene` is specified, only matching
variants are included in the computation. This enables direct comparison of
divergence at functional versus neutral sites. For example, elevated F_ST at
missense sites relative to synonymous sites in the same region is a signature
of divergent selection.

Multiple conditioning options can be combined. The `--min-af` filter applies to
both populations: a variant is excluded if its frequency is below the threshold
in either POP1 or POP2.

### Performance

FAST PATH complexity: O(V x K). Typical performance for a full chromosome
(~1M variants) is 1-5 seconds. Adding `--pop3` for PBS increases computation
by roughly 50% because three pairwise F_ST values must be computed.

### Persistence

This command does not write results to the graph. To persist divergence
statistics per window, use `graphpop genome-scan` with `--pop2`.

## Examples

```bash
# Basic: Fst and Dxy between EUR and EAS on chromosome 22
graphpop divergence chr22 1 50818468 EUR EAS

# PBS with African outgroup, output to file
graphpop divergence chr22 1 50818468 EAS EUR \
    --pop3 YRI \
    -o pbs_eas_chr22.tsv

# Conditioned on missense variants for functional divergence
graphpop divergence chr22 1 50818468 CEU YRI \
    --consequence missense_variant \
    -o fst_missense.tsv

# Pathway-conditioned divergence between rice subpopulations
graphpop divergence chr1 1 43270923 GJ-tmp XI-1A \
    --pathway "Starch and sucrose metabolism" \
    --format json \
    -o starch_divergence.json

# Gene-level divergence with common-variant filter
graphpop divergence chr22 23522552 23838780 CEU CHB \
    --gene EWSR1 \
    --min-af 0.05 \
    -o ewsr1_fst.tsv
```

## See Also

- `graphpop diversity` -- within-population diversity statistics (pi, theta_w, Tajima's D)
- `graphpop genome-scan` -- sliding-window scan including divergence statistics
- `graphpop joint-sfs` -- joint site frequency spectrum between two populations (the full frequency distribution underlying F_ST)
- `graphpop pop-summary` -- whole-chromosome population summary
