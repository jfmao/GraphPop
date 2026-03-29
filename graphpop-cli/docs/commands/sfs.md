# Compute the site frequency spectrum

## Description

`graphpop sfs` computes the site frequency spectrum (SFS) for a single
population over a genomic region. The SFS is a vector counting how many variant
sites have each possible allele count, from singletons through fixed derived
alleles. It is the fundamental summary of allele frequency variation in a
population and underlies most population genetics inference methods.

By default, the command returns the **folded** spectrum, which uses minor allele
counts and does not require ancestral state information. When `--unfolded` is
specified, the command returns the **unfolded** (polarized) spectrum based on
derived allele counts, which requires that ancestral allele annotations are
present in the graph. Only variants with known ancestral state contribute to the
unfolded SFS.

This command operates on the **FAST PATH**, reading allele count arrays directly
from Variant nodes with O(V x K) complexity.

Use `sfs` when you need the full frequency distribution for demographic
inference (e.g., input to dadi, moments, or fastsimcoal2), for goodness-of-fit
tests against neutral models, or for visualizing the shape of variation in a
population.

## Usage

```
graphpop sfs CHR START END POPULATION [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome name (e.g., `chr22`, `chr1`) |
| `START` | integer | yes | -- | Start position of the genomic region (1-based, inclusive) |
| `END` | integer | yes | -- | End position of the genomic region (1-based, inclusive) |
| `POPULATION` | string | yes | -- | Population identifier as stored in the graph |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--unfolded` | flag | false | Compute the unfolded (derived allele) SFS. Requires ancestral allele annotations in the graph. Variants without ancestral state are excluded. |
| `--consequence` | string | none | Filter variants by VEP consequence type |
| `--pathway` | string | none | Filter variants to those in genes belonging to the named pathway |
| `--gene` | string | none | Filter variants to those within the named gene |
| `--min-af` | float | none | Exclude variants with allele frequency below this threshold |
| `--max-af` | float | none | Exclude variants with allele frequency above this threshold |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json` |

## Value

Returns a single row with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `sfs` | list of integers | The site frequency spectrum as an array. For the folded SFS, entry i gives the count of sites where the minor allele count is i (i = 0, 1, ..., floor(n/2)). For the unfolded SFS, entry i gives the count of sites where the derived allele count is i (i = 0, 1, ..., n). |
| `n_variants` | integer | Total number of variant sites used in the computation |
| `max_ac` | integer | Maximum possible allele count (equal to 2N where N is the population sample size, i.e., the number of chromosomes sampled) |
| `n_polarized` | integer | Number of variants with known ancestral allele. For the folded SFS, this may exceed the number of variants actually used in the unfolded spectrum. For the unfolded SFS, this equals the number of sites that contributed. |

## Details

### Algorithm

The SFS is constructed by binning variants according to their allele count in
the target population:

- **Folded SFS**: For each variant, the minor allele count is
  min(ac, n - ac) where ac is the allele count and n is the allele number.
  The spectrum has floor(n/2) + 1 entries.

- **Unfolded SFS**: For each variant with known ancestral state, the derived
  allele count is used directly. If the reference allele is ancestral, the
  derived count equals ac; if the reference is derived, the derived count
  equals n - ac. Variants without ancestral annotation are excluded, and their
  count is reflected in the difference between `n_variants` and `n_polarized`.
  The spectrum has n + 1 entries.

The zeroth entry (sfs[0]) counts monomorphic sites in the region if present in
the graph; in practice, most graphs store only polymorphic sites, so sfs[0] may
be zero or absent.

### Conditioning behavior

Annotation conditioning (`--consequence`, `--pathway`, `--gene`) restricts the
SFS to a functional subset of variants. This is particularly useful for
comparing the shape of the SFS across functional categories. For example,
missense variants typically show a shift toward rare alleles compared to
synonymous variants, reflecting purifying selection.

The `--min-af` and `--max-af` options truncate the spectrum, which is useful for
focusing on common or rare variation but changes the interpretation of the SFS
shape.

### Performance

FAST PATH: O(V x K). Construction of the SFS for a full chromosome typically
completes in 1-3 seconds. The spectrum is a lightweight summary: even for very
large sample sizes (n = 6000 chromosomes), the output is a single vector of
at most 3001 entries.

### Persistence

This command does not write results to the graph. Output is returned to stdout
or to the specified file.

## Examples

```bash
# Basic: folded SFS for YRI on chromosome 22
graphpop sfs chr22 1 50818468 YRI

# Unfolded (polarized) SFS with output to file
graphpop sfs chr22 1 50818468 YRI \
    --unfolded \
    -o sfs_yri_unfolded.tsv

# SFS conditioned on missense variants (compare against synonymous)
graphpop sfs chr22 1 50818468 CEU \
    --consequence missense_variant \
    -o sfs_ceu_missense.tsv

graphpop sfs chr22 1 50818468 CEU \
    --consequence synonymous_variant \
    -o sfs_ceu_synonymous.tsv

# SFS for variants in a specific pathway, JSON format
graphpop sfs chr1 1 43270923 GJ-tmp \
    --pathway "Starch and sucrose metabolism" \
    --format json \
    -o sfs_starch.json

# Rare-variant SFS only (singletons through 5% frequency)
graphpop sfs chr22 1 50818468 EUR --max-af 0.05
```

## See Also

- `graphpop joint-sfs` -- two-dimensional joint SFS between two populations
- `graphpop diversity` -- summary statistics derived from the SFS (pi, theta_w, Tajima's D)
- `graphpop pop-summary` -- whole-chromosome summary including SFS-derived statistics
