# Compute number of segregating sites by length (nSL)

## Description

Computes the number of segregating sites by length (nSL) for every variant on a
chromosome within a single population. nSL is a haplotype-based selection
statistic that measures the extent of haplotype homozygosity around a focal
variant by counting the number of segregating sites over which haplotypes remain
identical, rather than integrating over physical distance. This makes nSL robust
to variation in local recombination rate, which can confound iHS.

This command implements the method of Ferrer-Admetlla et al. (2014). It is a
FULL PATH procedure: it traverses individual-level CARRIES relationships to
reconstruct haplotypes. Use nSL as a complement or alternative to iHS,
particularly in regions with poorly characterized recombination maps or
heterogeneous recombination rates.

## Usage

```
graphpop nsl CHR POPULATION [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome identifier (e.g., `chr22`, `chr1`). Must match the `chr` property on Variant nodes. |
| `POPULATION` | string | yes | -- | Population label (e.g., `EUR`, `GJ-tmp`). Must match a population defined in the graph. |
| `--min-af` | float | no | none | Minimum allele frequency. Variants with AF below this threshold are excluded. |
| `--persist` | flag | no | false | Write computed scores back to Variant nodes as properties `nsl_{pop}` (standardized) and `nsl_unstd_{pop}` (unstandardized). |
| `-o`, `--output` | path | no | stdout | Output file path. If omitted, results are printed to standard output. |
| `--format` | choice | no | `tsv` | Output format. One of: `tsv`, `csv`, `json`. |

## Value

Returns one row per variant that passes the allele frequency filter. Columns:

| Column | Type | Description |
|--------|------|-------------|
| `variantId` | string | Variant identifier in `chr:pos:ref:alt` format. |
| `pos` | integer | Genomic position (1-based). |
| `af` | float | Derived allele frequency in the target population. |
| `nsl_unstd` | float | Unstandardized nSL: `ln(SL_A / SL_D)`, where SL_A and SL_D are the mean number of segregating sites over which haplotypes remain identical for the ancestral and derived alleles, respectively. |
| `nsl` | float | Standardized nSL. Scores are grouped into allele frequency bins and z-score normalized within each bin, analogous to iHS standardization. |

## Details

### Algorithm

1. Variants on the target chromosome are loaded and sorted by position.
2. For each focal variant, the pairwise shared suffix length (SSL) is computed
   for all haplotype pairs carrying the same allele (ancestral or derived).
   SSL counts the number of consecutive variants over which two haplotypes
   are identical, extending outward from the focal site.
3. The mean SSL is computed separately for ancestral-allele haplotypes (SL_A)
   and derived-allele haplotypes (SL_D).
4. The unstandardized score is `ln(SL_A / SL_D)`.
5. Scores are standardized within allele frequency bins (default: 20 bins)
   using z-score normalization.

### Comparison with iHS

| Property | iHS | nSL |
|----------|-----|-----|
| Distance metric | Physical distance (bp) | Number of segregating sites |
| Recombination rate sensitivity | High | Low |
| Genetic map required | Recommended | Not needed |
| Best for | Uniform recombination regions | Variable recombination regions |
| Statistical power | Slightly higher when map is accurate | More robust overall |

### Conditioning on annotation

nSL does **not** support a `--consequence` option at compute time. As with all
haplotype-based statistics, the computation requires an unbroken chain of
variants to correctly measure haplotype extent. Removing variants by annotation
would introduce gaps that corrupt the SSL measurement.

To analyze nSL in a functional context, first persist the scores and then use
`graphpop filter` for post-hoc annotation filtering:

```
graphpop nsl chr22 EUR --persist
graphpop filter nsl chr22 EUR --consequence missense_variant
```

### Persistence

When `--persist` is set, two properties are written to each qualifying Variant
node:

- `nsl_{pop}` -- the standardized score
- `nsl_unstd_{pop}` -- the unstandardized score

For example, with population `GJ-tmp`, properties `nsl_GJ-tmp` and
`nsl_unstd_GJ-tmp` are created.

### Memory and performance

nSL uses pairwise SSL computation, which is O(n * h) where n is the number of
variants and h is the number of haplotypes. The chunked processing approach
keeps peak memory at approximately 200 MB. Because nSL counts segregating sites
rather than integrating over physical distance, it avoids the need to load a
genetic map, simplifying the computation.

## Examples

Compute nSL for a European population on chromosome 22:

```
graphpop nsl chr22 EUR -o nsl_eur_chr22.tsv
```

Persist scores with a minimum allele frequency filter:

```
graphpop nsl chr22 EUR --min-af 0.05 --persist -o nsl_eur.tsv
```

Persist and then filter for variants in a specific gene:

```
graphpop nsl chr1 GJ-tmp --persist
graphpop filter nsl chr1 GJ-tmp --gene GW5 --min-score 2.0
```

Output as CSV for import into a spreadsheet:

```
graphpop nsl chr22 AFR --format csv -o nsl_afr_chr22.csv
```

Compare nSL and iHS on the same population (run both, then analyze):

```
graphpop nsl chr22 EUR --persist
graphpop ihs chr22 EUR --persist
graphpop filter nsl chr22 EUR --consequence missense_variant -o nsl_missense.tsv
graphpop filter ihs chr22 EUR --consequence missense_variant -o ihs_missense.tsv
```

## See Also

- [graphpop ihs](ihs.md) -- integrated haplotype score (physical-distance based)
- [graphpop xpehh](xpehh.md) -- cross-population extended haplotype homozygosity
- [graphpop filter](filter.md) -- post-hoc annotation filtering of persisted statistics
- [graphpop garud-h](garud-h.md) -- haplotype homozygosity statistics for sweep detection

## References

Ferrer-Admetlla A, Liang M, Korneliussen T, Nielsen R. On detecting incomplete
soft or hard selective sweeps using haplotype structure. *Molecular Biology and
Evolution*. 2014;31(5):1275-1291.
