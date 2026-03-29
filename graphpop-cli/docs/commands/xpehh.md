# Compute cross-population extended haplotype homozygosity (XP-EHH)

## Description

Computes the cross-population extended haplotype homozygosity (XP-EHH) statistic
for every variant on a chromosome, comparing two populations. XP-EHH detects
selective sweeps that have reached or nearly reached fixation in one population
but remain polymorphic in the other. Unlike iHS, which detects ongoing sweeps
within a population, XP-EHH is powered to detect completed or near-completed
sweeps by contrasting haplotype lengths between populations.

This command implements the method of Sabeti et al. (2007). It is a FULL PATH
procedure: it traverses individual-level CARRIES relationships in both
populations to reconstruct haplotypes. Positive XP-EHH values indicate
selection in `POP1`; negative values indicate selection in `POP2`.

## Usage

```
graphpop xpehh CHR POP1 POP2 [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome identifier (e.g., `chr22`, `chr1`). Must match the `chr` property on Variant nodes. |
| `POP1` | string | yes | -- | First population label (e.g., `EUR`). Positive XP-EHH values indicate longer haplotypes (selection) in this population. |
| `POP2` | string | yes | -- | Second population label (e.g., `AFR`). Negative XP-EHH values indicate longer haplotypes (selection) in this population. |
| `--min-af` | float | no | none | Minimum allele frequency. Variants below this AF in either population are excluded. |
| `--persist` | flag | no | false | Write computed scores back to Variant nodes as properties `xpehh_{pop1}_{pop2}` (standardized) and `xpehh_unstd_{pop1}_{pop2}` (unstandardized). |
| `-o`, `--output` | path | no | stdout | Output file path. If omitted, results are printed to standard output. |
| `--format` | choice | no | `tsv` | Output format. One of: `tsv`, `csv`, `json`. |

## Value

Returns one row per variant that passes the allele frequency filter. Columns:

| Column | Type | Description |
|--------|------|-------------|
| `variantId` | string | Variant identifier in `chr:pos:ref:alt` format. |
| `pos` | integer | Genomic position (1-based). |
| `af_pop1` | float | Derived allele frequency in `POP1`. |
| `af_pop2` | float | Derived allele frequency in `POP2`. |
| `xpehh_unstd` | float | Unstandardized XP-EHH: `ln(iHH_pop1 / iHH_pop2)`. Positive values indicate longer haplotypes in POP1. |
| `xpehh` | float | Standardized XP-EHH. Genome-wide z-score normalization (mean 0, standard deviation 1) applied across all variants, without allele frequency binning. |

## Details

### Algorithm

1. Variants on the target chromosome are loaded and sorted by position.
2. For each focal variant, integrated EHH (iHH) is computed separately for
   each population by integrating EHH outward from the focal site until decay
   drops below 0.05.
3. The unstandardized score is `ln(iHH_pop1 / iHH_pop2)`.
4. Unlike iHS, XP-EHH uses genome-wide z-score standardization (not
   frequency-binned) because the comparison is between populations rather than
   between alleles at the same site.

### Conditioning on annotation

XP-EHH does **not** support a `--consequence` option at compute time. As with
iHS, EHH decay requires integrating over a contiguous chain of variants.
Removing variants by functional annotation would corrupt the haplotype
integration distances.

To analyze XP-EHH in a functional context, first persist the scores and then use
`graphpop filter` for post-hoc annotation filtering:

```
graphpop xpehh chr22 EUR AFR --persist
graphpop filter xpehh chr22 EUR --pop2 AFR --consequence missense_variant
```

### Persistence

When `--persist` is set, two properties are written to each qualifying Variant
node:

- `xpehh_{pop1}_{pop2}` -- the standardized score
- `xpehh_unstd_{pop1}_{pop2}` -- the unstandardized score

For example, with populations `EUR` and `AFR`, properties `xpehh_EUR_AFR` and
`xpehh_unstd_EUR_AFR` are created. Note that population order matters:
`xpehh_EUR_AFR` and `xpehh_AFR_EUR` are different properties (equal in
magnitude but opposite in sign).

### Interpretation

| XP-EHH value | Interpretation |
|--------------|----------------|
| > 2.0 | Strong evidence for selection in POP1 |
| < -2.0 | Strong evidence for selection in POP2 |
| -2.0 to 2.0 | No strong signal |

### Memory and performance

XP-EHH uses the same chunked processing strategy as iHS (5 Mb core + 2 Mb
margin). Because two populations must be traversed, runtime is roughly twice
that of a single-population iHS computation. Peak memory remains approximately
200 MB.

## Examples

Compute XP-EHH between Europeans and Africans on chromosome 22:

```
graphpop xpehh chr22 EUR AFR -o xpehh_eur_afr_chr22.tsv
```

Persist scores to the graph with a minimum allele frequency threshold:

```
graphpop xpehh chr22 EUR AFR --min-af 0.05 --persist -o xpehh_eur_afr.tsv
```

Persist and then filter for missense variants in a specific pathway:

```
graphpop xpehh chr22 EUR AFR --persist
graphpop filter xpehh chr22 EUR --pop2 AFR --pathway "Immune System"
```

Output JSON for downstream analysis:

```
graphpop xpehh chr1 GJ-tmp XI-1A --format json -o xpehh_rice.json
```

Compare two rice subpopulations with annotation filtering:

```
graphpop xpehh chr1 GJ-tmp XI-1A --min-af 0.05 --persist
graphpop filter xpehh chr1 GJ-tmp --pop2 XI-1A --gene GW5 --min-score 2.0
```

## See Also

- [graphpop ihs](ihs.md) -- within-population integrated haplotype score
- [graphpop nsl](nsl.md) -- number of segregating sites by length
- [graphpop filter](filter.md) -- post-hoc annotation filtering of persisted statistics
- [graphpop divergence](divergence.md) -- Fst, Dxy, PBS between populations (FAST PATH)

## References

Sabeti PC, Varilly P, Fry B, et al. Genome-wide detection and characterization
of positive selection in human populations. *Nature*. 2007;449(7164):913-918.
