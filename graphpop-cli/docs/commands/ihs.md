# Compute integrated haplotype score (iHS)

## Description

Computes the integrated haplotype score (iHS) for every variant on a chromosome
within a single population. iHS measures the decay of extended haplotype
homozygosity (EHH) around a focal variant, comparing the ancestral and derived
alleles. Variants under recent positive selection exhibit unusually long
haplotypes on the favored allele, producing extreme iHS values.

This command implements the method of Voight et al. (2006). It is a FULL PATH
procedure: it traverses individual-level CARRIES relationships to reconstruct
haplotypes, rather than operating on pre-aggregated allele count arrays. Use iHS
when screening a chromosome for ongoing or recent selective sweeps within a
population. For cross-population comparisons, see `graphpop xpehh`.

## Usage

```
graphpop ihs CHR POPULATION [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome identifier (e.g., `chr22`, `chr1`). Must match the `chr` property on Variant nodes. |
| `POPULATION` | string | yes | -- | Population label (e.g., `EUR`, `GJ-tmp`). Must match a population defined in the graph. |
| `--min-af` | float | no | none | Minimum allele frequency. Variants with AF below this threshold in the target population are excluded. Typical value: 0.05. |
| `--persist` | flag | no | false | Write computed scores back to Variant nodes as properties `ihs_{pop}` (standardized) and `ihs_unstd_{pop}` (unstandardized). |
| `-o`, `--output` | path | no | stdout | Output file path. If omitted, results are printed to standard output. |
| `--format` | choice | no | `tsv` | Output format. One of: `tsv`, `csv`, `json`. |

## Value

Returns one row per variant that passes the allele frequency filter. Columns:

| Column | Type | Description |
|--------|------|-------------|
| `variantId` | string | Variant identifier in `chr:pos:ref:alt` format. |
| `pos` | integer | Genomic position (1-based). |
| `af` | float | Derived allele frequency in the target population. |
| `ihs_unstd` | float | Unstandardized iHS: `ln(iHH_A / iHH_D)`, where iHH_A and iHH_D are the integrated EHH for the ancestral and derived alleles, respectively. |
| `ihs` | float | Standardized iHS. The unstandardized scores are grouped into allele frequency bins and z-score normalized within each bin, yielding a distribution with mean 0 and standard deviation 1 genome-wide. |

## Details

### Algorithm

1. Variants on the target chromosome are loaded and sorted by position.
2. The chromosome is processed in chunks of 5 Mb (core region) with a 2 Mb
   margin on each side to avoid edge effects in EHH decay.
3. For each focal variant, EHH is computed separately for haplotypes carrying
   the ancestral allele and for those carrying the derived allele, integrating
   outward until EHH decays below 0.05.
4. The unstandardized score is `ln(iHH_ancestral / iHH_derived)`.
5. Scores are standardized by allele frequency bin (default: 20 equal-width bins
   across the [0, 1] AF range). Within each bin, scores are z-score normalized.

### Conditioning on annotation

iHS does **not** support a `--consequence` option at compute time. This is by
design: EHH decay is measured across a contiguous chain of variants along the
chromosome. Removing variants by functional annotation would break the physical
linkage chain and produce incorrect haplotype integration distances. iHS must
always be computed on the full set of variants.

To analyze iHS in a functional context, first persist the scores and then use
`graphpop filter` for post-hoc annotation filtering:

```
graphpop ihs chr22 EUR --persist
graphpop filter ihs chr22 EUR --consequence missense_variant
```

### Persistence

When `--persist` is set, two properties are written to each qualifying Variant
node:

- `ihs_{pop}` -- the standardized score
- `ihs_unstd_{pop}` -- the unstandardized score

For example, with population `EUR`, properties `ihs_EUR` and `ihs_unstd_EUR`
are created. These properties are required for downstream `graphpop filter`
queries.

### Memory and performance

The chunked processing strategy (5 Mb core + 2 Mb margin) keeps peak memory
usage at approximately 200 MB regardless of chromosome length. A typical human
chromosome (e.g., chr22 with ~1 M variants and ~2,500 samples) completes in
2-5 minutes depending on hardware and Neo4j configuration.

## Examples

Compute iHS for Europeans on chromosome 22, writing results to a TSV file:

```
graphpop ihs chr22 EUR -o ihs_eur_chr22.tsv
```

Compute iHS with a minimum allele frequency filter and persist scores to the
graph for downstream filtering:

```
graphpop ihs chr22 EUR --min-af 0.05 --persist -o ihs_eur_chr22.tsv
```

Persist scores and then query only missense variants with extreme iHS:

```
graphpop ihs chr22 EUR --persist
graphpop filter ihs chr22 EUR --consequence missense_variant --min-score 2.0
```

Output as JSON for programmatic consumption:

```
graphpop ihs chr1 GJ-tmp --format json -o ihs_gj_chr1.json
```

Run iHS on rice indica population across chromosome 1:

```
graphpop ihs chr1 XI-1A --min-af 0.05 --persist
```

## See Also

- [graphpop xpehh](xpehh.md) -- cross-population extended haplotype homozygosity
- [graphpop nsl](nsl.md) -- number of segregating sites by length (recombination-rate robust alternative)
- [graphpop filter](filter.md) -- post-hoc annotation filtering of persisted statistics
- [graphpop genome-scan](genome-scan.md) -- sliding-window genome scan (FAST PATH)

## References

Voight BF, Kudaravalli S, Wen X, Pritchard JK. A map of recent positive
selection in the human genome. *PLoS Biology*. 2006;4(3):e72.
