# Compute pairwise linkage disequilibrium (r2 and D')

## Description

Computes pairwise linkage disequilibrium (LD) between all variant pairs within a
specified genomic region, reporting both r-squared (r2) and D-prime (D')
statistics. LD measures the non-random association of alleles at different loci
and is fundamental to understanding haplotype structure, mapping disease
associations, and designing tag-SNP panels.

This command is a FULL PATH procedure: it traverses individual-level CARRIES
relationships to reconstruct two-locus haplotypes for every variant pair within
the distance threshold. It can optionally persist LD relationships as edges in
the graph, enabling downstream graph-based analyses such as LD-clumping or
haplotype block detection via graph traversals.

## Usage

```
graphpop ld CHR START END POPULATION MAX_DIST THRESHOLD [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome identifier (e.g., `chr22`, `chr1`). Must match the `chr` property on Variant nodes. |
| `START` | integer | yes | -- | Start position of the region (1-based, inclusive). |
| `END` | integer | yes | -- | End position of the region (inclusive). |
| `POPULATION` | string | yes | -- | Population label (e.g., `EUR`, `GJ-tmp`). LD is population-specific; different populations may show different LD patterns. |
| `MAX_DIST` | integer | yes | -- | Maximum physical distance in base pairs between variant pairs. Pairs separated by more than this distance are not computed. Typical values: 100000 (100 kb) for LD decay analysis, 500000 (500 kb) for comprehensive mapping. |
| `THRESHOLD` | float | yes | -- | Minimum r2 value to include in output. Pairs with r2 below this threshold are discarded. Default recommendation: 0.2 for typical analyses, 0.0 to report all pairs (use with caution -- output can be very large). |
| `--persist` | flag | no | false | Write LD edges between Variant nodes in the graph. Creates `LD` relationships with `r2`, `dprime`, and `distance` properties. |
| `--min-af` | float | no | none | Minimum allele frequency. Variants with AF below this threshold are excluded from all pairwise comparisons. Rare variants tend to have inflated r2 with nearby common variants. |
| `-o`, `--output` | path | no | stdout | Output file path. If omitted, results are printed to standard output. |
| `--format` | choice | no | `tsv` | Output format. One of: `tsv`, `csv`, `json`. |

## Value

Returns one row per variant pair that meets the distance and r2 thresholds.
Columns:

| Column | Type | Description |
|--------|------|-------------|
| `variant1` | string | First variant identifier in `chr:pos:ref:alt` format (lower position). |
| `variant2` | string | Second variant identifier in `chr:pos:ref:alt` format (higher position). |
| `r2` | float | Squared Pearson correlation between allele dosages at the two loci. Ranges from 0 (no association) to 1 (perfect association). |
| `dprime` | float | Normalized coefficient of linkage disequilibrium (Lewontin's D'). Ranges from 0 to 1. D' = 1 indicates no historical recombination between the two sites (within the sample). |
| `distance` | integer | Physical distance in base pairs between the two variants. |

## Details

### Algorithm

1. All variants in the region [START, END] with AF above `--min-af` are loaded
   and sorted by position.
2. For each variant pair (i, j) where `pos_j - pos_i <= MAX_DIST`:
   a. Two-locus haplotypes are reconstructed from CARRIES relationships.
   b. Haplotype frequencies are computed: p_AB, p_Ab, p_aB, p_ab.
   c. D = p_AB - p_A * p_B (deviation from linkage equilibrium).
   d. D' = D / D_max, where D_max depends on the sign of D.
   e. r2 = D^2 / (p_A * p_a * p_B * p_b).
3. Pairs with r2 >= THRESHOLD are included in the output.

### Persistence

When `--persist` is set, LD relationships are written as edges between Variant
nodes:

```
(v1:Variant)-[:LD {r2: 0.85, dprime: 0.97, distance: 1234, population: 'EUR'}]->(v2:Variant)
```

These edges enable graph-native LD queries such as:
- Finding all variants in LD with a candidate SNP
- LD-clumping via graph traversal
- Haplotype block detection using connected components on LD edges

Persisted LD edges are population-specific (tagged with the `population`
property). Running `--persist` for different populations creates separate edges.

### Output size

LD computation is O(n^2) in the number of variants within the distance window.
For dense regions (e.g., > 10,000 variants in 500 kb), the number of pairs can
be very large. Use a reasonable `THRESHOLD` (>= 0.2) and `MAX_DIST` (<= 500 kb)
to keep output manageable. For genome-wide LD maps, consider processing
chromosome segments sequentially.

### Memory and performance

Memory usage scales with the number of variants in the region and the sample
size. For a region with 10,000 variants and 2,500 samples, peak memory is
approximately 300-500 MB. Larger regions or lower thresholds increase both
computation time and memory usage significantly.

## Examples

Compute LD in a 1 Mb region on chromosome 22, reporting pairs with r2 >= 0.2:

```
graphpop ld chr22 20000000 21000000 EUR 100000 0.2 -o ld_eur.tsv
```

Compute all LD pairs (no r2 threshold) within 50 kb for a small region:

```
graphpop ld chr22 20000000 20500000 EUR 50000 0.0 -o ld_all.tsv
```

Persist LD edges to the graph for downstream graph queries:

```
graphpop ld chr22 20000000 21000000 EUR 100000 0.2 --persist -o ld_eur.tsv
```

Filter by minimum allele frequency to focus on common variant LD:

```
graphpop ld chr22 20000000 21000000 EUR 200000 0.2 --min-af 0.05 -o ld_common.tsv
```

Compute LD for a rice population in a candidate gene region:

```
graphpop ld chr1 5000000 5500000 GJ-tmp 100000 0.1 --min-af 0.05 -o ld_rice.tsv
```

## See Also

- [graphpop ihs](ihs.md) -- integrated haplotype score (uses LD structure implicitly)
- [graphpop nsl](nsl.md) -- nSL (haplotype-based, related to LD extent)
- [graphpop genome-scan](genome-scan.md) -- sliding-window genome scan (FAST PATH)
- [graphpop filter](filter.md) -- post-hoc annotation filtering of persisted statistics

## References

Lewontin RC. The interaction of selection and linkage. I. General
considerations; heterotic models. *Genetics*. 1964;49(1):49-67.

Hill WG, Robertson A. Linkage disequilibrium in finite populations.
*Theoretical and Applied Genetics*. 1968;38(6):226-231.
