# Compute Garud's haplotype homozygosity statistics (H1, H12, H2/H1)

## Description

Computes Garud's H statistics (H1, H12, H2/H1) in sliding windows across a
chromosome for a single population. These statistics summarize the haplotype
frequency spectrum to detect recent selective sweeps and distinguish hard sweeps
(single beneficial mutation) from soft sweeps (multiple beneficial mutations or
standing variation).

This command implements the method of Garud et al. (2015). It is a FULL PATH
procedure: it traverses individual-level CARRIES relationships to reconstruct
haplotypes within each window, then computes haplotype frequencies. Use Garud's
H when you need a window-level sweep signal that can differentiate sweep types,
complementing per-variant statistics like iHS and nSL.

## Usage

```
graphpop garud-h CHR POPULATION WINDOW_SIZE STEP_SIZE [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome identifier (e.g., `chr22`, `chr1`). Must match the `chr` property on Variant nodes. |
| `POPULATION` | string | yes | -- | Population label (e.g., `EUR`, `GJ-tmp`). Must match a population defined in the graph. |
| `WINDOW_SIZE` | integer | yes | -- | Window size in base pairs. Typical values: 50000 (50 kb) for fine-grained scans, 200000 (200 kb) for broad sweeps. |
| `STEP_SIZE` | integer | yes | -- | Step size in base pairs. Controls overlap between consecutive windows. Set equal to `WINDOW_SIZE` for non-overlapping windows, or to a fraction (e.g., half) for overlapping scans. |
| `--min-af` | float | no | none | Minimum allele frequency. Variants below this AF are excluded from haplotype construction. |
| `-o`, `--output` | path | no | stdout | Output file path. If omitted, results are printed to standard output. |
| `--format` | choice | no | `tsv` | Output format. One of: `tsv`, `csv`, `json`. |

## Value

Returns one row per window. Columns:

| Column | Type | Description |
|--------|------|-------------|
| `chr` | string | Chromosome identifier. |
| `start` | integer | Window start position (1-based, inclusive). |
| `end` | integer | Window end position (inclusive). |
| `population` | string | Population label. |
| `h1` | float | H1: sum of squared haplotype frequencies. Equivalent to expected haplotype homozygosity. Ranges from 1/n (all haplotypes unique) to 1.0 (single haplotype fixed). |
| `h12` | float | H12: H1 computed after merging the two most frequent haplotype classes. Sensitive to sweeps because a sweep elevates one or two haplotype classes. |
| `h2_h1` | float | H2/H1: ratio of H2 (H1 minus the squared frequency of the most common haplotype) to H1. Distinguishes hard from soft sweeps. |
| `hap_diversity` | float | Haplotype diversity: `1 - H1`. The probability that two randomly chosen haplotypes differ. |
| `n_haplotypes` | integer | Number of distinct haplotypes observed in the window. |
| `n_variants` | integer | Number of variants in the window used for haplotype construction. |

## Details

### Algorithm

1. The chromosome is divided into windows of size `WINDOW_SIZE`, advancing by
   `STEP_SIZE`.
2. Within each window, haplotypes are reconstructed from CARRIES relationships
   for all samples in the population.
3. Haplotype frequencies p_1 >= p_2 >= ... >= p_k are computed.
4. Statistics are calculated:
   - H1 = sum(p_i^2)
   - H12 = (p_1 + p_2)^2 + sum(p_i^2 for i >= 3)
   - H2 = H1 - p_1^2
   - H2/H1 = H2 / H1

### Interpretation

| Signal | H12 | H2/H1 | Interpretation |
|--------|-----|-------|----------------|
| Hard sweep | > 0.3 | < 0.1 | Single haplotype at high frequency |
| Soft sweep | > 0.3 | > 0.1 | Multiple haplotypes at high frequency |
| Neutral | < 0.1 | variable | No sweep signal |
| Intermediate | 0.1 - 0.3 | -- | Weak or ancient sweep, or drift |

A threshold of H12 > 0.3 is commonly used to identify candidate sweep regions.
Among candidates, H2/H1 < 0.1 suggests a hard sweep (one dominant haplotype),
while H2/H1 > 0.1 suggests a soft sweep (multiple haplotypes swept up).

### Persistence

This command is **read-only** and does not persist results to the graph.
Garud's H statistics are window-level summaries. If you need to query these
results later, save them to a file and use standard filtering tools, or use
`graphpop filter h12` if the values have been stored on GenomicWindow nodes by
`graphpop genome-scan`.

### Window size considerations

- **Small windows (20-50 kb):** Higher resolution but noisier. Fewer variants
  per window may reduce haplotype diversity estimates. Best for fine-mapping
  within a candidate region.
- **Large windows (100-500 kb):** Smoother signal, better for genome-wide
  scans. May merge distinct signals in gene-dense regions.
- **Step size:** Using `STEP_SIZE = WINDOW_SIZE / 2` provides good balance
  between resolution and computation time.

### Memory and performance

Each window requires reconstructing haplotypes for all samples in the
population. Memory usage depends on window size and sample count. For typical
parameters (200 kb windows, 2,500 samples), peak memory is approximately
300 MB. Genome-wide scans on a human chromosome complete in 5-15 minutes.

## Examples

Compute Garud's H in 200 kb windows with 100 kb steps on chromosome 22:

```
graphpop garud-h chr22 EUR 200000 100000 -o garud_eur_chr22.tsv
```

Fine-grained scan with 50 kb non-overlapping windows:

```
graphpop garud-h chr22 EUR 50000 50000 -o garud_fine_chr22.tsv
```

Rice population scan with a minimum allele frequency filter:

```
graphpop garud-h chr1 GJ-tmp 100000 50000 --min-af 0.05 -o garud_gj_chr1.tsv
```

Output as JSON for programmatic downstream analysis:

```
graphpop garud-h chr22 AFR 200000 100000 --format json -o garud_afr.json
```

Broad scan to identify candidate regions, then refine with iHS:

```
graphpop garud-h chr22 EUR 200000 100000 -o garud_broad.tsv
# Inspect garud_broad.tsv for windows with H12 > 0.3
graphpop ihs chr22 EUR --persist
graphpop filter ihs chr22 EUR --min-score 2.0
```

## See Also

- [graphpop ihs](ihs.md) -- per-variant integrated haplotype score
- [graphpop nsl](nsl.md) -- per-variant nSL (recombination-rate robust)
- [graphpop genome-scan](genome-scan.md) -- sliding-window diversity and divergence (FAST PATH)
- [graphpop diversity](diversity.md) -- nucleotide diversity and Tajima's D (FAST PATH)
- [graphpop filter](filter.md) -- post-hoc annotation filtering of persisted statistics

## References

Garud NR, Messer PW, Buzbas EO, Petrov DA. Recent selective sweeps in North
American Drosophila melanogaster show signatures of soft sweeps. *PLoS
Genetics*. 2015;11(2):e1005004.
