# Detect runs of homozygosity (ROH)

## Description

Detects runs of homozygosity (ROH) for each sample in a population on a given
chromosome. ROH are long contiguous stretches of homozygous genotypes that arise
from consanguinity, population bottlenecks, or background selection. The
distribution of ROH lengths carries information about demographic history: long
ROH (> 1 Mb) indicate recent inbreeding, while short ROH (< 100 kb) reflect
ancient population size reduction.

This command implements two detection methods: a 2-state Hidden Markov Model
(HMM) following Narasimhan et al. (2016), and a simpler sliding-window
approach. It is a FULL PATH procedure: it traverses individual-level CARRIES
relationships to reconstruct per-sample genotype sequences. Output is per-sample
summary statistics; this command does not persist results to the graph.

## Usage

```
graphpop roh CHR POPULATION [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome identifier (e.g., `chr22`, `chr1`). Must match the `chr` property on Variant nodes. |
| `POPULATION` | string | yes | -- | Population label (e.g., `EUR`, `GJ-tmp`). Must match a population defined in the graph. |
| `--method` | choice | no | `hmm` | ROH detection method. One of: `hmm` (2-state Viterbi HMM, recommended), `window` (sliding window with heterozygosity threshold). |
| `--min-length` | integer | no | none | Minimum ROH length in base pairs. ROH shorter than this threshold are discarded before summary computation. Typical values: 100000 (100 kb) for ancient signals, 1000000 (1 Mb) for recent inbreeding. |
| `-o`, `--output` | path | no | stdout | Output file path. If omitted, results are printed to standard output. |
| `--format` | choice | no | `tsv` | Output format. One of: `tsv`, `csv`, `json`. |

## Value

Returns one row per sample in the population. Columns:

| Column | Type | Description |
|--------|------|-------------|
| `sampleId` | string | Sample identifier. |
| `n_roh` | integer | Number of ROH segments detected on this chromosome for this sample. |
| `total_length` | integer | Sum of all ROH segment lengths in base pairs. |
| `froh` | float | Fraction of the chromosome covered by ROH: `total_length / chromosome_length`. Values range from 0 to 1. FROH > 0.0625 is consistent with offspring of first-cousin mating. |
| `mean_length` | float | Mean ROH segment length in base pairs. |
| `max_length` | integer | Length of the longest ROH segment in base pairs. |

## Details

### HMM method (default)

The HMM method uses a 2-state Viterbi algorithm following Narasimhan et al.
(2016):

- **State 1 (non-ROH):** Expected heterozygosity matches the population average.
- **State 2 (ROH):** Expected heterozygosity is near zero (small emission
  probability for heterozygous sites, accounting for genotyping error).

Transition probabilities are estimated from inter-variant distances and a
configurable recombination rate. The Viterbi algorithm finds the most likely
sequence of states, and contiguous ROH-state runs are reported as segments.

### Window method

The sliding-window method scans the chromosome with a fixed-size window and
flags windows where the heterozygosity rate falls below a threshold.
Consecutive flagged windows are merged into ROH segments. This method is faster
but less accurate than the HMM, particularly at ROH boundaries.

### Persistence

This command is **read-only** and does not persist results to the graph. ROH
statistics are per-sample summaries that do not map naturally to Variant node
properties. Results are written to the output file or standard output only.

### FROH interpretation

| FROH range | Interpretation |
|------------|----------------|
| < 0.01 | Outbred, large ancestral population |
| 0.01 - 0.0625 | Moderate background inbreeding or bottleneck |
| 0.0625 - 0.125 | Consistent with second-cousin or first-cousin-once-removed mating |
| > 0.125 | Strong inbreeding signal (first-cousin or closer) |

### ROH length classes

| Length class | Origin | Approximate age |
|-------------|--------|-----------------|
| < 100 kb | Ancient Ne reduction | > 50 generations |
| 100 kb - 1 Mb | Intermediate bottleneck | 10 - 50 generations |
| > 1 Mb | Recent consanguinity | < 10 generations |

### Memory and performance

ROH detection loads genotype data for all samples in the population on the
target chromosome. Memory usage scales linearly with sample count and variant
count. For a typical dataset (3,000 samples, 1 M variants), expect
approximately 500 MB peak memory. The HMM method is computationally heavier
than the window method but remains tractable for chromosome-scale analysis.

## Examples

Detect ROH using the default HMM method for Europeans on chromosome 22:

```
graphpop roh chr22 EUR -o roh_eur_chr22.tsv
```

Use the sliding-window method with a minimum ROH length of 1 Mb:

```
graphpop roh chr22 EUR --method window --min-length 1000000 -o roh_eur_long.tsv
```

Detect ROH for a rice population, focusing on recent inbreeding (long ROH):

```
graphpop roh chr1 GJ-tmp --min-length 500000 -o roh_gj_chr1.tsv
```

Output as JSON for programmatic analysis:

```
graphpop roh chr22 AFR --format json -o roh_afr_chr22.json
```

Compare HMM and window methods on the same data:

```
graphpop roh chr22 EUR --method hmm -o roh_hmm.tsv
graphpop roh chr22 EUR --method window -o roh_window.tsv
```

## See Also

- [graphpop ihs](ihs.md) -- integrated haplotype score (detects sweeps, complementary signal)
- [graphpop garud-h](garud-h.md) -- haplotype homozygosity in windows (related haplotype analysis)
- [graphpop diversity](diversity.md) -- nucleotide diversity and heterozygosity (FAST PATH)
- [graphpop filter](filter.md) -- post-hoc annotation filtering (for persisted statistics)

## References

Narasimhan V, Danecek P, Scally A, Xue Y, Tyler-Smith C, Durbin R. BCFtools/RoH:
a hidden Markov model approach for detecting autozygosity from next-generation
sequencing data. *Bioinformatics*. 2016;32(11):1749-1751.
