# Compute whole-chromosome population summary statistics

## Description

`graphpop pop-summary` computes summary statistics for an entire chromosome in a
single population and **writes the results to the Population node** in the
graph. This provides a chromosome-level characterization of genetic variation:
nucleotide diversity, Watterson's theta, Tajima's D, Fay and Wu's H, mean
heterozygosity, and the inbreeding coefficient.

Unlike `graphpop diversity`, which computes statistics for an arbitrary region
and returns results without modifying the graph, `pop-summary` covers the full
chromosome and persists the results as properties on the Population node. This
makes the summary statistics available for subsequent Cypher queries without
recomputation.

This command operates on the **FAST PATH**, reading allele count arrays from
Variant nodes with O(V x K) complexity.

Use `pop-summary` to generate the definitive per-chromosome summary for each
population in your dataset, to populate Population nodes for downstream queries,
or as the first step in a full-genome characterization workflow.

## Usage

```
graphpop pop-summary CHR POPULATION [OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `CHR` | string | yes | -- | Chromosome name (e.g., `chr22`, `chr1`). Statistics are computed over all variants on this chromosome. |
| `POPULATION` | string | yes | -- | Population identifier as stored in the graph |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--consequence` | string | none | Filter variants by VEP consequence type (e.g., `missense_variant`). When specified, summary statistics reflect only the matching variant subset. |
| `--pathway` | string | none | Filter variants to those in genes belonging to the named pathway |
| `--gene` | string | none | Filter variants to those within the named gene |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json` |

## Value

Returns a single row with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `pi` | float | Nucleotide diversity across the entire chromosome |
| `theta_w` | float | Watterson's theta for the chromosome |
| `tajima_d` | float | Tajima's D for the chromosome (genome-wide values are typically near zero for neutral populations; significantly negative values across a whole chromosome suggest population expansion) |
| `fay_wu_h` | float | Fay and Wu's H for the chromosome (requires ancestral allele annotation) |
| `mean_he` | float | Mean expected heterozygosity across all sites on the chromosome |
| `mean_ho` | float | Mean observed heterozygosity across all sites on the chromosome |
| `mean_fis` | float | Mean inbreeding coefficient (1 - mean_ho / mean_he) across the chromosome |
| `n_variants` | integer | Total number of variant sites on the chromosome (after any filtering) |
| `n_segregating` | integer | Number of segregating sites in the target population |
| `n_polarized` | integer | Number of variants with known ancestral allele |

## Details

### Algorithm

`pop-summary` operates identically to `graphpop diversity` in its statistical
computations, but with two key differences:

1. **Scope**: The entire chromosome is treated as a single region. All Variant
   nodes on the specified chromosome are included (subject to any conditioning
   filters).

2. **Persistence**: Results are written as properties on the Population node in
   the graph. For a population named `EUR` on `chr22`, the following properties
   are set:

   ```
   Population(name='EUR') {
       chr22_pi: 0.000123,
       chr22_theta_w: 0.000118,
       chr22_tajima_d: -0.31,
       chr22_fay_wu_h: -0.045,
       chr22_mean_he: 0.124,
       chr22_mean_ho: 0.121,
       chr22_mean_fis: 0.024,
       chr22_n_variants: 1103547,
       chr22_n_segregating: 523412,
       chr22_n_polarized: 498211
   }
   ```

   These properties are queryable in Cypher:

   ```cypher
   MATCH (p:Population)
   WHERE p.name = 'EUR'
   RETURN p.chr22_pi, p.chr22_tajima_d
   ```

### Conditioning behavior

When `--consequence`, `--pathway`, or `--gene` is specified, the summary
statistics reflect only the matching variants. This is useful for comparing
genome-wide diversity at functional versus neutral sites. For example, running
`pop-summary` with `--consequence synonymous_variant` and again with
`--consequence missense_variant` produces a direct comparison of neutral versus
functional diversity.

Note: when conditioning is applied, the persisted property names may include the
conditioning label to avoid overwriting unconditional summaries. Check the
Population node properties after running conditioned summaries.

### Performance

FAST PATH: O(V x K) where V is the total number of variants on the chromosome.
For a chromosome with ~1M variants, computation typically takes 2-5 seconds.
The persistence step (writing properties to the Population node) is fast
(single node update).

### Relationship to other commands

| Command | Scope | Persistence | Use case |
|---------|-------|-------------|----------|
| `diversity` | Arbitrary region | None | Ad-hoc queries, scripting |
| `pop-summary` | Whole chromosome | Population node | Definitive chromosome summary |
| `genome-scan` | Sliding windows | GenomicWindow nodes | Spatial variation along chromosome |

For a complete analysis, run `pop-summary` for each population and chromosome,
then `genome-scan` for windowed resolution. The `graphpop run-all` command
automates this workflow across all populations and chromosomes.

## Examples

```bash
# Basic: chromosome-wide summary for EUR
graphpop pop-summary chr22 EUR

# Summary with output to file
graphpop pop-summary chr22 EUR -o pop_summary_eur_chr22.tsv

# Conditioned on missense variants
graphpop pop-summary chr22 CEU \
    --consequence missense_variant \
    -o pop_summary_missense.tsv

# Pathway-conditioned summary in JSON format
graphpop pop-summary chr1 GJ-tmp \
    --pathway "Starch and sucrose metabolism" \
    --format json \
    -o pop_summary_starch.json

# Run for all major 1000 Genomes populations (scripted)
for pop in EUR EAS AFR SAS AMR; do
    graphpop pop-summary chr22 "$pop" -o "pop_summary_${pop}_chr22.tsv"
done
```

## See Also

- `graphpop diversity` -- region-level diversity statistics (without persistence)
- `graphpop genome-scan` -- sliding-window scan with GenomicWindow node persistence
- `graphpop run-all` -- automated full-genome analysis across all populations and chromosomes
- `graphpop aggregate` -- generate summary tables from run-all results
