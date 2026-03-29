# Rank genes by composite selection score

## Description

`graphpop rank-genes` produces a ranked list of genes for a target population,
scored by a composite metric that combines multiple dimensions of selection
evidence: maximum |iHS|, maximum |XP-EHH|, maximum H12, mean Fst across gene
body, and the count of HIGH-impact variants. This provides a single prioritized
gene list for follow-up investigation.

The command reads persisted statistics from the graph. Each component score is
rank-normalized to [0, 1] across all genes and then combined with equal weights
(or user-specified weights) into a composite score. Genes are returned sorted by
descending composite score.

Use `rank-genes` after running selection scans and persisting scores to get a
quick overview of the strongest selection candidates in a population.

## Usage

```
graphpop rank-genes --pop POPULATION [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop` | string | required | Target population for within-population statistics (iHS, H12, HIGH-impact count). |
| `--pop2` | string | none | Second population for pairwise statistics (XP-EHH, Fst). If omitted, pairwise components are excluded from the composite score. |
| `--chr` | string | all | Restrict to a single chromosome. |
| `--top` | integer | 50 | Number of top-ranked genes to return. Use `0` for all genes. |
| `--weights` | string | `1,1,1,1,1` | Comma-separated weights for the five components: ihs, xpehh, h12, fst, high_impact. Components with weight 0 are excluded. |
| `--min-variants` | integer | 5 | Exclude genes with fewer than this many variants (avoids noise from very small genes). |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file. |

## Value

Returns one row per gene, sorted by descending composite score:

| Column | Type | Description |
|--------|------|-------------|
| `rank` | integer | Rank position (1 = strongest candidate). |
| `gene` | string | Gene symbol. |
| `chr` | string | Chromosome. |
| `start` | integer | Gene start position. |
| `end` | integer | Gene end position. |
| `n_variants` | integer | Number of variants in the gene. |
| `max_abs_ihs` | float | Maximum |iHS| within the gene. |
| `max_abs_xpehh` | float | Maximum |XP-EHH| within the gene (if `--pop2` given). |
| `max_h12` | float | Maximum H12 within the gene. |
| `mean_fst` | float | Mean Fst across gene variants (if `--pop2` given). |
| `n_high_impact` | integer | Count of HIGH-impact variants (stop_gained, frameshift, splice_donor, splice_acceptor). |
| `composite_score` | float | Weighted sum of rank-normalized component scores [0, 1]. |

## Details

### Scoring algorithm

1. For each gene, the command queries persisted variant-level statistics and
   computes the five component values (max |iHS|, max |XP-EHH|, max H12, mean
   Fst, HIGH-impact count).
2. Each component is rank-normalized across all qualifying genes: rank / N_genes,
   producing a uniform [0, 1] distribution. This avoids scale differences between
   statistics.
3. The composite score is the weighted average: sum(w_i * rank_i) / sum(w_i).
4. Genes are sorted by composite score descending. Ties are broken by max |iHS|.

### Missing components

If a component statistic has not been persisted (e.g., iHS was never computed),
that component is assigned rank 0 for all genes and effectively contributes
nothing. A warning is printed listing which components are missing. To exclude a
component explicitly, set its weight to 0.

### Performance

The query aggregates statistics across Gene nodes and their variant members.
For a chromosome with ~500 genes and ~1M scored variants, typical runtime is
3-8 seconds. Genome-wide ranking across all chromosomes may take 30-60 seconds.

## Examples

```bash
# Top 20 selection candidates in Europeans
graphpop rank-genes --pop EUR --pop2 AFR --top 20

# Rank all genes on chr22 with custom weights (emphasize iHS and Fst)
graphpop rank-genes --pop EUR --pop2 EAS --chr chr22 --top 0 \
    --weights 2,1,1,2,0.5 -o ranked_chr22.tsv

# Rice: rank genes for indica vs japonica comparison
graphpop rank-genes --pop GJ-tmp --pop2 XI-1A --top 50 \
    --format json -o rice_ranked.json

# Exclude pairwise stats, rank by within-population evidence only
graphpop rank-genes --pop YRI --top 30 --weights 1,0,1,0,1

# Strict filter: only genes with at least 20 variants
graphpop rank-genes --pop CEU --pop2 YRI --min-variants 20 --top 100 \
    -o top100_genes.tsv
```

## See Also

- `graphpop lookup gene` -- detailed annotation for a single gene
- `graphpop neighbors` -- explore graph neighborhood of a gene
- `graphpop converge` -- find regions where multiple statistics converge
- `graphpop filter` -- single-statistic filtering of persisted scores
- `graphpop ihs` -- compute and persist iHS
- `graphpop xpehh` -- compute and persist XP-EHH
