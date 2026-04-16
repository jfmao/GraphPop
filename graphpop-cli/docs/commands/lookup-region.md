# Look up genes and variant counts in a genomic region

## Description

`graphpop lookup region` queries the graph for all Gene nodes that overlap a
specified genomic interval and returns one row per gene with a summary of the
variants in that gene within the queried region. Variants that fall inside the
region but outside any annotated gene are aggregated into a single
`intergenic` row.

This command provides a rapid overview of what the graph contains for any
genomic window: how many genes are covered, how many variants each gene
contributes, and the positional span of those variants. It is a useful first
step before running a focused scan with `graphpop genome-scan` or inspecting
individual genes with `graphpop lookup gene`.

## Usage

```
graphpop lookup region CHR START END [OPTIONS]
```

## Arguments

| Name | Type | Required | Description |
|------|------|----------|-------------|
| `CHR` | string | yes | Chromosome name (must match the naming convention used at import, e.g., `chr22` or `Chr01`). |
| `START` | integer | yes | Start position of the region (1-based, inclusive). |
| `END` | integer | yes | End position of the region (1-based, inclusive). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout. |

## Value

Returns one row per gene (or intergenic block) overlapping the region:

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Gene symbol, or `intergenic` for variants with no gene annotation. |
| `gene_id` | string | Gene identifier (e.g., Ensembl ID). Empty for intergenic rows. |
| `gene_start` | integer | Gene start position. Empty for intergenic rows. |
| `gene_end` | integer | Gene end position. Empty for intergenic rows. |
| `variant_count` | integer | Number of variants in this gene within the queried region. |
| `min_pos` | integer | Minimum variant position within this gene/block in the queried region. |
| `max_pos` | integer | Maximum variant position within this gene/block in the queried region. |

## Details

### Region overlap semantics

A gene is included if its genomic span (`gene_start` to `gene_end`) overlaps
the queried interval `[START, END]`. Variants are counted only when their
`pos` falls within `[START, END]`, so a gene that extends beyond the region
boundary will have a `variant_count` that reflects only the variants inside
the query window, and `min_pos`/`max_pos` are likewise clipped to the window.

### Intergenic variants

Variants whose `pos` falls within `[START, END]` but that have no
`HAS_CONSEQUENCE` edge to any Gene node are aggregated into a single row with
`gene = 'intergenic'`. If all variants in the region are annotated to genes,
no intergenic row is emitted.

### Chromosome naming

The `CHR` argument must match the chromosome naming convention used when the
data was imported (e.g., `chr22` for 1000 Genomes data, `Chr01` for Rice 3K
data). Use `graphpop inventory` to list the chromosomes present in the graph
if you are unsure of the convention.

## Examples

```bash
# Summarize genes in a 600 kb window on chromosome 6 (rice GW5 region)
graphpop lookup region chr6 9000000 9600000

# Summarize a 1 Mb window on human chromosome 22 and save to file
graphpop lookup region chr22 16000000 17000000 -o region.tsv

# JSON output for programmatic use
graphpop lookup region chr22 16000000 17000000 --format json
```

## See Also

- `graphpop lookup gene` -- detailed per-variant annotation for a single gene
- `graphpop lookup pathway` -- show genes and variant counts in a pathway
- `graphpop lookup variant` -- full annotation for a single variant
- `graphpop lookup` -- parent command overview
- `graphpop neighbors` -- explore graph neighborhood of a gene
