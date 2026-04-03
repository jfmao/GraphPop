# Look up genes, pathways, variants, and genomic regions

## Description

`graphpop lookup` queries the graph for annotation information about a gene,
pathway, variant, or genomic region. It is the fastest way to answer "what does
the graph know about X?" without writing Cypher. Four subcommands are available:

| Subcommand | Query target |
|------------|-------------|
| `gene` | Gene node: per-variant detail with consequences, pathways, persisted iHS/XP-EHH scores |
| `pathway` | Pathway node: member genes with variant counts (substring match on pathway name) |
| `variant` | Single Variant node: full property dump including allele frequencies, gene, pathways |
| `region` | Genomic interval: genes overlapping the region with variant counts |

All subcommands operate as read-only queries against the existing graph.

## Usage

```
graphpop lookup gene GENE_NAME [OPTIONS]
graphpop lookup pathway PW_NAME [OPTIONS]
graphpop lookup variant VAR_ID [OPTIONS]
graphpop lookup region CHR START END [OPTIONS]
```

## Arguments

| Name | Type | Required | Description |
|------|------|----------|-------------|
| `GENE_NAME` | string | yes (gene) | Gene symbol (e.g., `KCNE1`) or gene ID (e.g., `ENSG00000180509`). |
| `PW_NAME` | string | yes (pathway) | Pathway name substring (e.g., `"starch"`). Matched with CONTAINS. |
| `VAR_ID` | string | yes (variant) | Variant identifier in `chr:pos:ref:alt` format (e.g., `chr22:16050075:A:G`). |
| `CHR` | string | yes (region) | Chromosome name. |
| `START` | integer | yes (region) | Start position (1-based, inclusive). |
| `END` | integer | yes (region) | End position (1-based, inclusive). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout. |

## Value

### graphpop lookup gene

Returns one row per variant in the gene, with gene-level and variant-level annotation:

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Gene symbol. |
| `gene_id` | string | Gene identifier (e.g., Ensembl ID). |
| `chr` | string | Chromosome. |
| `start` | integer | Gene start position. |
| `end` | integer | Gene end position. |
| `variant_id` | string | Variant identifier (chr:pos:ref:alt). |
| `pos` | integer | Variant genomic position. |
| `ref` | string | Reference allele. |
| `alt` | string | Alternate allele. |
| `pathways` | list | Pathway names the gene belongs to. |
| `ihs_scores` | list | Persisted iHS scores as `key=value` strings (e.g., `ihs_EUR=2.31`). |
| `xpehh_scores` | list | Persisted XP-EHH scores as `key=value` strings. |

### graphpop lookup pathway

Returns one row per gene in each matching pathway:

| Column | Type | Description |
|--------|------|-------------|
| `pathway` | string | Pathway name. |
| `pathway_id` | string | Pathway identifier. |
| `gene` | string | Gene symbol. |
| `gene_id` | string | Gene identifier. |
| `chr` | string | Chromosome. |
| `gene_start` | integer | Gene start position. |
| `gene_end` | integer | Gene end position. |
| `variant_count` | integer | Number of variants with consequences in this gene. |

### graphpop lookup variant

Returns all properties stored on the Variant node, flattened into columns.
Additional annotation columns are appended:

| Column | Type | Description |
|--------|------|-------------|
| *(all variant properties)* | various | Every property on the Variant node (variantId, pos, ref, alt, chr, pop_ids, ac, an, af, etc.). |
| `gene` | string | Gene symbol (if the variant has a HAS_CONSEQUENCE edge). |
| `gene_id` | string | Gene identifier. |
| `pathways` | list | Pathway names via the gene. |

### graphpop lookup region

Returns one row per gene (or intergenic block) in the region:

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Gene symbol, or `intergenic` if no gene annotation. |
| `gene_id` | string | Gene identifier. |
| `gene_start` | integer | Gene start position. |
| `gene_end` | integer | Gene end position. |
| `variant_count` | integer | Number of variants in this gene within the queried region. |
| `min_pos` | integer | Minimum variant position in the gene/region. |
| `max_pos` | integer | Maximum variant position in the gene/region. |

## Examples

```bash
# Look up a gene: per-variant detail with pathways and scores
graphpop lookup gene KCNE1
graphpop lookup gene GW5 -o gw5_info.tsv
graphpop lookup gene ENSG00000180509 --format json

# Look up a pathway by substring match
graphpop lookup pathway "Cardiac repolarization"
graphpop lookup pathway "starch" -o starch_pathway.tsv

# Look up a single variant with all annotations
graphpop lookup variant chr22:16050075:A:G --format json

# Look up all genes in a genomic region
graphpop lookup region chr6 9000000 9600000
graphpop lookup region chr22 16000000 17000000 -o region.tsv
```

## See Also

- `graphpop neighbors` -- explore graph neighborhood around a gene
- `graphpop rank-genes` -- rank genes by composite selection score
- `graphpop converge` -- find regions where multiple statistics exceed thresholds
- `graphpop filter` -- post-hoc filtering of persisted statistics
- `graphpop inventory` -- check what data and statistics are loaded in the graph
