# Look up genes, pathways, variants, and genomic regions

## Description

`graphpop lookup` queries the graph for annotation and summary information about
a gene, pathway, variant, or genomic region. It is the fastest way to answer
"what does the graph know about X?" without writing Cypher. Four subcommands
are available:

| Subcommand | Query target |
|------------|-------------|
| `gene` | Gene node: variant count, consequence breakdown, pathways, persisted selection statistics |
| `pathway` | Pathway node: member gene list, mean Fst across gene bodies, pathway-level statistics |
| `variant` | Single Variant node: full annotation including allele frequencies, consequences, gene membership, persisted scores |
| `region` | Genomic interval: genes overlapping the region, windowed statistics, variant density |

All subcommands operate as read-only queries against the existing graph. They do
not compute new statistics; they report what has already been persisted by prior
commands such as `graphpop genome-scan`, `graphpop ihs --persist`, or
`graphpop pop-summary`.

## Usage

```
graphpop lookup gene GENE_NAME [OPTIONS]
graphpop lookup pathway PATHWAY_NAME [OPTIONS]
graphpop lookup variant VARIANT_ID [OPTIONS]
graphpop lookup region CHR START END [OPTIONS]
```

## Arguments

| Name | Type | Required | Description |
|------|------|----------|-------------|
| `GENE_NAME` | string | yes (gene) | Gene symbol (e.g., `BRCA2`, `Os01g0100100`). Case-sensitive. |
| `PATHWAY_NAME` | string | yes (pathway) | Pathway name as stored on Pathway nodes (e.g., `"Starch and sucrose metabolism"`). |
| `VARIANT_ID` | string | yes (variant) | Variant identifier in `chr:pos:ref:alt` format (e.g., `chr22:17080794:G:A`). |
| `CHR` | string | yes (region) | Chromosome name. |
| `START` | integer | yes (region) | Start position (1-based, inclusive). |
| `END` | integer | yes (region) | End position (1-based, inclusive). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pop` | string | all | Restrict reported statistics to a specific population. |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout. |
| `--include-variants` | flag | false | (gene/pathway/region) Include per-variant detail rows in addition to the summary. |
| `--include-scores` | flag | false | Include all persisted selection scores (iHS, XP-EHH, nSL, H12) on each variant row. |

## Value

### graphpop lookup gene

Returns a summary block followed by optional variant detail:

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Gene symbol. |
| `chr` | string | Chromosome. |
| `start` | integer | Gene start position. |
| `end` | integer | Gene end position. |
| `n_variants` | integer | Total variants within gene boundaries. |
| `n_missense` | integer | Count of missense variants. |
| `n_synonymous` | integer | Count of synonymous variants. |
| `n_high_impact` | integer | Count of HIGH-impact variants (stop_gained, frameshift, splice). |
| `pathways` | string | Comma-separated list of pathway names. |
| `mean_fst_{pop_pair}` | float | Mean Fst across gene variants (one column per persisted population pair). |
| `max_abs_ihs_{pop}` | float | Maximum |iHS| within the gene (one column per population with persisted scores). |

### graphpop lookup pathway

| Column | Type | Description |
|--------|------|-------------|
| `pathway` | string | Pathway name. |
| `n_genes` | integer | Number of member genes. |
| `genes` | string | Comma-separated gene symbols. |
| `mean_fst_{pop_pair}` | float | Mean Fst across all pathway genes. |
| `mean_pi_{pop}` | float | Mean nucleotide diversity across pathway genes. |

### graphpop lookup variant

| Column | Type | Description |
|--------|------|-------------|
| `variantId` | string | Variant identifier. |
| `pos` | integer | Genomic position. |
| `ref` | string | Reference allele. |
| `alt` | string | Alternate allele. |
| `consequence` | string | VEP consequence type. |
| `gene` | string | Gene symbol (if in a gene). |
| `af_{pop}` | float | Allele frequency per population. |
| `ihs_{pop}` | float | Persisted iHS score (if available). |
| `xpehh_{pop_pair}` | float | Persisted XP-EHH score (if available). |

### graphpop lookup region

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Gene symbol. |
| `gene_start` | integer | Gene start. |
| `gene_end` | integer | Gene end. |
| `n_variants` | integer | Variants in this gene within the queried region. |
| `mean_fst` | float | Mean Fst (if persisted). |
| `max_abs_ihs` | float | Max |iHS| (if persisted). |

## Details

### Resolution order

For `lookup gene`, the command first matches by exact gene symbol on Gene nodes.
If no match is found, it attempts a case-insensitive search and reports the
closest match with a warning. For rice datasets, both MSU and RAP identifiers
are supported (e.g., `LOC_Os01g01010` or `Os01g0100100`).

### Performance

All lookups are index-backed. Gene and pathway lookups resolve in under 100 ms.
Region lookups scale linearly with the number of genes in the interval but
remain fast for typical query sizes (under 1 second for a 5 Mb region). The
`--include-variants` flag may increase response time for large genes with
thousands of variants.

## Examples

```bash
# Look up a gene: variant count, pathways, selection stats
graphpop lookup gene BRCA2 --pop EUR --format json

# Look up a rice gene with full variant detail
graphpop lookup gene Os01g0100100 --pop GJ-tmp --include-variants \
    --include-scores -o gene_detail.tsv

# Look up a pathway and its gene members
graphpop lookup pathway "Starch and sucrose metabolism" --pop GJ-tmp

# Look up a single variant with all annotations
graphpop lookup variant chr22:17080794:G:A --format json

# Look up all genes in a 2 Mb region
graphpop lookup region chr22 20000000 22000000 --pop CEU -o region_genes.tsv
```

## See Also

- `graphpop neighbors` -- explore graph neighborhood around a gene
- `graphpop rank-genes` -- rank genes by composite selection score
- `graphpop converge` -- find regions where multiple statistics exceed thresholds
- `graphpop filter` -- post-hoc filtering of persisted statistics
- `graphpop inventory` -- check what data and statistics are loaded in the graph
