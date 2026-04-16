# Look up genes and variant counts in a pathway

## Description

`graphpop lookup pathway` queries the graph for all Pathway nodes whose name
contains the given substring and returns one row per member gene per matching
pathway. Each row carries the pathway identifier, gene coordinates, and the
number of variants with consequences in that gene.

Pathway name matching is performed with a case-insensitive CONTAINS predicate,
so a short substring such as `"starch"` matches `"Starch biosynthesis"`,
`"Starch degradation"`, and any other pathway whose name contains that string.
Quoting the argument is recommended when it contains spaces.

Use `lookup pathway` to explore which genes belong to a biological pathway of
interest before running selection scans or `neighbors` traversals.

## Usage

```
graphpop lookup pathway PW_NAME [OPTIONS]
```

## Arguments

| Name | Type | Required | Description |
|------|------|----------|-------------|
| `PW_NAME` | string | yes | Pathway name or substring (matched with CONTAINS). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout. |

## Value

Returns one row per gene in each matching pathway:

| Column | Type | Description |
|--------|------|-------------|
| `pathway` | string | Full pathway name. |
| `pathway_id` | string | Pathway identifier. |
| `gene` | string | Gene symbol. |
| `gene_id` | string | Gene identifier (e.g., Ensembl ID). |
| `chr` | string | Chromosome. |
| `gene_start` | integer | Gene start position. |
| `gene_end` | integer | Gene end position. |
| `variant_count` | integer | Number of variants with consequences in this gene. |

## Details

### Substring matching

The CONTAINS match is case-insensitive and anchored to no particular position
in the pathway name. Providing a very short or common substring (e.g., `"bio"`)
may match a large number of pathways and produce a correspondingly large result
set. Use a more specific substring or pipe the output through a filter if the
result is unexpectedly large.

### Multiple matching pathways

When `PW_NAME` matches more than one pathway, the output contains rows from all
matching pathways. The `pathway` and `pathway_id` columns identify which pathway
each row belongs to.

### variant_count semantics

`variant_count` is the number of Variant nodes connected to the gene via
`HAS_CONSEQUENCE` edges. Variants with no consequence annotation on the gene
are not counted. A gene with `variant_count = 0` has been imported but has no
variant consequence edges in the current graph.

## Examples

```bash
# Find all genes in pathways related to cardiac repolarization
graphpop lookup pathway "Cardiac repolarization"

# Find starch-related pathway genes and save to file
graphpop lookup pathway "starch" -o starch_pathway.tsv

# JSON output for programmatic use
graphpop lookup pathway "insulin signaling" --format json
```

## See Also

- `graphpop lookup gene` -- detailed per-variant annotation for a single gene
- `graphpop lookup variant` -- full annotation for a single variant
- `graphpop lookup region` -- genes and variant counts in a genomic region
- `graphpop lookup` -- parent command overview
- `graphpop neighbors` -- explore graph neighborhood of a gene
