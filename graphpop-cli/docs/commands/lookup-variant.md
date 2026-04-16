# Look up full annotation for a single variant

## Description

`graphpop lookup variant` retrieves all information stored for a single Variant
node and returns it as a one-row table. Every property on the node is included
(allele frequencies, allele counts, persisted selection statistics), along with
gene annotation from `HAS_CONSEQUENCE` edges and pathway membership via the
annotated gene.

This is the fastest way to inspect a specific variant in full detail. It is
particularly useful for validating that allele frequency arrays were imported
correctly, checking which gene a variant has been annotated to, or verifying
that selection statistics were persisted after a scan.

## Usage

```
graphpop lookup variant VAR_ID [OPTIONS]
```

## Arguments

| Name | Type | Required | Description |
|------|------|----------|-------------|
| `VAR_ID` | string | yes | Variant identifier in `chr:pos:ref:alt` format (e.g., `chr22:16050075:A:G`). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout. |

## Value

Returns a single row with all Variant node properties flattened into columns,
followed by gene and pathway annotation:

| Column | Type | Description |
|--------|------|-------------|
| *(all variant properties)* | various | Every property stored on the Variant node: `variant_id`, `pos`, `ref`, `alt`, `chr`, `pop_ids`, `ac`, `an`, `af`, and any persisted statistics (e.g., `ihs_EUR`, `xpehh_EUR_YRI`). |
| `gene` | string | Gene symbol, if the variant has a `HAS_CONSEQUENCE` edge to a Gene node. |
| `gene_id` | string | Gene identifier. |
| `pathways` | list | Pathway names the annotated gene belongs to (semicolon-separated in TSV/CSV). |

## Details

### Property flattening

Neo4j list properties (`pop_ids`, `ac`, `an`, `af`) are rendered as
semicolon-separated strings in TSV and CSV output. In JSON output they are
arrays. The exact set of columns depends on which statistics have been persisted
in the current graph: a variant gains additional columns each time `graphpop ihs`,
`graphpop xpehh`, or other scan commands write scores back to Variant nodes.

### Variant ID format

The identifier must exactly match the `variant_id` property stored in the graph,
which follows the convention `chr:pos:ref:alt`. Chromosome names must use the
same prefix convention as the import (e.g., `chr22` vs `Chr01`). Use
`graphpop inventory` to confirm naming conventions used in the current graph.

### Missing annotation

If the variant has no `HAS_CONSEQUENCE` edge, `gene` and `gene_id` are empty
and `pathways` is an empty list. If no persisted selection statistics exist for
the variant, the corresponding columns are absent from the output.

## Examples

```bash
# Full annotation for a human variant on chromosome 22
graphpop lookup variant chr22:16050075:A:G

# Rice variant output as JSON (note Chr capitalization convention)
graphpop lookup variant Chr01:524650:C:T --format json

# Save variant details to a file
graphpop lookup variant chr6:88110059:G:A -o variant_detail.tsv
```

## See Also

- `graphpop lookup gene` -- detailed per-variant annotation for a single gene
- `graphpop lookup pathway` -- show genes and variant counts in a pathway
- `graphpop lookup region` -- genes and variant counts in a genomic region
- `graphpop lookup` -- parent command overview
- `graphpop neighbors` -- explore graph neighborhood of a gene
