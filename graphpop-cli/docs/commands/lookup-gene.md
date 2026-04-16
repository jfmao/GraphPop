# Look up all variants, pathways, and selection scores for a gene

## Description

`graphpop lookup gene` queries the graph for a single Gene node and returns
one row per variant in that gene. Each row carries gene-level annotation
(pathways, Ensembl ID) together with variant-level annotation (position,
alleles, persisted iHS and XP-EHH scores). This is the primary entry point
for inspecting a gene of interest in full detail without writing Cypher.

Both gene symbols (e.g., `KCNE1`, `GW5`) and Ensembl-style gene IDs (e.g.,
`ENSG00000180509`, `Os01g0100100`) are accepted. The lookup resolves the
identifier against the `gene_symbol` and `gene_id` properties on Gene nodes.

Use `lookup gene` after a genome scan or `rank-genes` run to inspect a
top-ranked candidate in context.

## Usage

```
graphpop lookup gene GENE_NAME [OPTIONS]
```

## Arguments

| Name | Type | Required | Description |
|------|------|----------|-------------|
| `GENE_NAME` | string | yes | Gene symbol (e.g., `KCNE1`) or Ensembl ID (e.g., `ENSG00000180509`). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file instead of stdout. |

## Value

Returns one row per variant in the gene:

| Column | Type | Description |
|--------|------|-------------|
| `gene` | string | Gene symbol. |
| `gene_id` | string | Gene identifier (e.g., Ensembl ID). |
| `chr` | string | Chromosome. |
| `start` | integer | Gene start position. |
| `end` | integer | Gene end position. |
| `variant_id` | string | Variant identifier in `chr:pos:ref:alt` format. |
| `pos` | integer | Variant genomic position. |
| `ref` | string | Reference allele. |
| `alt` | string | Alternate allele. |
| `pathways` | list | Pathway names the gene belongs to (semicolon-separated in TSV/CSV). |
| `ihs_scores` | list | Persisted iHS scores as `key=value` strings (e.g., `ihs_EUR=2.31`). |
| `xpehh_scores` | list | Persisted XP-EHH scores as `key=value` strings (e.g., `xpehh_EUR_YRI=1.87`). |

## Details

### Gene identifier resolution

The command first attempts to match `GENE_NAME` against the `gene_symbol`
property on Gene nodes (case-sensitive). If no match is found it falls back
to matching against `gene_id`. This means Ensembl IDs can be used directly
without any flag.

### Selection score columns

`ihs_scores` and `xpehh_scores` collect all properties on the Variant node
whose names start with `ihs_` or `xpehh_` respectively. These are persisted
by `graphpop ihs` and `graphpop xpehh` and follow the naming convention
`ihs_{pop}` and `xpehh_{pop1}_{pop2}`. If no scores have been persisted, the
columns are present but empty.

### List encoding

In TSV and CSV output, list columns (`pathways`, `ihs_scores`, `xpehh_scores`)
are encoded as semicolon-separated strings. In JSON output they are arrays.

## Examples

```bash
# Per-variant detail for a human cardiac gene
graphpop lookup gene KCNE1

# Save full rice gene annotation to a file
graphpop lookup gene GW5 -o gw5_info.tsv

# Look up by Ensembl ID and output as JSON
graphpop lookup gene ENSG00000180509 --format json
```

## See Also

- `graphpop lookup pathway` -- show genes and variant counts in a pathway
- `graphpop lookup variant` -- full annotation for a single variant
- `graphpop lookup region` -- genes and variant counts in a genomic region
- `graphpop lookup` -- parent command overview
- `graphpop neighbors` -- explore graph neighborhood of a gene
