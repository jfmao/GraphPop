# Extract Variant nodes with flexible annotation and frequency filters

## Description

`graphpop extract variants` queries Variant nodes from the graph and returns
tabular records with the fields you choose. Filters span chromosome, genomic
region, population-specific allele frequency, VEP consequence type, pathway
membership, and gene association. Any combination of filters can be applied
simultaneously; all conditions are ANDed.

The Cypher query is built dynamically: only the clauses required by the options
you supply are included. Annotation joins (consequence, pathway, gene) add graph
traversal steps automatically; no pre-join or preprocessing is required.

When `--pop` is given alongside `--min-af` or `--max-af`, allele frequency
filtering uses the population-specific value from the per-population `af[]`
array on each Variant node. Without `--pop`, AF filters are not available. When
`af` appears in `--fields` and `--pop` is specified, the output column is named
`af_{population}` to make the source population explicit.

This command does not compute new statistics. It extracts what is already stored
in the graph. For computing summary statistics over a filtered variant set, see
`graphpop diversity` or `graphpop sfs`.

## Usage

```
graphpop extract variants [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--chr` | string | none | Restrict to a single chromosome (e.g., `chr22`, `chr1`). |
| `--start` | integer | none | Start position of the region (1-based, inclusive). Requires `--chr`. |
| `--end` | integer | none | End position of the region (1-based, inclusive). Requires `--chr`. |
| `--pop` | string | none | Population name for AF lookup and AF-based filtering. |
| `--min-af` | float | none | Minimum allele frequency in `--pop`. Requires `--pop`. |
| `--max-af` | float | none | Maximum allele frequency in `--pop`. Requires `--pop`. |
| `--consequence` | string | none | Filter by VEP consequence type (e.g., `missense_variant`, `synonymous_variant`, `stop_gained`). |
| `--pathway` | string | none | Filter to variants in genes belonging to a biological pathway (substring match against pathway name). |
| `--gene` | string | none | Filter to variants within a specific gene (gene symbol or gene ID). |
| `--fields` | string | `variantId,pos,ref,alt,af` | Comma-separated list of Variant node properties to include in the output. |
| `--limit` | integer | `10000` | Maximum number of rows to return. |
| `-o`, `--output` | path | stdout | Output file path. |
| `--format` | choice | `tsv` | Output format: `tsv`, `csv`, or `json`. |

## Details

### Field resolution

The `--fields` option accepts any Variant node property name stored in the graph
(e.g., `variantId`, `pos`, `chr`, `ref`, `alt`, `is_polarized`,
`ancestral_allele`, `vep_consequence`). When `--pop` is specified and `af`
appears in the field list, the value is resolved to the population-specific
allele frequency from the indexed `af[]` array, and the output column is renamed
`af_{population}`. Properties not present on a given Variant node are returned
as empty (TSV/CSV) or `null` (JSON).

### Annotation joins

Specifying `--consequence` adds a traversal over the VEP annotation attached to
each Variant node. Specifying `--pathway` or `--gene` adds graph traversal steps
along `Variant→Gene→Pathway` edges. These joins are performed inside Neo4j and
benefit from relationship indices. Combining multiple annotation filters
restricts results to variants that satisfy all conditions simultaneously.

### Limit and pagination

The default `--limit 10000` guards against accidentally returning millions of
rows for an unfiltered query. Increase the limit explicitly when extracting
large regions, or combine with region and annotation filters to narrow the
result set before raising the limit.

### Performance

Queries filtered by chromosome and region use the positional index on Variant
nodes and typically complete in under 2 seconds for regions spanning 1 M
variants. Annotation joins (consequence, pathway, gene) add 1-5 seconds
depending on the number of matching variants and pathway traversal depth.
Unfiltered whole-genome queries against large datasets should include a
reasonable `--limit` or be replaced with `graphpop dump` for bulk export.

## Examples

```bash
# Missense variants on chr22 for EUR with AF > 0.01
graphpop extract variants --chr chr22 --pop EUR \
    --consequence missense_variant --min-af 0.01 \
    -o missense_chr22_eur.tsv

# Region extraction with custom field selection
graphpop extract variants --chr chr22 --start 16000000 --end 17000000 \
    --pop EUR --fields variantId,pos,ref,alt,af_EUR,vep_consequence \
    -o region_chr22_16-17Mb.tsv

# Gene-specific extraction in JSON format
graphpop extract variants --gene KCNE1 --pop EUR \
    --fields variantId,pos,ref,alt,af,vep_consequence \
    --format json -o kcne1_variants.json

# Rare variants in a pathway across all chromosomes
graphpop extract variants --pathway "DNA repair" \
    --pop AFR --max-af 0.01 --limit 50000 \
    -o dna_repair_rare_afr.tsv

# Rice: coding variants in a gene with frequency filter
graphpop extract variants --gene OsWRKY45 --pop GJ-tmp \
    --consequence missense_variant --min-af 0.05 \
    -o oswrky45_missense.tsv
```

## See Also

- `graphpop extract samples` -- export Sample nodes for a population
- `graphpop extract genotypes` -- export a sample x variant genotype matrix for a region
- `graphpop extract` -- overview of all extract subcommands
- `graphpop export-bed` -- export high-signal regions in BED format
- `graphpop filter` -- filter and display persisted statistics with threshold queries
- `graphpop lookup` -- look up a single variant or gene record
