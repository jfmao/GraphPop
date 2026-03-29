# Show graph contents and computed statistics

## Description

`graphpop inventory` reports what data is currently loaded in the active Neo4j
database: node and relationship counts by type, population definitions,
chromosomes present, annotation layers (Gene, Pathway, GOTerm), and which
statistics have been persisted to the graph. It is the equivalent of a quick
health check and data manifest for the graph.

Use `inventory` at the start of an analysis session to confirm that import and
computation steps completed, before running downstream commands that depend on
persisted results. It is also useful for documenting the state of a database
before sharing or archiving.

## Usage

```
graphpop inventory [OPTIONS]
```

## Arguments

This command takes no positional arguments. It queries the currently active
database (see `graphpop db-switch`).

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--section` | string | `all` | Report a specific section only: `nodes`, `relationships`, `populations`, `chromosomes`, `annotations`, `statistics`. |
| `--format` | choice | `text` | Output format: `text` (human-readable summary), `tsv`, `json`. |
| `-o`, `--output` | path | stdout | Write output to a file. |

## Value

When `--format text` (default), the output is a structured summary:

```
=== GraphPop Inventory ===

Nodes:
  Variant          1,103,547
  Sample               2,504
  Gene                   546
  Pathway                 82
  GOTerm               1,247
  GenomicWindow        4,126
  Population              26

Relationships:
  CARRIES         12,847,231
  NEXT             1,103,546
  IN_GENE            987,312
  IN_PATHWAY           2,148
  HAS_GO_TERM          3,871
  LD                 487,209

Populations:
  AFR, AMR, EAS, EUR, SAS (5 super-populations)
  ACB, ASW, BEB, ... (26 populations)

Chromosomes:
  chr22 (1,103,547 variants)

Annotations:
  VEP consequences: loaded (987,312 annotated variants)
  Pathways: 82 loaded (Reactome)
  GO terms: 1,247 loaded

Persisted statistics:
  GenomicWindow: pi, theta_w, tajima_d, fst_wc, dxy (4,126 windows)
  Variant-level: ihs_EUR, ihs_EAS, xpehh_EUR_EAS (1,047,832 scored variants)
  Population: pi, theta_w, tajima_d (26 populations)
```

When `--format json`, returns a structured JSON object with the same sections.

When `--format tsv`, returns one row per item with columns: `section`, `key`,
`value`.

## Details

### What is counted

The inventory queries Neo4j label counts (which are maintained as metadata and
are instantaneous) for node totals. Relationship counts use type-specific
aggregation. Persisted statistics are detected by sampling Variant and
GenomicWindow nodes for known property prefixes (`ihs_`, `xpehh_`, `nsl_`,
`h12_`, `fst_`, `pi_`).

### Performance

The inventory query completes in under 2 seconds regardless of database size.
Node label counts are O(1) in Neo4j. Statistic detection uses a LIMIT-based
sampling strategy.

## Examples

```bash
# Full inventory of the active database
graphpop inventory

# Check only what statistics have been persisted
graphpop inventory --section statistics

# Machine-readable inventory for scripting
graphpop inventory --format json -o inventory.json

# Check node counts in TSV format
graphpop inventory --section nodes --format tsv
```

## See Also

- `graphpop db-info` -- database connection and configuration details
- `graphpop db-list` -- list available databases
- `graphpop status` -- Neo4j server status
- `graphpop validate` -- verify data integrity and correctness
