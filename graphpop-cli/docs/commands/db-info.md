# graphpop db info

## Title

Show Detailed Database Information

## Description

Displays comprehensive information about the currently active database,
including node counts by label, relationship counts by type, and the list of
installed GraphPop stored procedures.

## Usage

```
graphpop db info
```

## Arguments

This command takes no positional arguments.

## Options

This command uses only the global connection options (`--uri`, `--user`,
`--password`, `--database`, `--config`).

## Value

Prints a multi-section report to stdout:

```
Database: rice3k

Node counts:
  Variant                    2,103,456
  Sample                         3,000
  Population                        12
  Chromosome                        12
  Gene                          28,236
  Pathway                         312
  GenomicWindow                 4,800

Relationship counts:
  CARRIES                  18,456,789
  NEXT                      2,103,455
  ON_CHROMOSOME             2,103,456
  IN_POPULATION                 3,000
  HAS_CONSEQUENCE           1,856,234
  IN_PATHWAY                   45,678

GraphPop procedures (12):
  graphpop.diversity
  graphpop.divergence
  graphpop.garud_h
  graphpop.genome_scan
  graphpop.ihs
  graphpop.joint_sfs
  graphpop.ld
  graphpop.nsl
  graphpop.pop_summary
  graphpop.roh
  graphpop.sfs
  graphpop.xpehh
```

## Details

Runs three Cypher queries against the active database:

1. **Node counts** -- Uses `db.labels()` to enumerate all node labels and
   counts nodes for each. Common labels include `Variant`, `Sample`,
   `Population`, `Chromosome`, `Gene`, `Pathway`, `GOTerm`, and
   `GenomicWindow` (created by `genome-scan --persist`).

2. **Relationship counts** -- Uses `db.relationshipTypes()` to enumerate all
   relationship types and counts edges for each. Key types include `CARRIES`
   (sample genotypes), `NEXT` (variant positional chain), `ON_CHROMOSOME`,
   `IN_POPULATION`, `HAS_CONSEQUENCE`, `IN_PATHWAY`, and `LD` (linkage
   disequilibrium edges).

3. **Installed procedures** -- Uses `SHOW PROCEDURES` filtered to names
   starting with `graphpop`. A fully operational installation should list 12
   procedures (6 fast-path + 6 full-path).

If no GraphPop procedures are found, the command prints `NONE INSTALLED` as a
warning. This means the `graphpop-procedures.jar` plugin has not been deployed.

## Examples

### Check the current database

```bash
graphpop db info
```

### Check a specific database

```bash
graphpop --database human_chr22 db info
```

### Verify import completeness

```bash
graphpop import --vcf data.vcf.gz --panel panel.txt --database mydb
graphpop start
graphpop db info  # confirm node/edge counts match expectations
```

## See Also

- `graphpop db list` -- List all databases.
- `graphpop validate` -- Run integrity checks on the database.
- `graphpop status` -- Check server status and plugin deployment.
