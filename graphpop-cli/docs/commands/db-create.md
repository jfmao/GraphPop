# graphpop db create

## Title

Create a New Neo4j Database

## Description

Creates a new named database in the Neo4j server. This command requires
Neo4j Enterprise Edition; on Community Edition it will fail with a descriptive
error message.

## Usage

```
graphpop db create NAME
```

## Arguments

| Argument | Type | Description |
|---|---|---|
| `NAME` | TEXT | The name for the new database. Must be a valid Neo4j database name (alphanumeric and hyphens, starting with a letter). |

## Options

This command uses only the global connection options (`--uri`, `--user`,
`--password`, `--database`, `--config`).

## Value

No return value. Prints a confirmation message and suggests switching to the
new database. Exits with code 1 if the creation fails.

## Details

Runs `CREATE DATABASE <name> IF NOT EXISTS` against the Neo4j `system` database.
This requires:

- **Neo4j Enterprise Edition** -- Community Edition does not support the
  `CREATE DATABASE` command. With Community, use `neo4j-admin` directly or
  create databases through the import pipeline (`graphpop import --database`).
- **Admin privileges** -- The connecting user must have database management
  permissions.

After creation, the database starts in an empty state. To populate it, use
`graphpop import --database <name>` or `graphpop load --database <name>`.

The `IF NOT EXISTS` clause means running this command on an already-existing
database is a no-op.

Note: `graphpop import` creates the database as part of its bulk-import step,
so you do not need to run `db create` before importing.

## Examples

### Create a database for a new project

```bash
graphpop db create human_chr22
graphpop db switch human_chr22
```

### Create and import in one workflow

```bash
graphpop db create rice3k
graphpop import --vcf rice.vcf.gz --panel panel.txt --database rice3k
```

## See Also

- `graphpop db list` -- List all databases.
- `graphpop db switch` -- Set the active database.
- `graphpop db drop` -- Drop a database.
- `graphpop import` -- Import VCF data (creates the database automatically).
