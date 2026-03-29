# graphpop db list

## Title

List All Neo4j Databases

## Description

Lists all databases in the Neo4j server with their status, size on disk, and
whether each is the currently active GraphPop database. The active database is
marked with an asterisk (`*`).

## Usage

```
graphpop db list
```

## Arguments

This command takes no positional arguments.

## Options

This command uses only the global connection options (`--uri`, `--user`,
`--password`, `--database`, `--config`).

## Value

Prints a formatted table to stdout:

```
Database                  Status       Size            Active
------------------------------------------------------------
neo4j                     online       4.2 GB
rice3k                    online       6.8 GB           *
system                    online       250.3 KB
```

The `Active` column marks the database that GraphPop commands will use by
default (as set in `~/.graphpop/config.yaml`).

## Details

Runs the Cypher command `SHOW DATABASES` to retrieve database metadata. On
Neo4j Community Edition, the `sizeOnDisk` field may not be available; the
command gracefully falls back to displaying name and status only.

The active database is read from the `database` key in
`~/.graphpop/config.yaml`. This is the database that all GraphPop procedure
calls will target unless overridden with `--database`.

Neo4j Community Edition supports only one user database (`neo4j`) plus the
`system` database. Neo4j Enterprise Edition supports multiple user databases.

## Examples

### List all databases

```bash
graphpop db list
```

### List databases on a remote server

```bash
graphpop --uri bolt://server:7687 --password secret db list
```

### Pipe to grep for a specific database

```bash
graphpop db list | grep rice
```

## See Also

- `graphpop db create` -- Create a new database (Enterprise only).
- `graphpop db switch` -- Change the active database.
- `graphpop db info` -- Show detailed node/edge counts for the current database.
