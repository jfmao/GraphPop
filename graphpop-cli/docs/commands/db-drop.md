# graphpop db drop

## Title

Drop a Neo4j Database

## Description

Permanently deletes a named database from the Neo4j server. Requires
confirmation before proceeding. This command requires Neo4j Enterprise Edition.

## Usage

```
graphpop db drop NAME [OPTIONS]
```

## Arguments

| Argument | Type | Description |
|---|---|---|
| `NAME` | TEXT | The name of the database to drop. |

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--force` | FLAG | `false` | Skip the interactive confirmation prompt. Use with caution. |

## Value

No return value. Prints a confirmation message on success. Exits with code 1 if
the drop fails or is cancelled.

## Details

Runs `DROP DATABASE <name> IF EXISTS` against the Neo4j `system` database.

**Safety guards:**

- The `neo4j` and `system` databases are protected and cannot be dropped. The
  command exits with an error if you attempt to drop either.
- An interactive confirmation prompt is shown unless `--force` is passed.
- The `IF EXISTS` clause means dropping a non-existent database is a no-op.

**Requirements:**

- Neo4j Enterprise Edition (Community does not support `DROP DATABASE`).
- Admin privileges on the Neo4j server.

**This operation is irreversible.** All data in the database is permanently
deleted. Use `graphpop dump` to create a backup before dropping.

If you are using Community Edition and need to remove a database, stop Neo4j
and delete the database directory manually from `<neo4j-home>/data/databases/`.

## Examples

### Drop with confirmation prompt

```bash
graphpop db drop old_analysis
# Drop database 'old_analysis'? This cannot be undone [y/N]: y
# Database 'old_analysis' dropped.
```

### Drop without confirmation (scripting)

```bash
graphpop db drop temp_project --force
```

### Dump before dropping

```bash
graphpop stop
graphpop dump --database old_analysis -o old_analysis_backup.dump
graphpop start
graphpop db drop old_analysis --force
```

## See Also

- `graphpop db list` -- List all databases.
- `graphpop db create` -- Create a new database.
- `graphpop dump` -- Back up a database before dropping.
