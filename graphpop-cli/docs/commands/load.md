# graphpop load

## Title

Restore a Database from a Dump File

## Description

Loads a previously dumped Neo4j database from a `.dump` file. The database can
be restored under a different name than the original. After loading, the
GraphPop configuration is updated to use the restored database.

## Usage

```
graphpop load [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--dump-file` | PATH | *(required)* | Path to the `.dump` file to restore. Must exist. |
| `--database` | TEXT | *(required)* | Name for the restored database. Can differ from the original name. |
| `--neo4j-home` | PATH | *(from config)* | Neo4j installation directory. |
| `--overwrite` | FLAG | `false` | Overwrite the database if it already exists. Without this flag, loading into an existing database will fail. |

## Value

No return value. Side effects:

- Database restored from the dump file into `<neo4j-home>/data/databases/<database>/`.
- `~/.graphpop/config.yaml` updated: `database` set to the restored database name.

## Details

Runs `neo4j-admin database load --from-path=<dir> <database>`.

**Prerequisites:**

- **Neo4j must be stopped** before loading. Run `graphpop stop` first.
- The `neo4j-admin` binary must exist at `<neo4j-home>/bin/neo4j-admin`.

If the target database already exists and `--overwrite` is not set, the command
will fail with an error suggesting the flag. With `--overwrite`, the existing
database is replaced entirely.

The database name does not need to match the original. This allows you to load
the same dump under multiple names for comparison or testing.

After loading, start Neo4j and verify the database:

```bash
graphpop start
graphpop db info
graphpop validate
```

## Examples

### Load a shared database

```bash
graphpop stop
graphpop load --dump-file rice3k_v1.dump --database rice3k
graphpop start
graphpop db info
```

### Load under a different name

```bash
graphpop stop
graphpop load --dump-file colleague_analysis.dump --database review_copy
graphpop start
graphpop db switch review_copy
```

### Overwrite an existing database

```bash
graphpop stop
graphpop load --dump-file rice3k_v2.dump --database rice3k --overwrite
graphpop start
```

### Load and immediately analyze

```bash
graphpop stop
graphpop load --dump-file human_chr22.dump --database human_chr22
graphpop start
graphpop diversity chr22 1 50818468 EUR
```

## See Also

- `graphpop dump` -- Create a dump file from an existing database.
- `graphpop stop` -- Stop Neo4j before loading (required).
- `graphpop start` -- Start Neo4j after loading.
- `graphpop db info` -- Verify the loaded database.
- `graphpop validate` -- Check integrity of the loaded database.
