# graphpop db switch

## Title

Set the Active Database in GraphPop Configuration

## Description

Changes the default database that all subsequent GraphPop commands will use.
Updates the `database` key in `~/.graphpop/config.yaml`. This is a local
configuration change only; it does not modify anything in Neo4j.

## Usage

```
graphpop db switch NAME
```

## Arguments

| Argument | Type | Description |
|---|---|---|
| `NAME` | TEXT | The name of the database to switch to. This should match an existing Neo4j database name. |

## Options

This command takes no additional options.

## Value

No return value. Prints the old and new database names:

```
Active database: neo4j -> rice3k
All graphpop commands will now use database 'rice3k'.
```

## Details

Reads `~/.graphpop/config.yaml`, updates the `database` field, and writes it
back. If the config file does not exist, it is created with just the database
setting.

This command does not verify that the named database actually exists in Neo4j.
If you switch to a non-existent database, subsequent commands will fail with a
connection error.

The database setting can also be overridden per-command using the `--database`
global option or the `GRAPHPOP_DATABASE` environment variable:

```bash
# Per-command override (does not change config)
graphpop --database human_chr22 diversity chr22 1 50000000 EUR

# Environment variable override
export GRAPHPOP_DATABASE=human_chr22
graphpop diversity chr22 1 50000000 EUR
```

## Examples

### Switch to a rice dataset

```bash
graphpop db switch rice3k
graphpop db info  # now shows rice3k node counts
```

### Switch back to the default

```bash
graphpop db switch neo4j
```

### Compare databases

```bash
graphpop db switch human_chr22
graphpop diversity chr22 1 50000000 EUR -o human.tsv

graphpop db switch rice3k
graphpop diversity chr1 1 43270923 GJ-tmp -o rice.tsv
```

## See Also

- `graphpop db list` -- List all databases with active marker.
- `graphpop db info` -- Show node/edge counts for the current database.
- `graphpop config show` -- Display the full configuration including active database.
