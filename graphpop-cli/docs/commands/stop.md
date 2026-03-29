# graphpop stop

## Title

Stop the Neo4j Database Server

## Description

Stops a running Neo4j database server. The Neo4j installation location is
resolved from the GraphPop config file or the `--neo4j-home` option. Equivalent
to running `neo4j stop` from the Neo4j `bin/` directory.

## Usage

```
graphpop stop [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--neo4j-home` | PATH | *(from config)* | Neo4j installation directory. Overrides the value in `~/.graphpop/config.yaml`. |

## Value

No return value. Prints the Neo4j shutdown message to stdout. Exits with code 0
on success.

## Details

Runs `<neo4j-home>/bin/neo4j stop`, which sends a graceful shutdown signal to
the Neo4j JVM process. The server flushes pending transactions and closes
cleanly.

You must stop Neo4j before running `graphpop dump`, `graphpop load`, or
`graphpop import` (Step 2: bulk import), because `neo4j-admin` requires
exclusive access to the data directory.

If Neo4j is not running, the command prints a message indicating this and exits
normally.

## Examples

### Stop the server

```bash
graphpop stop
```

### Stop before dumping a database

```bash
graphpop stop
graphpop dump --database rice3k -o rice3k_v1.dump
graphpop start
```

### Stop a specific installation

```bash
graphpop stop --neo4j-home /opt/neo4j
```

## See Also

- `graphpop start` -- Start the Neo4j server.
- `graphpop status` -- Check whether Neo4j is running.
- `graphpop dump` -- Dump a database (requires Neo4j to be stopped).
- `graphpop load` -- Load a dump file (requires Neo4j to be stopped).
