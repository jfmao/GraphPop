# graphpop start

## Title

Start the Neo4j Database Server

## Description

Starts the Neo4j database server as a background daemon. The Neo4j installation
location is read from `~/.graphpop/config.yaml` (the `neo4j_home` key) or can
be overridden with `--neo4j-home`. Equivalent to running `neo4j start` from the
Neo4j `bin/` directory.

## Usage

```
graphpop start [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--neo4j-home` | PATH | *(from config)* | Neo4j installation directory. Overrides the value in `~/.graphpop/config.yaml`. |

## Value

No return value. Prints the Neo4j start-up output (including PID and bolt port)
to stdout. Exits with code 0 on success.

## Details

The command locates the Neo4j binary using the following resolution order:

1. The `--neo4j-home` option, if provided.
2. The `neo4j_home` value in `~/.graphpop/config.yaml`.
3. Fallback candidates: `~/neo4j`, `/var/lib/neo4j`.

It then runs `<neo4j-home>/bin/neo4j start`, which launches the JVM as a
background process. Neo4j typically takes 5-15 seconds to become available on
the bolt port (7687). Use `graphpop status` to confirm readiness.

If Neo4j is already running, the command prints a message indicating this and
exits normally.

### Prerequisites

- Neo4j must be installed (run `graphpop setup` first).
- Java 17+ must be on the `PATH`.
- Port 7687 (bolt) and 7474 (HTTP browser) must be available.

## Examples

### Start with default configuration

```bash
graphpop start
```

### Start a specific Neo4j installation

```bash
graphpop start --neo4j-home /opt/neo4j
```

### Start and verify

```bash
graphpop start && graphpop status
```

## See Also

- `graphpop stop` -- Stop the Neo4j server.
- `graphpop status` -- Check whether Neo4j is running and show version info.
- `graphpop setup` -- Download, configure, and initialize Neo4j.
