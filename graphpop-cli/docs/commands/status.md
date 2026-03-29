# graphpop status

## Title

Check Neo4j Server Status, Version, Plugin, and Configuration

## Description

Reports whether Neo4j is running, shows the Neo4j version, displays the current
GraphPop configuration, and checks whether the GraphPop stored procedures plugin
is installed.

## Usage

```
graphpop status [OPTIONS]
```

## Arguments

This command takes no positional arguments.

## Options

| Option | Type | Default | Description |
|---|---|---|---|
| `--neo4j-home` | PATH | *(from config)* | Neo4j installation directory. Overrides the value in `~/.graphpop/config.yaml`. |

## Value

No return value. Prints a multi-section status report to stdout:

```
Neo4j is running at pid 12345
Version: neo4j 5.26.0

GraphPop config (~/.graphpop/config.yaml):
  URI:      bolt://localhost:7687
  Database: rice3k
  Neo4j:    /home/user/neo4j

GraphPop plugin: graphpop-procedures-1.0-SNAPSHOT.jar
```

## Details

The command performs four checks:

1. **Process status** -- Runs `neo4j status` to detect if the server is running
   and to report its PID.

2. **Version** -- Runs `neo4j version` to display the installed Neo4j version
   string.

3. **Configuration** -- Reads `~/.graphpop/config.yaml` and displays the
   connection URI, active database, and Neo4j home directory. Passwords are
   not shown.

4. **Plugin status** -- Scans `<neo4j-home>/plugins/` for any JAR file matching
   `graphpop*.jar`. If found, prints the filename. If not, prints
   `NOT INSTALLED` with instructions on how to build and deploy it.

This command does not require Neo4j to be running. It checks the filesystem
directly for the plugin and configuration.

## Examples

### Quick health check

```bash
graphpop status
```

### Check a specific installation

```bash
graphpop status --neo4j-home /opt/neo4j
```

### Scripting: check if Neo4j is running

```bash
if graphpop status 2>&1 | grep -q "is running"; then
    echo "Neo4j is up"
else
    graphpop start
fi
```

## See Also

- `graphpop start` -- Start the Neo4j server.
- `graphpop stop` -- Stop the Neo4j server.
- `graphpop config show` -- Display full configuration details.
- `graphpop db info` -- Show node/edge counts and installed procedures.
