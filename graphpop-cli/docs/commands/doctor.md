# Run a full health check on the GraphPop installation

## Description

`graphpop doctor` verifies that all components of the GraphPop stack are
correctly installed and operational. It checks Java, the Neo4j installation
directory, running process status, port reachability, plugin deployment,
configuration, and database connectivity.

Each check prints a color-coded `[OK]` or `[FAIL]` result with diagnostic
detail. The command exits with code 0 if all checks pass, or 1 if any check
fails, making it suitable for use in CI pipelines and automated validation.

No arguments or database connection are required — `doctor` reads the
configuration file at `~/.graphpop/config.yaml` (if present) and probes the
local system.

## Usage

```
graphpop doctor
```

## Arguments

This command takes no positional arguments.

## Options

This command takes no options.

## Checks performed

| # | Check | What it verifies |
|---|-------|-----------------|
| 1 | Java 21+ | `java -version` returns a JDK/JRE >= 21 |
| 2 | Config file | `~/.graphpop/config.yaml` exists and is readable |
| 3 | Neo4j home | The `neo4j_home` directory from config contains a `bin/neo4j` binary |
| 4 | Neo4j version | `neo4j version` command succeeds |
| 5 | GraphPop plugin | `graphpop-procedures*.jar` exists in the Neo4j plugins directory |
| 6 | Neo4j running | `neo4j status` reports the process is running |
| 7 | Bolt port | TCP connection to the configured Bolt port (default 7687) succeeds |
| 8 | HTTP port | TCP connection to the configured HTTP port (default 7474) succeeds |
| 9 | Bolt connection | Authenticates via the Neo4j Python driver and retrieves server info |

Checks 7-9 are only meaningful when Neo4j is running. Check 9 requires the
`neo4j` Python driver to be installed; if missing, it reports the import error
rather than failing silently.

## Details

### Configuration discovery

`doctor` reads `~/.graphpop/config.yaml` for `neo4j_home`, `uri`, and
`password`. If the file does not exist, it falls back to `~/neo4j` as the
Neo4j home directory and default ports (7687 Bolt, 7474 HTTP).

### Exit codes

- **0** — all checks passed
- **1** — one or more checks failed

### Typical failure scenarios

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Java FAIL | JDK not installed or < 21 | `conda install -c conda-forge openjdk=21` |
| Config FAIL | Never ran setup | `graphpop setup --password mypass` |
| Neo4j home FAIL | Neo4j removed or path changed | Re-run `graphpop setup` or fix `config.yaml` |
| Plugin FAIL | JAR not deployed | `graphpop setup --skip-download --password mypass` |
| Neo4j running FAIL | Server not started | `graphpop start` |
| Bolt port FAIL | Wrong port or firewall | Check `--bolt-port` in setup, or firewall rules |
| Bolt connection FAIL | Wrong password | Update password in `~/.graphpop/config.yaml` |

## Examples

```bash
# Basic health check
graphpop doctor

# Use in a CI script (exits non-zero on failure)
graphpop doctor && echo "Ready" || echo "Fix issues above"

# After initial setup
graphpop setup --password mypass
graphpop start
graphpop doctor
```

## See Also

- [graphpop setup](setup.md) — install and configure Neo4j
- [graphpop start](start.md) — start the Neo4j server
- [graphpop status](status.md) — quick status check (less thorough than doctor)
- [graphpop validate](validate.md) — validate graph database structure (not installation)
