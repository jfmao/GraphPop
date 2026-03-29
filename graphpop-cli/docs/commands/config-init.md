# graphpop config init

## Title

Create a GraphPop Configuration File Interactively

## Description

Launches an interactive prompt to create the GraphPop configuration file at
`~/.graphpop/config.yaml`. Asks for the Neo4j connection URI, username,
password, default database, and Neo4j home directory. Sets file permissions
to owner-only (mode 600) to protect the stored password.

## Usage

```
graphpop config init
```

## Arguments

This command takes no positional arguments.

## Options

This command takes no options. All values are prompted interactively.

## Value

No return value. Creates `~/.graphpop/config.yaml` with the following
structure:

```yaml
database: rice3k
neo4j_home: /home/user/neo4j
password: mySecurePass123
uri: bolt://localhost:7687
user: neo4j
```

File permissions are set to `0600` (owner read/write only).

## Details

The interactive prompts and their defaults:

| Prompt | Default | Notes |
|---|---|---|
| Neo4j URI | `bolt://localhost:7687` | Use `bolt+s://` for encrypted connections. |
| Neo4j user | `neo4j` | The default Neo4j admin user. |
| Neo4j password | *(none)* | Input is hidden. No confirmation prompt. |
| Default database | `neo4j` | The database used by all GraphPop commands. |
| Neo4j home directory | `~/neo4j` | Used by `start`, `stop`, `status`, `dump`, `load`. |

If `~/.graphpop/config.yaml` already exists, you will be asked to confirm
before overwriting. Declining preserves the existing file.

The configuration file is also created automatically by `graphpop setup`. Use
`config init` when you need to reconfigure without re-running the full setup
(e.g., connecting to a remote server or changing the password).

### Configuration Precedence

GraphPop resolves configuration values in this order (highest wins):

1. Command-line options (`--uri`, `--user`, `--password`, `--database`)
2. Environment variables (`GRAPHPOP_URI`, `GRAPHPOP_USER`, `GRAPHPOP_PASSWORD`, `GRAPHPOP_DATABASE`)
3. Config file (`~/.graphpop/config.yaml`)
4. Built-in defaults

## Examples

### First-time configuration

```bash
graphpop config init
# Neo4j URI [bolt://localhost:7687]:
# Neo4j user [neo4j]:
# Neo4j password:
# Default database [neo4j]: rice3k
# Neo4j home directory [/home/user/neo4j]:
#
# Config written to /home/user/.graphpop/config.yaml
# Permissions set to owner-only (600).
```

### Reconfigure for a remote server

```bash
graphpop config init
# Neo4j URI [bolt://localhost:7687]: bolt://server.lab.org:7687
# Neo4j user [neo4j]: admin
# Neo4j password:
# Default database [neo4j]: production
# Neo4j home directory [/home/user/neo4j]: /opt/neo4j
```

## See Also

- `graphpop config show` -- Display the current configuration.
- `graphpop config set` -- Change a single configuration value.
- `graphpop config path` -- Print the config file path.
- `graphpop setup` -- Full setup including Neo4j download and config creation.
