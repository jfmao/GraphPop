# graphpop config

## Title

Manage GraphPop Configuration (show, set, path)

## Description

The `graphpop config` group contains three subcommands for viewing and modifying
the GraphPop configuration file at `~/.graphpop/config.yaml`. This page
documents `config show`, `config set`, and `config path`.

For `config init` (interactive creation), see `config-init.md`.

---

## graphpop config show

### Usage

```
graphpop config show
```

### Description

Displays the current GraphPop configuration. Passwords are masked as `****`.
Also shows any active environment variable overrides.

### Value

Prints the configuration to stdout:

```
Config: /home/user/.graphpop/config.yaml

  uri: bolt://localhost:7687
  user: neo4j
  password: ****
  database: rice3k
  neo4j_home: /home/user/neo4j

Environment overrides:
  GRAPHPOP_DATABASE=human_chr22
```

If no config file exists, prints a message suggesting `graphpop config init`.

### Details

Reads `~/.graphpop/config.yaml` and displays all key-value pairs. The password
is never shown in plaintext; it is displayed as `****` if set or `(not set)` if
empty.

Environment variable overrides are checked for `GRAPHPOP_URI`, `GRAPHPOP_USER`,
`GRAPHPOP_PASSWORD`, and `GRAPHPOP_DATABASE`. If any are set, they are listed
in a separate section (passwords still masked).

### Examples

```bash
graphpop config show
```

```bash
# Check which database is active
graphpop config show | grep database
```

---

## graphpop config set

### Usage

```
graphpop config set KEY VALUE
```

### Arguments

| Argument | Type | Description |
|---|---|---|
| `KEY` | TEXT | Configuration key to set (e.g., `database`, `uri`, `neo4j_home`, `pagecache`). |
| `VALUE` | TEXT | New value for the key. |

### Description

Sets a single configuration value in `~/.graphpop/config.yaml`. Creates the
file and parent directory if they do not exist. Prints the old and new values.

### Value

```
database: neo4j -> rice3k
```

### Details

Any key can be set, not just the standard connection keys. This allows storing
custom settings (e.g., `pagecache`, `heap`) that other tools or scripts can
read from the same config file.

Standard keys recognized by GraphPop commands:

| Key | Used By | Description |
|---|---|---|
| `uri` | All commands | Neo4j bolt connection URI |
| `user` | All commands | Neo4j username |
| `password` | All commands | Neo4j password |
| `database` | All commands | Default database name |
| `neo4j_home` | `start`, `stop`, `status`, `setup`, `dump`, `load`, `import` | Neo4j installation directory |

### Examples

```bash
# Switch the default database
graphpop config set database rice3k

# Point to a different Neo4j installation
graphpop config set neo4j_home /opt/neo4j-enterprise

# Store a custom setting
graphpop config set pagecache 20g
```

---

## graphpop config path

### Usage

```
graphpop config path
```

### Description

Prints the absolute path to the GraphPop configuration file. Useful in scripts
and for debugging.

### Value

```
/home/user/.graphpop/config.yaml
```

### Details

Always prints `~/.graphpop/config.yaml` (expanded to the full home directory
path). This is the default config path used by all GraphPop commands unless
overridden with `--config`.

### Examples

```bash
# View the raw config file
cat $(graphpop config path)

# Edit the config manually
nano $(graphpop config path)

# Back up the config
cp $(graphpop config path) ~/config_backup.yaml
```

---

## Configuration Precedence

GraphPop resolves each setting from the following sources, in order of
decreasing priority:

1. **Command-line options** (`--uri`, `--user`, `--password`, `--database`)
2. **Environment variables** (`GRAPHPOP_URI`, `GRAPHPOP_USER`, `GRAPHPOP_PASSWORD`, `GRAPHPOP_DATABASE`)
3. **Config file** (`~/.graphpop/config.yaml`)
4. **Built-in defaults** (`bolt://localhost:7687`, `neo4j`, `neo4j`, `neo4j`)

## See Also

- `graphpop config init` -- Create a new config file interactively.
- `graphpop setup` -- Full setup that also creates the config file.
- `graphpop db switch` -- Shortcut to change the active database.
