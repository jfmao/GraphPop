# Set a configuration value

## Description

`graphpop config set` updates a single key in the GraphPop configuration file
(`~/.graphpop/config.yaml`). If the file does not exist, it is created.

## Usage

```
graphpop config set KEY VALUE
```

## Arguments

| Name | Type | Required | Description |
|------|------|----------|-------------|
| `KEY` | string | yes | Configuration key to set (e.g., `database`, `uri`, `neo4j_home`). |
| `VALUE` | string | yes | New value for the key. |

## Value

Prints the old and new value: `key: old → new`.

## Details

Common configuration keys:

| Key | Description | Example |
|-----|-------------|---------|
| `uri` | Neo4j bolt URI | `bolt://localhost:7687` |
| `user` | Neo4j username | `neo4j` |
| `password` | Neo4j password | *(prompted interactively via `config init`)* |
| `database` | Default database name | `rice3k` |
| `neo4j_home` | Neo4j installation directory | `/home/user/neo4j` |
| `pagecache` | Neo4j page cache size | `20g` |

Any arbitrary key can be set — the file is a plain YAML dictionary.

## Examples

```bash
# Switch active database
graphpop config set database rice3k

# Set Neo4j home directory
graphpop config set neo4j_home /opt/neo4j

# Set page cache size
graphpop config set pagecache 20g
```

## See Also

- `graphpop config show` -- display current configuration
- `graphpop config init` -- create config file interactively
- `graphpop config path` -- print config file path
- `graphpop db switch` -- switch active database (also updates config)
