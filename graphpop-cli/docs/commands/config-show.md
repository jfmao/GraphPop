# Display current configuration

## Description

`graphpop config show` prints the current GraphPop configuration from
`~/.graphpop/config.yaml`, with the password masked. It also shows any
environment variable overrides that are currently active.

## Usage

```
graphpop config show
```

## Value

Prints each configuration key and its value. The `password` field is always
masked as `****`. If environment variables (`GRAPHPOP_URI`, `GRAPHPOP_USER`,
`GRAPHPOP_PASSWORD`, `GRAPHPOP_DATABASE`) are set, they are listed separately
as overrides.

## Examples

```bash
# Show current config
graphpop config show

# Example output:
# Config: /home/user/.graphpop/config.yaml
#
#   uri: bolt://localhost:7687
#   user: neo4j
#   password: ****
#   database: rice3k
#   neo4j_home: /home/user/neo4j
#
# Environment overrides:
#   GRAPHPOP_DATABASE=human1000g
```

## See Also

- `graphpop config init` -- create config file interactively
- `graphpop config set` -- set a configuration value
- `graphpop config path` -- print config file path
