# Print config file path

## Description

`graphpop config path` prints the path to the GraphPop configuration file.
The default location is `~/.graphpop/config.yaml`.

## Usage

```
graphpop config path
```

## Value

Prints the absolute path to the config file (e.g., `/home/user/.graphpop/config.yaml`).
The file may or may not exist — use `graphpop config init` to create it.

## Examples

```bash
# Print the config path
graphpop config path

# Use in a script
CONFIG=$(graphpop config path)
cat "$CONFIG"
```

## See Also

- `graphpop config init` -- create config file interactively
- `graphpop config show` -- display current configuration
- `graphpop config set` -- set a configuration value
