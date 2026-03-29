"""graphpop config — manage GraphPop configuration."""
from __future__ import annotations

from pathlib import Path

import click
import yaml


CONFIG_PATH = Path.home() / ".graphpop" / "config.yaml"

DEFAULTS = {
    "uri": "bolt://localhost:7687",
    "user": "neo4j",
    "password": "",
    "database": "neo4j",
    "neo4j_home": str(Path.home() / "neo4j"),
}


@click.group()
def config():
    """Manage GraphPop configuration."""
    pass


@config.command()
def init():
    """Create a new GraphPop config file interactively."""
    if CONFIG_PATH.exists():
        if not click.confirm(f"{CONFIG_PATH} already exists. Overwrite?"):
            return

    cfg = {}
    cfg["uri"] = click.prompt("Neo4j URI", default=DEFAULTS["uri"])
    cfg["user"] = click.prompt("Neo4j user", default=DEFAULTS["user"])
    cfg["password"] = click.prompt("Neo4j password", hide_input=True)
    cfg["database"] = click.prompt("Default database", default=DEFAULTS["database"])
    cfg["neo4j_home"] = click.prompt("Neo4j home directory",
                                      default=DEFAULTS["neo4j_home"])

    CONFIG_PATH.parent.mkdir(exist_ok=True)
    with open(CONFIG_PATH, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False)
    CONFIG_PATH.chmod(0o600)  # Restrict permissions (contains password)

    click.echo(f"\nConfig written to {CONFIG_PATH}")
    click.echo("Permissions set to owner-only (600).")


@config.command()
def show():
    """Display the current configuration."""
    if not CONFIG_PATH.exists():
        click.echo(f"No config file found at {CONFIG_PATH}")
        click.echo("Run 'graphpop config init' to create one.")
        return

    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f) or {}

    click.echo(f"Config: {CONFIG_PATH}\n")
    for key, value in cfg.items():
        if key == "password":
            display = "****" if value else "(not set)"
        else:
            display = value
        click.echo(f"  {key}: {display}")

    # Show env var overrides
    import os
    overrides = []
    if os.environ.get("GRAPHPOP_URI"):
        overrides.append(f"  GRAPHPOP_URI={os.environ['GRAPHPOP_URI']}")
    if os.environ.get("GRAPHPOP_USER"):
        overrides.append(f"  GRAPHPOP_USER={os.environ['GRAPHPOP_USER']}")
    if os.environ.get("GRAPHPOP_PASSWORD"):
        overrides.append("  GRAPHPOP_PASSWORD=****")
    if os.environ.get("GRAPHPOP_DATABASE"):
        overrides.append(f"  GRAPHPOP_DATABASE={os.environ['GRAPHPOP_DATABASE']}")
    if overrides:
        click.echo("\nEnvironment overrides:")
        for o in overrides:
            click.echo(o)


@config.command()
@click.argument("key")
@click.argument("value")
def set(key, value):
    """Set a configuration value.

    \b
    Examples:
      graphpop config set database rice3k
      graphpop config set pagecache 20g
      graphpop config set neo4j_home /opt/neo4j
    """
    cfg = {}
    if CONFIG_PATH.exists():
        with open(CONFIG_PATH) as f:
            cfg = yaml.safe_load(f) or {}

    old = cfg.get(key, "(not set)")
    cfg[key] = value

    CONFIG_PATH.parent.mkdir(exist_ok=True)
    with open(CONFIG_PATH, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False)

    click.echo(f"{key}: {old} → {value}")


@config.command()
def path():
    """Print the config file path."""
    click.echo(str(CONFIG_PATH))
