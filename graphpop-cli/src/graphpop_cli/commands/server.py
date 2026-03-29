"""graphpop start/stop/status — Neo4j server lifecycle management."""
from __future__ import annotations

import subprocess
from pathlib import Path

import click
import yaml


def _get_neo4j_home() -> Path:
    """Get Neo4j home from config or default."""
    config_path = Path.home() / ".graphpop" / "config.yaml"
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
        if "neo4j_home" in cfg:
            return Path(cfg["neo4j_home"])
    # Fallbacks
    for candidate in [Path.home() / "neo4j", Path("/var/lib/neo4j")]:
        if (candidate / "bin" / "neo4j").exists():
            return candidate
    return Path.home() / "neo4j"


def _run_neo4j_cmd(command: str, neo4j_home: Path | None = None) -> tuple[int, str]:
    """Run a neo4j command and return (returncode, output)."""
    home = neo4j_home or _get_neo4j_home()
    neo4j_bin = home / "bin" / "neo4j"
    if not neo4j_bin.exists():
        return 1, f"Neo4j not found at {home}. Run 'graphpop setup' first."
    result = subprocess.run(
        [str(neo4j_bin), command],
        capture_output=True, text=True,
    )
    output = (result.stdout + result.stderr).strip()
    return result.returncode, output


@click.command()
@click.option("--neo4j-home", type=click.Path(), help="Neo4j installation directory")
def start(neo4j_home):
    """Start the Neo4j database server."""
    home = Path(neo4j_home) if neo4j_home else None
    click.echo("Starting Neo4j...")
    rc, output = _run_neo4j_cmd("start", home)
    click.echo(output)
    if rc == 0:
        click.echo("\nNeo4j started. Use 'graphpop status' to verify.")


@click.command()
@click.option("--neo4j-home", type=click.Path(), help="Neo4j installation directory")
def stop(neo4j_home):
    """Stop the Neo4j database server."""
    home = Path(neo4j_home) if neo4j_home else None
    click.echo("Stopping Neo4j...")
    rc, output = _run_neo4j_cmd("stop", home)
    click.echo(output)


@click.command()
@click.option("--neo4j-home", type=click.Path(), help="Neo4j installation directory")
def status(neo4j_home):
    """Check whether Neo4j is running and show database info."""
    home = Path(neo4j_home) if neo4j_home else _get_neo4j_home()

    # Check Neo4j process
    rc, output = _run_neo4j_cmd("status", home)
    click.echo(output)

    # Show version
    neo4j_bin = home / "bin" / "neo4j"
    if neo4j_bin.exists():
        result = subprocess.run([str(neo4j_bin), "version"],
                                capture_output=True, text=True)
        if result.returncode == 0:
            click.echo(f"Version: {result.stdout.strip()}")

    # Show config
    config_path = Path.home() / ".graphpop" / "config.yaml"
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
        click.echo(f"\nGraphPop config ({config_path}):")
        click.echo(f"  URI:      {cfg.get('uri', 'not set')}")
        click.echo(f"  Database: {cfg.get('database', 'not set')}")
        click.echo(f"  Neo4j:    {cfg.get('neo4j_home', 'not set')}")

    # Show plugin status
    plugins_dir = home / "plugins"
    jar_files = list(plugins_dir.glob("graphpop*.jar")) if plugins_dir.exists() else []
    if jar_files:
        click.echo(f"\nGraphPop plugin: {jar_files[0].name}")
    else:
        click.echo("\nGraphPop plugin: NOT INSTALLED")
        click.echo("  Build with: cd graphpop-procedures && mvn package")
        click.echo("  Deploy with: graphpop setup --deploy-plugin target/graphpop-procedures-*.jar")
