"""graphpop db — database management (list, create, switch, drop, info)."""
from __future__ import annotations

from pathlib import Path

import click
import yaml

from ..cli import pass_ctx


@click.group()
def db():
    """Manage Neo4j databases for GraphPop."""
    pass


@db.command()
@pass_ctx
def list(ctx):
    """List all databases with sizes and status."""
    cypher = "SHOW DATABASES YIELD name, currentStatus, sizeOnDisk ORDER BY name"
    try:
        records = ctx.run(cypher)
    except Exception:
        # Fallback for Neo4j Community (SHOW DATABASES may not return sizeOnDisk)
        try:
            records = ctx.run("SHOW DATABASES YIELD name, currentStatus ORDER BY name")
        except Exception as e:
            click.echo(f"Error: {e}", err=True)
            raise SystemExit(1)

    if not records:
        click.echo("No databases found.")
        return

    # Show current active database
    config_path = Path.home() / ".graphpop" / "config.yaml"
    active_db = "neo4j"
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
        active_db = cfg.get("database", "neo4j")

    click.echo(f"{'Database':<25} {'Status':<12} {'Size':<15} {'Active'}")
    click.echo("-" * 60)
    for rec in records:
        name = rec.get("name", "?")
        status = rec.get("currentStatus", "?")
        size = rec.get("sizeOnDisk", "")
        if isinstance(size, (int, float)) and size > 0:
            size = _format_size(size)
        active = " *" if name == active_db else ""
        click.echo(f"{name:<25} {status:<12} {str(size):<15}{active}")


@db.command()
@click.argument("name")
@pass_ctx
def create(ctx, name):
    """Create a new database."""
    click.echo(f"Creating database '{name}'...")
    try:
        # Must run against system database
        from neo4j import GraphDatabase
        driver = GraphDatabase.driver(ctx.cfg["uri"],
                                       auth=(ctx.cfg["user"], ctx.cfg["password"]))
        with driver.session(database="system") as session:
            session.run(f"CREATE DATABASE `{name}` IF NOT EXISTS")
        driver.close()
        click.echo(f"Database '{name}' created.")
        click.echo(f"Switch to it with: graphpop db switch {name}")
    except Exception as e:
        if "Unsupported" in str(e) or "Enterprise" in str(e):
            click.echo(
                "Error: CREATE DATABASE requires Neo4j Enterprise Edition.\n"
                "With Community Edition, use 'neo4j' as the default database\n"
                "or create databases via neo4j-admin.",
                err=True,
            )
        else:
            click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)


@db.command()
@click.argument("name")
def switch(name):
    """Set the active database in GraphPop config."""
    config_path = Path.home() / ".graphpop" / "config.yaml"
    cfg = {}
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
    old = cfg.get("database", "neo4j")
    cfg["database"] = name
    config_path.parent.mkdir(exist_ok=True)
    with open(config_path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False)
    click.echo(f"Active database: {old} → {name}")
    click.echo(f"All graphpop commands will now use database '{name}'.")


@db.command()
@click.argument("name")
@click.option("--force", is_flag=True, help="Skip confirmation prompt")
@pass_ctx
def drop(ctx, name, force):
    """Drop a database (requires confirmation)."""
    if name in ("neo4j", "system"):
        click.echo(f"Error: Cannot drop the '{name}' system database.", err=True)
        raise SystemExit(1)

    if not force:
        click.confirm(f"Drop database '{name}'? This cannot be undone", abort=True)

    try:
        from neo4j import GraphDatabase
        driver = GraphDatabase.driver(ctx.cfg["uri"],
                                       auth=(ctx.cfg["user"], ctx.cfg["password"]))
        with driver.session(database="system") as session:
            session.run(f"DROP DATABASE `{name}` IF EXISTS")
        driver.close()
        click.echo(f"Database '{name}' dropped.")
    except Exception as e:
        if "Unsupported" in str(e) or "Enterprise" in str(e):
            click.echo(
                "Error: DROP DATABASE requires Neo4j Enterprise Edition.\n"
                "With Community Edition, use neo4j-admin to manage databases.",
                err=True,
            )
        else:
            click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)


@db.command()
@pass_ctx
def info(ctx):
    """Show detailed information about the current database."""
    click.echo(f"Database: {ctx.database}\n")

    # Node counts
    try:
        records = ctx.run(
            "CALL db.labels() YIELD label "
            "CALL { WITH label MATCH (n) WHERE label IN labels(n) "
            "RETURN count(n) AS cnt } RETURN label, cnt ORDER BY cnt DESC"
        )
        if records:
            click.echo("Node counts:")
            for rec in records:
                click.echo(f"  {rec['label']:<20} {rec['cnt']:>12,}")

        # Relationship counts
        records = ctx.run(
            "CALL db.relationshipTypes() YIELD relationshipType AS type "
            "CALL { WITH type MATCH ()-[r]->() WHERE type(r) = type "
            "RETURN count(r) AS cnt } RETURN type, cnt ORDER BY cnt DESC"
        )
        if records:
            click.echo("\nRelationship counts:")
            for rec in records:
                click.echo(f"  {rec['type']:<25} {rec['cnt']:>12,}")

        # Check GraphPop procedures
        records = ctx.run(
            "SHOW PROCEDURES YIELD name WHERE name STARTS WITH 'graphpop' "
            "RETURN name ORDER BY name"
        )
        if records:
            click.echo(f"\nGraphPop procedures ({len(records)}):")
            for rec in records:
                click.echo(f"  {rec['name']}")
        else:
            click.echo("\nGraphPop procedures: NONE INSTALLED")

    except Exception as e:
        click.echo(f"Error querying database: {e}", err=True)


def _format_size(size_bytes: int | float) -> str:
    """Format bytes as human-readable size."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(size_bytes) < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"
