"""graphpop dump/load — database dump and restore for sharing."""
from __future__ import annotations

import json
import subprocess
from datetime import datetime
from pathlib import Path

import click
import yaml


def _get_neo4j_home() -> Path:
    """Get Neo4j home from config."""
    config_path = Path.home() / ".graphpop" / "config.yaml"
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
        if "neo4j_home" in cfg:
            return Path(cfg["neo4j_home"])
    return Path.home() / "neo4j"


@click.command()
@click.option("--database", required=True, help="Database name to dump")
@click.option("-o", "--output", "output_path", type=click.Path(),
              help="Output dump file path (default: <database>_<date>.dump)")
@click.option("--neo4j-home", type=click.Path(), help="Neo4j installation directory")
@click.option("--manifest/--no-manifest", default=True,
              help="Generate a JSON manifest with database metadata (default: yes)")
def dump(database, output_path, neo4j_home, manifest):
    """Dump a Neo4j database to a file for sharing.

    Creates a neo4j-admin dump file that can be shared and loaded on another
    machine. Optionally generates a JSON manifest with node/edge counts,
    populations, chromosomes, and computed statistics metadata.

    \b
    Examples:
      graphpop dump --database rice3k
      graphpop dump --database rice3k -o rice3k_v1.dump
      graphpop dump --database neo4j --no-manifest
    """
    home = Path(neo4j_home) if neo4j_home else _get_neo4j_home()
    admin_bin = home / "bin" / "neo4j-admin"

    if not admin_bin.exists():
        click.echo(f"Error: neo4j-admin not found at {admin_bin}", err=True)
        raise SystemExit(1)

    # Default output path
    if not output_path:
        date_str = datetime.now().strftime("%Y%m%d")
        output_path = f"graphpop_{database}_{date_str}.dump"

    output_file = Path(output_path)

    click.echo(f"Dumping database '{database}' to {output_file}...")

    # Run neo4j-admin dump
    cmd = [
        str(admin_bin), "database", "dump",
        f"--to-path={output_file.parent}",
        database,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        click.echo(f"Error: {result.stderr.strip()}", err=True)
        click.echo("\nNote: Neo4j must be stopped before dumping.", err=True)
        click.echo("Run 'graphpop stop' first, then retry.", err=True)
        raise SystemExit(1)

    # Rename if neo4j-admin used its own naming
    expected_dump = output_file.parent / f"{database}.dump"
    if expected_dump.exists() and expected_dump != output_file:
        expected_dump.rename(output_file)

    size = output_file.stat().st_size if output_file.exists() else 0
    click.echo(f"Dump complete: {output_file} ({_format_size(size)})")

    # Generate manifest
    if manifest:
        manifest_path = output_file.with_suffix(".manifest.json")
        _generate_manifest(home, database, manifest_path, output_file)
        click.echo(f"Manifest: {manifest_path}")


@click.command()
@click.option("--dump-file", required=True, type=click.Path(exists=True),
              help="Path to the .dump file")
@click.option("--database", required=True,
              help="Name for the restored database")
@click.option("--neo4j-home", type=click.Path(), help="Neo4j installation directory")
@click.option("--overwrite", is_flag=True, help="Overwrite existing database")
def load(dump_file, database, neo4j_home, overwrite):
    """Load a database from a dump file.

    Restores a previously dumped Neo4j database. The database name can be
    different from the original.

    \b
    Examples:
      graphpop load --dump-file rice3k_v1.dump --database rice3k
      graphpop load --dump-file shared_db.dump --database myanalysis --overwrite
    """
    home = Path(neo4j_home) if neo4j_home else _get_neo4j_home()
    admin_bin = home / "bin" / "neo4j-admin"

    if not admin_bin.exists():
        click.echo(f"Error: neo4j-admin not found at {admin_bin}", err=True)
        raise SystemExit(1)

    dump_path = Path(dump_file)
    click.echo(f"Loading database '{database}' from {dump_path}...")

    cmd = [
        str(admin_bin), "database", "load",
        f"--from-path={dump_path.parent}",
    ]
    if overwrite:
        cmd.append("--overwrite-destination=true")
    cmd.append(database)

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        click.echo(f"Error: {result.stderr.strip()}", err=True)
        if "already exists" in result.stderr:
            click.echo("Use --overwrite to replace the existing database.", err=True)
        elif "running" in result.stderr.lower():
            click.echo("Stop Neo4j first: graphpop stop", err=True)
        raise SystemExit(1)

    # Update config
    config_path = Path.home() / ".graphpop" / "config.yaml"
    cfg = {}
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
    cfg["database"] = database
    with open(config_path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False)

    click.echo(f"Database '{database}' loaded successfully.")
    click.echo(f"Config updated to use database '{database}'.")
    click.echo(f"\nNext: graphpop start && graphpop db info")


def _generate_manifest(neo4j_home: Path, database: str,
                       manifest_path: Path, dump_path: Path):
    """Generate a JSON manifest with database metadata."""
    manifest = {
        "database": database,
        "date": datetime.now().isoformat(),
        "dump_file": str(dump_path),
        "dump_size_bytes": dump_path.stat().st_size if dump_path.exists() else 0,
    }

    # Try to query the database for node/edge counts
    # (only works if Neo4j is running — otherwise just save basic info)
    try:
        from ..connection import load_config, get_driver
        cfg = load_config()
        driver = get_driver(cfg)
        with driver.session(database=database) as session:
            # Node counts
            result = session.run(
                "CALL db.labels() YIELD label "
                "CALL { WITH label MATCH (n) WHERE label IN labels(n) "
                "RETURN count(n) AS cnt } RETURN label, cnt"
            )
            manifest["node_counts"] = {r["label"]: r["cnt"] for r in result}

            # Relationship counts
            result = session.run(
                "CALL db.relationshipTypes() YIELD relationshipType AS type "
                "CALL { WITH type MATCH ()-[r]->() WHERE type(r) = type "
                "RETURN count(r) AS cnt } RETURN type, cnt"
            )
            manifest["edge_counts"] = {r["type"]: r["cnt"] for r in result}

            # Populations
            result = session.run(
                "MATCH (p:Population) RETURN p.populationId AS pop, "
                "p.n_samples AS n ORDER BY n DESC"
            )
            manifest["populations"] = {r["pop"]: r["n"] for r in result}

            # Chromosomes
            result = session.run(
                "MATCH (c:Chromosome) RETURN c.chromosomeId AS chr, "
                "c.length AS len ORDER BY chr"
            )
            manifest["chromosomes"] = {r["chr"]: r["len"] for r in result}

        driver.close()
    except Exception:
        manifest["note"] = "Neo4j not running; metadata not available"

    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, default=str)


def _format_size(size_bytes: int) -> str:
    """Format bytes as human-readable."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(size_bytes) < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"
