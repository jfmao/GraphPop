"""graphpop setup — download, configure, and initialize Neo4j for GraphPop."""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
import tarfile
from pathlib import Path

import click
import yaml

from ..cli import pass_ctx


DEFAULT_NEO4J_HOME = Path.home() / "neo4j"
NEO4J_VERSION = "5.26.0"
NEO4J_DOWNLOAD_URL = (
    f"https://dist.neo4j.org/neo4j-community-{NEO4J_VERSION}-unix.tar.gz"
)


@click.command()
@click.option("--neo4j-home", type=click.Path(), default=str(DEFAULT_NEO4J_HOME),
              help=f"Neo4j installation directory (default: {DEFAULT_NEO4J_HOME})")
@click.option("--pagecache", default="16g",
              help="Neo4j page cache size (default: 16g)")
@click.option("--heap", default="4g",
              help="Neo4j JVM heap size (default: 4g)")
@click.option("--password", prompt=True, hide_input=True,
              confirmation_prompt=True,
              help="Neo4j password for the 'neo4j' user")
@click.option("--skip-download", is_flag=True,
              help="Skip downloading Neo4j (use existing installation)")
@click.option("--deploy-plugin", type=click.Path(exists=True),
              help="Path to graphpop-procedures.jar to deploy")
def setup(neo4j_home, pagecache, heap, password, skip_download, deploy_plugin):
    """Set up Neo4j for GraphPop.

    Downloads Neo4j Community Edition, configures memory settings, sets the
    initial password, deploys the GraphPop procedures plugin, and creates
    the GraphPop config file.

    \b
    Examples:
      graphpop setup --password mypass
      graphpop setup --neo4j-home /opt/neo4j --pagecache 20g --heap 8g
      graphpop setup --skip-download --deploy-plugin target/graphpop-procedures.jar
    """
    neo4j_path = Path(neo4j_home)

    # Step 1: Download Neo4j
    if not skip_download:
        if neo4j_path.exists() and (neo4j_path / "bin" / "neo4j").exists():
            click.echo(f"Neo4j already installed at {neo4j_path}")
            if not click.confirm("Re-install?"):
                skip_download = True

    if not skip_download:
        _download_neo4j(neo4j_path)

    # Verify installation
    neo4j_bin = neo4j_path / "bin" / "neo4j"
    if not neo4j_bin.exists():
        click.echo(f"Error: Neo4j not found at {neo4j_path}", err=True)
        click.echo("Use --neo4j-home to specify the installation directory.", err=True)
        raise SystemExit(1)

    # Step 2: Configure Neo4j
    click.echo("\nConfiguring Neo4j...")
    _configure_neo4j(neo4j_path, pagecache, heap)

    # Step 3: Set initial password
    click.echo("Setting Neo4j password...")
    _set_password(neo4j_path, password)

    # Step 4: Deploy GraphPop plugin
    if deploy_plugin:
        click.echo(f"Deploying GraphPop plugin from {deploy_plugin}...")
        plugin_dest = neo4j_path / "plugins" / "graphpop-procedures.jar"
        shutil.copy2(deploy_plugin, plugin_dest)
        click.echo(f"  Deployed to {plugin_dest}")

    # Step 5: Create GraphPop config
    config_dir = Path.home() / ".graphpop"
    config_dir.mkdir(exist_ok=True)
    config_path = config_dir / "config.yaml"

    config = {
        "uri": "bolt://localhost:7687",
        "user": "neo4j",
        "password": password,
        "database": "neo4j",
        "neo4j_home": str(neo4j_path),
    }
    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False)
    click.echo(f"\nGraphPop config written to {config_path}")

    # Step 6: Summary
    click.echo(f"""
Setup complete!

  Neo4j home:    {neo4j_path}
  Page cache:    {pagecache}
  Heap:          {heap}
  Config:        {config_path}
  Plugin:        {'deployed' if deploy_plugin else 'not deployed (use --deploy-plugin)'}

Next steps:
  graphpop start                         # Start Neo4j
  graphpop import --vcf data.vcf.gz \\
    --panel panel.txt --database mydb    # Import data
  graphpop diversity chr1 1 50000000 POP # Run analysis
""")


def _download_neo4j(dest: Path):
    """Download and extract Neo4j Community Edition."""
    import urllib.request

    tarball = Path(f"/tmp/neo4j-community-{NEO4J_VERSION}-unix.tar.gz")
    if tarball.exists():
        click.echo(f"Using cached download: {tarball}")
    else:
        click.echo(f"Downloading Neo4j {NEO4J_VERSION}...")
        click.echo(f"  URL: {NEO4J_DOWNLOAD_URL}")
        urllib.request.urlretrieve(NEO4J_DOWNLOAD_URL, tarball)
        click.echo(f"  Downloaded to {tarball}")

    click.echo(f"Extracting to {dest}...")
    if dest.exists():
        shutil.rmtree(dest)

    with tarfile.open(tarball) as tf:
        tf.extractall(dest.parent)

    # The tarball extracts to neo4j-community-X.Y.Z/
    extracted = dest.parent / f"neo4j-community-{NEO4J_VERSION}"
    if extracted.exists() and extracted != dest:
        extracted.rename(dest)
    click.echo(f"  Installed to {dest}")


def _configure_neo4j(neo4j_home: Path, pagecache: str, heap: str):
    """Configure Neo4j memory and settings."""
    conf_path = neo4j_home / "conf" / "neo4j.conf"

    # Read existing config
    lines = conf_path.read_text().splitlines() if conf_path.exists() else []

    # Settings to apply
    settings = {
        "server.memory.pagecache.size": pagecache,
        "server.memory.heap.initial_size": heap,
        "server.memory.heap.max_size": heap,
        "server.directories.import": "import",
        "db.tx_log.rotation.retention_policy": "2 days 2G",
        "dbms.security.procedures.unrestricted": "graphpop.*",
    }

    # Update or append settings
    updated_keys = set()
    new_lines = []
    for line in lines:
        key = line.split("=")[0].strip() if "=" in line and not line.startswith("#") else None
        if key and key in settings:
            new_lines.append(f"{key}={settings[key]}")
            updated_keys.add(key)
        else:
            new_lines.append(line)

    # Append settings not yet in config
    for key, value in settings.items():
        if key not in updated_keys:
            new_lines.append(f"{key}={value}")

    conf_path.write_text("\n".join(new_lines) + "\n")
    for k, v in settings.items():
        click.echo(f"  {k}={v}")


def _set_password(neo4j_home: Path, password: str):
    """Set the initial Neo4j password."""
    admin_bin = neo4j_home / "bin" / "neo4j-admin"
    try:
        result = subprocess.run(
            [str(admin_bin), "dbms", "set-initial-password", password],
            capture_output=True, text=True,
        )
        if result.returncode == 0:
            click.echo("  Password set successfully")
        else:
            # May already be set
            if "already" in result.stderr.lower() or "already" in result.stdout.lower():
                click.echo("  Password already set (use Neo4j browser to change)")
            else:
                click.echo(f"  Warning: {result.stderr.strip()}")
    except FileNotFoundError:
        click.echo("  Warning: neo4j-admin not found, skipping password setup")
