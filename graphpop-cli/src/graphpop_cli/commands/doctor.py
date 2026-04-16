"""graphpop doctor — installation health check."""
from __future__ import annotations

import re
import socket
import subprocess
from pathlib import Path

import click
import yaml

from .setup import DEFAULT_BOLT_PORT, DEFAULT_HTTP_PORT


def _check(label: str, ok: bool, detail: str = "") -> bool:
    """Print a check result and return whether it passed."""
    mark = click.style("OK", fg="green") if ok else click.style("FAIL", fg="red")
    msg = f"  [{mark}] {label}"
    if detail:
        msg += f" — {detail}"
    click.echo(msg)
    return ok


@click.command()
def doctor():
    """Run a full health check on the GraphPop installation.

    Verifies Java, Neo4j home directory, running process, port reachability,
    plugin deployment, config file, and password connectivity.
    """
    click.echo("GraphPop Doctor\n")
    all_ok = True

    # 1. Java
    click.echo("Checking Java...")
    java_ok, java_detail = _check_java_health()
    all_ok &= _check("Java 21+", java_ok, java_detail)

    # 2. Config file
    click.echo("\nChecking configuration...")
    config_path = Path.home() / ".graphpop" / "config.yaml"
    cfg = {}
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
        all_ok &= _check("Config file", True, str(config_path))
    else:
        all_ok &= _check("Config file", False,
                          f"{config_path} not found — run 'graphpop setup'")

    # 3. Neo4j home
    click.echo("\nChecking Neo4j installation...")
    neo4j_home = Path(cfg.get("neo4j_home", Path.home() / "neo4j"))
    neo4j_bin = neo4j_home / "bin" / "neo4j"
    if neo4j_bin.exists():
        all_ok &= _check("Neo4j home", True, str(neo4j_home))
    else:
        all_ok &= _check("Neo4j home", False,
                          f"neo4j binary not found at {neo4j_home}")

    # 4. Neo4j version
    if neo4j_bin.exists():
        result = subprocess.run(
            [str(neo4j_bin), "version"], capture_output=True, text=True,
        )
        version_str = result.stdout.strip() if result.returncode == 0 else "unknown"
        all_ok &= _check("Neo4j version", result.returncode == 0, version_str)

    # 5. Plugin deployment
    click.echo("\nChecking plugin...")
    plugins_dir = neo4j_home / "plugins"
    jar_files = list(plugins_dir.glob("graphpop*.jar")) if plugins_dir.exists() else []
    if jar_files:
        all_ok &= _check("GraphPop plugin", True, jar_files[0].name)
    else:
        all_ok &= _check("GraphPop plugin", False,
                          "not found in plugins/ — run 'graphpop setup'")

    # 6. Neo4j process
    click.echo("\nChecking Neo4j process...")
    if neo4j_bin.exists():
        result = subprocess.run(
            [str(neo4j_bin), "status"], capture_output=True, text=True,
        )
        output = (result.stdout + result.stderr).strip()
        running = result.returncode == 0 and "running" in output.lower()
        all_ok &= _check("Neo4j running", running,
                          output.splitlines()[0] if output else "no output")
    else:
        all_ok &= _check("Neo4j running", False, "neo4j binary not found")

    # 7. Port reachability
    click.echo("\nChecking ports...")
    uri = cfg.get("uri", f"bolt://localhost:{DEFAULT_BOLT_PORT}")
    # Parse port from URI
    port_match = re.search(r":(\d+)$", uri)
    bolt_port = int(port_match.group(1)) if port_match else DEFAULT_BOLT_PORT

    bolt_ok = _is_port_open("127.0.0.1", bolt_port)
    all_ok &= _check(f"Bolt port {bolt_port}", bolt_ok,
                      "reachable" if bolt_ok else "not reachable")

    # Check HTTP port from config
    http_port = cfg.get("http_port", DEFAULT_HTTP_PORT)
    http_ok = _is_port_open("127.0.0.1", http_port)
    all_ok &= _check(f"HTTP port {http_port}", http_ok,
                      "reachable" if http_ok else "not reachable")

    # 8. Bolt connectivity (if neo4j driver available)
    click.echo("\nChecking database connectivity...")
    password = cfg.get("password")
    if password and bolt_ok:
        conn_ok, conn_detail = _check_bolt_connectivity(uri, password)
        all_ok &= _check("Bolt connection", conn_ok, conn_detail)
    elif not bolt_ok:
        all_ok &= _check("Bolt connection", False,
                          "skipped — port not reachable")
    else:
        all_ok &= _check("Bolt connection", False,
                          "skipped — no password in config")

    # Summary
    click.echo("")
    if all_ok:
        click.echo(click.style("All checks passed.", fg="green"))
    else:
        click.echo(click.style(
            "Some checks failed. Review the output above.", fg="yellow"))
    raise SystemExit(0 if all_ok else 1)


def _check_java_health() -> tuple[bool, str]:
    """Return (ok, detail) for Java version check."""
    try:
        result = subprocess.run(
            ["java", "-version"], capture_output=True, text=True,
        )
        output = result.stderr + result.stdout
        first_line = output.splitlines()[0].strip() if output else "unknown"
        m = re.search(r'"(\d+)', output)
        if m and int(m.group(1)) >= 21:
            return True, first_line
        elif m:
            return False, f"{first_line} (need 21+, found {m.group(1)})"
        return False, f"{first_line} (could not parse version)"
    except FileNotFoundError:
        return False, "java not found — install via: conda install -c conda-forge openjdk=21"


def _is_port_open(host: str, port: int, timeout: float = 2.0) -> bool:
    """Check if a TCP port is accepting connections."""
    try:
        with socket.create_connection((host, port), timeout=timeout):
            return True
    except (OSError, ConnectionRefusedError):
        return False


def _check_bolt_connectivity(uri: str, password: str) -> tuple[bool, str]:
    """Try a Bolt connection and return (ok, detail)."""
    try:
        from neo4j import GraphDatabase
        driver = GraphDatabase.driver(uri, auth=("neo4j", password))
        driver.verify_connectivity()
        info = driver.get_server_info()
        driver.close()
        return True, f"connected to {info.agent}"
    except ImportError:
        return False, "neo4j-driver not installed — pip install neo4j"
    except Exception as e:
        return False, str(e)
