"""Tests for graphpop_cli.connection."""
from __future__ import annotations

from graphpop_cli.connection import load_config


# ---------------------------------------------------------------------------
# load_config
# ---------------------------------------------------------------------------

def test_load_config_defaults(tmp_path, monkeypatch):
    # Point to a non-existent config file so no YAML is loaded
    fake_path = tmp_path / "nonexistent.yaml"
    # Clear any env vars that could interfere
    for var in ("GRAPHPOP_URI", "GRAPHPOP_USER",
                "GRAPHPOP_PASSWORD", "GRAPHPOP_DATABASE"):
        monkeypatch.delenv(var, raising=False)

    cfg = load_config(config_path=fake_path)
    assert cfg["uri"] == "bolt://localhost:7687"
    assert cfg["user"] == "neo4j"
    assert cfg["password"] == "neo4j"
    assert cfg["database"] == "neo4j"


def test_load_config_env_override(tmp_path, monkeypatch):
    fake_path = tmp_path / "nonexistent.yaml"
    monkeypatch.setenv("GRAPHPOP_URI", "bolt://remote:7687")
    monkeypatch.setenv("GRAPHPOP_PASSWORD", "s3cret")
    # Ensure other vars don't leak in
    monkeypatch.delenv("GRAPHPOP_USER", raising=False)
    monkeypatch.delenv("GRAPHPOP_DATABASE", raising=False)

    cfg = load_config(config_path=fake_path)
    assert cfg["uri"] == "bolt://remote:7687"
    assert cfg["password"] == "s3cret"
    # Defaults still apply for unset keys
    assert cfg["user"] == "neo4j"
    assert cfg["database"] == "neo4j"


def test_load_config_yaml_file(tmp_path, monkeypatch):
    config_file = tmp_path / "config.yaml"
    config_file.write_text(
        "uri: bolt://yaml-host:7687\n"
        "user: yaml_user\n"
        "password: yaml_pass\n"
        "database: yaml_db\n"
    )
    for var in ("GRAPHPOP_URI", "GRAPHPOP_USER",
                "GRAPHPOP_PASSWORD", "GRAPHPOP_DATABASE"):
        monkeypatch.delenv(var, raising=False)

    cfg = load_config(config_path=config_file)
    assert cfg["uri"] == "bolt://yaml-host:7687"
    assert cfg["user"] == "yaml_user"
    assert cfg["password"] == "yaml_pass"
    assert cfg["database"] == "yaml_db"


def test_load_config_env_beats_yaml(tmp_path, monkeypatch):
    config_file = tmp_path / "config.yaml"
    config_file.write_text(
        "uri: bolt://yaml-host:7687\n"
        "password: yaml_pass\n"
    )
    monkeypatch.setenv("GRAPHPOP_URI", "bolt://env-host:7687")
    monkeypatch.delenv("GRAPHPOP_USER", raising=False)
    monkeypatch.delenv("GRAPHPOP_PASSWORD", raising=False)
    monkeypatch.delenv("GRAPHPOP_DATABASE", raising=False)

    cfg = load_config(config_path=config_file)
    # Env var wins over YAML
    assert cfg["uri"] == "bolt://env-host:7687"
    # YAML value still used when no env override
    assert cfg["password"] == "yaml_pass"
