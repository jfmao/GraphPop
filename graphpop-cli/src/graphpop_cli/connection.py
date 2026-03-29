"""Neo4j connection management for GraphPop CLI."""
from __future__ import annotations

import os
from pathlib import Path

import yaml
from neo4j import GraphDatabase


_DEFAULT_CONFIG_PATH = Path.home() / ".graphpop" / "config.yaml"


def load_config(config_path: Path | None = None) -> dict:
    """Load connection config from file, env vars, or defaults."""
    cfg = {
        "uri": "bolt://localhost:7687",
        "user": "neo4j",
        "password": "neo4j",
        "database": "neo4j",
    }

    # Config file
    path = config_path or _DEFAULT_CONFIG_PATH
    if path.exists():
        with open(path) as f:
            file_cfg = yaml.safe_load(f) or {}
        cfg.update({k: v for k, v in file_cfg.items() if v is not None})

    # Env vars override
    if v := os.environ.get("GRAPHPOP_URI"):
        cfg["uri"] = v
    if v := os.environ.get("GRAPHPOP_USER"):
        cfg["user"] = v
    if v := os.environ.get("GRAPHPOP_PASSWORD"):
        cfg["password"] = v
    if v := os.environ.get("GRAPHPOP_DATABASE"):
        cfg["database"] = v

    return cfg


def get_driver(cfg: dict):
    """Create a Neo4j driver from config dict."""
    return GraphDatabase.driver(cfg["uri"], auth=(cfg["user"], cfg["password"]))


def run_procedure(driver, database: str, cypher: str, **params) -> list[dict]:
    """Run a Cypher procedure call and return list of record dicts."""
    with driver.session(database=database) as session:
        result = session.run(cypher, **params)
        return [record.data() for record in result]
