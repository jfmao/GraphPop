"""Shared validation utilities for GraphPop CLI commands."""
from __future__ import annotations

import re

import click

_IDENT_RE = re.compile(r'^[A-Za-z0-9_-]+$')


def validate_identifier(value: str, label: str = "identifier") -> str:
    """Validate that a value is safe for use as a Cypher property name.

    Only alphanumeric characters, hyphens, and underscores are allowed.
    Raises click.BadParameter if the value contains unsafe characters.
    """
    if not _IDENT_RE.match(value):
        raise click.BadParameter(
            f"Invalid {label}: {value!r}. Only alphanumeric, hyphen, "
            "and underscore characters are allowed."
        )
    return value
