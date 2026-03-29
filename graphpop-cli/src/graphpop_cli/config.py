"""Shared options and Cypher builders for GraphPop CLI."""
from __future__ import annotations


def build_options_map(consequence: str | None = None,
                      pathway: str | None = None,
                      gene: str | None = None,
                      min_af: float | None = None,
                      max_af: float | None = None,
                      variant_type: str | None = None,
                      **extra) -> dict:
    """Build the options map passed to GraphPop procedures."""
    opts = {}
    if consequence:
        opts["consequence"] = consequence
    if pathway:
        opts["pathway"] = pathway
    if gene:
        opts["gene"] = gene
    if min_af is not None:
        opts["min_af"] = min_af
    if max_af is not None:
        opts["max_af"] = max_af
    if variant_type:
        opts["variant_type"] = variant_type
    opts.update({k: v for k, v in extra.items() if v is not None})
    return opts


def build_cypher(procedure: str, positional: list[str],
                 options: dict | None = None,
                 yield_cols: list[str] | None = None) -> str:
    """Build a CALL ... YIELD Cypher statement."""
    args = ", ".join(positional)
    if options:
        opts_str = ", ".join(
            f"{k}: {_cypher_literal(v)}" for k, v in options.items()
        )
        args += f", {{{opts_str}}}"
    cypher = f"CALL {procedure}({args})"
    if yield_cols:
        cypher += " YIELD " + ", ".join(yield_cols)
    return cypher


def _cypher_literal(v) -> str:
    """Convert Python value to Cypher literal."""
    if isinstance(v, str):
        return f"'{v}'"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, list):
        inner = ", ".join(_cypher_literal(x) for x in v)
        return f"[{inner}]"
    return f"'{v}'"
