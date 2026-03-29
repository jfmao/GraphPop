"""Output formatters for GraphPop CLI."""
from __future__ import annotations

import csv
import io
import json
import sys
from typing import TextIO


def write_header_comment(out: TextIO, command: str, args: dict):
    """Write a comment header with command info."""
    parts = [f"# graphpop {command}"]
    for k, v in args.items():
        if v is not None and v is not False:
            parts.append(f"# {k}: {v}")
    out.write("\n".join(parts) + "\n")


def write_tsv(records: list[dict], out: TextIO, header: bool = True):
    """Write records as TSV."""
    if not records:
        return
    keys = list(records[0].keys())
    if header:
        out.write("\t".join(keys) + "\n")
    for rec in records:
        vals = []
        for k in keys:
            v = rec[k]
            if isinstance(v, float):
                vals.append(f"{v:.6g}")
            elif isinstance(v, list):
                vals.append(",".join(str(x) for x in v))
            elif v is None:
                vals.append("NA")
            else:
                vals.append(str(v))
        out.write("\t".join(vals) + "\n")


def write_csv(records: list[dict], out: TextIO):
    """Write records as CSV."""
    if not records:
        return
    writer = csv.DictWriter(out, fieldnames=records[0].keys())
    writer.writeheader()
    for rec in records:
        writer.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v)
                         for k, v in rec.items()})


def write_json(records: list[dict], out: TextIO):
    """Write records as JSON."""
    json.dump(records, out, indent=2, default=str)
    out.write("\n")


def get_output(output_path: str | None) -> TextIO:
    """Get output file handle (stdout if None)."""
    if output_path:
        return open(output_path, "w")
    return sys.stdout


def format_output(records: list[dict], output_path: str | None,
                  fmt: str = "tsv", command: str = "", args: dict | None = None):
    """Format and write records to output."""
    out = get_output(output_path)
    try:
        if args and output_path:
            write_header_comment(out, command, args)
        if fmt == "tsv":
            write_tsv(records, out)
        elif fmt == "csv":
            write_csv(records, out)
        elif fmt == "json":
            write_json(records, out)
    finally:
        if output_path and out is not sys.stdout:
            out.close()
