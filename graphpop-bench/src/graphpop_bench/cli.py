"""graphpop-bench CLI entry point.

Initial command surface (Step G):

  graphpop-bench profile <output_dir> -- <cmd> <args...>

Subsequent steps (Step H1+, Step I) add:

  graphpop-bench run <competitor> ...
  graphpop-bench run-figure <fig_id> --phase 1
"""
from __future__ import annotations

import json
from pathlib import Path

import click

from .profiling import profile_command, write_receipt


@click.group()
def main():
    """Cross-paper benchmarking + comparison harness for GraphPop."""


@main.command()
@click.argument("output_dir", type=click.Path(file_okay=False))
@click.option("--tool", default="graphpop", show_default=True,
              help="Tool label recorded in the receipt")
@click.option("--tool-version", default=None,
              help="Tool version string recorded in the receipt")
@click.option("--seed", type=int, default=None,
              help="Seed to record in the receipt")
@click.option("--graphpop-commit", default=None,
              help="GraphPop commit SHA to record in the receipt")
@click.argument("cmd", nargs=-1, required=True, type=click.UNPROCESSED)
def profile(output_dir, tool, tool_version, seed, graphpop_commit, cmd):
    """Profile a subprocess and write a JSON receipt.

    Example:

      graphpop-bench profile ./out --tool graphpop -- sleep 0.5
    """
    cmd_list = list(cmd)
    if not cmd_list:
        raise click.UsageError("no command to profile (separate with --)")
    result = profile_command(cmd_list, output_dir)
    receipt = result.as_receipt(
        tool=tool,
        tool_version=tool_version,
        seed=seed,
        graphpop_commit=graphpop_commit,
    )
    path = write_receipt(receipt, output_dir)
    click.echo(
        f"profile exit={result.exit_code} "
        f"wall={result.wall_clock_s:.3f}s "
        f"rss_peak={result.rss_peak_mb:.1f}MB "
        f"backend={result.backend} "
        f"receipt={path}",
        err=True,
    )
    if result.stdout:
        click.echo(result.stdout, nl=False)
    raise SystemExit(result.exit_code)


if __name__ == "__main__":
    main()
