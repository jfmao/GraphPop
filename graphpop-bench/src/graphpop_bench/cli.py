"""graphpop-bench CLI entry point.

Command surface:

  graphpop-bench profile <output_dir> -- <cmd> <args...>
  graphpop-bench run plink_grm --input <vcf-or-bfile-prefix> \\
                                 --output <dir>
  (subsequent H2+ wrappers add more `run <competitor>` subcommands)

Subsequent steps (Step I) add:

  graphpop-bench run-figure <fig_id> --phase 1
"""
from __future__ import annotations

import json
from pathlib import Path

import click

from .competitors import KingRunner, PlinkGrmRunner
from .profiling import profile_command, write_receipt


@click.group()
def main():
    """Cross-paper benchmarking + comparison harness for GraphPop."""


@main.group()
def run():
    """Run a registered competitor wrapper.

    Subcommands:

      plink_grm  PLINK 2.0 `--make-grm-bin` on a VCF or BED prefix.
      king       KING-robust `--kinship` (or `--related`).
    """


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


@run.command("plink_grm")
@click.option("--input", "input_path", required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="Input .vcf(.gz) or BED-prefix (no extension)")
@click.option("--output", "output_dir", required=True,
              type=click.Path(file_okay=False),
              help="Output directory (created if missing)")
@click.option("--input-kind", default="auto", show_default=True,
              type=click.Choice(["auto", "vcf", "bfile"]))
@click.option("--seed", type=int, default=None,
              help="Seed recorded in the receipt")
@click.option("--graphpop-commit", default=None,
              help="GraphPop commit SHA recorded in the receipt")
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def run_plink_grm(input_path, output_dir, input_kind, seed,
                   graphpop_commit, extra_args):
    """PLINK 2.0 `--make-grm-bin` wrapper.

    Example:

      graphpop-bench run plink_grm \\
          --input cohort.vcf.gz --output ./out \\
          -- --maf 0.05
    """
    if not PlinkGrmRunner.is_available():
        raise click.ClickException(
            "PLINK binary not on PATH (looked for plink2, plink); "
            "install plink2 first")
    runner = PlinkGrmRunner()
    result = runner.run(
        Path(input_path), Path(output_dir),
        input_kind=input_kind,
        extra_args=list(extra_args),
        seed=seed,
        graphpop_commit=graphpop_commit,
    )
    click.echo(
        f"plink_grm n_samples={len(result.sample_ids)} "
        f"wall={result.profiling.wall_clock_s:.3f}s "
        f"rss_peak={result.profiling.rss_peak_mb:.1f}MB "
        f"tsv={result.normalised_tsv}",
        err=True,
    )


@run.command("king")
@click.option("--input", "input_path", required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="Input .vcf(.gz) or BED-prefix (no extension)")
@click.option("--output", "output_dir", required=True,
              type=click.Path(file_okay=False),
              help="Output directory (created if missing)")
@click.option("--input-kind", default="auto", show_default=True,
              type=click.Choice(["auto", "vcf", "bfile"]))
@click.option("--mode", default="kinship", show_default=True,
              type=click.Choice(["kinship", "related"]),
              help="KING analysis mode")
@click.option("--related-degree", type=int, default=3, show_default=True,
              help="Degree threshold when --mode=related")
@click.option("--seed", type=int, default=None,
              help="Seed recorded in the receipt")
@click.option("--graphpop-commit", default=None,
              help="GraphPop commit SHA recorded in the receipt")
@click.argument("extra_args", nargs=-1, type=click.UNPROCESSED)
def run_king(input_path, output_dir, input_kind, mode, related_degree,
              seed, graphpop_commit, extra_args):
    """KING-robust kinship inference wrapper (Manichaikul 2010).

    KING outputs the kinship coefficient phi (≈ 0.25 for parent-
    child, ≈ 0 for unrelated); not directly comparable with
    PLINK GRM entries.

    Example:

      graphpop-bench run king --input cohort.vcf.gz --output ./out
    """
    if not KingRunner.is_available():
        raise click.ClickException(
            "KING binary not on PATH; install KING first")
    runner = KingRunner()
    result = runner.run(
        Path(input_path), Path(output_dir),
        input_kind=input_kind,
        mode=mode,
        related_degree=related_degree,
        extra_args=list(extra_args),
        seed=seed,
        graphpop_commit=graphpop_commit,
    )
    click.echo(
        f"king n_samples={len(result.sample_ids)} "
        f"n_pairs={len(result.pair_kinship)} "
        f"wall={result.profiling.wall_clock_s:.3f}s "
        f"rss_peak={result.profiling.rss_peak_mb:.1f}MB "
        f"tsv={result.normalised_tsv}",
        err=True,
    )


if __name__ == "__main__":
    main()
