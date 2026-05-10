"""graphpop sim — simulation integration (M12, Phase 6).

Thin facade over graphpop-sim. Same pattern as `graphpop gnn`
and `graphpop pedigree`.
"""
from __future__ import annotations

import subprocess
import sys

import click

from ..cli import pass_ctx


@click.group()
def sim():
    """Simulate populations and run ABC inference (Phase 6).

    Subcommands:

      msprime  Run msprime from a YAML config; optional :ARGRun ingest.
      slim     Run a curated SLiM template; optional :ARGRun ingest.
      abc      Rejection-ABC: simulate, summarise, accept top-ε.
    """


@sim.command("msprime")
@click.argument("config_path",
                type=click.Path(exists=True, dir_okay=False))
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True))
@click.option("--run-id", default=None,
              help="ARGRun id when --ingest is set")
@click.option("--ingest", is_flag=True,
              help="Ingest the resulting tree sequence")
@pass_ctx
def msprime_cmd(ctx, config_path, output, run_id, ingest):
    """Delegate to graphpop-sim msprime."""
    args = [
        sys.executable, "-m", "graphpop_sim.cli", "msprime",
        config_path, "--output", output,
    ]
    if run_id:
        args += ["--run-id", run_id]
    if ingest:
        args += ["--ingest"]
    raise SystemExit(subprocess.call(args))


@sim.command("slim")
@click.option("--template", required=True,
              type=click.Choice(["neutral", "sweep_genic", "bottleneck"]))
@click.option("--params", default="",
              help='Comma-separated key=value overrides')
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True))
@click.option("--run-id", default=None)
@click.option("--ingest", is_flag=True)
@pass_ctx
def slim_cmd(ctx, template, params, output, run_id, ingest):
    """Delegate to graphpop-sim slim."""
    args = [
        sys.executable, "-m", "graphpop_sim.cli", "slim",
        "--template", template,
        "--params", params,
        "--output", output,
    ]
    if run_id:
        args += ["--run-id", run_id]
    if ingest:
        args += ["--ingest"]
    raise SystemExit(subprocess.call(args))


@sim.command("abc")
@click.option("--prior", "prior_path", required=True,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--observed", "observed_path", required=True,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--base-config", "base_config_path", default=None,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--n-sim", type=int, default=1000, show_default=True)
@click.option("--epsilon", type=float, default=0.05, show_default=True)
@click.option("--summary-keys", default="pi,theta_w,tajimas_d",
              show_default=True)
@click.option("--seed", type=int, default=42, show_default=True)
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True))
@pass_ctx
def abc_cmd(ctx, prior_path, observed_path, base_config_path,
             n_sim, epsilon, summary_keys, seed, output):
    """Delegate to graphpop-sim abc."""
    args = [
        sys.executable, "-m", "graphpop_sim.cli", "abc",
        "--prior", prior_path,
        "--observed", observed_path,
        "--n-sim", str(n_sim),
        "--epsilon", str(epsilon),
        "--summary-keys", summary_keys,
        "--seed", str(seed),
        "--output", output,
    ]
    if base_config_path:
        args += ["--base-config", base_config_path]
    raise SystemExit(subprocess.call(args))
