"""graphpop-sim CLI (M12).

Three subcommands: msprime / slim / abc.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import click

from .abc import load_priors, posterior_to_tsv, run_abc
from .msprime_runner import MsprimeConfig, MsprimeRunner
from .slim_runner import SlimRunner
from .summary_stats import SummaryStats, compute_summary_stats


@click.group()
def main():
    """Simulation integration for GraphPop (Phase 6)."""


def _build_driver():
    from neo4j import GraphDatabase

    uri = os.environ.get("GRAPHPOP_URI", "bolt://localhost:7687")
    user = os.environ.get("GRAPHPOP_USER", "neo4j")
    password = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
    return GraphDatabase.driver(uri, auth=(user, password))


@main.command()
@click.argument("config_path",
                type=click.Path(exists=True, dir_okay=False))
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True),
              help="Output .trees file path")
@click.option("--run-id", default=None,
              help="ARGRun id when --ingest is set")
@click.option("--ingest", is_flag=True,
              help="Ingest the resulting tree sequence as a :ARGRun")
def msprime(config_path: str, output: str, run_id: str | None,
             ingest: bool) -> None:
    """Run msprime from a YAML config; dump .trees + optional ingest."""
    cfg = MsprimeConfig.from_yaml(config_path)
    ts = MsprimeRunner(cfg).simulate_and_dump(output)
    click.echo(
        f"Simulated {ts.num_samples} samples, {ts.num_trees} trees, "
        f"{ts.num_mutations} mutations → {output}",
        err=True,
    )
    if ingest:
        if run_id is None:
            raise click.UsageError("--run-id is required with --ingest")
        try:
            from graphpop_import.arg_importer import ARGIngester
        except ImportError as exc:
            raise click.ClickException(
                "graphpop-import not installed; "
                "`pip install graphpop-sim[ingest]` + the import package."
            ) from exc
        driver = _build_driver()
        try:
            ing = ARGIngester(driver)
            summary = ing.ingest(ts, run_id=run_id, source="msprime",
                                  params=cfg.to_dict())
            click.echo(
                f"Ingested ARGRun {summary.run_id}: {summary.n_nodes} "
                f"TreeNodes, {summary.n_edges} edges.",
                err=True,
            )
        finally:
            driver.close()


@main.command()
@click.option("--template", required=True,
              type=click.Choice(SlimRunner.list_templates()),
              help="Curated SLiM template")
@click.option("--params", default="",
              help='Comma-separated key=value overrides '
                   '(e.g. "Ne=5000,sweep_s=0.01")')
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True),
              help="Output .trees file path")
@click.option("--run-id", default=None,
              help="ARGRun id when --ingest is set")
@click.option("--ingest", is_flag=True,
              help="Ingest the resulting tree sequence as a :ARGRun")
def slim(template: str, params: str, output: str,
          run_id: str | None, ingest: bool) -> None:
    """Run a SLiM template; dump .trees + optional ingest."""
    if not SlimRunner.is_available():
        raise click.ClickException(
            "SLiM binary not found on PATH; install SLiM first.")
    parsed: dict[str, object] = {}
    for pair in params.split(",") if params else []:
        if not pair.strip():
            continue
        k, _, v = pair.partition("=")
        try:
            parsed[k.strip()] = int(v)
        except ValueError:
            try:
                parsed[k.strip()] = float(v)
            except ValueError:
                parsed[k.strip()] = v.strip()
    runner = SlimRunner()
    result = runner.run(template=template, params=parsed, output_path=output)
    if result.return_code != 0:
        raise click.ClickException(
            f"SLiM exited with code {result.return_code}")
    click.echo(f"SLiM run complete → {output}", err=True)

    if ingest:
        if run_id is None:
            raise click.UsageError("--run-id is required with --ingest")
        try:
            import tskit
            from graphpop_import.arg_importer import ARGIngester
        except ImportError as exc:
            raise click.ClickException(
                "graphpop-import + tskit required for --ingest."
            ) from exc
        ts = tskit.load(output)
        driver = _build_driver()
        try:
            ing = ARGIngester(driver)
            summary = ing.ingest(ts, run_id=run_id, source="slim",
                                  params={"template": template,
                                          "params": parsed})
            click.echo(
                f"Ingested ARGRun {summary.run_id}: {summary.n_nodes} "
                f"TreeNodes, {summary.n_edges} edges.",
                err=True,
            )
        finally:
            driver.close()


@main.command()
@click.option("--prior", "prior_path", required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="YAML file specifying per-parameter prior distributions")
@click.option("--observed", "observed_path", required=True,
              type=click.Path(exists=True, dir_okay=False),
              help="TSV with observed summary stats "
                   "(headers: pi, theta_w, tajimas_d, ...)")
@click.option("--base-config", "base_config_path", default=None,
              type=click.Path(exists=True, dir_okay=False),
              help="YAML msprime config used as base for sampled draws")
@click.option("--n-sim", type=int, default=1000, show_default=True)
@click.option("--epsilon", type=float, default=0.05, show_default=True,
              help="Fraction of closest simulations to accept")
@click.option("--summary-keys", default="pi,theta_w,tajimas_d",
              show_default=True,
              help="Comma-separated stats to use as the distance vector")
@click.option("--seed", type=int, default=42, show_default=True)
@click.option("--output", "-o", required=True,
              type=click.Path(dir_okay=False, writable=True),
              help="Posterior TSV path")
def abc(prior_path: str, observed_path: str,
         base_config_path: str | None, n_sim: int,
         epsilon: float, summary_keys: str, seed: int,
         output: str) -> None:
    """Run rejection-ABC."""
    priors = load_priors(prior_path)
    obs = _load_observed(observed_path)
    keys = [s.strip() for s in summary_keys.split(",") if s.strip()]
    base = (MsprimeConfig.from_yaml(base_config_path)
            if base_config_path else None)
    result = run_abc(
        priors=priors, observed=obs, n_sim=n_sim, epsilon=epsilon,
        summary_keys=keys, base_config=base, seed=seed,
    )
    n = posterior_to_tsv(result, output)
    click.echo(
        f"ABC: ran {result.n_total} sims, accepted {result.n_accepted} "
        f"(ε = {result.epsilon}). Wrote {n} posterior rows → {output}",
        err=True,
    )


def _load_observed(path: str) -> SummaryStats:
    """Parse a TSV with header row → SummaryStats dataclass."""
    with open(path) as fh:
        lines = [line.rstrip() for line in fh if line.strip()]
    if not lines:
        raise click.ClickException(
            f"observed stats file is empty: {path}")
    header = lines[0].split("\t")
    row = lines[1].split("\t") if len(lines) > 1 else []
    d = dict(zip(header, row))

    def to_float(key: str, default: float = 0.0) -> float:
        v = d.get(key)
        return float(v) if v not in (None, "") else default

    fst_raw = d.get("mean_fst")
    return SummaryStats(
        pi=to_float("pi"),
        theta_w=to_float("theta_w"),
        tajimas_d=to_float("tajimas_d"),
        mean_fst=(float(fst_raw) if fst_raw not in (None, "") else None),
        n_samples=int(to_float("n_samples", 0)),
        n_segregating_sites=int(to_float("n_segregating_sites", 0)),
        sequence_length=int(to_float("sequence_length", 0)),
    )


if __name__ == "__main__":
    main()
