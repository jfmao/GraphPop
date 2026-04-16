"""GraphPop CLI — command-line interface for graph-native population genomics.

Usage:
    graphpop diversity chr22 1 50000000 EUR -o diversity.tsv
    graphpop ihs chr22 EUR --min-af 0.05 --persist -o ihs.tsv
    graphpop genome-scan chr22 EUR 100000 50000 --persist -o scan.tsv
"""
from __future__ import annotations

import sys
from pathlib import Path

import click

from .connection import load_config, get_driver


class GraphPopContext:
    """Shared context passed to all commands."""
    def __init__(self, uri=None, user=None, password=None, database=None,
                 config_path=None):
        cfg = load_config(Path(config_path) if config_path else None)
        if uri:
            cfg["uri"] = uri
        if user:
            cfg["user"] = user
        if password:
            cfg["password"] = password
        if database:
            cfg["database"] = database
        self.cfg = cfg
        self._driver = None

    @property
    def driver(self):
        if self._driver is None:
            self._driver = get_driver(self.cfg)
        return self._driver

    @property
    def database(self):
        return self.cfg["database"]

    def run(self, cypher: str, parameters: dict | None = None) -> list[dict]:
        """Run Cypher and return records as list of dicts."""
        try:
            with self.driver.session(database=self.database) as session:
                return [rec.data() for rec in session.run(cypher, parameters)]
        except Exception as e:
            err_msg = str(e)
            if "Connection refused" in err_msg or "Failed to establish" in err_msg:
                click.echo(
                    "Error: Cannot connect to Neo4j at "
                    f"{self.cfg['uri']}.\n"
                    "Is Neo4j running? Check connection with:\n"
                    f"  export GRAPHPOP_URI={self.cfg['uri']}\n"
                    "  or create ~/.graphpop/config.yaml",
                    err=True,
                )
            else:
                click.echo(f"Error: {e}", err=True)
            raise SystemExit(1)

    def close(self):
        if self._driver:
            self._driver.close()


pass_ctx = click.make_pass_decorator(GraphPopContext, ensure=True)


@click.group()
@click.option("--uri", envvar="GRAPHPOP_URI", help="Neo4j bolt URI")
@click.option("--user", envvar="GRAPHPOP_USER", help="Neo4j username")
@click.option("--password", envvar="GRAPHPOP_PASSWORD", help="Neo4j password")
@click.option("--database", envvar="GRAPHPOP_DATABASE", help="Neo4j database name")
@click.option("--config", "config_path", type=click.Path(),
              help="Config file path (default: ~/.graphpop/config.yaml)")
@click.version_option(package_name="graphpop-cli")
@click.pass_context
def main(ctx, uri, user, password, database, config_path):
    """GraphPop — graph-native population genomics from the command line.

    Compute population genetics statistics via Neo4j stored procedures with
    default TSV output. Use --persist to write results to graph nodes.
    """
    ctx.ensure_object(dict)
    ctx.obj = GraphPopContext(uri=uri, user=user, password=password,
                              database=database, config_path=config_path)


# Import all command modules
from .commands import (  # noqa: E402
    diversity, divergence, sfs, joint_sfs,
    genome_scan, pop_summary,
    ld, ihs, xpehh, nsl, roh, garud_h,
    query, run_all, aggregate, export_windows,
    setup, server, doctor, db, import_data, dump,
    config_cmd, validate, filter_results, plot,
    lookup, converge, inventory, rank_genes,
    extract, export_bed, batch, compare,
    report, neighbors,
)

# Individual procedures (12)
main.add_command(diversity.diversity)
main.add_command(divergence.divergence)
main.add_command(sfs.sfs)
main.add_command(joint_sfs.joint_sfs)
main.add_command(genome_scan.genome_scan)
main.add_command(pop_summary.pop_summary)
main.add_command(ld.ld)
main.add_command(ihs.ihs)
main.add_command(xpehh.xpehh)
main.add_command(nsl.nsl)
main.add_command(roh.roh)
main.add_command(garud_h.garud_h)

# Orchestration and export
main.add_command(run_all.run_all)
main.add_command(aggregate.aggregate)
main.add_command(export_windows.export_windows)
main.add_command(query.query)
main.add_command(filter_results.filter_results)

# Setup and server management
main.add_command(setup.setup)
main.add_command(server.start)
main.add_command(server.stop)
main.add_command(server.status)
main.add_command(doctor.doctor)

# Database management
main.add_command(db.db)
main.add_command(import_data.import_data)
main.add_command(dump.dump)
main.add_command(dump.load)

# Configuration and validation
main.add_command(config_cmd.config)
main.add_command(validate.validate)
main.add_command(plot.plot)

# Phase 1 high-priority commands
main.add_command(lookup.lookup)
main.add_command(converge.converge)
main.add_command(inventory.inventory)
main.add_command(rank_genes.rank_genes)

# Phase 2 commands
main.add_command(extract.extract)
main.add_command(export_bed.export_bed)
main.add_command(batch.batch)
main.add_command(compare.compare)

# Phase 3 commands
main.add_command(report.report)
main.add_command(neighbors.neighbors)


if __name__ == "__main__":
    main()
