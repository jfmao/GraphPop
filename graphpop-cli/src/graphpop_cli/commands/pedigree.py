"""graphpop pedigree — PRIMUS-style reconstruction (M9, graphpop-pedigree)."""
from __future__ import annotations

from pathlib import Path

import click

from ..cli import pass_ctx


@click.group()
def pedigree():
    """Pedigree reconstruction from :RELATIVE labels (M9).

    Subcommands:

      reconstruct  Build a PED-format pedigree for one or all families.
    """


@pedigree.command("reconstruct")
@click.option("--family-id", default=None,
              help="Restrict output to a single family (default: all)")
@click.option("--source", default="king", show_default=True,
              help=":RELATIVE.source to pull labels from")
@click.option("-o", "--output", "output_path", required=True,
              type=click.Path(dir_okay=False, writable=True),
              help="Output PED file path")
@pass_ctx
def reconstruct(ctx, family_id, source, output_path):
    """Reconstruct PED-format pedigree from :RELATIVE labels (M9).

    Pulls family clusters (:IN_FAMILY) and per-pair :RELATIVE labels
    from Neo4j, runs PRIMUS-lite reconstruction, and writes a PED
    file. Requires `pip install graphpop-pedigree`.
    """
    try:
        from graphpop_pedigree import (  # noqa: WPS433
            PedExporter,
            PedigreeReconstructor,
        )
    except ImportError as exc:
        raise click.ClickException(
            "graphpop-pedigree not installed; run "
            "`pip install graphpop-pedigree` (or `pip install -e ./graphpop-"
            "pedigree` from the repo root)."
        ) from exc

    # Pull samples + their family + sex.
    sample_query = (
        "MATCH (s:Sample)-[:IN_FAMILY {source: $src}]->(f:Family) "
        "RETURN f.family_id AS family_id, s.sampleId AS sid, "
        "s.sex AS sex"
    )
    sample_rows = ctx.run(sample_query, {"src": source})
    if family_id:
        sample_rows = [r for r in sample_rows if r["family_id"] == family_id]
    if not sample_rows:
        raise click.ClickException(
            f"No samples found for source='{source}'"
            + (f", family_id='{family_id}'" if family_id else "")
            + ". Have you run graphpop relate families?"
        )

    edge_query = (
        "MATCH (a:Sample)-[r:RELATIVE {source: $src}]->(b:Sample) "
        "RETURN a.sampleId AS sa, b.sampleId AS sb, "
        "r.relationship AS rel, r.degree AS deg"
    )
    edge_rows = ctx.run(edge_query, {"src": source})

    samples = [
        (r["family_id"], r["sid"], r.get("sex"))
        for r in sample_rows
    ]
    edges = [
        (r["sa"], r["sb"], r["rel"], r["deg"])
        for r in edge_rows
    ]
    rec = PedigreeReconstructor(samples, edges, source=source)
    rows = (rec.reconstruct_family(family_id)
            if family_id else rec.reconstruct_all())
    n = PedExporter.write(rows, output_path)
    click.echo(
        f"Wrote {n} pedigree rows to {output_path} "
        f"(source={source}"
        + (f", family_id={family_id}" if family_id else "")
        + ").",
        err=True,
    )
