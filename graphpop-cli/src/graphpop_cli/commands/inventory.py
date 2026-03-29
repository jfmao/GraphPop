"""graphpop inventory — comprehensive database inventory."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..formatters import format_output


@click.command("inventory")
@click.option("-o", "--output", "output_path", help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def inventory(ctx, output_path, fmt):
    """Show comprehensive database inventory.

    Reports node/relationship counts, populations, chromosomes, loaded
    annotations, and persisted statistics. No arguments needed.

    \b
    Examples:
      graphpop inventory
      graphpop inventory --format json -o db_inventory.json
    """
    sections = []

    # --- 1. Node label counts ---
    click.echo("Querying node counts...", err=True)
    labels = ["Variant", "Sample", "Population", "Gene", "Pathway",
              "GOTerm", "GenomicWindow"]
    for label in labels:
        recs = ctx.run(f"MATCH (n:{label}) RETURN count(n) AS count")
        count = recs[0]["count"] if recs else 0
        sections.append({"section": "nodes", "item": label, "value": str(count)})

    # --- 2. Relationship type counts ---
    click.echo("Querying relationship counts...", err=True)
    rel_types = ["CARRIES", "HAS_CONSEQUENCE", "IN_PATHWAY", "HAS_GO_TERM",
                 "NEXT", "LD", "BELONGS_TO"]
    for rel in rel_types:
        recs = ctx.run(f"MATCH ()-[r:{rel}]->() RETURN count(r) AS count")
        count = recs[0]["count"] if recs else 0
        sections.append({"section": "relationships", "item": rel, "value": str(count)})

    # --- 3. Populations and sample counts ---
    click.echo("Querying populations...", err=True)
    recs = ctx.run(
        "MATCH (p:Population) "
        "OPTIONAL MATCH (s:Sample)-[:BELONGS_TO]->(p) "
        "RETURN p.popId AS population, count(s) AS sample_count "
        "ORDER BY p.popId"
    )
    for rec in recs:
        sections.append({
            "section": "populations",
            "item": rec["population"],
            "value": str(rec["sample_count"]),
        })

    # --- 4. Chromosomes and variant counts ---
    click.echo("Querying chromosomes...", err=True)
    recs = ctx.run(
        "MATCH (v:Variant) "
        "RETURN v.chr AS chr, count(v) AS variant_count "
        "ORDER BY v.chr"
    )
    for rec in recs:
        sections.append({
            "section": "chromosomes",
            "item": rec["chr"],
            "value": str(rec["variant_count"]),
        })

    # --- 5. Annotation coverage ---
    click.echo("Querying annotations...", err=True)
    # HAS_CONSEQUENCE edges
    recs = ctx.run("MATCH ()-[r:HAS_CONSEQUENCE]->() RETURN count(r) AS count")
    has_conseq = recs[0]["count"] if recs else 0
    sections.append({"section": "annotations", "item": "HAS_CONSEQUENCE edges",
                     "value": str(has_conseq)})

    # IN_PATHWAY edges
    recs = ctx.run("MATCH ()-[r:IN_PATHWAY]->() RETURN count(r) AS count")
    has_pw = recs[0]["count"] if recs else 0
    sections.append({"section": "annotations", "item": "IN_PATHWAY edges",
                     "value": str(has_pw)})

    # HAS_GO_TERM edges
    recs = ctx.run("MATCH ()-[r:HAS_GO_TERM]->() RETURN count(r) AS count")
    has_go = recs[0]["count"] if recs else 0
    sections.append({"section": "annotations", "item": "HAS_GO_TERM edges",
                     "value": str(has_go)})

    # Ancestral allele coverage
    recs = ctx.run(
        "MATCH (v:Variant) WHERE v.ancestral_allele IS NOT NULL "
        "RETURN count(v) AS count"
    )
    aa_count = recs[0]["count"] if recs else 0
    sections.append({"section": "annotations", "item": "variants_with_ancestral_allele",
                     "value": str(aa_count)})

    # --- 6. Persisted statistics ---
    click.echo("Querying persisted statistics...", err=True)
    # Check for ihs/xpehh/nsl properties on Variant nodes (sample a few)
    for stat_prefix in ["ihs_", "xpehh_", "nsl_"]:
        recs = ctx.run(
            f"MATCH (v:Variant) "
            f"WITH v LIMIT 1 "
            f"UNWIND keys(v) AS k "
            f"WITH k WHERE k STARTS WITH '{stat_prefix}' "
            f"RETURN COLLECT(DISTINCT k) AS props"
        )
        props = recs[0]["props"] if recs and recs[0]["props"] else []
        if props:
            for p in props:
                sections.append({"section": "persisted_stats", "item": p,
                                 "value": "on Variant nodes"})
        else:
            sections.append({"section": "persisted_stats",
                             "item": f"{stat_prefix}*",
                             "value": "none found"})

    # GenomicWindow statistics
    recs = ctx.run("MATCH (w:GenomicWindow) RETURN count(w) AS count")
    gw_count = recs[0]["count"] if recs else 0
    sections.append({"section": "persisted_stats", "item": "GenomicWindow count",
                     "value": str(gw_count)})

    # Fst properties on GenomicWindow
    recs = ctx.run(
        "MATCH (w:GenomicWindow) "
        "WITH w LIMIT 1 "
        "UNWIND keys(w) AS k "
        "WITH k WHERE k STARTS WITH 'fst_' "
        "RETURN COLLECT(DISTINCT k) AS props"
    )
    fst_props = recs[0]["props"] if recs and recs[0]["props"] else []
    for p in fst_props:
        sections.append({"section": "persisted_stats", "item": p,
                         "value": "on GenomicWindow nodes"})

    # --- Print summary ---
    if fmt == "tsv" and not output_path:
        _print_inventory(sections)
    else:
        format_output(sections, output_path, fmt, "inventory", {})


def _print_inventory(sections: list[dict]):
    """Pretty-print inventory to stderr/stdout."""
    current_section = None
    for row in sections:
        sec = row["section"]
        if sec != current_section:
            current_section = sec
            click.echo(f"\n=== {sec.upper().replace('_', ' ')} ===")
        click.echo(f"  {row['item']:40s} {row['value']}")
    click.echo()
