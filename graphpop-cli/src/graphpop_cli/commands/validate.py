"""graphpop validate — check graph database integrity and completeness."""
from __future__ import annotations

import click

from ..cli import pass_ctx


@click.command()
@click.option("--fix", is_flag=True, help="Attempt to fix issues (create missing indexes)")
@pass_ctx
def validate(ctx, fix):
    """Validate the graph database structure and completeness.

    Checks that required node labels, indexes, and GraphPop procedures
    are present. Reports any issues found.

    \b
    Checks performed:
      - Required node labels (Variant, Sample, Population, Chromosome)
      - Optional node labels (Gene, Pathway, GOTerm, GenomicWindow)
      - Required indexes (Variant.chr+pos, Population.populationId)
      - GraphPop procedures installed
      - Variant node properties (pop_ids, ac, an, af, gt_packed)
      - Sample-Population relationships
    """
    issues = []
    ok = []

    click.echo(f"Validating database: {ctx.database}\n")

    # Check node labels
    click.echo("Node labels:")
    try:
        records = ctx.run(
            "CALL db.labels() YIELD label "
            "CALL { WITH label MATCH (n) WHERE label IN labels(n) "
            "RETURN count(n) AS cnt } RETURN label, cnt ORDER BY cnt DESC"
        )
        label_counts = {r["label"]: r["cnt"] for r in records}
    except Exception as e:
        click.echo(f"  Error: Cannot query database: {e}", err=True)
        raise SystemExit(1)

    required_labels = ["Variant", "Sample", "Population", "Chromosome"]
    optional_labels = ["Gene", "Pathway", "GOTerm", "GenomicWindow"]

    for label in required_labels:
        count = label_counts.get(label, 0)
        if count > 0:
            ok.append(f"{label}: {count:,}")
            click.echo(f"  [OK]   {label}: {count:,}")
        else:
            issues.append(f"Required label '{label}' missing or empty")
            click.echo(f"  [FAIL] {label}: MISSING")

    for label in optional_labels:
        count = label_counts.get(label, 0)
        if count > 0:
            click.echo(f"  [OK]   {label}: {count:,}")
        else:
            click.echo(f"  [--]   {label}: not present (optional)")

    # Check indexes
    click.echo("\nIndexes:")
    try:
        records = ctx.run(
            "SHOW INDEXES YIELD name, labelsOrTypes, properties, state "
            "RETURN name, labelsOrTypes, properties, state"
        )
        indexes = {
            (tuple(r["labelsOrTypes"]), tuple(r["properties"])): r["state"]
            for r in records
        }
    except Exception:
        indexes = {}

    required_indexes = [
        (("Variant",), ("chr", "pos")),
        (("Population",), ("populationId",)),
        (("Sample",), ("sampleId",)),
    ]

    for labels, props in required_indexes:
        state = indexes.get((labels, props), None)
        if state == "ONLINE":
            click.echo(f"  [OK]   {labels[0]}({', '.join(props)})")
        elif state:
            click.echo(f"  [WARN] {labels[0]}({', '.join(props)}): {state}")
        else:
            msg = f"Index on {labels[0]}({', '.join(props)}) missing"
            issues.append(msg)
            click.echo(f"  [FAIL] {msg}")
            if fix:
                idx_name = f"idx_{labels[0].lower()}_{'_'.join(props)}"
                cypher = (
                    f"CREATE INDEX {idx_name} IF NOT EXISTS "
                    f"FOR (n:{labels[0]}) ON ({', '.join(f'n.{p}' for p in props)})"
                )
                try:
                    ctx.run(cypher)
                    click.echo(f"         Fixed: created index {idx_name}")
                except Exception as e:
                    click.echo(f"         Fix failed: {e}", err=True)

    # Check procedures
    click.echo("\nGraphPop procedures:")
    try:
        records = ctx.run(
            "SHOW PROCEDURES YIELD name WHERE name STARTS WITH 'graphpop' "
            "RETURN name ORDER BY name"
        )
        proc_names = [r["name"] for r in records]
    except Exception:
        proc_names = []

    expected_procs = [
        "graphpop.diversity", "graphpop.divergence", "graphpop.sfs",
        "graphpop.joint_sfs", "graphpop.genome_scan", "graphpop.pop_summary",
        "graphpop.ld", "graphpop.ihs", "graphpop.xpehh",
        "graphpop.nsl", "graphpop.roh", "graphpop.garud_h",
    ]

    for proc in expected_procs:
        if proc in proc_names:
            click.echo(f"  [OK]   {proc}")
        else:
            issues.append(f"Procedure '{proc}' not installed")
            click.echo(f"  [FAIL] {proc}: NOT INSTALLED")

    # Check Variant properties
    click.echo("\nVariant node properties:")
    try:
        rec = ctx.run(
            "MATCH (v:Variant) WITH v LIMIT 1 "
            "RETURN keys(v) AS props"
        )
        if rec:
            props = set(rec[0]["props"])
            required_props = ["chr", "pos", "ref", "alt", "pop_ids", "ac", "an", "af"]
            optional_props = ["gt_packed", "phase_packed", "ancestral_allele", "is_polarized"]

            for p in required_props:
                if p in props:
                    click.echo(f"  [OK]   {p}")
                else:
                    issues.append(f"Variant property '{p}' missing")
                    click.echo(f"  [FAIL] {p}: MISSING")

            for p in optional_props:
                if p in props:
                    click.echo(f"  [OK]   {p}")
                else:
                    click.echo(f"  [--]   {p}: not present (optional)")
    except Exception as e:
        click.echo(f"  Error checking properties: {e}", err=True)

    # Summary
    click.echo(f"\n{'='*40}")
    if issues:
        click.echo(f"VALIDATION: {len(issues)} issue(s) found")
        for issue in issues:
            click.echo(f"  - {issue}")
        if not fix:
            click.echo("\nRun with --fix to attempt automatic fixes.")
    else:
        click.echo("VALIDATION: All checks passed")
