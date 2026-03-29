"""graphpop import — import VCF data into a Neo4j graph database."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import click
import yaml


def _get_neo4j_home() -> Path:
    """Get Neo4j home from config."""
    config_path = Path.home() / ".graphpop" / "config.yaml"
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
        if "neo4j_home" in cfg:
            return Path(cfg["neo4j_home"])
    return Path.home() / "neo4j"


@click.command("import")
@click.option("--vcf", required=True, type=click.Path(exists=True),
              help="Input VCF file (bgzipped recommended)")
@click.option("--panel", required=True, type=click.Path(exists=True),
              help="Population panel file (TSV: sample_id, population)")
@click.option("--database", required=True,
              help="Name for the Neo4j database")
@click.option("--vep", type=click.Path(exists=True),
              help="VEP/SnpEff annotation file")
@click.option("--pathways", type=click.Path(exists=True),
              help="Reactome/Plant Reactome pathway file")
@click.option("--go-terms", type=click.Path(exists=True),
              help="GO term annotation file (UniProt GOA format)")
@click.option("--ancestral", type=click.Path(exists=True),
              help="Ancestral allele FASTA (Ensembl EPO)")
@click.option("--csv-dir", type=click.Path(),
              help="Directory for intermediate CSV files (default: temp)")
@click.option("--neo4j-home", type=click.Path(),
              help="Neo4j installation directory")
@click.option("--threads", type=int, default=4,
              help="Import threads (default: 4)")
@click.option("--skip-csv", is_flag=True,
              help="Skip CSV generation (reuse existing CSVs)")
@click.option("--skip-import", is_flag=True,
              help="Skip neo4j-admin import (CSVs only)")
@click.option("--skip-annotations", is_flag=True,
              help="Skip annotation loading")
def import_data(vcf, panel, database, vep, pathways, go_terms, ancestral,
                csv_dir, neo4j_home, threads, skip_csv, skip_import,
                skip_annotations):
    """Import VCF data into a Neo4j graph database.

    This command orchestrates the full import pipeline:

    \b
    1. Parse VCF + panel → generate CSV files (Variant, Sample, Population, etc.)
    2. Run neo4j-admin database import to bulk-load CSVs
    3. Load functional annotations (VEP, pathways, GO terms, ancestral alleles)

    The database name is user-specified and stored in the GraphPop config.

    \b
    Examples:
      graphpop import --vcf data.vcf.gz --panel panel.txt --database myproject
      graphpop import --vcf rice.vcf.gz --panel rice_panel.txt \\
        --database rice3k --vep rice_vep.vcf --pathways plant_reactome.tsv
    """
    neo4j_path = Path(neo4j_home) if neo4j_home else _get_neo4j_home()
    csv_path = Path(csv_dir) if csv_dir else Path(f"/tmp/graphpop_csv_{database}")
    csv_path.mkdir(parents=True, exist_ok=True)

    click.echo(f"GraphPop Import Pipeline")
    click.echo(f"  VCF:       {vcf}")
    click.echo(f"  Panel:     {panel}")
    click.echo(f"  Database:  {database}")
    click.echo(f"  Neo4j:     {neo4j_path}")
    click.echo(f"  CSV dir:   {csv_path}")
    click.echo()

    # Step 1: Generate CSVs
    if not skip_csv:
        click.echo("Step 1/3: Generating CSV files from VCF...")
        _run_csv_generation(vcf, panel, csv_path, threads)
    else:
        click.echo("Step 1/3: Skipping CSV generation (--skip-csv)")

    # Step 2: neo4j-admin import
    if not skip_import:
        click.echo("\nStep 2/3: Running neo4j-admin bulk import...")
        _run_bulk_import(neo4j_path, csv_path, database)
    else:
        click.echo("\nStep 2/3: Skipping bulk import (--skip-import)")

    # Step 3: Load annotations
    if not skip_annotations:
        click.echo("\nStep 3/3: Loading annotations...")
        _load_annotations(neo4j_path, database, vep, pathways, go_terms, ancestral)
    else:
        click.echo("\nStep 3/3: Skipping annotations (--skip-annotations)")

    # Update config with new database
    _update_config(database)

    click.echo(f"""
Import complete!

  Database:  {database}
  GraphPop config updated to use database '{database}'.

Next steps:
  graphpop start                                     # Start Neo4j (if not running)
  graphpop db info                                   # Verify node/edge counts
  graphpop diversity chr1 1 50000000 YOUR_POP        # Run first analysis
  graphpop run-all --database {database} -d results/ # Full-genome analysis
""")


def _run_csv_generation(vcf: str, panel: str, csv_dir: Path, threads: int):
    """Run the graphpop-import CSV generation."""
    try:
        import importlib
        # Try importing graphpop_import directly
        spec = importlib.util.find_spec("graphpop_import")
        if spec:
            click.echo("  Using graphpop-import Python package...")
            from graphpop_import.vcf_parser import VCFParser
            from graphpop_import.csv_emitter import CSVEmitter
            parser = VCFParser(vcf, panel)
            emitter = CSVEmitter(str(csv_dir))
            parser.parse(emitter)
            click.echo(f"  CSVs written to {csv_dir}")
            return
    except ImportError:
        pass

    # Fallback: run as subprocess
    click.echo("  Running graphpop-import as subprocess...")
    scripts = [
        Path("graphpop-import/src/graphpop_import/vcf_parser.py"),
        Path("scripts/rice_csv_parallel.py"),
    ]
    for script in scripts:
        if script.exists():
            result = subprocess.run(
                [sys.executable, str(script),
                 "--vcf", vcf, "--panel", panel, "--output", str(csv_dir),
                 "--threads", str(threads)],
                capture_output=True, text=True,
            )
            if result.returncode == 0:
                click.echo(f"  CSVs written to {csv_dir}")
                return
            else:
                click.echo(f"  Warning: {result.stderr[:200]}", err=True)

    click.echo(
        "  Error: graphpop-import not found.\n"
        "  Install with: pip install -e graphpop-import/\n"
        "  Or generate CSVs manually and use --skip-csv",
        err=True,
    )
    raise SystemExit(1)


def _run_bulk_import(neo4j_home: Path, csv_dir: Path, database: str):
    """Run neo4j-admin database import."""
    admin_bin = neo4j_home / "bin" / "neo4j-admin"
    if not admin_bin.exists():
        click.echo(f"  Error: neo4j-admin not found at {admin_bin}", err=True)
        raise SystemExit(1)

    # Check if database already exists
    db_dir = neo4j_home / "data" / "databases" / database
    if db_dir.exists():
        if not click.confirm(f"  Database '{database}' already exists. Overwrite?"):
            click.echo("  Import cancelled.")
            raise SystemExit(0)

    # Build neo4j-admin import command
    cmd = [
        str(admin_bin), "database", "import", "full",
        f"--nodes=Variant={csv_dir}/variant_header.csv,{csv_dir}/variants_*.csv",
        f"--nodes=Sample={csv_dir}/sample_header.csv,{csv_dir}/samples.csv",
        f"--nodes=Population={csv_dir}/population_header.csv,{csv_dir}/populations.csv",
        f"--nodes=Chromosome={csv_dir}/chromosome_header.csv,{csv_dir}/chromosomes.csv",
        f"--relationships=NEXT={csv_dir}/next_header.csv,{csv_dir}/next_*.csv",
        f"--relationships=ON_CHROMOSOME={csv_dir}/on_chromosome_header.csv,{csv_dir}/on_chromosome_*.csv",
        f"--relationships=IN_POPULATION={csv_dir}/in_population_header.csv,{csv_dir}/in_population.csv",
        "--overwrite-destination=true",
        database,
    ]

    click.echo(f"  Running: neo4j-admin database import {database}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        click.echo(f"  Import failed: {result.stderr[:500]}", err=True)
        click.echo("  You may need to stop Neo4j first: graphpop stop", err=True)
        raise SystemExit(1)
    click.echo("  Bulk import complete.")


def _load_annotations(neo4j_home: Path, database: str,
                      vep: str | None, pathways: str | None,
                      go_terms: str | None, ancestral: str | None):
    """Load functional annotations via Cypher transactions."""
    if not any([vep, pathways, go_terms, ancestral]):
        click.echo("  No annotations specified, skipping.")
        return

    # Load annotations by running the appropriate Python scripts
    scripts_dir = Path("scripts")
    annotation_scripts = []

    if vep:
        click.echo(f"  Loading VEP annotations from {vep}...")
        annotation_scripts.append(("load_annotations", ["--vep", vep]))

    if pathways:
        click.echo(f"  Loading pathway annotations from {pathways}...")
        annotation_scripts.append(("load_annotations", ["--pathways", pathways]))

    if go_terms:
        click.echo(f"  Loading GO term annotations from {go_terms}...")
        annotation_scripts.append(("load_annotations", ["--go", go_terms]))

    if ancestral:
        click.echo(f"  Loading ancestral alleles from {ancestral}...")
        annotation_scripts.append(("load_annotations", ["--ancestral", ancestral]))

    for script_name, args in annotation_scripts:
        # Try to find the annotation loading script
        candidates = [
            scripts_dir / f"{script_name}.py",
            scripts_dir / "load_rice_annotations.py",
            Path(f"graphpop-import/src/graphpop_import/{script_name}.py"),
        ]
        for script in candidates:
            if script.exists():
                result = subprocess.run(
                    [sys.executable, str(script), "--database", database] + args,
                    capture_output=True, text=True,
                )
                if result.returncode == 0:
                    click.echo(f"  Loaded: {script_name}")
                    break
                else:
                    click.echo(f"  Warning: {result.stderr[:200]}", err=True)
        else:
            click.echo(f"  Annotation script not found for: {script_name}")
            click.echo("  You can load annotations manually after import.")


def _update_config(database: str):
    """Update GraphPop config to use the new database."""
    config_path = Path.home() / ".graphpop" / "config.yaml"
    cfg = {}
    if config_path.exists():
        with open(config_path) as f:
            cfg = yaml.safe_load(f) or {}
    cfg["database"] = database
    config_path.parent.mkdir(exist_ok=True)
    with open(config_path, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False)
    click.echo(f"  Config updated: database = {database}")
