"""Orchestrate the full VCF-to-Neo4j bulk import pipeline."""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

from .csv_emitter import CSVEmitter
from .vcf_parser import VCFParser

logger = logging.getLogger(__name__)


def run_bulk_import(
    vcf_path: Path,
    panel_path: Path,
    out_dir: Path,
    *,
    stratify_by: str = "superpopulation",
    region: str | None = None,
    chunk_size: int = 100_000,
    database: str = "graphpop",
    skip_csv: bool = False,
    skip_import: bool = False,
) -> None:
    """Run the end-to-end import: parse VCF, emit CSVs, invoke neo4j-admin import."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    # -- Step 1: Generate CSVs --
    if not skip_csv:
        logger.info("Parsing VCF and emitting CSVs to %s", out_dir)
        parser = VCFParser(
            vcf_path, panel_path, stratify_by=stratify_by, region=region
        )
        emitter = CSVEmitter.run(parser, out_dir, chunk_size=chunk_size)
        logger.info(
            "CSV generation complete: %d variants, %d carries edges",
            emitter.n_variants,
            emitter.n_carries,
        )
    else:
        logger.info("Skipping CSV generation (--skip-csv)")

    # -- Step 2: Run neo4j-admin import --
    if not skip_import:
        import_script = Path(__file__).resolve().parents[3] / "scripts" / "neo4j-import.sh"
        if not import_script.exists():
            logger.error("Import script not found: %s", import_script)
            logger.info("Run manually: sudo bash scripts/neo4j-import.sh %s %s", out_dir, database)
            return

        logger.info("Running neo4j-admin import via %s", import_script)
        result = subprocess.run(
            ["sudo", "bash", str(import_script), str(out_dir), database],
            check=False,
        )
        if result.returncode != 0:
            logger.error("neo4j-admin import failed with exit code %d", result.returncode)
            sys.exit(result.returncode)
    else:
        logger.info("Skipping neo4j-admin import (--skip-import)")
        logger.info("Run manually: sudo bash scripts/neo4j-import.sh %s %s", out_dir, database)


def main() -> None:
    """CLI entry point for graphpop-import."""
    parser = argparse.ArgumentParser(
        description="GraphPop: VCF to Neo4j bulk import pipeline",
    )
    parser.add_argument("vcf", type=Path, help="Path to VCF/BCF file")
    parser.add_argument("panel", type=Path, help="Path to panel/PED file")
    parser.add_argument(
        "-o", "--out-dir", type=Path, default=Path("csv_out"),
        help="Output directory for CSV files (default: csv_out)",
    )
    parser.add_argument(
        "--stratify-by", choices=["population", "superpopulation"],
        default="superpopulation", help="Population stratification level",
    )
    parser.add_argument("--region", help="Genomic region to parse (e.g., chr22)")
    parser.add_argument(
        "--chunk-size", type=int, default=100_000,
        help="Variants per processing chunk (default: 100000)",
    )
    parser.add_argument(
        "--database", default="graphpop",
        help="Neo4j database name (default: graphpop)",
    )
    parser.add_argument(
        "--skip-csv", action="store_true",
        help="Skip CSV generation (use existing CSVs)",
    )
    parser.add_argument(
        "--skip-import", action="store_true",
        help="Skip neo4j-admin import (generate CSVs only)",
    )

    args = parser.parse_args()

    run_bulk_import(
        vcf_path=args.vcf,
        panel_path=args.panel,
        out_dir=args.out_dir,
        stratify_by=args.stratify_by,
        region=args.region,
        chunk_size=args.chunk_size,
        database=args.database,
        skip_csv=args.skip_csv,
        skip_import=args.skip_import,
    )
