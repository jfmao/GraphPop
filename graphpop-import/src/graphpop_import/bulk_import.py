"""Orchestrate the full VCF-to-Neo4j bulk import pipeline."""


def run_bulk_import() -> None:
    """Run the end-to-end import: parse VCF, emit CSVs, invoke neo4j-admin import.

    TODO: Implement pipeline orchestration.
    """


def main() -> None:
    """CLI entry point for graphpop-import."""
    run_bulk_import()
