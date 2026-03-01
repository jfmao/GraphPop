#!/usr/bin/env python3
"""Annotate Variant nodes with ancestral allele from Ensembl EPO FASTA.

Reads the Ensembl EPO 6-primate ancestral allele FASTA for a chromosome,
compares each variant's REF/ALT to the ancestral base, and writes
`ancestral_allele` ('REF' or 'ALT') and `is_polarized` (boolean) properties.

Usage:
    conda run -n graphevo python scripts/annotate_ancestral.py

Requirements:
    - neo4j Python driver (`pip install neo4j`)
    - Ancestral FASTA at data/raw/ensembl_ancestor/homo_sapiens_ancestor_22.fa
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-5s %(message)s",
)
logger = logging.getLogger(__name__)


def load_ancestral_fasta(fasta_path: Path) -> str:
    """Load an Ensembl ancestral FASTA file into a single string.

    The FASTA uses uppercase for high-confidence EPO alignment and
    lowercase for low-confidence. '.' or 'N' indicates unknown.
    Returns the sequence as-is (preserving case) with 0-based indexing.
    """
    logger.info("Loading ancestral FASTA: %s", fasta_path)
    lines = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            lines.append(line.strip())
    seq = "".join(lines)
    logger.info("Loaded %d positions", len(seq))
    return seq


def annotate_chromosome(
    uri: str,
    user: str,
    password: str,
    chrom: str,
    ancestral_seq: str,
    batch_size: int = 10000,
) -> dict[str, int]:
    """Annotate all Variant nodes on a chromosome with ancestral allele.

    Returns counts: {annotated, skipped_unknown, skipped_no_match, total}.
    """
    from neo4j import GraphDatabase

    driver = GraphDatabase.driver(uri, auth=(user, password))

    stats = {"total": 0, "annotated": 0, "skipped_unknown": 0, "skipped_no_match": 0}

    with driver.session() as session:
        # Fetch all variants on the chromosome
        result = session.run(
            "MATCH (v:Variant) WHERE v.chr = $chr "
            "RETURN v.variantId AS vid, v.pos AS pos, v.ref AS ref, v.alt AS alt "
            "ORDER BY v.pos",
            chr=chrom,
        )

        batch = []
        for record in result:
            stats["total"] += 1
            vid = record["vid"]
            pos = record["pos"]
            ref_allele = record["ref"].upper()
            alt_allele = record["alt"].upper()

            # FASTA is 0-indexed, VCF positions are 1-indexed
            fasta_idx = pos - 1
            if fasta_idx < 0 or fasta_idx >= len(ancestral_seq):
                stats["skipped_unknown"] += 1
                continue

            ancestral_char = ancestral_seq[fasta_idx]
            ancestral_upper = ancestral_char.upper()

            # Unknown ancestral allele
            if ancestral_upper in (".", "N", "-"):
                stats["skipped_unknown"] += 1
                continue

            # Determine polarization
            is_high_confidence = ancestral_char.isupper()

            if ancestral_upper == ref_allele[0]:  # Compare first base for indels
                batch.append({
                    "vid": vid,
                    "ancestral_allele": "REF",
                    "is_polarized": is_high_confidence,
                })
                stats["annotated"] += 1
            elif ancestral_upper == alt_allele[0]:
                batch.append({
                    "vid": vid,
                    "ancestral_allele": "ALT",
                    "is_polarized": is_high_confidence,
                })
                stats["annotated"] += 1
            else:
                stats["skipped_no_match"] += 1

            # Write in batches
            if len(batch) >= batch_size:
                _write_batch(session, batch)
                batch = []

        if batch:
            _write_batch(session, batch)

    driver.close()
    return stats


def _write_batch(session, batch: list[dict]) -> None:
    """Write a batch of ancestral allele annotations to Neo4j."""
    session.run(
        """
        UNWIND $batch AS item
        MATCH (v:Variant {variantId: item.vid})
        SET v.ancestral_allele = item.ancestral_allele,
            v.is_polarized = item.is_polarized
        """,
        batch=batch,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Annotate Variant nodes with ancestral allele from Ensembl EPO FASTA"
    )
    parser.add_argument(
        "--chr", default="chr22",
        help="Chromosome name (must match Variant.chr in the database)",
    )
    parser.add_argument(
        "--fasta",
        default="data/raw/ensembl_ancestor/homo_sapiens_ancestor_22.fa",
        help="Path to ancestral FASTA file",
    )
    parser.add_argument("--uri", default="bolt://localhost:7687")
    parser.add_argument("--user", default="neo4j")
    parser.add_argument("--password", default="graphpop")
    parser.add_argument("--batch-size", type=int, default=10000)
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        logger.error("FASTA file not found: %s", fasta_path)
        sys.exit(1)

    ancestral_seq = load_ancestral_fasta(fasta_path)

    logger.info("Annotating variants on %s...", args.chr)
    stats = annotate_chromosome(
        uri=args.uri,
        user=args.user,
        password=args.password,
        chrom=args.chr,
        ancestral_seq=ancestral_seq,
        batch_size=args.batch_size,
    )

    logger.info("Annotation complete:")
    logger.info("  Total variants: %d", stats["total"])
    logger.info("  Annotated: %d (%.1f%%)", stats["annotated"],
                100 * stats["annotated"] / max(stats["total"], 1))
    logger.info("  Skipped (unknown ancestor): %d", stats["skipped_unknown"])
    logger.info("  Skipped (no REF/ALT match): %d", stats["skipped_no_match"])


if __name__ == "__main__":
    main()
