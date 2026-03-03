#!/usr/bin/env python3
"""Fix CARRIES edge phase values from the imported data.

The original import stored gt_phases (boolean: phased=True/False) as the phase
value, meaning all phased het genotypes got phase=1 regardless of which haplotype
carries the ALT allele. This script reads the VCF to determine the correct phase
(0 = ALT on first haplotype, 1 = ALT on second haplotype) and updates the Neo4j
CARRIES edges per-variant using Cypher.

Strategy: for each VCF variant, collect the sample IDs whose phase should be 0,
then update all their CARRIES edges in a single Cypher query per variant.

Usage:
    conda run -n graphevo python scripts/fix_carries_phase.py
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

import time
import numpy as np
from cyvcf2 import VCF
from neo4j import GraphDatabase

VCF_PATH = "data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz"
CHR = "chr22"
NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"

# How many variants to batch in one transaction
TX_BATCH = 100


def main():
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    vcf = VCF(VCF_PATH)
    sample_ids = list(vcf.samples)
    print(f"VCF samples: {len(sample_ids)}")
    print(f"Processing {CHR}...")

    t0 = time.time()
    n_variants = 0
    n_hets_total = 0
    n_to_fix = 0
    n_variants_with_fixes = 0

    # Batch: list of (variant_id, [sample_ids_to_fix])
    batch = []

    for v in vcf(CHR):
        alts = v.ALT
        if not alts or alts[0] in (".", "<*>", "*"):
            continue
        if len(alts) > 1:
            continue

        pos = v.POS
        ref = v.REF
        alt = alts[0]
        variant_id = f"{CHR}:{pos}:{ref}:{alt}"

        gt_types = v.gt_types
        het_indices = np.flatnonzero(gt_types == 1)

        if len(het_indices) == 0:
            n_variants += 1
            continue

        n_hets_total += len(het_indices)

        genotypes = v.genotypes
        fix_sids = []
        for i in het_indices:
            a0, a1, _ = genotypes[i]
            if a0 == 1:  # 1|0 → phase should be 0 (currently all are 1)
                fix_sids.append(sample_ids[i])

        if fix_sids:
            batch.append((variant_id, fix_sids))
            n_to_fix += len(fix_sids)
            n_variants_with_fixes += 1

        n_variants += 1

        if len(batch) >= TX_BATCH:
            _flush_batch(driver, batch)
            batch = []

        if n_variants % 100000 == 0:
            elapsed = time.time() - t0
            print(f"  {n_variants:>8,d} variants | {n_hets_total:>10,d} hets | "
                  f"{n_to_fix:>10,d} fixes ({n_variants_with_fixes:,d} variants) | "
                  f"{elapsed:.0f}s")

    if batch:
        _flush_batch(driver, batch)

    vcf.close()

    elapsed = time.time() - t0
    print(f"\nDone: {n_variants:,} variants, {n_hets_total:,} total hets, "
          f"{n_to_fix:,} phase corrections across {n_variants_with_fixes:,} variants "
          f"in {elapsed:.1f}s")

    # Quick verification on a small sample
    print("\nVerifying a few edges...")
    with driver.session() as session:
        result = session.run(
            "MATCH (s:Sample)-[c:CARRIES]->(v:Variant) "
            "WHERE c.gt = 1 "
            "WITH c.phase AS phase, count(*) AS cnt "
            "RETURN phase, cnt ORDER BY phase LIMIT 5"
        )
        for record in result:
            print(f"  phase={record['phase']}: {record['cnt']:,}")

    driver.close()


def _flush_batch(driver, batch):
    """Update CARRIES edges for a batch of variants."""
    with driver.session() as session:
        for variant_id, fix_sids in batch:
            session.execute_write(_update_variant, variant_id, fix_sids)


def _update_variant(tx, variant_id, fix_sids):
    """Set phase=0 on CARRIES edges for specific samples at one variant."""
    query = """
    MATCH (s:Sample)-[c:CARRIES]->(v:Variant {variantId: $vid})
    WHERE c.gt = 1 AND s.sampleId IN $sids
    SET c.phase = 0
    """
    tx.run(query, vid=variant_id, sids=fix_sids)


if __name__ == "__main__":
    main()
