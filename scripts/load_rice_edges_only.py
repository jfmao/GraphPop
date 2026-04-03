#!/usr/bin/env python3
"""Load rice HAS_CONSEQUENCE edges only (genes already loaded).

Streams the SnpEff VCF and loads edges in batches. Requires variant_id index.
"""

import os
import sys
import time

import cyvcf2
from neo4j import GraphDatabase

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")

SNPEFF_VCF = "data/raw/3kRG_data/NB_bialSNP_pseudo_canonical_ALL.vcf.gz"
CHR_MAP = {str(i): f"Chr{i}" for i in range(1, 13)}
FUNCTIONAL_IMPACTS = {"HIGH", "MODERATE", "LOW"}
BATCH_SIZE = 5000


def parse_ann(ann_str):
    for entry in ann_str.split(","):
        parts = entry.split("|")
        if len(parts) < 7:
            continue
        impact = parts[2]
        if impact not in FUNCTIONAL_IMPACTS:
            continue
        yield parts[3], impact, parts[1], parts[5], parts[6]


def main():
    t0 = time.time()
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    print("==> Streaming SnpEff VCF and loading HAS_CONSEQUENCE edges...")
    print(f"    Batch size: {BATCH_SIZE}")
    sys.stdout.flush()

    vcf = cyvcf2.VCF(SNPEFF_VCF)
    batch = []
    total_loaded = 0
    n_variants = 0

    for v in vcf:
        chrom = CHR_MAP.get(v.CHROM)
        if not chrom:
            continue

        ann = v.INFO.get("ANN")
        if not ann:
            continue

        variant_id = f"{chrom}:{v.POS}:{v.REF}:{v.ALT[0]}"
        seen_genes = set()

        for gene_id, impact, consequence, feature_type, feature in parse_ann(ann):
            if not gene_id or gene_id in seen_genes:
                continue
            seen_genes.add(gene_id)
            batch.append({
                "variantId": variant_id,
                "geneId": gene_id,
                "impact": impact,
                "consequence": consequence,
                "feature_type": feature_type,
                "feature": feature,
            })

        n_variants += 1

        if len(batch) >= BATCH_SIZE:
            with driver.session() as session:
                session.run("""
                    UNWIND $edges AS e
                    MATCH (v:Variant {variantId: e.variantId})
                    MATCH (g:Gene {geneId: e.geneId})
                    MERGE (v)-[r:HAS_CONSEQUENCE]->(g)
                    SET r.impact = e.impact, r.consequence = e.consequence,
                        r.feature_type = e.feature_type, r.feature = e.feature
                """, edges=batch)
            total_loaded += len(batch)
            batch = []

            if total_loaded % 50_000 == 0:
                elapsed = time.time() - t0
                rate = total_loaded / elapsed
                print(f"    {total_loaded:,} edges loaded ({rate:.0f}/s, {n_variants:,} variants scanned)")
                sys.stdout.flush()

    # Final batch
    if batch:
        with driver.session() as session:
            session.run("""
                UNWIND $edges AS e
                MATCH (v:Variant {variantId: e.variantId})
                MATCH (g:Gene {geneId: e.geneId})
                MERGE (v)-[r:HAS_CONSEQUENCE]->(g)
                SET r.impact = e.impact, r.consequence = e.consequence,
                    r.feature_type = e.feature_type, r.feature = e.feature
            """, edges=batch)
        total_loaded += len(batch)

    vcf.close()
    driver.close()

    total_time = time.time() - t0
    print(f"\n=== Complete: {total_loaded:,} edges loaded in {total_time:.1f}s ({total_time/60:.1f} min) ===")
    print(f"    Variants scanned: {n_variants:,}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
