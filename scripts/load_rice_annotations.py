#!/usr/bin/env python3
"""Load rice SnpEff annotations into Neo4j.

Reads the SnpEff-annotated VCF (NB_bialSNP_pseudo_canonical_ALL.vcf.gz),
extracts Gene nodes and HAS_CONSEQUENCE edges, and loads them into Neo4j.

Only loads functional annotations (HIGH, MODERATE, LOW impact).
MODIFIER annotations (intergenic, upstream, intronic) are skipped.

Usage:
    conda run -n graphevo python scripts/load_rice_annotations.py
"""

import os
import sys
import time
from collections import defaultdict

import cyvcf2
from neo4j import GraphDatabase

NEO4J_URI = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.environ.get("NEO4J_USER", "neo4j")
NEO4J_PASS = os.environ.get("NEO4J_PASS", "graphpop")

SNPEFF_VCF = "data/raw/3kRG_data/NB_bialSNP_pseudo_canonical_ALL.vcf.gz"

# Rice chromosomes: SnpEff VCF uses "1","2",...,"12" while our DB uses "Chr1","Chr2",...
CHR_MAP = {str(i): f"Chr{i}" for i in range(1, 13)}

FUNCTIONAL_IMPACTS = {"HIGH", "MODERATE", "LOW"}
BATCH_SIZE = 5000
GENE_BATCH = 10000


def parse_ann(ann_str):
    """Parse SnpEff ANN field, yield (gene_id, impact, consequence, feature_type, feature)."""
    for entry in ann_str.split(","):
        parts = entry.split("|")
        if len(parts) < 5:
            continue
        impact = parts[2]
        if impact not in FUNCTIONAL_IMPACTS:
            continue
        gene_id = parts[3]
        consequence = parts[1]
        feature_type = parts[5] if len(parts) > 5 else ""
        feature = parts[6] if len(parts) > 6 else ""
        yield gene_id, impact, consequence, feature_type, feature


def main():
    t0 = time.time()
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    # Pass 1: Collect all genes and edges
    print("==> Pass 1: Scanning SnpEff VCF for genes and consequences...")
    vcf = cyvcf2.VCF(SNPEFF_VCF)

    genes = {}  # gene_id -> set of biotypes
    edges = []  # (variant_id, gene_id, impact, consequence, feature_type, feature)
    n_variants = 0
    n_ann = 0

    for v in vcf:
        chrom_raw = v.CHROM
        chrom = CHR_MAP.get(chrom_raw)
        if not chrom:
            continue

        variant_id = f"{chrom}:{v.POS}:{v.REF}:{v.ALT[0]}"

        ann = v.INFO.get("ANN")
        if not ann:
            continue

        seen_genes = set()  # deduplicate per variant
        for gene_id, impact, consequence, feature_type, feature in parse_ann(ann):
            if not gene_id or gene_id in seen_genes:
                continue
            seen_genes.add(gene_id)

            if gene_id not in genes:
                genes[gene_id] = "protein_coding"  # default; SnpEff doesn't always give biotype

            edges.append((variant_id, gene_id, impact, consequence, feature_type, feature))
            n_ann += 1

        n_variants += 1
        if n_variants % 1_000_000 == 0:
            print(f"    {n_variants:,} variants scanned, {n_ann:,} edges, {len(genes):,} genes")

    vcf.close()
    scan_time = time.time() - t0
    print(f"    Scan complete: {n_variants:,} variants, {n_ann:,} edges, {len(genes):,} genes in {scan_time:.1f}s")

    # Pass 2: Load Gene nodes
    print(f"\n==> Loading {len(genes):,} Gene nodes...")
    t1 = time.time()
    gene_list = [{"geneId": gid, "symbol": gid, "biotype": bt} for gid, bt in genes.items()]

    with driver.session() as session:
        # Create index first
        session.run("CREATE INDEX gene_id_idx IF NOT EXISTS FOR (g:Gene) ON (g.geneId)")
        session.run("CREATE INDEX gene_symbol_idx IF NOT EXISTS FOR (g:Gene) ON (g.symbol)")

        for i in range(0, len(gene_list), GENE_BATCH):
            batch = gene_list[i:i + GENE_BATCH]
            session.run("""
                UNWIND $genes AS g
                MERGE (n:Gene {geneId: g.geneId})
                SET n.symbol = g.symbol, n.biotype = g.biotype
            """, genes=batch)
            if (i + GENE_BATCH) % 50000 == 0 or i + GENE_BATCH >= len(gene_list):
                print(f"    {min(i + GENE_BATCH, len(gene_list)):,} / {len(gene_list):,} genes loaded")

    gene_time = time.time() - t1
    print(f"    Gene nodes loaded in {gene_time:.1f}s")

    # Pass 3: Load HAS_CONSEQUENCE edges
    print(f"\n==> Loading {len(edges):,} HAS_CONSEQUENCE edges...")
    t2 = time.time()
    loaded = 0

    with driver.session() as session:
        for i in range(0, len(edges), BATCH_SIZE):
            batch = [
                {
                    "variantId": e[0],
                    "geneId": e[1],
                    "impact": e[2],
                    "consequence": e[3],
                    "feature_type": e[4],
                    "feature": e[5],
                }
                for e in edges[i:i + BATCH_SIZE]
            ]
            session.run("""
                UNWIND $edges AS e
                MATCH (v:Variant {variantId: e.variantId})
                MATCH (g:Gene {geneId: e.geneId})
                MERGE (v)-[r:HAS_CONSEQUENCE]->(g)
                SET r.impact = e.impact,
                    r.consequence = e.consequence,
                    r.feature_type = e.feature_type,
                    r.feature = e.feature
            """, edges=batch)
            loaded += len(batch)
            if loaded % 100_000 == 0 or loaded >= len(edges):
                elapsed = time.time() - t2
                rate = loaded / elapsed if elapsed > 0 else 0
                print(f"    {loaded:,} / {len(edges):,} edges loaded ({rate:.0f} edges/s)")

    edge_time = time.time() - t2
    total_time = time.time() - t0

    print(f"\n=== Rice Annotation Loading Complete ===")
    print(f"  Genes:  {len(genes):,}")
    print(f"  Edges:  {len(edges):,}")
    print(f"  VCF scan:   {scan_time:.1f}s")
    print(f"  Gene load:  {gene_time:.1f}s")
    print(f"  Edge load:  {edge_time:.1f}s")
    print(f"  Total:      {total_time:.1f}s")

    driver.close()


if __name__ == "__main__":
    main()
