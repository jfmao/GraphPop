#!/usr/bin/env python3
"""Load gene_nodes.csv and has_consequence_edges.csv into Neo4j via Python driver.

Bypasses LOAD CSV (which requires files in Neo4j import dir) by reading CSVs
directly and using UNWIND batching through the neo4j Python driver.

Usage:
    python scripts/load_annotations_driver.py [CSV_DIR]
    python scripts/load_annotations_driver.py data/raw/1000g/full_genome_annotations/
"""

import csv
import os
import sys
import time

from neo4j import GraphDatabase

NEO4J_URI = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.environ.get("NEO4J_USER", "neo4j")
NEO4J_PASS = os.environ.get("NEO4J_PASS", "graphpop")

BATCH_SIZE = 10000
EDGE_BATCH_SIZE = 5000  # Smaller for edges (MATCH lookups)

# Only load functional annotations (HIGH, MODERATE, LOW impact)
# MODIFIER = 89M edges (intergenic, upstream, intronic) — skip unless --all-impacts
FUNCTIONAL_IMPACTS = {"HIGH", "MODERATE", "LOW"}


def load_genes(driver, csv_path):
    """Load Gene nodes from gene_nodes.csv."""
    print(f"==> Loading Gene nodes from {csv_path}")
    genes = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            genes.append({
                "geneId": row["geneId:ID(Gene)"],
                "symbol": row["symbol"],
                "biotype": row["biotype"],
            })

    print(f"    Read {len(genes):,} genes")
    t0 = time.time()
    loaded = 0

    with driver.session() as session:
        for i in range(0, len(genes), BATCH_SIZE):
            batch = genes[i:i + BATCH_SIZE]
            session.run(
                """
                UNWIND $genes AS g
                MERGE (n:Gene {geneId: g.geneId})
                SET n.symbol = g.symbol, n.biotype = g.biotype
                """,
                genes=batch,
            )
            loaded += len(batch)
            if loaded % 50000 == 0 or loaded == len(genes):
                print(f"    {loaded:,} / {len(genes):,} genes loaded")

    elapsed = time.time() - t0
    print(f"    Done in {elapsed:.1f}s ({len(genes):,} genes)")
    return len(genes)


def load_edges(driver, csv_path, all_impacts=False):
    """Load HAS_CONSEQUENCE edges from has_consequence_edges.csv."""
    impact_filter = "all" if all_impacts else "HIGH/MODERATE/LOW only"
    print(f"==> Loading HAS_CONSEQUENCE edges ({impact_filter})")
    print(f"    Source: {csv_path}")
    t0 = time.time()
    loaded = 0
    skipped_match = 0
    skipped_impact = 0
    batch = []
    read_count = 0

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            read_count += 1
            impact = row["impact"]

            # Filter by impact unless --all-impacts
            if not all_impacts and impact not in FUNCTIONAL_IMPACTS:
                skipped_impact += 1
                continue

            edge = {
                "variantId": row[":START_ID(Variant)"],
                "geneId": row[":END_ID(Gene)"],
                "consequence": row["consequence"],
                "impact": impact,
                "feature": row["feature"],
                "feature_type": row["feature_type"],
            }
            batch.append(edge)

            if len(batch) >= EDGE_BATCH_SIZE:
                result = _flush_edges(driver, batch)
                loaded += result
                skipped_match += len(batch) - result
                batch = []

                if loaded % 100000 < EDGE_BATCH_SIZE:
                    elapsed = time.time() - t0
                    rate = loaded / elapsed if elapsed > 0 else 0
                    print(
                        f"    {loaded:,} loaded, {skipped_match:,} no-match, "
                        f"{skipped_impact:,} MODIFIER-skipped "
                        f"({read_count:,} rows read) [{rate:.0f} edges/s]"
                    )

    # Flush remaining
    if batch:
        result = _flush_edges(driver, batch)
        loaded += result
        skipped_match += len(batch) - result

    elapsed = time.time() - t0
    print(f"    Done in {elapsed:.1f}s")
    print(f"    {loaded:,} edges loaded, {skipped_match:,} no-match, "
          f"{skipped_impact:,} MODIFIER-skipped ({read_count:,} total rows)")
    return loaded


def _flush_edges(driver, batch):
    """Write a batch of edges. Returns number of edges created."""
    with driver.session() as session:
        result = session.run(
            """
            UNWIND $edges AS e
            MATCH (v:Variant {variantId: e.variantId})
            MATCH (g:Gene {geneId: e.geneId})
            CREATE (v)-[:HAS_CONSEQUENCE {
                consequence: e.consequence,
                impact: e.impact,
                feature: e.feature,
                feature_type: e.feature_type
            }]->(g)
            RETURN count(*) AS created
            """,
            edges=batch,
        )
        record = result.single()
        return record["created"] if record else 0


def create_indexes(driver):
    """Create indexes for Gene nodes and HAS_CONSEQUENCE edges."""
    print("==> Creating indexes...")
    with driver.session() as session:
        session.run("CREATE CONSTRAINT gene_id IF NOT EXISTS FOR (g:Gene) REQUIRE g.geneId IS UNIQUE")
        session.run("CREATE INDEX gene_symbol IF NOT EXISTS FOR (g:Gene) ON (g.symbol)")
        session.run("CREATE INDEX consequence_impact IF NOT EXISTS FOR ()-[r:HAS_CONSEQUENCE]-() ON (r.impact)")
        session.run("CREATE INDEX consequence_type IF NOT EXISTS FOR ()-[r:HAS_CONSEQUENCE]-() ON (r.consequence)")
    print("    Done.")


def verify(driver):
    """Run verification queries."""
    print("\n==> Verification:")
    with driver.session() as session:
        # Gene count
        result = session.run("MATCH (g:Gene) RETURN count(g) AS cnt")
        print(f"    Genes: {result.single()['cnt']:,}")

        # Edge count
        result = session.run("MATCH ()-[r:HAS_CONSEQUENCE]->() RETURN count(r) AS cnt")
        print(f"    HAS_CONSEQUENCE edges: {result.single()['cnt']:,}")

        # Sample edges by impact
        result = session.run("""
            MATCH ()-[r:HAS_CONSEQUENCE]->()
            RETURN r.impact AS impact, count(*) AS cnt
            ORDER BY cnt DESC
        """)
        print("    By impact:")
        for record in result:
            print(f"      {record['impact']}: {record['cnt']:,}")

        # Sample HIGH impact
        result = session.run("""
            MATCH (v:Variant)-[r:HAS_CONSEQUENCE {impact: 'HIGH'}]->(g:Gene)
            RETURN v.variantId AS variant, g.symbol AS gene, r.consequence AS consequence
            LIMIT 5
        """)
        print("    HIGH impact examples:")
        for record in result:
            print(f"      {record['variant']} → {record['gene']} ({record['consequence']})")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Load annotations into Neo4j")
    parser.add_argument("csv_dir", nargs="?", default="data/raw/1000g/full_genome_annotations")
    parser.add_argument("--all-impacts", action="store_true",
                        help="Load ALL impacts including MODIFIER (89M+ edges, slow)")
    args = parser.parse_args()

    csv_dir = os.path.abspath(args.csv_dir)

    gene_path = os.path.join(csv_dir, "gene_nodes.csv")
    edge_path = os.path.join(csv_dir, "has_consequence_edges.csv")

    for p in [gene_path, edge_path]:
        if not os.path.exists(p):
            print(f"ERROR: {p} not found")
            sys.exit(1)

    print(f"CSV directory: {csv_dir}")
    print(f"Neo4j:         {NEO4J_URI}")
    print()

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    # Check connectivity
    driver.verify_connectivity()
    print("Connected to Neo4j.\n")

    # Check if annotations already loaded
    with driver.session() as session:
        result = session.run("MATCH (g:Gene) RETURN count(g) AS cnt")
        existing = result.single()["cnt"]
        if existing > 0:
            print(f"WARNING: {existing:,} Gene nodes already exist.")
            result = session.run("MATCH ()-[r:HAS_CONSEQUENCE]->() RETURN count(r) AS cnt")
            existing_edges = result.single()["cnt"]
            print(f"WARNING: {existing_edges:,} HAS_CONSEQUENCE edges already exist.")
            print("Skipping load — annotations already present.")
            print("To reload, first run: MATCH (g:Gene) DETACH DELETE g")
            verify(driver)
            driver.close()
            return

    create_indexes(driver)
    load_genes(driver, gene_path)
    load_edges(driver, edge_path, all_impacts=args.all_impacts)
    verify(driver)

    driver.close()
    print("\n==> Annotation loading complete.")


if __name__ == "__main__":
    main()
