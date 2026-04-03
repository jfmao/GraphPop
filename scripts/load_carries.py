#!/usr/bin/env python3
"""Stream CARRIES edges from CSV into Neo4j. Designed to run overnight.

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/load_carries.py
"""
import csv
import time
from pathlib import Path
from neo4j import GraphDatabase
import os

NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"
CSV_PATH   = "data/raw/1000g/csv_out/carries_edges.csv"
BATCH_SIZE = 20_000

CYPHER = """
UNWIND $batch AS row
MATCH (s:Sample  {sampleId:  row.sample})
MATCH (v:Variant {variantId: row.variant})
MERGE (s)-[:CARRIES {gt: row.gt, phase: row.phase}]->(v)
"""

def main():
    path = Path(CSV_PATH)
    total_lines = sum(1 for _ in open(path)) - 1  # subtract header
    print(f"CARRIES CSV: {path} ({total_lines:,} rows)", flush=True)

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()
    print("Neo4j connected.", flush=True)

    # Create index to speed up MATCH on Sample if not exists
    with driver.session(database=NEO4J_DB) as s:
        s.run("CREATE INDEX carries_sample IF NOT EXISTS FOR (s:Sample) ON (s.sampleId)")
        print("Sample index ensured.", flush=True)

    batch = []
    total_created = 0
    t0 = time.time()
    n_batches = 0

    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            batch.append({
                "sample":  row[":START_ID(Sample)"],
                "variant": row[":END_ID(Variant)"],
                "gt":      int(row["gt:int"]),
                "phase":   int(row["phase:int"]),
            })
            if len(batch) >= BATCH_SIZE:
                with driver.session(database=NEO4J_DB) as s:
                    r = s.run(CYPHER, batch=batch).consume()
                    total_created += r.counters.relationships_created
                n_batches += 1
                elapsed = time.time() - t0
                rate = (n_batches * BATCH_SIZE) / elapsed
                pct  = (n_batches * BATCH_SIZE) / total_lines * 100
                eta  = (total_lines - n_batches * BATCH_SIZE) / rate / 3600
                print(
                    f"  batch {n_batches}: {total_created:,} edges created | "
                    f"{pct:.1f}% | {rate:.0f} rows/s | ETA {eta:.1f}h",
                    flush=True
                )
                batch = []

    # Final partial batch
    if batch:
        with driver.session(database=NEO4J_DB) as s:
            r = s.run(CYPHER, batch=batch).consume()
            total_created += r.counters.relationships_created

    driver.close()
    elapsed = time.time() - t0
    print(f"\nDone. {total_created:,} CARRIES edges created in {elapsed/3600:.2f}h", flush=True)

if __name__ == "__main__":
    main()
