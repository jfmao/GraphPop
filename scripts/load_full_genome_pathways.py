#!/usr/bin/env python3
"""Load full-genome Pathway and GO term annotations into Neo4j.

Steps:
  1. Parse cached Reactome (Ensembl2Reactome.txt) + BioMart GO TSV
  2. Collect all unique ENSG IDs (Homo sapiens only)
  3. Batch-query Ensembl REST /lookup/id to get ENSG → symbol mapping
  4. Translate edges to use gene symbols (matching Neo4j Gene.geneId)
  5. Clear old IN_PATHWAY / HAS_GO_TERM edges and load new ones

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/load_full_genome_pathways.py

Requires cached files in data/raw/1000g/csv_out/_cache/:
  - Ensembl2Reactome.txt
  - biomart_go.tsv
"""

import json
import time
import urllib.request
import urllib.error
from collections import defaultdict
from pathlib import Path

from neo4j import GraphDatabase

# ── Config ─────────────────────────────────────────────────────────────────────
CACHE_DIR   = Path("data/raw/1000g/csv_out/_cache")
REACTOME    = CACHE_DIR / "Ensembl2Reactome.txt"
BIOMART_GO  = CACHE_DIR / "biomart_go.tsv"

NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB   = "neo4j"
BATCH_SIZE = 5000

ENSEMBL_REST = "https://rest.ensembl.org"
LOOKUP_BATCH = 1000   # Ensembl REST max IDs per POST
RETRY_WAIT   = 5      # seconds between retries


# ── Step 1: Parse Reactome ─────────────────────────────────────────────────────

def parse_reactome(path: Path):
    """Returns (pathways dict, edges list) with ENSG IDs."""
    pathways = {}   # pathwayId → name
    edges = set()   # (ensg, pathwayId, evidence)
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            gene_id, reactome_id, _url, pathway_name, evidence, species = parts[:6]
            if species != "Homo sapiens":
                continue
            pathways[reactome_id] = pathway_name
            edges.add((gene_id, reactome_id, evidence))
    print(f"  Reactome: {len(pathways):,} pathways, {len(edges):,} gene-pathway edges", flush=True)
    return pathways, edges


# ── Step 2: Parse BioMart GO ────────────────────────────────────────────────────

def parse_biomart_go(path: Path):
    """Returns (goterms dict, edges set) with ENSG IDs."""
    import csv
    goterms = {}    # goId → (name, aspect)
    edges = set()   # (ensg, goId)
    aspect_map = {
        "biological_process": "biological_process",
        "molecular_function": "molecular_function",
        "cellular_component": "cellular_component",
    }
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # header
        for row in reader:
            if len(row) < 4:
                continue
            gene_id, go_id, go_name, go_namespace = row[:4]
            if not go_id:
                continue
            goterms[go_id] = (go_name, aspect_map.get(go_namespace, go_namespace))
            edges.add((gene_id, go_id))
    print(f"  BioMart GO: {len(goterms):,} GO terms, {len(edges):,} gene-GO edges", flush=True)
    return goterms, edges


# ── Step 3: Ensembl REST lookup ─────────────────────────────────────────────────

def ensembl_lookup_batch(ids: list[str]) -> dict[str, str]:
    """POST /lookup/id with up to LOOKUP_BATCH IDs. Returns {id: symbol}."""
    payload = json.dumps({"ids": ids}).encode()
    req = urllib.request.Request(
        f"{ENSEMBL_REST}/lookup/id",
        data=payload,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        method="POST",
    )
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.loads(resp.read())
                result = {}
                for eid, info in data.items():
                    if info and isinstance(info, dict):
                        sym = info.get("display_name") or info.get("id")
                        if sym:
                            result[eid] = sym
                return result
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = int(e.headers.get("Retry-After", RETRY_WAIT))
                print(f"    Rate limited, waiting {wait}s ...", flush=True)
                time.sleep(wait)
            else:
                print(f"    HTTP error {e.code}, attempt {attempt+1}/3", flush=True)
                time.sleep(RETRY_WAIT)
        except Exception as exc:
            print(f"    Error: {exc}, attempt {attempt+1}/3", flush=True)
            time.sleep(RETRY_WAIT)
    return {}


def build_ensg_to_symbol(all_ensg: set[str]) -> dict[str, str]:
    """Batch-query Ensembl REST to map ENSG → symbol."""
    ids = sorted(all_ensg)
    mapping = {}
    n_batches = (len(ids) + LOOKUP_BATCH - 1) // LOOKUP_BATCH
    print(f"  Querying Ensembl REST for {len(ids):,} ENSG IDs ({n_batches} batches) ...", flush=True)
    for i in range(0, len(ids), LOOKUP_BATCH):
        batch = ids[i:i + LOOKUP_BATCH]
        bn = i // LOOKUP_BATCH + 1
        result = ensembl_lookup_batch(batch)
        mapping.update(result)
        print(f"    batch {bn}/{n_batches}: {len(result)} mapped ({len(mapping):,} total)", flush=True)
        time.sleep(0.2)   # be polite
    print(f"  ENSG→symbol: {len(mapping):,} / {len(ids):,} mapped", flush=True)
    return mapping


# ── Step 4: Translate edges ─────────────────────────────────────────────────────

def translate_edges(edges, ensg_to_symbol):
    """Replace ENSG IDs with gene symbols. Skip if no mapping."""
    translated = []
    skipped = 0
    for edge in edges:
        ensg = edge[0]
        sym = ensg_to_symbol.get(ensg)
        if sym:
            translated.append((sym,) + edge[1:])
        else:
            skipped += 1
    return translated, skipped


# ── Step 5: Load into Neo4j ─────────────────────────────────────────────────────

def neo4j_flush(session, cypher, batch, label):
    total = 0
    n_batches = (len(batch) + BATCH_SIZE - 1) // BATCH_SIZE
    for i in range(0, len(batch), BATCH_SIZE):
        chunk = batch[i:i + BATCH_SIZE]
        bn = i // BATCH_SIZE + 1
        r = session.run(cypher, batch=chunk).consume()
        total += r.counters.relationships_created + r.counters.properties_set
        print(f"  {label}: batch {bn}/{n_batches} — {total:,} total", flush=True)
    return total


IN_PATHWAY_CYPHER = """
UNWIND $batch AS row
MATCH (g:Gene   {geneId:    row.gene})
MATCH (p:Pathway {pathwayId: row.pathway})
MERGE (g)-[:IN_PATHWAY {evidence: row.evidence}]->(p)
"""

HAS_GO_TERM_CYPHER = """
UNWIND $batch AS row
MATCH (g:Gene   {geneId: row.gene})
MATCH (t:GOTerm {goId:   row.go})
MERGE (g)-[:HAS_GO_TERM]->(t)
"""

PATHWAY_NODE_CYPHER = """
UNWIND $batch AS row
MERGE (p:Pathway {pathwayId: row.id})
SET p.name = row.name, p.source = 'Reactome'
"""

GOTERM_NODE_CYPHER = """
UNWIND $batch AS row
MERGE (t:GOTerm {goId: row.id})
SET t.name = row.name, t.aspect = row.aspect
"""


def load_pathway_nodes(session, pathways: dict):
    batch = [{"id": pid, "name": name} for pid, name in pathways.items()]
    neo4j_flush(session, PATHWAY_NODE_CYPHER, batch, "Pathway nodes")


def load_goterm_nodes(session, goterms: dict):
    batch = [{"id": gid, "name": name, "aspect": aspect}
             for gid, (name, aspect) in goterms.items()]
    neo4j_flush(session, GOTERM_NODE_CYPHER, batch, "GOTerm nodes")


def load_in_pathway(session, edges):
    batch = [{"gene": sym, "pathway": pid, "evidence": ev}
             for sym, pid, ev in edges]
    neo4j_flush(session, IN_PATHWAY_CYPHER, batch, "IN_PATHWAY edges")


def load_has_go_term(session, edges):
    batch = [{"gene": sym, "go": gid} for sym, gid in edges]
    neo4j_flush(session, HAS_GO_TERM_CYPHER, batch, "HAS_GO_TERM edges")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("=== Full-Genome Pathway + GO Term Loading ===", flush=True)

    # Parse source files
    print("\n[1/5] Parsing Reactome ...", flush=True)
    pathways, reactome_edges = parse_reactome(REACTOME)

    print("\n[2/5] Parsing BioMart GO ...", flush=True)
    goterms, go_edges = parse_biomart_go(BIOMART_GO)

    # Collect all ENSG IDs
    all_ensg = set()
    for ensg, *_ in reactome_edges:
        all_ensg.add(ensg)
    for ensg, _ in go_edges:
        all_ensg.add(ensg)
    print(f"\n  Total unique ENSG IDs: {len(all_ensg):,}", flush=True)

    # Map ENSG → symbol
    print("\n[3/5] Building ENSG → symbol mapping ...", flush=True)
    ensg_to_symbol = build_ensg_to_symbol(all_ensg)

    # Translate edges
    print("\n[4/5] Translating edges ...", flush=True)
    pathway_edges_sym, skip_p = translate_edges(reactome_edges, ensg_to_symbol)
    go_edges_sym, skip_g = translate_edges(go_edges, ensg_to_symbol)
    print(f"  IN_PATHWAY  : {len(pathway_edges_sym):,} edges (skipped {skip_p:,} unmapped ENSGs)", flush=True)
    print(f"  HAS_GO_TERM : {len(go_edges_sym):,} edges (skipped {skip_g:,} unmapped ENSGs)", flush=True)

    # Connect to Neo4j and load
    print("\n[5/5] Loading into Neo4j ...", flush=True)
    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()
    print("  Neo4j connected.", flush=True)

    with driver.session(database=NEO4J_DB) as session:
        # Ensure indexes
        session.run("CREATE INDEX pathway_id IF NOT EXISTS FOR (p:Pathway) ON (p.pathwayId)")
        session.run("CREATE INDEX goterm_id  IF NOT EXISTS FOR (t:GOTerm)  ON (t.goId)")
        print("  Indexes ensured.", flush=True)

        # Clear old edges (they were chr22-only)
        r = session.run("MATCH ()-[r:IN_PATHWAY]->()  DELETE r RETURN count(r) AS n").single()
        print(f"  Cleared {r['n']:,} old IN_PATHWAY edges.", flush=True)
        r = session.run("MATCH ()-[r:HAS_GO_TERM]->() DELETE r RETURN count(r) AS n").single()
        print(f"  Cleared {r['n']:,} old HAS_GO_TERM edges.", flush=True)

        # Load nodes
        print("\n  Loading Pathway nodes ...", flush=True)
        load_pathway_nodes(session, pathways)

        print("\n  Loading GOTerm nodes ...", flush=True)
        load_goterm_nodes(session, goterms)

        # Load edges
        print("\n  Loading IN_PATHWAY edges ...", flush=True)
        load_in_pathway(session, pathway_edges_sym)

        print("\n  Loading HAS_GO_TERM edges ...", flush=True)
        load_has_go_term(session, go_edges_sym)

    driver.close()
    print("\n=== Done ===", flush=True)

    # Summary query
    driver2 = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    with driver2.session(database=NEO4J_DB) as session:
        counts = {
            "Pathway nodes": session.run("MATCH (p:Pathway) RETURN count(p) AS n").single()["n"],
            "GOTerm nodes":  session.run("MATCH (t:GOTerm)  RETURN count(t) AS n").single()["n"],
            "IN_PATHWAY":    session.run("MATCH ()-[:IN_PATHWAY]->()  RETURN count(*) AS n").single()["n"],
            "HAS_GO_TERM":   session.run("MATCH ()-[:HAS_GO_TERM]->() RETURN count(*) AS n").single()["n"],
        }
    driver2.close()
    print("\nFinal counts:", flush=True)
    for k, v in counts.items():
        print(f"  {k:<20}: {v:>10,}", flush=True)


if __name__ == "__main__":
    main()
