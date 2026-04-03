#!/usr/bin/env python3
"""Download and process rice pathway data from Plant Reactome.

Steps:
1. Read Plant Reactome gene_ids_by_pathway_and_species.tab (already downloaded)
2. Extract Oryza sativa UniProt IDs → pathway mappings
3. Map UniProt IDs → LOC_Os gene IDs via UniProt REST API
4. Save final LOC_Os → pathway mapping
5. Load INTO Neo4j as Pathway nodes + IN_PATHWAY edges

Usage:
    conda run -n graphevo python scripts/download_rice_pathways.py
"""

import json
import sys
import time
import urllib.request
import urllib.parse
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PLANT_REACTOME_DIR = ROOT / "data" / "raw" / "plant_reactome"
PATHWAY_FILE = PLANT_REACTOME_DIR / "gene_ids_by_pathway_and_species.tab"
OUTPUT_MAPPING = PLANT_REACTOME_DIR / "rice_loc_pathway_mapping.tsv"

NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")


def load_rice_pathways():
    """Load Oryza sativa entries from Plant Reactome."""
    print("==> Loading Plant Reactome rice pathways...")
    pathways = {}  # pathway_id -> pathway_name
    uniprot_to_pathways = defaultdict(set)  # uniprot_id -> set of pathway_ids

    with open(PATHWAY_FILE) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            pathway_id, pathway_name, species, gene_id = parts[0], parts[1], parts[2], parts[3]
            if "Oryza sativa" not in species:
                continue
            pathways[pathway_id] = pathway_name
            uniprot_to_pathways[gene_id].add(pathway_id)

    print(f"  {len(pathways)} pathways, {len(uniprot_to_pathways)} unique gene IDs")
    return pathways, uniprot_to_pathways


def map_uniprot_to_loc(uniprot_ids):
    """Map UniProt accessions to LOC_Os gene IDs via UniProt REST API."""
    print(f"==> Mapping {len(uniprot_ids)} UniProt IDs to LOC_Os via UniProt API...")
    sys.stdout.flush()

    uniprot_to_loc = {}
    batch_size = 200
    ids_list = sorted(uniprot_ids)

    for i in range(0, len(ids_list), batch_size):
        batch = ids_list[i:i + batch_size]
        query = " OR ".join(f"accession:{uid}" for uid in batch)
        url = (
            f"https://rest.uniprot.org/uniprotkb/search?"
            f"query=({urllib.parse.quote(query)})&"
            f"fields=accession,gene_oln,gene_names&"
            f"format=tsv&size=500"
        )

        for attempt in range(3):
            try:
                req = urllib.request.Request(url)
                req.add_header("User-Agent", "GraphPop/1.0 (jfmao@research)")
                resp = urllib.request.urlopen(req, timeout=60)
                data = resp.read().decode()
                break
            except Exception as e:
                if attempt < 2:
                    time.sleep(5 * (attempt + 1))
                    continue
                print(f"  ERROR batch {i}: {e}")
                data = ""

        for line in data.strip().split("\n")[1:]:  # skip header
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            accession = parts[0]
            gene_oln = parts[1]  # Ordered Locus Names (LOC_Os...)
            gene_names = parts[2]

            # Extract LOC_Os ID from ordered locus names
            loc_id = None
            for name in gene_oln.split():
                if name.startswith("LOC_Os"):
                    loc_id = name
                    break
            if not loc_id:
                for name in gene_names.split():
                    if name.startswith("LOC_Os"):
                        loc_id = name
                        break

            if loc_id and accession in uniprot_ids:
                uniprot_to_loc[accession] = loc_id

        if (i + batch_size) % 1000 == 0 or i + batch_size >= len(ids_list):
            print(f"  {min(i + batch_size, len(ids_list))}/{len(ids_list)} queried, "
                  f"{len(uniprot_to_loc)} mapped")
            sys.stdout.flush()

        time.sleep(0.5)  # rate limit

    print(f"  Mapped {len(uniprot_to_loc)}/{len(uniprot_ids)} UniProt → LOC_Os")
    return uniprot_to_loc


def load_into_neo4j(loc_to_pathways, pathways):
    """Load Pathway nodes and IN_PATHWAY edges into Neo4j."""
    from neo4j import GraphDatabase
import os

    print(f"==> Loading {len(pathways)} pathways and {sum(len(v) for v in loc_to_pathways.values())} edges into Neo4j...")
    sys.stdout.flush()

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    # Create Pathway nodes
    pathway_list = [{"pathwayId": pid, "name": pname, "source": "PlantReactome"}
                    for pid, pname in pathways.items()]

    with driver.session() as session:
        session.run("CREATE INDEX pathway_id_idx IF NOT EXISTS FOR (p:Pathway) ON (p.pathwayId)")

        for i in range(0, len(pathway_list), 500):
            batch = pathway_list[i:i + 500]
            session.run("""
                UNWIND $pathways AS p
                MERGE (pw:Pathway {pathwayId: p.pathwayId})
                SET pw.name = p.name, pw.source = p.source
            """, pathways=batch)
        print(f"  {len(pathway_list)} Pathway nodes created")

    # Create IN_PATHWAY edges (Gene → Pathway)
    edges = []
    for loc_id, pathway_ids in loc_to_pathways.items():
        for pid in pathway_ids:
            edges.append({"geneId": loc_id, "pathwayId": pid})

    loaded = 0
    with driver.session() as session:
        for i in range(0, len(edges), 2000):
            batch = edges[i:i + 2000]
            session.run("""
                UNWIND $edges AS e
                MATCH (g:Gene {geneId: e.geneId})
                MATCH (pw:Pathway {pathwayId: e.pathwayId})
                MERGE (g)-[:IN_PATHWAY]->(pw)
            """, edges=batch)
            loaded += len(batch)
            if loaded % 10000 == 0 or loaded >= len(edges):
                print(f"  {loaded}/{len(edges)} IN_PATHWAY edges loaded")
                sys.stdout.flush()

    driver.close()
    print(f"  Done: {len(pathway_list)} pathways, {len(edges)} edges")


def main():
    t0 = time.time()

    # Step 1: Load Plant Reactome rice data
    pathways, uniprot_to_pathways = load_rice_pathways()

    # Step 2: Map UniProt → LOC_Os
    all_uniprots = set(uniprot_to_pathways.keys())
    uniprot_to_loc = map_uniprot_to_loc(all_uniprots)

    # Step 3: Build LOC_Os → pathway mapping
    loc_to_pathways = defaultdict(set)
    unmapped = 0
    for uniprot_id, pathway_ids in uniprot_to_pathways.items():
        loc_id = uniprot_to_loc.get(uniprot_id)
        if loc_id:
            for pid in pathway_ids:
                loc_to_pathways[loc_id].add(pid)
        else:
            unmapped += 1

    print(f"\n==> Mapping summary:")
    print(f"  LOC_Os genes with pathways: {len(loc_to_pathways)}")
    print(f"  Total gene-pathway pairs: {sum(len(v) for v in loc_to_pathways.values())}")
    print(f"  Unmapped UniProt IDs: {unmapped}")

    # Step 4: Save mapping file
    with open(OUTPUT_MAPPING, "w") as f:
        f.write("loc_os_id\tpathway_id\tpathway_name\n")
        for loc_id in sorted(loc_to_pathways.keys()):
            for pid in sorted(loc_to_pathways[loc_id]):
                f.write(f"{loc_id}\t{pid}\t{pathways[pid]}\n")
    print(f"  Saved to {OUTPUT_MAPPING}")

    # Step 5: Load into Neo4j
    load_into_neo4j(loc_to_pathways, pathways)

    total_time = time.time() - t0
    print(f"\n=== Complete in {total_time:.1f}s ===")


if __name__ == "__main__":
    main()
