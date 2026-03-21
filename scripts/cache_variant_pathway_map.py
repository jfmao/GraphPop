#!/usr/bin/env python3
"""One-time cache of variant → pathway mapping from Neo4j.

Streams all HIGH/MODERATE impact variant-pathway edges from the graph and
saves a compact JSON cache so that Analysis 5 (pathway trajectory embedding)
can run without Neo4j.

Cache structure (data/results/variant_pathway_cache.json):
  {
    "pathway_meta": {
      "R-HSA-xxxx": {"name": "...", "mean_fst": 0.123},
      ...
    },
    "variant_pathways": {
      "chr1": {"12345": ["R-HSA-xxx", ...], ...},
      "chr2": {...},
      ...
    },
    "n_variants":   int,   # unique (chr, pos) entries
    "n_rows":       int,   # total variant-pathway rows streamed
    "n_pathways":   int,   # unique pathway IDs
    "rare_af":      float  # AF threshold used downstream
  }

Runtime: ~60–90 seconds (same query as Inv.18).
Output:  ~50–100 MB JSON (can be loaded in <5s for Analysis 5).

Usage:
  conda run -n graphevo python -u scripts/cache_variant_pathway_map.py
"""

import json
import time
from collections import defaultdict
from pathlib import Path

RESULTS  = Path("data/results")
OUT_PATH = RESULTS / "variant_pathway_cache.json"


def main():
    t0 = time.time()
    print("=== Cache: Variant → Pathway Mapping ===\n")

    # Load pathway FST metadata from Inv.01 results
    print("Loading pathway FST metadata...", flush=True)
    pws = json.load(open(RESULTS / "pathway_fst.json"))
    pathway_meta = {
        p["pathway_id"]: {
            "name":     p["pathway_name"],
            "mean_fst": round(p["mean_fst"], 6) if p.get("mean_fst") is not None else None,
        }
        for p in pws
    }
    print(f"  {len(pathway_meta):,} pathways with FST data", flush=True)

    # Stream Neo4j query
    print("\nStreaming Neo4j (all chromosomes, HIGH/MODERATE impact)...", flush=True)
    from neo4j import GraphDatabase
    driver = GraphDatabase.driver("bolt://localhost:7687", auth=("neo4j", "graphpop"))

    # variant_pathways[chr][pos_str] = [pathway_id, ...]
    variant_pathways = defaultdict(lambda: defaultdict(list))
    n_rows = 0
    pathway_ids_seen = set()

    try:
        with driver.session(database="neo4j") as s:
            result = s.run(
                "MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway) "
                "WHERE c.impact IN ['HIGH', 'MODERATE'] "
                "RETURN v.chr AS chr, v.pos AS pos, p.pathwayId AS pathway_id"
            )
            for rec in result:
                pid = rec["pathway_id"]
                if pid:
                    chrom = rec["chr"]
                    pos   = str(int(rec["pos"]))
                    variant_pathways[chrom][pos].append(pid)
                    pathway_ids_seen.add(pid)
                n_rows += 1
                if n_rows % 250_000 == 0:
                    elapsed = time.time() - t0
                    print(f"  {n_rows:,} rows  |  "
                          f"{sum(len(v) for v in variant_pathways.values()):,} positions  |  "
                          f"{elapsed:.0f}s", flush=True)
    finally:
        driver.close()

    n_variants = sum(len(v) for v in variant_pathways.values())
    print(f"\nDone: {n_rows:,} rows streamed  →  {n_variants:,} unique positions  "
          f"|  {len(pathway_ids_seen):,} pathways  ({time.time()-t0:.0f}s)", flush=True)

    # Deduplicate pathway lists (same pathway can appear multiple times via multiple genes)
    print("Deduplicating pathway lists...", flush=True)
    for chrom in variant_pathways:
        for pos in variant_pathways[chrom]:
            variant_pathways[chrom][pos] = list(set(variant_pathways[chrom][pos]))

    # Build output
    out = {
        "pathway_meta":     pathway_meta,
        "variant_pathways": {c: dict(v) for c, v in variant_pathways.items()},
        "n_variants":       n_variants,
        "n_rows":           n_rows,
        "n_pathways":       len(pathway_ids_seen),
        "rare_af":          0.01,
        "query": (
            "MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway) "
            "WHERE c.impact IN ['HIGH', 'MODERATE'] "
            "RETURN v.chr, v.pos, p.pathwayId"
        ),
    }

    # Write JSON
    print(f"Writing {OUT_PATH} ...", flush=True)
    with open(OUT_PATH, "w") as f:
        json.dump(out, f, separators=(",", ":"))   # compact, no spaces

    size_mb = OUT_PATH.stat().st_size / 1e6
    print(f"  Done: {size_mb:.1f} MB  ({time.time()-t0:.0f}s total)")
    print(f"\nCache ready at: {OUT_PATH}")
    print("Analysis 5 can now run without Neo4j.")


if __name__ == "__main__":
    main()
