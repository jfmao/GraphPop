#!/usr/bin/env python3
"""Write back computed population genetics results to Neo4j nodes.

Operations:
  roh      - ROH summaries → Sample node properties (roh_{chr}_*)
  garud_h  - Garud's H sweep windows → GenomicWindow node properties
  xpehh    - XP-EHH top peak scores → Variant node properties
  all      - run all three (default)

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/import_popgen_to_neo4j.py --op all --dataset human

Note: Rice write-backs are skipped (--dataset human only). The Neo4j database
currently contains only human (1000 Genomes) data. Rice has no Sample, Variant,
or GenomicWindow nodes in Neo4j, so --dataset rice or --dataset both will find
no matching nodes and set 0 properties.
"""
import argparse
import json

from neo4j import GraphDatabase
import os

NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = os.environ.get("GRAPHPOP_USER", "neo4j")
NEO4J_PASS = os.environ.get("GRAPHPOP_PASSWORD", "graphpop")
NEO4J_DB   = "neo4j"
BATCH_SIZE = 5000

HUMAN_RESULTS = "human_interpretation_results.json"
RICE_RESULTS  = "results/rice/rice_interpretation_results.json"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def flush(session, cypher, batch, label):
    """Run a batched UNWIND cypher and print a summary."""
    if not batch:
        print(f"  {label}: nothing to write", flush=True)
        return 0
    total_props = 0
    n_batches = (len(batch) + BATCH_SIZE - 1) // BATCH_SIZE
    for i in range(0, len(batch), BATCH_SIZE):
        chunk = batch[i : i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        print(f"  {label}: batch {batch_num}/{n_batches} ({len(chunk)} items)...", flush=True)
        r = session.run(cypher, batch=chunk).consume()
        total_props += r.counters.properties_set
    print(f"  {label}: DONE — {total_props} properties set ({len(batch)} items)", flush=True)
    return total_props


# ---------------------------------------------------------------------------
# Priority 1: ROH → Sample nodes
# ---------------------------------------------------------------------------

ROH_CYPHER = """
UNWIND $batch AS item
MATCH (s:Sample {sampleId: item.sampleId})
SET s += item.props
"""


def import_roh(session, data):
    roh_hmm = data.get("roh_hmm", {})
    if not roh_hmm:
        print("  ROH: no roh_hmm key found — skipping")
        return

    batch = []
    skipped = 0
    for key, entry in roh_hmm.items():
        per_sample = entry.get("per_sample")
        if not per_sample:
            skipped += 1
            continue
        # key format: "POP|chrN"  e.g. "YRI|chr22" or "XI-adm|Chr1"
        parts = key.split("|", 1)
        chrom = parts[1] if len(parts) == 2 else key
        for s in per_sample:
            sid = s.get("sampleId")
            if not sid:
                continue
            batch.append({
                "sampleId": sid,
                "props": {
                    f"roh_{chrom}_n_roh":    int(s["n_roh"]),
                    f"roh_{chrom}_froh":     float(s["froh"]),
                    f"roh_{chrom}_total_kb": float(s["total_length"]) / 1000.0,
                },
            })

    if skipped:
        print(f"  ROH: {skipped} pop|chr entries had no per_sample data (skipped)", flush=True)
    flush(session, ROH_CYPHER, batch, "ROH → Sample")


# ---------------------------------------------------------------------------
# Priority 2: Garud's H → GenomicWindow nodes
# ---------------------------------------------------------------------------

GH_CYPHER = """
UNWIND $batch AS item
MATCH (w:GenomicWindow {chr: item.chr, start: item.start,
                        end: item.end, population: item.population})
SET w += item.props
"""


def import_garud_h(session, data):
    ann = data.get("annotate", {})
    sweep_annotations = ann.get("sweep_annotations", {})
    if not sweep_annotations:
        print("  Garud's H: no sweep_annotations found — skipping")
        return

    batch = []
    for pop, pop_data in sweep_annotations.items():
        for sw in pop_data.get("sweeps", []):
            batch.append({
                "chr":        sw["chr"],
                "start":      int(sw["start"]),
                "end":        int(sw["end"]),
                "population": pop,
                "props": {
                    "h1":            float(sw["h1"]),
                    "h12":           float(sw["h12"]),
                    "h2_h1":         float(sw["h2_h1"]),
                    "hap_diversity": float(sw["hap_diversity"]) if sw.get("hap_diversity") is not None else None,
                    "sweep_type":    sw.get("sweep_type", ""),
                },
            })

    flush(session, GH_CYPHER, batch, "Garud's H → GenomicWindow")


# ---------------------------------------------------------------------------
# Priority 3: XP-EHH → Variant nodes
# ---------------------------------------------------------------------------

XP_CYPHER = """
UNWIND $batch AS item
MATCH (v:Variant {variantId: item.variantId})
SET v += item.props
"""


def import_xpehh(session, data):
    ann = data.get("annotate", {})
    xpehh_annotations = ann.get("xpehh_annotations", {})
    if not xpehh_annotations:
        print("  XP-EHH: no xpehh_annotations found — skipping")
        return

    batch = []
    for pair, pair_data in xpehh_annotations.items():
        pair_tag = pair.replace("-", "_").replace(" ", "_")
        prop_key = f"xpehh_{pair_tag}"
        for peak in pair_data.get("top_peaks", []):
            vid   = peak.get("peak_variantId")
            score = peak.get("peak_xpehh")
            if vid and score is not None:
                batch.append({
                    "variantId": vid,
                    "props": {prop_key: float(score)},
                })

    flush(session, XP_CYPHER, batch, "XP-EHH → Variant")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--op",      choices=["roh", "garud_h", "xpehh", "all"],
                        default="all", help="which write-back to run")
    parser.add_argument("--dataset", choices=["human", "rice", "both"],
                        default="both", help="which dataset to process")
    args = parser.parse_args()

    datasets = {"human": HUMAN_RESULTS, "rice": RICE_RESULTS}
    targets  = ["human", "rice"] if args.dataset == "both" else [args.dataset]
    ops      = ["roh", "garud_h", "xpehh"] if args.op == "all" else [args.op]

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))
    driver.verify_connectivity()
    print(f"Neo4j connected: {NEO4J_URI}", flush=True)

    try:
        for ds in targets:
            path = datasets[ds]
            print(f"\n=== {ds.upper()} ({path}) ===", flush=True)
            print(f"  Loading JSON...", flush=True)
            with open(path) as f:
                data = json.load(f)
            print(f"  JSON loaded.", flush=True)
            with driver.session(database=NEO4J_DB) as session:
                if "roh"     in ops: import_roh(session, data)
                if "garud_h" in ops: import_garud_h(session, data)
                if "xpehh"   in ops: import_xpehh(session, data)
    finally:
        driver.close()

    print("\nDone.")


if __name__ == "__main__":
    main()
