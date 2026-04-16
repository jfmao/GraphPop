#!/usr/bin/env bash
# snapshot_database.sh — Create a versioned, publication-ready snapshot of a GraphPop Neo4j database.
#
# Creates:
#   <SNAPSHOT_DIR>/graphpop_<LABEL>_<DATE>.dump          — Neo4j dump file
#   <SNAPSHOT_DIR>/graphpop_<LABEL>_<DATE>.manifest.json — metadata for reproducibility
#
# Usage:
#   sudo bash scripts/snapshot_database.sh [--label human] [--db neo4j] [--out data/snapshots]
#
# Examples:
#   sudo bash scripts/snapshot_database.sh --label human_1000g_v1
#   sudo bash scripts/snapshot_database.sh --label rice_3k_v1 --db rice
#
# Notes:
#   - Neo4j is stopped before dumping and restarted afterwards.
#   - Run as root (or via sudo) — neo4j-admin needs write access to /var/lib/neo4j/.
#   - Requires graphevo conda env for the metadata collection step.
#   - Dump files are suitable for upload to Zenodo/Figshare for publication.

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
LABEL="human"
DATABASE="neo4j"
SNAPSHOT_DIR="data/snapshots"
GRAPHEVO_PYTHON="python"
NEO4J_URI="bolt://localhost:7687"
NEO4J_USER="${GRAPHPOP_USER:-neo4j}"
NEO4J_PASS="${GRAPHPOP_PASSWORD:-graphpop}"

# ── Parse args ────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --label) LABEL="$2"; shift 2 ;;
        --db)    DATABASE="$2"; shift 2 ;;
        --out)   SNAPSHOT_DIR="$2"; shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

DATE=$(date +%Y%m%d)
DUMP_NAME="graphpop_${LABEL}_${DATE}.dump"
MANIFEST_NAME="graphpop_${LABEL}_${DATE}.manifest.json"

# Resolve to absolute path (snapshot dir may be relative to project root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
SNAPSHOT_ABS="$PROJECT_ROOT/$SNAPSHOT_DIR"

mkdir -p "$SNAPSHOT_ABS"

echo "================================================================"
echo "  GraphPop Database Snapshot"
echo "  Label    : $LABEL"
echo "  Database : $DATABASE"
echo "  Output   : $SNAPSHOT_ABS/$DUMP_NAME"
echo "================================================================"

# ── Step 1: Collect metadata while Neo4j is running ──────────────────────────
echo ""
echo "[1/4] Collecting database metadata..."

MANIFEST_TMP="$SNAPSHOT_ABS/${MANIFEST_NAME}.tmp"

"$GRAPHEVO_PYTHON" - <<PYEOF
import json, subprocess, datetime, platform, socket
from neo4j import GraphDatabase

uri, user, pw, db = "$NEO4J_URI", "$NEO4J_USER", "$NEO4J_PASS", "$DATABASE"
driver = GraphDatabase.driver(uri, auth=(user, pw))

meta = {
    "graphpop_version": "1.0",
    "label": "$LABEL",
    "database": db,
    "snapshot_date": datetime.date.today().isoformat(),
    "neo4j_version": None,
    "host": socket.gethostname(),
    "node_counts": {},
    "relationship_counts": {},
    "embedded_statistics": {},
    "schema": {},
    "restore_instructions": (
        "1. Install Neo4j 5.x (Community or Enterprise)\\n"
        "2. Place the .dump file in an accessible directory\\n"
        "3. Stop Neo4j: sudo neo4j stop\\n"
        "4. Load: sudo neo4j-admin database load graphpop --from-path=<dir> --overwrite-destination=true\\n"
        "5. Start Neo4j: sudo neo4j start\\n"
        "6. Activate: cypher-shell -u neo4j -p <pass> -d system 'CREATE DATABASE graphpop IF NOT EXISTS;'\\n"
        "7. See README: https://github.com/your-org/GraphPop"
    ),
}

with driver.session(database=db) as s:
    # Neo4j version
    try:
        r = s.run("CALL dbms.components() YIELD name, versions RETURN name, versions").single()
        meta["neo4j_version"] = r["versions"][0] if r else None
    except Exception:
        pass

    # Node counts by label
    for label in ["Variant", "Sample", "Population", "Chromosome", "Gene", "Pathway", "GOTerm", "GenomicWindow"]:
        try:
            n = s.run(f"MATCH (n:{label}) RETURN count(n) AS c").single()["c"]
            meta["node_counts"][label] = n
            print(f"  {label}: {n:,}")
        except Exception:
            meta["node_counts"][label] = None

    # Relationship counts
    for rel in ["CARRIES", "NEXT", "IN_GENE", "IN_PATHWAY", "LD"]:
        try:
            n = s.run(f"MATCH ()-[r:{rel}]->() RETURN count(r) AS c").single()["c"]
            meta["relationship_counts"][rel] = n
            print(f"  [{rel}]: {n:,}")
        except Exception:
            meta["relationship_counts"][rel] = None

    # Check which popgen statistics are embedded
    stats = {}
    try:
        r = s.run(
            "MATCH (s:Sample) WHERE any(k IN keys(s) WHERE k STARTS WITH 'roh_') "
            "RETURN count(s) AS n LIMIT 1"
        ).single()
        stats["roh_on_samples"] = r["n"] if r else 0
    except Exception:
        stats["roh_on_samples"] = None

    try:
        r = s.run(
            "MATCH (w:GenomicWindow) WHERE w.h12 IS NOT NULL RETURN count(w) AS n"
        ).single()
        stats["garud_h_sweep_windows"] = r["n"] if r else 0
    except Exception:
        stats["garud_h_sweep_windows"] = None

    try:
        r = s.run(
            "MATCH (v:Variant) WHERE any(k IN keys(v) WHERE k STARTS WITH 'xpehh_') "
            "RETURN count(v) AS n LIMIT 1"
        ).single()
        stats["xpehh_peak_variants"] = r["n"] if r else 0
    except Exception:
        stats["xpehh_peak_variants"] = None

    try:
        r = s.run(
            "MATCH (w:GenomicWindow) WHERE w.pi IS NOT NULL RETURN count(w) AS n"
        ).single()
        stats["genomic_windows_with_diversity"] = r["n"] if r else 0
    except Exception:
        stats["genomic_windows_with_diversity"] = None

    meta["embedded_statistics"] = stats
    for k, v in stats.items():
        print(f"  stat/{k}: {v:,}" if isinstance(v, int) else f"  stat/{k}: {v}")

    # Indexes
    try:
        indexes = []
        for rec in s.run("SHOW INDEXES YIELD name, labelsOrTypes, properties, state"):
            indexes.append({
                "name": rec["name"],
                "labels": rec["labelsOrTypes"],
                "properties": rec["properties"],
                "state": rec["state"],
            })
        meta["schema"]["indexes"] = indexes
    except Exception:
        pass

driver.close()

with open("$MANIFEST_TMP", "w") as f:
    json.dump(meta, f, indent=2)
print("  Metadata collected.")
PYEOF

echo "[1/4] Done."

# ── Step 2: Stop Neo4j ────────────────────────────────────────────────────────
echo ""
echo "[2/4] Stopping Neo4j..."
neo4j stop 2>/dev/null || true
# Wait until fully stopped
for i in $(seq 1 30); do
    sleep 2
    if ! neo4j status 2>/dev/null | grep -q "running"; then
        echo "  Neo4j stopped."
        break
    fi
    echo "  Waiting for Neo4j to stop... ($i)"
done

# ── Step 3: Dump the database ─────────────────────────────────────────────────
echo ""
echo "[3/4] Dumping database '$DATABASE' to $SNAPSHOT_ABS/$DUMP_NAME ..."
neo4j-admin database dump \
    --overwrite-destination=true \
    --to-path="$SNAPSHOT_ABS" \
    "$DATABASE"

# neo4j-admin names the file <database>.dump; rename to our convention
if [[ -f "$SNAPSHOT_ABS/${DATABASE}.dump" && "${DATABASE}.dump" != "$DUMP_NAME" ]]; then
    mv "$SNAPSHOT_ABS/${DATABASE}.dump" "$SNAPSHOT_ABS/$DUMP_NAME"
fi

DUMP_SIZE=$(du -sh "$SNAPSHOT_ABS/$DUMP_NAME" | cut -f1)
echo "  Dump complete: $DUMP_SIZE"

# Finalise manifest with dump size
"$GRAPHEVO_PYTHON" - <<PYEOF2
import json, os
with open("$MANIFEST_TMP") as f:
    meta = json.load(f)
meta["dump_file"]  = "$DUMP_NAME"
meta["dump_size_bytes"] = os.path.getsize("$SNAPSHOT_ABS/$DUMP_NAME")
with open("$SNAPSHOT_ABS/$MANIFEST_NAME", "w") as f:
    json.dump(meta, f, indent=2)
os.remove("$MANIFEST_TMP")
print("  Manifest written.")
PYEOF2

echo "[3/4] Done."

# ── Step 4: Restart Neo4j ─────────────────────────────────────────────────────
echo ""
echo "[4/4] Restarting Neo4j..."
neo4j start
# Wait until available
for i in $(seq 1 30); do
    sleep 3
    if neo4j status 2>/dev/null | grep -q "running"; then
        echo "  Neo4j running."
        break
    fi
    echo "  Waiting for Neo4j to start... ($i)"
done

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "================================================================"
echo "  Snapshot complete."
echo ""
echo "  Dump     : $SNAPSHOT_ABS/$DUMP_NAME  ($DUMP_SIZE)"
echo "  Manifest : $SNAPSHOT_ABS/$MANIFEST_NAME"
echo ""
echo "  To share for publication:"
echo "    Upload both files to Zenodo/Figshare/OSF."
echo "    Reference the manifest for node counts and embedded statistics."
echo "================================================================"
