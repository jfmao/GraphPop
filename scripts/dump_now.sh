#!/usr/bin/env bash
# One-shot dump of the GraphPop human_1000g_v1 database.
# Run as: sudo bash scripts/dump_now.sh
set -euo pipefail

LABEL="human_1000g_v1"
DATABASE="neo4j"
SNAP_DIR="/mnt/data/GraphPop/data/snapshots"
DUMP_FILE="$SNAP_DIR/graphpop_${LABEL}_$(date +%Y%m%d).dump"

mkdir -p "$SNAP_DIR"

echo "=== GraphPop Database Dump ==="
echo "  Target : $DUMP_FILE"

echo "[1/3] Stopping Neo4j ..."
neo4j stop
for i in $(seq 1 30); do
    sleep 2
    neo4j status 2>/dev/null | grep -q "running" || { echo "  Stopped."; break; }
    echo "  Waiting ($i/30) ..."
done

echo "[2/3] Dumping '$DATABASE' ..."
neo4j-admin database dump \
    --overwrite-destination=true \
    --to-path="$SNAP_DIR" \
    "$DATABASE"

# neo4j-admin names it <database>.dump — rename to our convention
[[ -f "$SNAP_DIR/${DATABASE}.dump" ]] && mv "$SNAP_DIR/${DATABASE}.dump" "$DUMP_FILE"
echo "  Dump size: $(du -sh "$DUMP_FILE" | cut -f1)"

# Patch dump filename into manifest
/home/jfmao/miniconda3/envs/graphevo/bin/python - <<PYEOF
import json, os
mf = "$SNAP_DIR/graphpop_${LABEL}_manifest.json"
with open(mf) as f: m = json.load(f)
m["dump_file"]       = os.path.basename("$DUMP_FILE")
m["dump_size_bytes"] = os.path.getsize("$DUMP_FILE")
with open(mf, "w") as f: json.dump(m, f, indent=2)
print(f"  Manifest updated: {mf}")
PYEOF

echo "[3/3] Restarting Neo4j ..."
neo4j start
for i in $(seq 1 30); do
    sleep 3
    neo4j status 2>/dev/null | grep -q "running" && { echo "  Running."; break; }
    echo "  Waiting ($i/30) ..."
done

echo ""
echo "=== Done ==="
echo "  $DUMP_FILE"
echo "  $SNAP_DIR/graphpop_${LABEL}_manifest.json"
echo ""
echo "To restore on another machine:"
echo "  sudo neo4j-admin database load neo4j --from-path=$SNAP_DIR --overwrite-destination=true"
