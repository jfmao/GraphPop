#!/usr/bin/env bash
# restore_database.sh — Load a GraphPop snapshot dump into a local Neo4j instance.
#
# Intended for researchers who download a GraphPop database from a public repository
# (Zenodo, Figshare, OSF, etc.) and want to reproduce analyses.
#
# Usage:
#   sudo bash scripts/restore_database.sh <path/to/graphpop_*.dump> [TARGET_DB]
#
# Examples:
#   sudo bash scripts/restore_database.sh data/snapshots/graphpop_human_1000g_v1_20260320.dump
#   sudo bash scripts/restore_database.sh graphpop_human.dump graphpop_human
#
# After restoring, connect with:
#   Neo4j Browser : http://localhost:7474
#   Bolt URI      : bolt://localhost:7687
#   Credentials   : neo4j / <your-password>
#   Database      : <TARGET_DB>

set -euo pipefail

# ── Args ──────────────────────────────────────────────────────────────────────
DUMP_FILE="${1:-}"
TARGET_DB="${2:-}"

if [[ -z "$DUMP_FILE" ]]; then
    echo "Usage: sudo bash scripts/restore_database.sh <dump_file> [target_db_name]"
    echo ""
    echo "  dump_file      — path to a .dump file created by snapshot_database.sh"
    echo "  target_db_name — name for the restored database (default: derived from filename)"
    exit 1
fi

if [[ ! -f "$DUMP_FILE" ]]; then
    echo "ERROR: Dump file not found: $DUMP_FILE"
    exit 1
fi

# Derive target DB name from filename if not given
# graphpop_human_1000g_v1_20260320.dump → graphpop_human_1000g_v1
if [[ -z "$TARGET_DB" ]]; then
    BASENAME=$(basename "$DUMP_FILE" .dump)
    # Strip trailing _YYYYMMDD date stamp
    TARGET_DB=$(echo "$BASENAME" | sed 's/_[0-9]\{8\}$//')
fi

DUMP_ABS="$(cd "$(dirname "$DUMP_FILE")" && pwd)/$(basename "$DUMP_FILE")"
DUMP_DIR="$(dirname "$DUMP_ABS")"

# ── Check for companion manifest ──────────────────────────────────────────────
MANIFEST="${DUMP_ABS%.dump}.manifest.json"
if [[ -f "$MANIFEST" ]]; then
    echo "================================================================"
    echo "  GraphPop Database Restore"
    echo "  Manifest found: $(basename "$MANIFEST")"
    echo ""
    python3 -c "
import json
with open('$MANIFEST') as f:
    m = json.load(f)
print(f'  Label          : {m.get(\"label\",\"?\")}')
print(f'  Snapshot date  : {m.get(\"snapshot_date\",\"?\")}')
print(f'  Neo4j version  : {m.get(\"neo4j_version\",\"?\")}')
print(f'  Dump size      : {m.get(\"dump_size_bytes\",0)/1e9:.1f} GB')
print()
nc = m.get('node_counts', {})
for k, v in nc.items():
    if v: print(f'  {k:<20}: {v:>15,}')
print()
st = m.get('embedded_statistics', {})
for k, v in st.items():
    if v: print(f'  {k:<30}: {v:>12,}')
" 2>/dev/null || true
    echo "================================================================"
else
    echo "================================================================"
    echo "  GraphPop Database Restore"
    echo "  Dump : $(basename "$DUMP_ABS")"
    echo "  Target database: $TARGET_DB"
    echo "================================================================"
fi

echo ""
read -r -p "Restore dump as database '$TARGET_DB'? This cannot be undone. [y/N] " CONFIRM
if [[ "$CONFIRM" != "y" && "$CONFIRM" != "Y" ]]; then
    echo "Aborted."
    exit 0
fi

# ── Stop Neo4j ────────────────────────────────────────────────────────────────
echo ""
echo "[1/4] Stopping Neo4j..."
neo4j stop 2>/dev/null || true
for i in $(seq 1 30); do
    sleep 2
    if ! neo4j status 2>/dev/null | grep -q "running"; then
        echo "  Neo4j stopped."
        break
    fi
    echo "  Waiting... ($i)"
done

# ── Load dump ─────────────────────────────────────────────────────────────────
echo ""
echo "[2/4] Loading dump into database '$TARGET_DB'..."
echo "  Source: $DUMP_ABS"
neo4j-admin database load \
    --overwrite-destination=true \
    --from-path="$DUMP_DIR" \
    "$TARGET_DB"
echo "  Load complete."

# ── Start Neo4j ───────────────────────────────────────────────────────────────
echo ""
echo "[3/4] Starting Neo4j..."
neo4j start
for i in $(seq 1 30); do
    sleep 3
    if neo4j status 2>/dev/null | grep -q "running"; then
        echo "  Neo4j running."
        break
    fi
    echo "  Waiting... ($i)"
done

# ── Activate database ─────────────────────────────────────────────────────────
echo ""
echo "[4/4] Activating database '$TARGET_DB'..."
# Try cypher-shell if available
if command -v cypher-shell &>/dev/null; then
    cypher-shell -d system \
        "CREATE DATABASE \`$TARGET_DB\` IF NOT EXISTS;" \
        2>/dev/null && echo "  Database activated via cypher-shell." || \
        echo "  (Database may already be active or needs manual activation — see below.)"
else
    echo "  cypher-shell not found. Activate manually in Neo4j Browser:"
    echo "    :use system"
    echo "    CREATE DATABASE \`$TARGET_DB\` IF NOT EXISTS;"
fi

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "================================================================"
echo "  Restore complete."
echo ""
echo "  Connect in Neo4j Browser: http://localhost:7474"
echo "    :use $TARGET_DB"
echo ""
echo "  Bolt URI : bolt://localhost:7687"
echo "  Database : $TARGET_DB"
echo ""
echo "  Quick verification:"
echo "    MATCH (v:Variant) RETURN count(v) AS variants;"
echo "    MATCH (s:Sample)  RETURN count(s) AS samples;"
echo "    MATCH (w:GenomicWindow) WHERE w.h12 IS NOT NULL"
echo "      RETURN w.chr, w.h12, w.sweep_type ORDER BY w.h12 DESC LIMIT 5;"
echo "================================================================"
