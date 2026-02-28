#!/usr/bin/env bash
# Load VEP annotation CSVs (Gene nodes + HAS_CONSEQUENCE edges) into Neo4j
# using cypher-shell LOAD CSV.
#
# Usage:
#   bash scripts/load-annotations.sh [CSV_DIR]
#
# Defaults:
#   CSV_DIR = data/raw/1000g/csv_out
#
# Prerequisites:
#   - Neo4j must be RUNNING
#   - gene_nodes.csv and has_consequence_edges.csv must exist in CSV_DIR
#   - CSVs must be copied to Neo4j's import directory (or use file:/// with
#     dbms.directories.import configured)
#
# The script copies CSVs to /var/lib/neo4j/import/ so LOAD CSV can read them.
set -euo pipefail

CSV_DIR="${1:-data/raw/1000g/csv_out}"
CSV_DIR="$(cd "$CSV_DIR" && pwd)"

NEO4J_USER="${NEO4J_USER:-neo4j}"
NEO4J_PASS="${NEO4J_PASS:-graphpop}"
NEO4J_DB="${NEO4J_DB:-neo4j}"

IMPORT_DIR="/var/lib/neo4j/import"

echo "==> Load VEP Annotations into Neo4j"
echo "    CSV directory : $CSV_DIR"
echo "    Database      : $NEO4J_DB"
echo ""

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------

for f in gene_nodes.csv has_consequence_edges.csv; do
    if [[ ! -f "$CSV_DIR/$f" ]]; then
        echo "ERROR: Missing required file: $CSV_DIR/$f" >&2
        echo "       Run VEPParser first to generate these files." >&2
        exit 1
    fi
done

echo "==> CSV files found:"
ls -lh "$CSV_DIR"/gene_nodes.csv "$CSV_DIR"/has_consequence_edges.csv | \
    awk '{printf "    %-40s %s\n", $NF, $5}'
echo ""

# ---------------------------------------------------------------------------
# Copy CSVs to Neo4j import directory
# ---------------------------------------------------------------------------

echo "==> Copying CSVs to $IMPORT_DIR ..."
sudo mkdir -p "$IMPORT_DIR"
sudo cp "$CSV_DIR/gene_nodes.csv" "$IMPORT_DIR/"
sudo cp "$CSV_DIR/has_consequence_edges.csv" "$IMPORT_DIR/"
sudo chown neo4j:neo4j "$IMPORT_DIR"/gene_nodes.csv "$IMPORT_DIR"/has_consequence_edges.csv
echo "    Done."
echo ""

# ---------------------------------------------------------------------------
# Create indexes first (idempotent)
# ---------------------------------------------------------------------------

echo "==> Creating Gene constraints and indexes..."

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
CREATE CONSTRAINT gene_id IF NOT EXISTS FOR (g:Gene) REQUIRE g.geneId IS UNIQUE;
CYPHER

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
CREATE INDEX gene_symbol IF NOT EXISTS FOR (g:Gene) ON (g.symbol);
CYPHER

echo "    Done."
echo ""

# ---------------------------------------------------------------------------
# Load Gene nodes
# ---------------------------------------------------------------------------

echo "==> Loading Gene nodes..."

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
LOAD CSV WITH HEADERS FROM 'file:///gene_nodes.csv' AS row
MERGE (g:Gene {geneId: row.`geneId:ID(Gene)`})
SET g.symbol = row.symbol, g.biotype = row.biotype;
CYPHER

echo "    Done."

# Count
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
    "MATCH (g:Gene) RETURN count(g) AS gene_count;"
echo ""

# ---------------------------------------------------------------------------
# Load HAS_CONSEQUENCE edges (batched with periodic commit)
# ---------------------------------------------------------------------------

echo "==> Loading HAS_CONSEQUENCE edges..."

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
LOAD CSV WITH HEADERS FROM 'file:///has_consequence_edges.csv' AS row
CALL {
  WITH row
  MATCH (v:Variant {variantId: row.`:START_ID(Variant)`})
  MATCH (g:Gene {geneId: row.`:END_ID(Gene)`})
  CREATE (v)-[:HAS_CONSEQUENCE {
    consequence: row.consequence,
    impact: row.impact,
    feature: row.feature,
    feature_type: row.feature_type,
    sift_score: toFloatOrNull(row.`sift_score:float`),
    sift_pred: row.sift_pred,
    polyphen_score: toFloatOrNull(row.`polyphen_score:float`),
    polyphen_pred: row.polyphen_pred,
    cadd_phred: toFloatOrNull(row.`cadd_phred:float`),
    revel: toFloatOrNull(row.`revel:float`)
  }]->(g)
} IN TRANSACTIONS OF 10000 ROWS;
CYPHER

echo "    Done."

# Count
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
    "MATCH ()-[r:HAS_CONSEQUENCE]->() RETURN count(r) AS edge_count;"
echo ""

# ---------------------------------------------------------------------------
# Create relationship index
# ---------------------------------------------------------------------------

echo "==> Creating HAS_CONSEQUENCE index on impact..."

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
CREATE INDEX consequence_impact IF NOT EXISTS FOR ()-[r:HAS_CONSEQUENCE]-() ON (r.impact);
CYPHER

echo "    Done."
echo ""

# ---------------------------------------------------------------------------
# Verification queries
# ---------------------------------------------------------------------------

echo "==> Verification:"
echo ""

echo "--- Gene count ---"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
    "MATCH (g:Gene) RETURN count(g) AS genes;"

echo ""
echo "--- HAS_CONSEQUENCE count ---"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
    "MATCH ()-[r:HAS_CONSEQUENCE]->() RETURN count(r) AS edges;"

echo ""
echo "--- HIGH impact variants (top 10) ---"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
MATCH (v:Variant)-[r:HAS_CONSEQUENCE {impact: 'HIGH'}]->(g:Gene)
RETURN v.variantId AS variant, g.symbol AS gene, r.consequence AS consequence
LIMIT 10;
CYPHER

echo ""
echo "--- Top 10 genes by variant count ---"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
MATCH (g:Gene)<-[:HAS_CONSEQUENCE]-(v:Variant)
RETURN g.symbol AS gene, count(v) AS variants
ORDER BY variants DESC
LIMIT 10;
CYPHER

echo ""
echo "==> Annotation loading complete."
