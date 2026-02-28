#!/usr/bin/env bash
# Load Pathway/GOTerm annotation CSVs into Neo4j using cypher-shell LOAD CSV.
#
# Usage:
#   bash scripts/load-pathways.sh [CSV_DIR]
#
# Defaults:
#   CSV_DIR = data/raw/1000g/csv_out
#
# Prerequisites:
#   - Neo4j must be RUNNING
#   - pathway_nodes.csv, goterm_nodes.csv, in_pathway_edges.csv,
#     has_go_term_edges.csv must exist in CSV_DIR
#   - Run PathwayParser first to generate these files
set -euo pipefail

CSV_DIR="${1:-data/raw/1000g/csv_out}"
CSV_DIR="$(cd "$CSV_DIR" && pwd)"

NEO4J_USER="${NEO4J_USER:-neo4j}"
NEO4J_PASS="${NEO4J_PASS:-graphpop}"
NEO4J_DB="${NEO4J_DB:-neo4j}"

IMPORT_DIR="/var/lib/neo4j/import"

echo "==> Load Pathway/GOTerm Annotations into Neo4j"
echo "    CSV directory : $CSV_DIR"
echo "    Database      : $NEO4J_DB"
echo ""

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------

REQUIRED_FILES=(pathway_nodes.csv in_pathway_edges.csv)
OPTIONAL_FILES=(goterm_nodes.csv has_go_term_edges.csv)

for f in "${REQUIRED_FILES[@]}"; do
    if [[ ! -f "$CSV_DIR/$f" ]]; then
        echo "ERROR: Missing required file: $CSV_DIR/$f" >&2
        echo "       Run PathwayParser first to generate these files." >&2
        exit 1
    fi
done

HAS_GO=true
for f in "${OPTIONAL_FILES[@]}"; do
    if [[ ! -f "$CSV_DIR/$f" ]]; then
        echo "WARNING: Missing GO file: $CSV_DIR/$f (BioMart may have failed)" >&2
        HAS_GO=false
        break
    fi
done

echo "==> CSV files found:"
ls -lh "$CSV_DIR"/pathway_nodes.csv "$CSV_DIR"/in_pathway_edges.csv | \
    awk '{printf "    %-40s %s\n", $NF, $5}'
if $HAS_GO; then
    ls -lh "$CSV_DIR"/goterm_nodes.csv "$CSV_DIR"/has_go_term_edges.csv | \
        awk '{printf "    %-40s %s\n", $NF, $5}'
fi
echo ""

# ---------------------------------------------------------------------------
# Copy CSVs to Neo4j import directory
# ---------------------------------------------------------------------------

echo "==> Copying CSVs to $IMPORT_DIR ..."
sudo mkdir -p "$IMPORT_DIR"
sudo cp "$CSV_DIR/pathway_nodes.csv" "$IMPORT_DIR/"
sudo cp "$CSV_DIR/in_pathway_edges.csv" "$IMPORT_DIR/"
sudo chown neo4j:neo4j "$IMPORT_DIR"/pathway_nodes.csv "$IMPORT_DIR"/in_pathway_edges.csv

if $HAS_GO; then
    sudo cp "$CSV_DIR/goterm_nodes.csv" "$IMPORT_DIR/"
    sudo cp "$CSV_DIR/has_go_term_edges.csv" "$IMPORT_DIR/"
    sudo chown neo4j:neo4j "$IMPORT_DIR"/goterm_nodes.csv "$IMPORT_DIR"/has_go_term_edges.csv
fi
echo "    Done."
echo ""

# ---------------------------------------------------------------------------
# Create constraints and indexes
# ---------------------------------------------------------------------------

echo "==> Creating Pathway constraints and indexes..."

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
CREATE CONSTRAINT pathway_id IF NOT EXISTS FOR (p:Pathway) REQUIRE p.pathwayId IS UNIQUE;
CYPHER

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
CREATE INDEX pathway_name IF NOT EXISTS FOR (p:Pathway) ON (p.name);
CYPHER

if $HAS_GO; then
    echo "==> Creating GOTerm constraints and indexes..."

    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
CREATE CONSTRAINT goterm_id IF NOT EXISTS FOR (t:GOTerm) REQUIRE t.goId IS UNIQUE;
CYPHER

    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
CREATE INDEX goterm_aspect IF NOT EXISTS FOR (t:GOTerm) ON (t.aspect);
CYPHER
fi

echo "    Done."
echo ""

# ---------------------------------------------------------------------------
# Load Pathway nodes
# ---------------------------------------------------------------------------

echo "==> Loading Pathway nodes..."

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
LOAD CSV WITH HEADERS FROM 'file:///pathway_nodes.csv' AS row
MERGE (p:Pathway {pathwayId: row.`pathwayId:ID(Pathway)`})
SET p.name = row.name, p.source = row.source;
CYPHER

echo "    Done."
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
    "MATCH (p:Pathway) RETURN count(p) AS pathway_count;"
echo ""

# ---------------------------------------------------------------------------
# Load IN_PATHWAY edges
# ---------------------------------------------------------------------------

echo "==> Loading IN_PATHWAY edges..."

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
LOAD CSV WITH HEADERS FROM 'file:///in_pathway_edges.csv' AS row
CALL {
  WITH row
  MATCH (g:Gene {geneId: row.`:START_ID(Gene)`})
  MATCH (p:Pathway {pathwayId: row.`:END_ID(Pathway)`})
  CREATE (g)-[:IN_PATHWAY {evidence: row.evidence}]->(p)
} IN TRANSACTIONS OF 5000 ROWS;
CYPHER

echo "    Done."
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
    "MATCH ()-[r:IN_PATHWAY]->() RETURN count(r) AS in_pathway_count;"
echo ""

# ---------------------------------------------------------------------------
# Load GOTerm nodes (if available)
# ---------------------------------------------------------------------------

if $HAS_GO; then
    echo "==> Loading GOTerm nodes..."

    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
LOAD CSV WITH HEADERS FROM 'file:///goterm_nodes.csv' AS row
MERGE (t:GOTerm {goId: row.`goId:ID(GOTerm)`})
SET t.name = row.name, t.aspect = row.aspect;
CYPHER

    echo "    Done."
    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
        "MATCH (t:GOTerm) RETURN count(t) AS goterm_count;"
    echo ""

    # -----------------------------------------------------------------------
    # Load HAS_GO_TERM edges
    # -----------------------------------------------------------------------

    echo "==> Loading HAS_GO_TERM edges..."

    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
LOAD CSV WITH HEADERS FROM 'file:///has_go_term_edges.csv' AS row
CALL {
  WITH row
  MATCH (g:Gene {geneId: row.`:START_ID(Gene)`})
  MATCH (t:GOTerm {goId: row.`:END_ID(GOTerm)`})
  CREATE (g)-[:HAS_GO_TERM]->(t)
} IN TRANSACTIONS OF 5000 ROWS;
CYPHER

    echo "    Done."
    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
        "MATCH ()-[r:HAS_GO_TERM]->() RETURN count(r) AS has_go_term_count;"
    echo ""
fi

# ---------------------------------------------------------------------------
# Verification queries
# ---------------------------------------------------------------------------

echo "==> Verification:"
echo ""

echo "--- Pathway count ---"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
    "MATCH (p:Pathway) RETURN count(p) AS pathways;"

echo ""
echo "--- Top 10 pathways by gene count ---"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
MATCH (g:Gene)-[:IN_PATHWAY]->(p:Pathway)
RETURN p.name AS pathway, count(g) AS genes
ORDER BY genes DESC
LIMIT 10;
CYPHER

if $HAS_GO; then
    echo ""
    echo "--- GOTerm count ---"
    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" \
        "MATCH (t:GOTerm) RETURN count(t) AS goterms;"

    echo ""
    echo "--- Top 10 GO terms by gene count ---"
    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
MATCH (g:Gene)-[:HAS_GO_TERM]->(t:GOTerm)
RETURN t.name AS term, t.aspect AS aspect, count(g) AS genes
ORDER BY genes DESC
LIMIT 10;
CYPHER

    echo ""
    echo "--- GO aspect distribution ---"
    cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" -d "$NEO4J_DB" <<'CYPHER'
MATCH (t:GOTerm)
RETURN t.aspect AS aspect, count(t) AS count
ORDER BY count DESC;
CYPHER
fi

echo ""
echo "==> Pathway/GOTerm loading complete."
