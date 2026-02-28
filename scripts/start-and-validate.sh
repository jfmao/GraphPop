#!/usr/bin/env bash
# Run as: sudo bash /mnt/e/GraphPop/scripts/start-and-validate.sh
set -euo pipefail

echo "==> Checking memory config..."
grep -n '^server\.memory\.' /etc/neo4j/neo4j.conf

echo ""
echo "==> Running neo4j-admin import..."
CSV_DIR="/mnt/e/GraphPop/data/raw/1000g/csv_out"

neo4j-admin database import full \
    --overwrite-destination=true \
    --array-delimiter=";" \
    --id-type=string \
    --skip-bad-relationships=true \
    --skip-duplicate-nodes=true \
    --bad-tolerance=1000 \
    --threads="$(nproc)" \
    --high-parallel-io=on \
    --nodes=Variant="$CSV_DIR/variant_nodes.csv" \
    --nodes=Sample="$CSV_DIR/sample_nodes.csv" \
    --nodes=Population="$CSV_DIR/population_nodes.csv" \
    --nodes=Chromosome="$CSV_DIR/chromosome_nodes.csv" \
    --relationships=CARRIES="$CSV_DIR/carries_edges.csv" \
    --relationships=NEXT="$CSV_DIR/next_edges.csv" \
    --relationships=ON_CHROMOSOME="$CSV_DIR/on_chromosome_edges.csv" \
    --relationships=IN_POPULATION="$CSV_DIR/in_population_edges.csv" \
    --report-file="$CSV_DIR/import.report" \
    graphpop

echo ""
echo "==> Starting Neo4j..."
neo4j start

echo "==> Waiting for Neo4j to be ready..."
for i in $(seq 1 30); do
    if cypher-shell -d system "RETURN 1;" >/dev/null 2>&1; then
        echo "    Neo4j is ready after ${i}s"
        break
    fi
    sleep 1
done

echo ""
echo "==> Creating database..."
cypher-shell -d system "CREATE DATABASE graphpop IF NOT EXISTS;"
sleep 3

echo ""
echo "==> Applying indexes..."
cypher-shell -d graphpop << 'CYPHER'
CREATE CONSTRAINT variant_id IF NOT EXISTS FOR (v:Variant) REQUIRE v.`variantId` IS UNIQUE;
CREATE CONSTRAINT sample_id IF NOT EXISTS FOR (s:Sample) REQUIRE s.`sampleId` IS UNIQUE;
CREATE CONSTRAINT population_id IF NOT EXISTS FOR (p:Population) REQUIRE p.`populationId` IS UNIQUE;
CREATE CONSTRAINT chromosome_id IF NOT EXISTS FOR (c:Chromosome) REQUIRE c.`chromosomeId` IS UNIQUE;
CREATE INDEX variant_chr_pos IF NOT EXISTS FOR (v:Variant) ON (v.chr, v.`pos`);
CREATE INDEX variant_type IF NOT EXISTS FOR (v:Variant) ON (v.variant_type);
CREATE INDEX variant_af_total IF NOT EXISTS FOR (v:Variant) ON (v.af_total);
CYPHER

echo ""
echo "==> Running validation queries..."
cypher-shell -d graphpop << 'CYPHER'
// Node counts
MATCH (v:Variant) RETURN 'Variant' AS label, count(v) AS count
UNION ALL
MATCH (s:Sample) RETURN 'Sample' AS label, count(s) AS count
UNION ALL
MATCH (p:Population) RETURN 'Population' AS label, count(p) AS count
UNION ALL
MATCH (c:Chromosome) RETURN 'Chromosome' AS label, count(c) AS count;
CYPHER

echo ""
cypher-shell -d graphpop << 'CYPHER'
// Relationship counts
MATCH ()-[r:CARRIES]->() RETURN 'CARRIES' AS type, count(r) AS count
UNION ALL
MATCH ()-[r:NEXT]->() RETURN 'NEXT' AS type, count(r) AS count
UNION ALL
MATCH ()-[r:ON_CHROMOSOME]->() RETURN 'ON_CHROMOSOME' AS type, count(r) AS count
UNION ALL
MATCH ()-[r:IN_POPULATION]->() RETURN 'IN_POPULATION' AS type, count(r) AS count;
CYPHER

echo ""
echo "==> Done! Neo4j is running with the graphpop database."
