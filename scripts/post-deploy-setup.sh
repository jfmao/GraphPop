#!/usr/bin/env bash
# Post-deployment: create indexes and run smoke tests on graphpop procedures.
#
# Usage:
#   bash scripts/post-deploy-setup.sh
#
# Prerequisites:
#   - Neo4j must be running with graphpop-procedures.jar deployed
#   - Data must be imported (variant nodes, population nodes, etc.)
set -euo pipefail

NEO4J_USER="${NEO4J_USER:-neo4j}"
NEO4J_PASS="${NEO4J_PASS:-graphpop}"

echo "==> Creating indexes for Milestone 1.2..."

# Composite range index on (chr, pos) — critical for windowed queries
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CREATE INDEX variant_chr_pos IF NOT EXISTS FOR (v:Variant) ON (v.chr, v.pos);
CYPHER

# GenomicWindow indexes
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CREATE CONSTRAINT gw_id IF NOT EXISTS FOR (gw:GenomicWindow) REQUIRE gw.windowId IS UNIQUE;
CYPHER

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CREATE INDEX gw_run IF NOT EXISTS FOR (gw:GenomicWindow) ON (gw.run_id);
CYPHER

cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CREATE INDEX gw_range IF NOT EXISTS FOR (gw:GenomicWindow) ON (gw.chr, gw.start, gw.end);
CYPHER

# Full-text index on GOTerm.name
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CREATE FULLTEXT INDEX goterm_name IF NOT EXISTS FOR (t:GOTerm) ON EACH [t.name];
CYPHER

echo "    Indexes created."
echo ""

# -----------------------------------------------------------------------
# Verify procedures are registered
# -----------------------------------------------------------------------

echo "==> Registered graphpop procedures:"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" \
    "SHOW PROCEDURES YIELD name, mode WHERE name STARTS WITH 'graphpop' RETURN name, mode ORDER BY name;"
echo ""

# -----------------------------------------------------------------------
# Smoke tests
# -----------------------------------------------------------------------

echo "==> Smoke test: graphpop.diversity on chr22 16M-17M for AFR"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CALL graphpop.diversity('chr22', 16000000, 17000000, 'AFR')
YIELD pi, theta_w, tajima_d, het_exp, het_obs, fis, n_variants, n_segregating
RETURN pi, theta_w, tajima_d, het_exp, het_obs, fis, n_variants, n_segregating;
CYPHER
echo ""

echo "==> Smoke test: graphpop.divergence AFR vs EUR"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CALL graphpop.divergence('chr22', 16000000, 17000000, 'AFR', 'EUR')
YIELD fst_hudson, dxy, da, n_variants
RETURN fst_hudson, dxy, da, n_variants;
CYPHER
echo ""

echo "==> Smoke test: graphpop.sfs (folded)"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CALL graphpop.sfs('chr22', 16000000, 17000000, 'EUR', true)
YIELD n_variants, max_ac
RETURN n_variants, max_ac;
CYPHER
echo ""

echo "==> Smoke test: graphpop.genome_scan (1 Mb window, 500 kb step, first 5 windows)"
cypher-shell -u "$NEO4J_USER" -p "$NEO4J_PASS" <<'CYPHER'
CALL graphpop.genome_scan('chr22', 'AFR', 1000000, 500000, {})
YIELD window_id, start, end, pi, tajima_d, n_variants
RETURN window_id, start, end, pi, tajima_d, n_variants
LIMIT 5;
CYPHER
echo ""

echo "==> All smoke tests complete."
