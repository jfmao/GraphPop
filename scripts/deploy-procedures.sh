#!/usr/bin/env bash
# Deploy graphpop-procedures JAR to Neo4j plugins directory.
#
# Usage:
#   sudo bash scripts/deploy-procedures.sh
#
# What it does:
#   1. Copies the shaded JAR to /var/lib/neo4j/plugins/
#   2. Adds the Vector API incubator module to Neo4j JVM args
#   3. Whitelists graphpop.* procedures in neo4j.conf
#   4. Restarts Neo4j
set -euo pipefail

JAR="graphpop-procedures/target/graphpop-procedures-0.1.0-SNAPSHOT.jar"
PLUGINS_DIR="/var/lib/neo4j/plugins"
NEO4J_CONF="/etc/neo4j/neo4j.conf"

if [[ ! -f "$JAR" ]]; then
    echo "ERROR: JAR not found at $JAR" >&2
    echo "       Run: cd graphpop-procedures && JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64 ./mvnw package" >&2
    exit 1
fi

echo "==> Deploying graphpop-procedures to Neo4j"

# 1. Copy JAR
echo "    Copying JAR to $PLUGINS_DIR/ ..."
mkdir -p "$PLUGINS_DIR"
cp "$JAR" "$PLUGINS_DIR/graphpop-procedures.jar"
chown neo4j:neo4j "$PLUGINS_DIR/graphpop-procedures.jar"
echo "    Done: $(ls -lh "$PLUGINS_DIR/graphpop-procedures.jar" | awk '{print $5, $NF}')"

# 2. Add Vector API module — uncomment if present, or add new line
if grep -q '^# *server.jvm.additional=--add-modules=jdk.incubator.vector' "$NEO4J_CONF" 2>/dev/null; then
    echo "    Uncommenting jdk.incubator.vector in neo4j.conf..."
    sed -i 's/^# *server.jvm.additional=--add-modules=jdk.incubator.vector/server.jvm.additional=--add-modules=jdk.incubator.vector/' "$NEO4J_CONF"
elif ! grep -q '^server.jvm.additional=--add-modules=jdk.incubator.vector' "$NEO4J_CONF" 2>/dev/null; then
    echo "    Adding jdk.incubator.vector to JVM args..."
    echo "" >> "$NEO4J_CONF"
    echo "# GraphPop: enable Java Vector API for SIMD procedures" >> "$NEO4J_CONF"
    echo "server.jvm.additional=--add-modules=jdk.incubator.vector" >> "$NEO4J_CONF"
fi

# 3. Whitelist graphpop.* procedures if not already present
if ! grep -q 'graphpop' "$NEO4J_CONF" 2>/dev/null; then
    echo "    Whitelisting graphpop.* procedures..."
    echo "" >> "$NEO4J_CONF"
    echo "# GraphPop: allow custom procedures" >> "$NEO4J_CONF"
    echo "dbms.security.procedures.unrestricted=graphpop.*" >> "$NEO4J_CONF"
fi

# 4. Restart Neo4j
echo "    Restarting Neo4j..."
neo4j restart 2>/dev/null || systemctl restart neo4j 2>/dev/null || service neo4j restart 2>/dev/null
echo "    Waiting for Neo4j to come online..."
sleep 5

# 5. Verify
for i in 1 2 3 4 5 6; do
    if cypher-shell -u neo4j -p graphpop "RETURN 1" >/dev/null 2>&1; then
        echo "    Neo4j is online."
        break
    fi
    if [ $i -eq 6 ]; then
        echo "    WARNING: Neo4j may not be fully started yet. Check: neo4j status"
    fi
    sleep 5
done

echo ""
echo "==> Checking registered procedures..."
cypher-shell -u neo4j -p graphpop "SHOW PROCEDURES YIELD name WHERE name STARTS WITH 'graphpop' RETURN name ORDER BY name;" 2>/dev/null || \
    echo "    (Could not verify — try: cypher-shell -u neo4j -p graphpop \"SHOW PROCEDURES YIELD name WHERE name STARTS WITH 'graphpop' RETURN name;\")"

echo ""
echo "==> Deployment complete."
