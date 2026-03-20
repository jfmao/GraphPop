#!/usr/bin/env bash
# Configure Neo4j for GraphPop development on WSL2.
# Run with: sudo bash /mnt/workspace/GraphPop/scripts/configure-neo4j.sh
set -euo pipefail

CONF="/etc/neo4j/neo4j.conf"

echo "==> Backing up $CONF to ${CONF}.bak"
cp "$CONF" "${CONF}.bak"

echo "==> Applying GraphPop configuration..."

# 1. Right-size memory for 31 GB WSL2
#    Heap 4g, page cache 20g → leaves ~7 GB for OS
sed -i 's/^server\.memory\.heap\.initial_size=.*/server.memory.heap.initial_size=4g/' "$CONF"
sed -i 's/^server\.memory\.heap\.max_size=.*/server.memory.heap.max_size=4g/' "$CONF"
sed -i 's/^server\.memory\.pagecache\.size=.*/server.memory.pagecache.size=20g/' "$CONF"

# 2. Enable Vector API incubator module for graphpop-procedures SIMD
sed -i 's|^# server\.jvm\.additional=--add-modules=jdk\.incubator\.vector|server.jvm.additional=--add-modules=jdk.incubator.vector|' "$CONF"

# 3. Set default database to graphpop
sed -i 's/^#initial\.dbms\.default_database=.*/initial.dbms.default_database=graphpop/' "$CONF"

# 4. Allow graphpop procedures unrestricted access
#    Append after the existing (commented) unrestricted line
if ! grep -q 'graphpop' "$CONF"; then
    sed -i '/^#dbms\.security\.procedures\.unrestricted=/a dbms.security.procedures.unrestricted=graphpop.*' "$CONF"
fi

# 5. Set plugin directory to also look at our project output (symlink approach)
echo ""
echo "==> Creating symlink for graphpop-procedures plugin"
ln -sf /mnt/workspace/GraphPop/graphpop-procedures/target/graphpop-procedures-0.1.0-SNAPSHOT.jar \
       /var/lib/neo4j/plugins/graphpop-procedures.jar 2>/dev/null || true

echo "==> Configuration applied. Key settings:"
grep -n 'heap\.\(initial_size\|max_size\)\|pagecache\.size\|default_database\|unrestricted\|incubator\.vector' "$CONF" | grep -v '^#' || true

echo ""
echo "==> Start Neo4j with:  sudo neo4j start"
echo "==> Or as service:     sudo systemctl start neo4j"
