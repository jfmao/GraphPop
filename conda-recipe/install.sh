#!/bin/bash
set -euo pipefail

# --- GraphPop conda installer ---
# This script is run by conda-build during package installation.
# It builds the Java procedures, installs all Python packages,
# and bundles the pre-compiled JAR with the CLI.

echo "=== Building GraphPop procedures (Java) ==="
cd graphpop-procedures
if [ -f "./mvnw" ]; then
    chmod +x ./mvnw
    ./mvnw package -DskipTests -q
else
    mvn package -DskipTests -q
fi
JAR_PATH=$(find target -name "graphpop-procedures-*.jar" -not -name "*original*" | head -1)
echo "  Built: ${JAR_PATH}"

# Install the JAR to a well-known location in the conda prefix
PLUGIN_DIR="${PREFIX}/share/graphpop/plugins"
mkdir -p "${PLUGIN_DIR}"
cp "${JAR_PATH}" "${PLUGIN_DIR}/graphpop-procedures.jar"
echo "  Installed JAR to: ${PLUGIN_DIR}/graphpop-procedures.jar"
cd ..

echo "=== Installing GraphPop CLI (Python) ==="
cd graphpop-cli
${PYTHON} -m pip install . --no-deps --no-build-isolation -q
cd ..

echo "=== Installing GraphPop Import Pipeline (Python) ==="
cd graphpop-import
${PYTHON} -m pip install . --no-deps --no-build-isolation -q
cd ..

echo "=== Installing GraphPop MCP Server (Python) ==="
cd graphpop-mcp
${PYTHON} -m pip install . --no-deps --no-build-isolation -q
cd ..

# Create a wrapper script that knows where the JAR is
WRAPPER="${PREFIX}/bin/graphpop-deploy-plugin"
cat > "${WRAPPER}" << 'SCRIPT'
#!/bin/bash
# Deploy the bundled GraphPop procedures plugin to a Neo4j installation.
# Usage: graphpop-deploy-plugin [NEO4J_HOME]
NEO4J_HOME="${1:-${HOME}/neo4j}"
PLUGIN_SRC="$(dirname $(dirname $(readlink -f $0)))/share/graphpop/plugins/graphpop-procedures.jar"
PLUGIN_DST="${NEO4J_HOME}/plugins/graphpop-procedures.jar"

if [ ! -f "${PLUGIN_SRC}" ]; then
    echo "Error: Plugin JAR not found at ${PLUGIN_SRC}" >&2
    exit 1
fi

mkdir -p "$(dirname ${PLUGIN_DST})"
cp "${PLUGIN_SRC}" "${PLUGIN_DST}"
echo "Deployed GraphPop plugin to ${PLUGIN_DST}"
echo "Restart Neo4j to activate: graphpop stop && graphpop start"
SCRIPT
chmod +x "${WRAPPER}"

echo "=== GraphPop installation complete ==="
echo ""
echo "Quick start:"
echo "  graphpop setup --password mypass    # Download & configure Neo4j"
echo "  graphpop start                      # Start the database"
echo "  graphpop --help                     # See all commands"
