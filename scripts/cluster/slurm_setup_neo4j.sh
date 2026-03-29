#!/bin/bash
#SBATCH --job-name=graphpop-setup
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=graphpop_setup_%j.log

# GraphPop — One-time Neo4j setup on a cluster node
#
# Downloads Neo4j, configures memory, sets password, deploys GraphPop plugin,
# and creates the GraphPop config file. Run once on the node that will host
# the database.
#
# Prerequisites:
#   - Java 21+ (module load java/21)
#   - graphpop CLI installed (pip install -e graphpop-cli/)
#   - GraphPop procedures JAR built (cd graphpop-procedures && mvn package)
#
# Usage:
#   sbatch slurm_setup_neo4j.sh
#   # Or interactively:
#   srun --cpus-per-task=4 --mem=64G --time=2:00:00 --pty bash
#   bash slurm_setup_neo4j.sh

set -euo pipefail

module load java/21 2>/dev/null || true
source activate graphpop 2>/dev/null || conda activate graphpop 2>/dev/null || true

# Configuration
NEO4J_HOME="${NEO4J_HOME:-$HOME/neo4j}"
PAGECACHE="${PAGECACHE:-20g}"
HEAP="${HEAP:-4g}"
PASSWORD="${GRAPHPOP_PASSWORD:-graphpop}"
PLUGIN_JAR="${PLUGIN_JAR:-graphpop-procedures/target/graphpop-procedures-0.1.0-SNAPSHOT.jar}"

# Use local SSD if available
if [ -d "/local/$USER" ]; then
    DATA_DIR="/local/$USER/neo4j_data"
elif [ -d "/scratch/$USER" ]; then
    DATA_DIR="/scratch/$USER/neo4j_data"
else
    DATA_DIR="$NEO4J_HOME/data"
fi

echo "=== GraphPop Neo4j Setup ==="
echo "Neo4j home:  $NEO4J_HOME"
echo "Data dir:    $DATA_DIR"
echo "Pagecache:   $PAGECACHE"
echo "Heap:        $HEAP"
echo "Node:        $(hostname)"
echo ""

# Setup Neo4j
if [ -f "$PLUGIN_JAR" ]; then
    graphpop setup --neo4j-home "$NEO4J_HOME" \
        --pagecache "$PAGECACHE" --heap "$HEAP" \
        --password "$PASSWORD" \
        --deploy-plugin "$PLUGIN_JAR"
else
    echo "WARNING: Plugin JAR not found at $PLUGIN_JAR"
    echo "Run: cd graphpop-procedures && mvn package -DskipTests"
    graphpop setup --neo4j-home "$NEO4J_HOME" \
        --pagecache "$PAGECACHE" --heap "$HEAP" \
        --password "$PASSWORD"
fi

# Symlink data directory to fast storage
if [ "$DATA_DIR" != "$NEO4J_HOME/data" ]; then
    mkdir -p "$DATA_DIR"
    if [ -d "$NEO4J_HOME/data" ] && [ ! -L "$NEO4J_HOME/data" ]; then
        mv "$NEO4J_HOME/data"/* "$DATA_DIR/" 2>/dev/null || true
        rm -rf "$NEO4J_HOME/data"
    fi
    ln -sf "$DATA_DIR" "$NEO4J_HOME/data"
    echo "Data directory symlinked to $DATA_DIR"
fi

# Configure Neo4j to listen on all interfaces (for remote connections)
CONF="$NEO4J_HOME/conf/neo4j.conf"
if ! grep -q "server.default_listen_address=0.0.0.0" "$CONF"; then
    echo "server.default_listen_address=0.0.0.0" >> "$CONF"
    echo "Configured Neo4j to accept remote connections"
fi

# Start and verify
graphpop start
sleep 5
graphpop status

echo ""
echo "=== Setup Complete ==="
echo "Database node: $(hostname)"
echo "Connect from other nodes with:"
echo "  export GRAPHPOP_URI=bolt://$(hostname):7687"
echo "  export GRAPHPOP_PASSWORD=$PASSWORD"
echo ""
echo "Save this for later:"
echo "  echo 'export GRAPHPOP_URI=bolt://$(hostname):7687' >> ~/.bashrc"
