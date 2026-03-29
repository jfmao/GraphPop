#!/bin/bash
# GraphPop — Interactive session setup for HPC clusters
#
# This script is sourced (not submitted) to set up an interactive session
# with Neo4j running on a compute node.
#
# Usage:
#   srun --nodes=1 --cpus-per-task=16 --mem=128G --time=8:00:00 --pty bash
#   source slurm_interactive.sh [neo4j_home] [data_dir]
#
# After sourcing, Neo4j is running and you can run GraphPop analyses directly.
# When done, run: python -c "from graphpop_import.cluster import stop_neo4j; stop_neo4j('$NEO4J_HOME')"

NEO4J_HOME="${1:-$HOME/neo4j}"
NEO4J_DATA="${2:-/scratch/$USER/graphpop_db}"

echo "=== GraphPop Interactive Session ==="
echo "Neo4j home:  $NEO4J_HOME"
echo "Data dir:    $NEO4J_DATA"
echo "Node:        $(hostname)"

# Check filesystem
python -c "
from graphpop_import.cluster import check_neo4j_data_dir
result = check_neo4j_data_dir('$NEO4J_DATA')
if result['is_network']:
    print('WARNING: ' + result['warning'])
else:
    print(f\"Filesystem OK: {result['fs_type']}\")
"

# Start Neo4j
python -c "
from graphpop_import.cluster import start_neo4j
result = start_neo4j('$NEO4J_HOME', data_dir='$NEO4J_DATA', wait=True, timeout=120)
print(f\"Neo4j status: {result['status']} on bolt://localhost:{result['bolt_port']}\")
"

echo ""
echo "Neo4j is running. You can now run GraphPop analyses:"
echo "  python scripts/human_full_analysis.py --phase all"
echo "  python scripts/human_interpret.py pinsps"
echo ""
echo "When done:"
echo "  python -c \"from graphpop_import.cluster import stop_neo4j; stop_neo4j('$NEO4J_HOME')\""
echo ""
