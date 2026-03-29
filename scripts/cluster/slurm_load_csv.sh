#!/bin/bash
#SBATCH --job-name=graphpop-load
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=graphpop_load_%j.log

# GraphPop — Step 2: Load pre-generated CSVs into Neo4j
#
# This requires Neo4j to be available. Run on a node with local SSD
# for the Neo4j data directory. Uses neo4j-admin import for bulk loading.
#
# Usage:
#   sbatch slurm_load_csv.sh /scratch/user/graphpop_csv /home/user/neo4j

set -euo pipefail

CSV_DIR="${1:?Usage: sbatch slurm_load_csv.sh <csv_dir> <neo4j_home> [data_dir]}"
NEO4J_HOME="${2:?Usage: sbatch slurm_load_csv.sh <csv_dir> <neo4j_home> [data_dir]}"
NEO4J_DATA="${3:-/scratch/$USER/graphpop_db}"

echo "=== GraphPop load-csv ==="
echo "CSV dir:        $CSV_DIR"
echo "Neo4j home:     $NEO4J_HOME"
echo "Neo4j data:     $NEO4J_DATA"
echo "Started:        $(date)"

# Activate environment
# module load java/21
# source ~/graphpop-env/bin/activate

# Check filesystem before proceeding
python -c "
from graphpop_import.cluster import check_neo4j_data_dir
result = check_neo4j_data_dir('$NEO4J_DATA')
if result['is_network']:
    import sys
    print('ERROR: ' + result['warning'], file=sys.stderr)
    sys.exit(1)
print(f\"Filesystem OK: {result['fs_type']} at {result['path']}\")
"

# Ensure data directory exists
mkdir -p "$NEO4J_DATA"

# Run neo4j-admin import
NEO4J_CONF="$NEO4J_HOME/conf" "$NEO4J_HOME/bin/neo4j-admin" database import full neo4j \
    --overwrite-destination \
    --nodes=Variant="$CSV_DIR/variant_header.csv,$CSV_DIR/variants.csv" \
    --nodes=Sample="$CSV_DIR/sample_header.csv,$CSV_DIR/samples.csv" \
    --nodes=Population="$CSV_DIR/population_header.csv,$CSV_DIR/populations.csv" \
    --nodes=Chromosome="$CSV_DIR/chromosome_header.csv,$CSV_DIR/chromosomes.csv" \
    --nodes=GenomicWindow="$CSV_DIR/window_header.csv,$CSV_DIR/windows.csv" \
    --relationships=NEXT="$CSV_DIR/next_header.csv,$CSV_DIR/next_edges.csv" \
    --relationships=BELONGS_TO="$CSV_DIR/belongs_to_header.csv,$CSV_DIR/belongs_to_edges.csv" \
    --delimiter=',' \
    --array-delimiter=';'

echo "=== Completed: $(date) ==="
echo "Database loaded at: $NEO4J_DATA"
