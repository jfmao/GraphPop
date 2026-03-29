#!/bin/bash
#SBATCH --job-name=graphpop-ingest
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=graphpop_ingest_%j.log

# GraphPop — Combined import in one job (CSV generation + Neo4j load)
#
# This combines both steps. Neo4j is started and stopped automatically.
# Requires Neo4j to be installed in user space (see setup instructions).
#
# Usage:
#   sbatch slurm_ingest_single.sh /path/to/data.vcf.gz /path/to/pops.tsv

set -euo pipefail

VCF_INPUT="${1:?Usage: sbatch slurm_ingest_single.sh <vcf> <population_map>}"
POP_MAP="${2:?Usage: sbatch slurm_ingest_single.sh <vcf> <population_map>}"
NEO4J_HOME="${3:-$HOME/neo4j}"
NEO4J_DATA="${4:-/scratch/$USER/graphpop_db}"
CSV_DIR="${5:-/scratch/$USER/graphpop_csv}"

echo "=== GraphPop ingest (single job) ==="
echo "Input VCF:      $VCF_INPUT"
echo "Population map: $POP_MAP"
echo "Neo4j home:     $NEO4J_HOME"
echo "Neo4j data:     $NEO4J_DATA"
echo "CSV dir:        $CSV_DIR"
echo "Threads:        $SLURM_CPUS_PER_TASK"
echo "Started:        $(date)"

# Activate environment
# module load java/21
# source ~/graphpop-env/bin/activate

# Check filesystem
python -c "
from graphpop_import.cluster import check_neo4j_data_dir
result = check_neo4j_data_dir('$NEO4J_DATA')
if result['is_network']:
    import sys
    print('ERROR: ' + result['warning'], file=sys.stderr)
    sys.exit(1)
print(f\"Filesystem OK: {result['fs_type']} at {result['path']}\")
"

# Step 1: Generate CSVs
echo "--- Step 1: Generating CSVs ---"
python -m graphpop_import.csv_emitter \
    --vcf "$VCF_INPUT" \
    --population-map "$POP_MAP" \
    --output-dir "$CSV_DIR"

# Step 2: Load into Neo4j
echo "--- Step 2: Loading into Neo4j ---"
mkdir -p "$NEO4J_DATA"

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

# Step 3: Deploy procedures and start Neo4j
echo "--- Step 3: Starting Neo4j ---"
GRAPHPOP_JAR=$(find "$(dirname "$0")/../../graphpop-procedures/target" -name "graphpop-procedures-*.jar" -not -name "*-sources*" 2>/dev/null | head -1)
if [ -n "$GRAPHPOP_JAR" ]; then
    cp "$GRAPHPOP_JAR" "$NEO4J_HOME/plugins/"
    echo "Deployed procedures JAR: $GRAPHPOP_JAR"
fi

python -c "
from graphpop_import.cluster import start_neo4j
result = start_neo4j('$NEO4J_HOME', data_dir='$NEO4J_DATA', wait=True, timeout=120)
print(f\"Neo4j status: {result['status']}\")
"

echo "=== Neo4j is running. You can now run analyses. ==="
echo "=== Completed: $(date) ==="
