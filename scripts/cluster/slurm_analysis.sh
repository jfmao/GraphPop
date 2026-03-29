#!/bin/bash
#SBATCH --job-name=graphpop-analysis
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=graphpop_analysis_%j.log

# GraphPop — Run population genetics analyses via stored procedures
#
# Neo4j must be running with the GraphPop database loaded and procedures
# deployed. Start it first with neo4j-start, or use an interactive session.
#
# Usage:
#   sbatch slurm_analysis.sh full_analysis.py --phase all
#   sbatch slurm_analysis.sh human_interpret.py pinsps
#   sbatch slurm_analysis.sh rice_full_analysis.py --phase all

set -euo pipefail

SCRIPT="${1:?Usage: sbatch slurm_analysis.sh <analysis_script> [args...]}"
shift

echo "=== GraphPop analysis ==="
echo "Script:         $SCRIPT"
echo "Arguments:      $*"
echo "Threads:        $SLURM_CPUS_PER_TASK"
echo "Started:        $(date)"

# Activate environment
# module load java/21
# source ~/graphpop-env/bin/activate

python "$SCRIPT" "$@"

echo "=== Completed: $(date) ==="
