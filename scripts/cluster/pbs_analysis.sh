#!/bin/bash
#PBS -N graphpop-analysis
#PBS -l ncpus=8
#PBS -l mem=64GB
#PBS -l walltime=12:00:00
#PBS -o graphpop_analysis.log
#PBS -j oe

# GraphPop — Run analyses (PBS/Torque)
#
# Neo4j must be running with GraphPop database loaded.
#
# Usage:
#   qsub -v SCRIPT=scripts/human_full_analysis.py,ARGS="--phase all" pbs_analysis.sh

set -euo pipefail

SCRIPT="${SCRIPT:?Set SCRIPT via qsub -v}"
ARGS="${ARGS:-}"

cd "$PBS_O_WORKDIR" 2>/dev/null || true

echo "=== GraphPop analysis (PBS) ==="
echo "Script:         $SCRIPT"
echo "Arguments:      $ARGS"
echo "Started:        $(date)"

# Activate environment
# module load java/21
# source ~/graphpop-env/bin/activate

python $SCRIPT $ARGS

echo "=== Completed: $(date) ==="
