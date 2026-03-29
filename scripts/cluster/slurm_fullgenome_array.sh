#!/bin/bash
#SBATCH --job-name=graphpop-fullgenome
#SBATCH --array=1-12                # Adjust: 1-22 for human, 1-12 for rice
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=graphpop_%A_%a.log

# GraphPop — Full-genome analysis via array jobs (one per chromosome)
#
# Runs all per-population procedures for one chromosome per array task.
# Neo4j must be running on the database node.
#
# Setup:
#   export GRAPHPOP_URI=bolt://db-node:7687
#   export GRAPHPOP_PASSWORD=mypassword
#   export GRAPHPOP_DATABASE=rice3k
#
# Usage:
#   sbatch --array=1-12 slurm_fullgenome_array.sh          # rice (12 chr)
#   sbatch --array=1-22 slurm_fullgenome_array.sh          # human (22 chr)
#   sbatch --array=1-12 slurm_fullgenome_array.sh GJ-tmp   # single population

set -euo pipefail

# Load environment
module load java/21 2>/dev/null || true
source activate graphpop 2>/dev/null || conda activate graphpop 2>/dev/null || true

# Chromosome from array index
# Adjust prefix for your dataset: "chr" for human, "Chr" for rice
CHR_PREFIX="${CHR_PREFIX:-Chr}"
CHR="${CHR_PREFIX}${SLURM_ARRAY_TASK_ID}"

# Population list (override with first argument or env var)
if [ -n "${1:-}" ]; then
    POPULATIONS="$1"
elif [ -n "${GRAPHPOP_POPULATIONS:-}" ]; then
    POPULATIONS="$GRAPHPOP_POPULATIONS"
else
    # Default: auto-detect from database (requires graphpop inventory)
    POPULATIONS=$(graphpop query "MATCH (p:Population) WHERE p.n_samples > 1 RETURN p.populationId ORDER BY p.populationId" --format tsv 2>/dev/null | tail -n +2 | tr '\n' ',' | sed 's/,$//')
    if [ -z "$POPULATIONS" ]; then
        echo "ERROR: Could not detect populations. Set GRAPHPOP_POPULATIONS or pass as argument." >&2
        exit 1
    fi
fi

OUTDIR="${GRAPHPOP_OUTDIR:-results}"
mkdir -p "$OUTDIR"/{diversity,ihs,nsl,garud_h,roh}

echo "$(date): Starting chr=$CHR populations=$POPULATIONS"

IFS=',' read -ra POPS <<< "$POPULATIONS"
for POP in "${POPS[@]}"; do
    POP=$(echo "$POP" | tr -d '[:space:]')
    echo "$(date): $CHR / $POP"

    # Diversity (FAST PATH)
    graphpop diversity "$CHR" 1 300000000 "$POP" \
        -o "$OUTDIR/diversity/${POP}_${CHR}.tsv" 2>&1 || echo "WARN: diversity $POP $CHR failed"

    # iHS (FULL PATH, persist to graph)
    graphpop ihs "$CHR" "$POP" --persist \
        -o "$OUTDIR/ihs/${POP}_${CHR}.tsv" 2>&1 || echo "WARN: ihs $POP $CHR failed"

    # nSL (FULL PATH, persist)
    graphpop nsl "$CHR" "$POP" --persist \
        -o "$OUTDIR/nsl/${POP}_${CHR}.tsv" 2>&1 || echo "WARN: nsl $POP $CHR failed"

    # Garud's H (FULL PATH)
    graphpop garud-h "$CHR" "$POP" 100000 50000 \
        -o "$OUTDIR/garud_h/${POP}_${CHR}.tsv" 2>&1 || echo "WARN: garud_h $POP $CHR failed"

    # ROH (FULL PATH)
    graphpop roh "$CHR" "$POP" \
        -o "$OUTDIR/roh/${POP}_${CHR}.tsv" 2>&1 || echo "WARN: roh $POP $CHR failed"
done

echo "$(date): Completed chr=$CHR"
