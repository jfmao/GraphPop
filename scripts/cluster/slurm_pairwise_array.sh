#!/bin/bash
#SBATCH --job-name=graphpop-pairwise
#SBATCH --array=1-12                # One per chromosome
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=graphpop_pair_%A_%a.log

# GraphPop — Pairwise statistics (XP-EHH, divergence) via array jobs
#
# Runs pairwise XP-EHH and divergence for specified population pairs
# on one chromosome per array task.
#
# Usage:
#   export GRAPHPOP_URI=bolt://db-node:7687
#   export GRAPHPOP_PAIRS="GJ-tmp:GJ-trp,GJ-tmp:XI-1A,XI-1A:cA-Aus"
#   sbatch --array=1-12 slurm_pairwise_array.sh

set -euo pipefail

module load java/21 2>/dev/null || true
source activate graphpop 2>/dev/null || conda activate graphpop 2>/dev/null || true

CHR_PREFIX="${CHR_PREFIX:-Chr}"
CHR="${CHR_PREFIX}${SLURM_ARRAY_TASK_ID}"

# Population pairs (colon-separated within pair, comma-separated between pairs)
PAIRS="${GRAPHPOP_PAIRS:?Set GRAPHPOP_PAIRS=POP1:POP2,POP3:POP4,...}"

OUTDIR="${GRAPHPOP_OUTDIR:-results}"
mkdir -p "$OUTDIR"/{xpehh,divergence}

echo "$(date): Starting pairwise analysis chr=$CHR"

IFS=',' read -ra PAIR_LIST <<< "$PAIRS"
for PAIR in "${PAIR_LIST[@]}"; do
    POP1=$(echo "$PAIR" | cut -d: -f1 | tr -d '[:space:]')
    POP2=$(echo "$PAIR" | cut -d: -f2 | tr -d '[:space:]')
    PAIR_ID="${POP1}_vs_${POP2}"

    echo "$(date): $CHR / $PAIR_ID"

    # XP-EHH (FULL PATH, persist)
    graphpop xpehh "$CHR" "$POP1" "$POP2" --persist \
        -o "$OUTDIR/xpehh/${PAIR_ID}_${CHR}.tsv" 2>&1 || echo "WARN: xpehh $PAIR_ID $CHR failed"

    # Divergence (FAST PATH)
    graphpop divergence "$CHR" 1 300000000 "$POP1" "$POP2" \
        -o "$OUTDIR/divergence/${PAIR_ID}_${CHR}.tsv" 2>&1 || echo "WARN: divergence $PAIR_ID $CHR failed"
done

echo "$(date): Completed pairwise chr=$CHR"
