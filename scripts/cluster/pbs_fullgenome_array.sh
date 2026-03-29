#!/bin/bash
#PBS -N graphpop-fullgenome
#PBS -J 1-12
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l walltime=4:00:00
#PBS -o graphpop_fullgenome.log
#PBS -j oe

# GraphPop — Full-genome per-population analysis (PBS array job)
#
# Usage:
#   qsub -v GRAPHPOP_URI=bolt://db-node:7687,GRAPHPOP_PASSWORD=mypass pbs_fullgenome_array.sh
#   qsub -v GRAPHPOP_URI=bolt://db-node:7687,GRAPHPOP_POPULATIONS=GJ-tmp,GJ-trp pbs_fullgenome_array.sh

set -euo pipefail

module load java/21 2>/dev/null || true
source activate graphpop 2>/dev/null || conda activate graphpop 2>/dev/null || true

CHR_PREFIX="${CHR_PREFIX:-Chr}"
CHR="${CHR_PREFIX}${PBS_ARRAY_INDEX}"
OUTDIR="${GRAPHPOP_OUTDIR:-results}"

POPULATIONS="${GRAPHPOP_POPULATIONS:-GJ-tmp,GJ-trp,XI-1A,XI-1B,XI-2,XI-3,XI-adm,cA-Aus,cB-Bas,GJ-adm,GJ-sbtrp,admix}"

mkdir -p "$OUTDIR"/{diversity,ihs,nsl,garud_h,roh}

IFS=',' read -ra POPS <<< "$POPULATIONS"
for POP in "${POPS[@]}"; do
    POP=$(echo "$POP" | tr -d '[:space:]')
    echo "$(date): $CHR / $POP"
    graphpop diversity "$CHR" 1 300000000 "$POP" -o "$OUTDIR/diversity/${POP}_${CHR}.tsv" 2>&1 || true
    graphpop ihs "$CHR" "$POP" --persist -o "$OUTDIR/ihs/${POP}_${CHR}.tsv" 2>&1 || true
    graphpop nsl "$CHR" "$POP" --persist -o "$OUTDIR/nsl/${POP}_${CHR}.tsv" 2>&1 || true
    graphpop garud-h "$CHR" "$POP" 100000 50000 -o "$OUTDIR/garud_h/${POP}_${CHR}.tsv" 2>&1 || true
    graphpop roh "$CHR" "$POP" -o "$OUTDIR/roh/${POP}_${CHR}.tsv" 2>&1 || true
done

echo "$(date): Completed $CHR"
