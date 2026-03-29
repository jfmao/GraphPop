#!/bin/bash
#SBATCH --job-name=graphpop-csv
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=graphpop_csv_%j.log

# GraphPop — Step 1: Generate CSV files from VCF (no Neo4j needed)
#
# This runs on any compute node. Embarrassingly parallel by chromosome.
# Produces neo4j-admin import CSVs for variants, samples, populations,
# annotations, and packed genotypes.
#
# Usage:
#   sbatch slurm_prepare_csv.sh /path/to/data.vcf.gz /path/to/pops.tsv
#   sbatch slurm_prepare_csv.sh data.vcf.gz pops.tsv /scratch/$USER/csv_out

set -euo pipefail

VCF_INPUT="${1:?Usage: sbatch slurm_prepare_csv.sh <vcf> <population_map> [output_dir]}"
POP_MAP="${2:?Usage: sbatch slurm_prepare_csv.sh <vcf> <population_map> [output_dir]}"
OUTPUT_DIR="${3:-/scratch/$USER/graphpop_csv}"

echo "=== GraphPop prepare-csv ==="
echo "Input VCF:      $VCF_INPUT"
echo "Population map: $POP_MAP"
echo "Output dir:     $OUTPUT_DIR"
echo "Threads:        $SLURM_CPUS_PER_TASK"
echo "Started:        $(date)"

# Activate Python environment (adjust to your cluster setup)
# module load java/21
# source ~/graphpop-env/bin/activate

python -m graphpop_import.csv_emitter \
    --vcf "$VCF_INPUT" \
    --population-map "$POP_MAP" \
    --output-dir "$OUTPUT_DIR"

echo "=== Completed: $(date) ==="
echo "CSV files written to: $OUTPUT_DIR"
ls -lh "$OUTPUT_DIR"
