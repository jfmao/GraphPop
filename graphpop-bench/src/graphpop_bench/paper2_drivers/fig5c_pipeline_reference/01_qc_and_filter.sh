#!/usr/bin/env bash
# Stage 1 — QC + filter the cohort VCF down to common variants.
# Equivalent in Cypher: the implicit Variant-node filtering done
# by the `restrict_to_pathway` predicate at query time.

set -euo pipefail

INPUT_VCF="${1:?usage: $0 <input.vcf.gz>}"
OUT_PREFIX="${2:?usage: $0 <input.vcf.gz> <out_prefix>}"

plink2 \
    --vcf "${INPUT_VCF}" \
    --maf 0.01 \
    --geno 0.05 \
    --hwe 1e-6 \
    --make-bed \
    --out "${OUT_PREFIX}_qc"

# Sanity: log the SNP count.
N_SNPS=$(wc -l < "${OUT_PREFIX}_qc.bim")
echo "[01_qc] SNPs after QC: ${N_SNPS}"
