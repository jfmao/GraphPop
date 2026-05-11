#!/usr/bin/env bash
# Stage 2 — Genotype GRM via PLINK 2.0.
# Equivalent in Cypher: `graphpop.kinship.branch_grm` directly.

set -euo pipefail

BFILE_PREFIX="${1:?usage: $0 <bfile_prefix>}"
OUT_PREFIX="${2:?usage: $0 <bfile_prefix> <out_prefix>}"

plink2 \
    --bfile "${BFILE_PREFIX}" \
    --make-grm-bin \
    --out "${OUT_PREFIX}_grm"
