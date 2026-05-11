#!/usr/bin/env bash
# Stage 3 — Robust kinship via KING (Manichaikul 2010).
# Equivalent in Cypher: composed predicate
# `graphpop.relate.classify` on the branch_grm output.

set -euo pipefail

BFILE_PREFIX="${1:?usage: $0 <bfile_prefix>}"
OUT_PREFIX="${2:?usage: $0 <bfile_prefix> <out_prefix>}"

king \
    -b "${BFILE_PREFIX}.bed" \
    --kinship \
    --degree 3 \
    --prefix "${OUT_PREFIX}_kin"

echo "[03_kin] kinship pair table at ${OUT_PREFIX}_kin.kin0"
