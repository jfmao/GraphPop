#!/usr/bin/env bash
# Stage 4 — Global ancestry via ADMIXTURE (Alexander 2009).
# Equivalent in Cypher: `:HAS_ANCESTRY` relationship traversal
# inside the same query (no second tool).

set -euo pipefail

BFILE_PREFIX="${1:?usage: $0 <bfile_prefix>}"
K="${2:-5}"

admixture \
    "${BFILE_PREFIX}.bed" \
    "${K}" \
    --cv

# Output goes to <bfile_prefix>.K.{Q,P} in the cwd.
echo "[04_admixture] Q matrix: ${BFILE_PREFIX}.${K}.Q"
