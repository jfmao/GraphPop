#!/usr/bin/env bash
# Stage 5 — Local-ancestry calls via RFMix (Maples 2013).
# Equivalent in Cypher: `ancestry: "EUR"` predicate against the
# `:Sample`-attached painting that GraphPop ingests once.

set -euo pipefail

QUERY_VCF="${1:?usage: $0 <query.vcf.gz>}"
REFERENCE_VCF="${2:?usage: $0 <query.vcf.gz> <reference.vcf.gz>}"
SAMPLE_MAP="${3:?usage: $0 <query> <reference> <sample_map.tsv>}"
GENETIC_MAP="${4:?usage: $0 <query> <reference> <sample_map> <genetic_map.tsv>}"
OUT_PREFIX="${5:?usage: $0 ... <out_prefix>}"

rfmix \
    -f "${QUERY_VCF}" \
    -r "${REFERENCE_VCF}" \
    -m "${SAMPLE_MAP}" \
    -g "${GENETIC_MAP}" \
    -o "${OUT_PREFIX}" \
    --chromosome=22
