#!/usr/bin/env bash
# prepare-vcf.sh — Decompose multi-allelic sites before GraphPop import.
#
# GraphPop's stored procedures assume biallelic variants. Rather than
# discarding multi-allelic sites entirely, we decompose them into
# separate biallelic records using bcftools norm.  This preserves all
# variant information while keeping downstream statistics valid and
# comparable to scikit-allel, PLINK2, vcftools, etc.
#
# Usage:
#   bash scripts/prepare-vcf.sh INPUT.vcf.gz OUTPUT.vcf.gz [REFERENCE.fa]
#
# Requirements:
#   - bcftools (https://samtools.github.io/bcftools/)
#   - tabix / htslib for indexing
#
# What it does:
#   1. Decompose multi-allelic sites into biallelic records (-m-)
#   2. Optionally left-align and normalise indels (if reference provided)
#   3. Remove exact duplicate records that decomposition can create
#   4. Re-index the output
#
# Example:
#   bash scripts/prepare-vcf.sh \
#       data/raw/1000g/chr22.vcf.gz \
#       data/raw/1000g/chr22.biallelic.vcf.gz \
#       data/raw/reference/GRCh38.fa

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 INPUT.vcf.gz OUTPUT.vcf.gz [REFERENCE.fa]" >&2
    exit 1
fi

INPUT="$1"
OUTPUT="$2"
REFERENCE="${3:-}"

# Check dependencies
if ! command -v bcftools &>/dev/null; then
    echo "ERROR: bcftools not found. Install with: conda install -c bioconda bcftools" >&2
    exit 1
fi

if ! command -v tabix &>/dev/null; then
    echo "ERROR: tabix not found. Install with: conda install -c bioconda htslib" >&2
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file not found: $INPUT" >&2
    exit 1
fi

echo "=== GraphPop VCF Preprocessing ==="
echo "Input:  $INPUT"
echo "Output: $OUTPUT"

# Count multi-allelic sites before decomposition
N_MULTI=$(bcftools view -H "$INPUT" | awk -F'\t' '{n=split($5,a,","); if(n>1) c++} END{print c+0}')
N_TOTAL=$(bcftools view -H "$INPUT" | wc -l)
echo "Total variants:        $N_TOTAL"
echo "Multi-allelic sites:   $N_MULTI ($(awk "BEGIN{printf \"%.1f\", $N_MULTI/$N_TOTAL*100}")%)"

# Build bcftools norm command
NORM_ARGS=("-m-" "-Oz" "-o" "$OUTPUT")

if [ -n "$REFERENCE" ]; then
    if [ ! -f "$REFERENCE" ]; then
        echo "ERROR: Reference file not found: $REFERENCE" >&2
        exit 1
    fi
    echo "Reference: $REFERENCE (left-aligning indels)"
    NORM_ARGS+=("-f" "$REFERENCE")
fi

echo ""
echo "Running: bcftools norm ${NORM_ARGS[*]} $INPUT"
bcftools norm "${NORM_ARGS[@]}" "$INPUT"

# Remove exact duplicates that decomposition can create
echo "Removing duplicate records..."
bcftools norm --rm-dup exact -Oz -o "${OUTPUT}.tmp" "$OUTPUT"
mv "${OUTPUT}.tmp" "$OUTPUT"

# Index
echo "Indexing output..."
tabix -p vcf "$OUTPUT"

# Report
N_OUTPUT=$(bcftools view -H "$OUTPUT" | wc -l)
N_NEW=$((N_OUTPUT - N_TOTAL))
echo ""
echo "=== Done ==="
echo "Output variants: $N_OUTPUT (+$N_NEW from decomposition)"
echo "Output file:     $OUTPUT"
echo "Index file:      ${OUTPUT}.tbi"
