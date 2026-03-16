#!/usr/bin/env bash
# =============================================================================
# download_3kRG_snpeff.sh
# Download 3K Rice Genomes SnpEff-annotated VCF + supporting files
# =============================================================================
#
# Dataset: 3K RG 29M biallelic SNPs with SnpEff functional annotations
#          against Rice Genome Annotation Project rel 7 (MSU7) gene models
# Source:  IRRI OryzaSNP / 3K RGP data repository
# License: Toronto Statement permissive license (Nature 461, 168-170, 2009)
#          By downloading, you agree to abide by the Toronto Statement.
# Samples: 3,024 rice accessions
# Variants: 29,635,224 biallelic SNPs
# Reference: Nipponbare MSU7 / IRGSP1.0
#
# Usage:
#   chmod +x download_3kRG_snpeff.sh
#   ./download_3kRG_snpeff.sh [DOWNLOAD_DIR]
#
# Default download directory: ./3kRG_data
# =============================================================================

set -euo pipefail

# --- Configuration ---
DOWNLOAD_DIR="${1:-./3kRG_data}"
MAX_RETRIES=3
RETRY_DELAY=10

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# --- File definitions ---
# Primary: SnpEff-annotated VCF (THE key file for GraphPop)
SNPEFF_VCF_URL="https://s3.amazonaws.com/3kricegenome/snpseek-dl/NB_bialSNP_pseudo_canonical_ALL.vcf.gz"
SNPEFF_VCF_FILE="NB_bialSNP_pseudo_canonical_ALL.vcf.gz"
SNPEFF_VCF_SIZE_BYTES=643561279  # ~614 MB compressed

# Supporting files
LICENSE_URL="https://s3.amazonaws.com/3kricegenome/README-3kRG-SNPs-Permissive-License.txt"
SNP_README_URL="https://s3.amazonaws.com/3kricegenome/reduced/README-3kRG-full-SNP-v1.txt"
PHENOTYPE_URL="https://s3-ap-southeast-1.amazonaws.com/oryzasnp-atcg-irri-org/3kRG-phenotypes/3kRG_PhenotypeData_v20170411.xlsx"

# --- Functions ---
log_info()  { echo -e "${GREEN}[INFO]${NC}  $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step()  { echo -e "${CYAN}[STEP]${NC}  $*"; }

check_disk_space() {
    local dir="$1"
    local needed_mb=2000  # ~2 GB to be safe (614 MB compressed + room for decompression)
    local available_mb
    available_mb=$(df -m "$dir" | awk 'NR==2 {print $4}')
    if [ "$available_mb" -lt "$needed_mb" ]; then
        log_error "Insufficient disk space: ${available_mb} MB available, need ${needed_mb} MB"
        log_error "The VCF is ~614 MB compressed; decompressed it will be several GB"
        exit 1
    fi
    log_info "Disk space OK: ${available_mb} MB available"
}

download_with_retry() {
    local url="$1"
    local output="$2"
    local description="$3"
    local attempt=1

    log_step "Downloading: ${description}"
    log_info "  URL: ${url}"
    log_info "  Output: ${output}"

    while [ $attempt -le $MAX_RETRIES ]; do
        if [ -f "$output" ]; then
            # Resume support: use -C - for partial downloads
            if curl -L -C - -o "$output" --progress-bar \
                    --connect-timeout 30 \
                    --retry 2 \
                    "$url"; then
                log_info "  ✓ Download complete (attempt ${attempt})"
                return 0
            fi
        else
            if curl -L -o "$output" --progress-bar \
                    --connect-timeout 30 \
                    --retry 2 \
                    "$url"; then
                log_info "  ✓ Download complete (attempt ${attempt})"
                return 0
            fi
        fi

        log_warn "  Attempt ${attempt}/${MAX_RETRIES} failed. Retrying in ${RETRY_DELAY}s..."
        attempt=$((attempt + 1))
        sleep $RETRY_DELAY
    done

    log_error "  ✗ Failed after ${MAX_RETRIES} attempts: ${description}"
    return 1
}

verify_vcf_gz() {
    local file="$1"
    log_step "Verifying VCF integrity..."

    # Check file size
    local actual_size
    actual_size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo 0)
    if [ "$actual_size" -ne "$SNPEFF_VCF_SIZE_BYTES" ]; then
        log_warn "File size mismatch: expected ${SNPEFF_VCF_SIZE_BYTES}, got ${actual_size}"
        log_warn "File may be incomplete or the server version has changed"
    else
        log_info "  ✓ File size matches expected (${actual_size} bytes)"
    fi

    # Check gzip integrity
    if command -v gzip &>/dev/null; then
        if gzip -t "$file" 2>/dev/null; then
            log_info "  ✓ Gzip integrity check passed"
        else
            log_error "  ✗ Gzip integrity check FAILED — file may be corrupt"
            return 1
        fi
    fi

    # Peek at VCF header to confirm it's a real SnpEff-annotated VCF
    if command -v zcat &>/dev/null; then
        local header_lines
        header_lines=$(zcat "$file" 2>/dev/null | head -50 | grep -c "^##" || true)
        if [ "$header_lines" -gt 5 ]; then
            log_info "  ✓ VCF header detected (${header_lines} header lines in first 50)"
        fi

        # Check for SnpEff annotation
        if zcat "$file" 2>/dev/null | head -100 | grep -qi "snpeff\|ANN=\|EFF=\|CSQ="; then
            log_info "  ✓ SnpEff/functional annotations detected in header"
        else
            log_warn "  SnpEff annotation tags not found in first 100 lines (may appear later)"
        fi
    fi

    return 0
}

print_summary() {
    local dir="$1"
    echo ""
    echo "============================================================"
    echo -e "${GREEN}  Download Complete${NC}"
    echo "============================================================"
    echo ""
    echo "  Directory: ${dir}"
    echo ""
    echo "  Files:"
    ls -lh "$dir"/ 2>/dev/null | grep -v "^total" | awk '{printf "    %-50s %s\n", $NF, $5}'
    echo ""
    echo "============================================================"
    echo "  Quick inspection commands:"
    echo "============================================================"
    echo ""
    echo "  # View VCF header (annotation fields, contigs, samples):"
    echo "  zcat ${dir}/${SNPEFF_VCF_FILE} | head -100"
    echo ""
    echo "  # Count variants per chromosome:"
    echo "  zcat ${dir}/${SNPEFF_VCF_FILE} | grep -v '^#' | cut -f1 | sort | uniq -c"
    echo ""
    echo "  # View first 5 data lines:"
    echo "  zcat ${dir}/${SNPEFF_VCF_FILE} | grep -v '^#' | head -5"
    echo ""
    echo "  # Extract sample names:"
    echo "  zcat ${dir}/${SNPEFF_VCF_FILE} | grep '^#CHROM' | cut -f10- | tr '\t' '\n' | head -20"
    echo ""
    echo "  # For GraphPop import: use bcftools to extract chr1"
    echo "  # bcftools view -r Chr1 ${dir}/${SNPEFF_VCF_FILE} -Oz -o ${dir}/3kRG_snpeff_chr1.vcf.gz"
    echo "  # bcftools index ${dir}/3kRG_snpeff_chr1.vcf.gz"
    echo ""
    echo "============================================================"
    echo -e "  ${CYAN}GraphPop Note:${NC} This VCF contains SnpEff annotations"
    echo "  against MSU7 gene models. During GraphPop import, parse the"
    echo "  ANN/EFF/CSQ field to create:"
    echo "    (:Variant)-[:HAS_CONSEQUENCE]->(:Gene)"
    echo "  edges directly — no separate annotation step needed."
    echo "============================================================"
}

# --- Main ---
echo ""
echo "============================================================"
echo "  3K Rice Genomes — SnpEff VCF Download Script"
echo "============================================================"
echo "  29,635,224 biallelic SNPs × 3,024 samples"
echo "  Functional annotations vs MSU7/IRGSP1.0 gene models"
echo "  Compressed size: ~614 MB"
echo "============================================================"
echo ""

# Check prerequisites
for cmd in curl; do
    if ! command -v "$cmd" &>/dev/null; then
        log_error "Required command not found: ${cmd}"
        exit 1
    fi
done

# Create download directory
mkdir -p "$DOWNLOAD_DIR"
log_info "Download directory: ${DOWNLOAD_DIR}"

# Check disk space
check_disk_space "$DOWNLOAD_DIR"

# 1. License
download_with_retry "$LICENSE_URL" \
    "${DOWNLOAD_DIR}/LICENSE_Toronto_Statement.txt" \
    "Data usage license (Toronto Statement)"

# 2. README
download_with_retry "$SNP_README_URL" \
    "${DOWNLOAD_DIR}/README_29M_SNPs.txt" \
    "29M SNP dataset README"

# 3. Phenotype data (small, useful for Sample nodes)
download_with_retry "$PHENOTYPE_URL" \
    "${DOWNLOAD_DIR}/3kRG_PhenotypeData.xlsx" \
    "Phenotype data (morpho-agronomic traits)"

# 4. THE MAIN FILE — SnpEff-annotated VCF (~614 MB)
echo ""
log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
log_info "Downloading main SnpEff VCF (~614 MB). This may take a while..."
log_info "Download supports resume (-C -) if interrupted."
log_info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

download_with_retry "$SNPEFF_VCF_URL" \
    "${DOWNLOAD_DIR}/${SNPEFF_VCF_FILE}" \
    "SnpEff-annotated VCF (29M biallelic SNPs, ~614 MB)"

# 5. Verify
verify_vcf_gz "${DOWNLOAD_DIR}/${SNPEFF_VCF_FILE}"

# 6. Summary
print_summary "$DOWNLOAD_DIR"

exit 0