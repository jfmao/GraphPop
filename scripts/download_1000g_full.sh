#!/usr/bin/env bash
# Download 1000 Genomes NYGC 30x phased VCFs (chr1-22) + tabix indices.
#
# Source: 1000G high-coverage 2022 phased release
# File pattern: 1kGP_high_coverage_Illumina.chrN.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
# Total: ~30 GB for chr1-22
#
# Usage:
#   bash scripts/download_1000g_full.sh [--dry-run]

set -euo pipefail

BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
OUT_DIR="data/raw/1000g/vcf"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()  { echo -e "${GREEN}[INFO]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
log_step()  { echo -e "${BLUE}[STEP]${NC} $*"; }

DRY_RUN=false
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=true
    log_warn "DRY RUN — will only print URLs, not download"
fi

mkdir -p "$OUT_DIR"

download_with_retry() {
    local url="$1"
    local dest="$2"
    local max_retries=3
    local retry_delay=10

    for attempt in $(seq 1 $max_retries); do
        if [[ -f "$dest" ]]; then
            # Resume partial download
            if wget -c -q --show-progress -O "$dest" "$url"; then
                return 0
            fi
        else
            if wget -q --show-progress -O "$dest" "$url"; then
                return 0
            fi
        fi

        if [[ $attempt -lt $max_retries ]]; then
            log_warn "  Attempt $attempt failed, retrying in ${retry_delay}s..."
            sleep $retry_delay
        fi
    done

    log_error "  Failed after $max_retries attempts: $url"
    return 1
}

# ---------------------------------------------------------------------------
# Phase 1: Download 1000G NYGC 30x VCFs (chr1-22)
# ---------------------------------------------------------------------------

log_step "Phase 1: Downloading 1000G NYGC 30x phased VCFs (chr1-22)"
echo ""

TOTAL_FILES=0
SKIPPED=0
DOWNLOADED=0
FAILED=0

for CHR_NUM in $(seq 1 22); do
    VCF_FILE="1kGP_high_coverage_Illumina.chr${CHR_NUM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    TBI_FILE="${VCF_FILE}.tbi"
    VCF_URL="${BASE_URL}/${VCF_FILE}"
    TBI_URL="${BASE_URL}/${TBI_FILE}"

    for FILE in "$VCF_FILE" "$TBI_FILE"; do
        URL="${BASE_URL}/${FILE}"
        DEST="${OUT_DIR}/${FILE}"
        TOTAL_FILES=$((TOTAL_FILES + 1))

        if $DRY_RUN; then
            echo "  $URL -> $DEST"
            continue
        fi

        if [[ -f "$DEST" ]] && [[ -s "$DEST" ]]; then
            log_info "  chr${CHR_NUM}: ${FILE##*.} exists, skipping"
            SKIPPED=$((SKIPPED + 1))
            continue
        fi

        log_info "  chr${CHR_NUM}: downloading $FILE..."
        if download_with_retry "$URL" "$DEST"; then
            DOWNLOADED=$((DOWNLOADED + 1))
        else
            FAILED=$((FAILED + 1))
        fi
    done
done

echo ""
log_step "Summary: ${TOTAL_FILES} total, ${DOWNLOADED} downloaded, ${SKIPPED} skipped, ${FAILED} failed"

# ---------------------------------------------------------------------------
# Phase 2: Verify downloads
# ---------------------------------------------------------------------------

if ! $DRY_RUN; then
    echo ""
    log_step "Phase 2: Verifying downloads"

    MISSING=0
    for CHR_NUM in $(seq 1 22); do
        VCF="${OUT_DIR}/1kGP_high_coverage_Illumina.chr${CHR_NUM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        TBI="${VCF}.tbi"
        if [[ ! -f "$VCF" ]]; then
            log_error "  Missing: chr${CHR_NUM} VCF"
            MISSING=$((MISSING + 1))
        fi
        if [[ ! -f "$TBI" ]]; then
            log_error "  Missing: chr${CHR_NUM} TBI"
            MISSING=$((MISSING + 1))
        fi
    done

    if [[ $MISSING -eq 0 ]]; then
        log_info "All 22 VCF + TBI files present"
    else
        log_error "$MISSING files missing"
    fi

    echo ""
    log_step "Disk usage:"
    du -sh "$OUT_DIR"

    echo ""
    log_step "File list:"
    ls -lhS "$OUT_DIR"/*.vcf.gz 2>/dev/null | awk '{printf "  %-80s %s\n", $NF, $5}'
fi

echo ""
log_step "Next steps:"
echo "  1. Extract ancestral FASTAs: cd data/raw/ensembl_ancestor && tar xzf homo_sapiens_ancestor_GRCh38.tar.gz"
echo "  2. Run import:  conda run -n graphevo python scripts/import_full_1000g.py"
