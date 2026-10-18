#!/usr/bin/env bash
# ==============================================================================
# utils/download_large_files.sh
#
# Download and prepare large reference and variant files not tracked in git:
# 1. Human Reference Genome GRCh38 (FASTA + .fai index) from NCBI
# 2. 1000 Genomes N=2,504 regional VCFs (Chr 4 - CXCL11) from EBI FTP
#
# Usage:
#   bash utils/download_large_files.sh [--all-autosomes]
#
# Options:
#   --all-autosomes   Extract Chr 1-22 VCFs instead of target chromosome (4)
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"
VCF_DIR="${DATA_DIR}/vcfs/n2504"

mkdir -p "${DATA_DIR}"
mkdir -p "${VCF_DIR}"

CHROM_MODE="target"
if [[ "${1:-}" == "--all-autosomes" ]]; then
    CHROM_MODE="autosomes"
fi

echo "================================================================="
echo " IFN646 Project: Large File Download & Preparation Script"
echo " Project root: ${PROJECT_ROOT}"
echo " Data directory: ${DATA_DIR}"
echo "================================================================="

# Ensure conda environment binaries are in PATH if available
if [[ -d "/home/jimmy/.conda/envs/ifn646/bin" && ":${PATH}:" != *":/home/jimmy/.conda/envs/ifn646/bin:"* ]]; then
    export PATH="/home/jimmy/.conda/envs/ifn646/bin:${PATH}"
fi

# ------------------------------------------------------------------------------
# 1. GRCh38 Reference Genome (NCBI GCF_000001405.26)
# ------------------------------------------------------------------------------
FASTA_PATH="${DATA_DIR}/GRCh38.fa"
FAI_PATH="${DATA_DIR}/GRCh38.fa.fai"
ARCHIVE_PATH="${DATA_DIR}/GCF_000001405.26_GRCh38_genomic.fna.gz"
FNA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GCF_000001405.26_GRCh38_genomic.fna.gz"

echo ""
echo "[1/2] Checking GRCh38 Reference Genome..."

# 1a. Compressed archive: GCF_000001405.26_GRCh38_genomic.fna.gz
if [[ -s "${ARCHIVE_PATH}" ]]; then
    echo "  [SKIP] Compressed archive ${ARCHIVE_PATH##*/} already present ($(du -h "${ARCHIVE_PATH}" | cut -f1))."
else
    echo "  Downloading ${ARCHIVE_PATH##*/} from NCBI..."
    if command -v curl &>/dev/null; then
        curl -C - -L --retry 3 -o "${ARCHIVE_PATH}" "${FNA_URL}"
    elif command -v wget &>/dev/null; then
        wget -c -O "${ARCHIVE_PATH}" "${FNA_URL}"
    else
        python3 -c "
import urllib.request
print('Downloading via urllib...')
urllib.request.urlretrieve('${FNA_URL}', '${ARCHIVE_PATH}')
"
    fi
fi

# 1b. Decompressed FASTA: GRCh38.fa
if [[ -s "${FASTA_PATH}" ]]; then
    echo "  [SKIP] GRCh38.fa already present ($(du -h "${FASTA_PATH}" | cut -f1))."
else
    echo "  Decompressing ${ARCHIVE_PATH##*/} -> ${FASTA_PATH}..."
    gzip -dc "${ARCHIVE_PATH}" > "${FASTA_PATH}.tmp"
    mv "${FASTA_PATH}.tmp" "${FASTA_PATH}"
    echo "  Decompression complete: ${FASTA_PATH}"
fi

# 1c. FASTA index: GRCh38.fa.fai
if [[ -s "${FAI_PATH}" ]]; then
    echo "  [SKIP] GRCh38.fa.fai index already present."
else
    echo "  Building FASTA index using samtools..."
    if command -v samtools &>/dev/null; then
        samtools faidx "${FASTA_PATH}"
    else
        echo "  WARNING: samtools not found in PATH; skipping index generation."
    fi
fi

# ------------------------------------------------------------------------------
# 2. 1000 Genomes N=2,504 Regional VCFs
# ------------------------------------------------------------------------------
echo ""
echo "[2/2] Checking 1000 Genomes N=2,504 Regional VCFs (${CHROM_MODE})..."

CHR4_VCF="${VCF_DIR}/sampled_chr4.norm.vcf.gz"

NEED_PREPARE=false

if [[ "${CHROM_MODE}" == "target" ]]; then
    if [[ -s "${CHR4_VCF}" && -s "${CHR4_VCF}.n2504_regions_complete" ]]; then
        echo "  [SKIP] Target regional VCF (Chr 4 - CXCL11) already present:"
        echo "         - ${CHR4_VCF}"
    else
        NEED_PREPARE=true
    fi
else
    # autosomes mode: check if all 22 chromosomes exist
    for c in {1..22}; do
        if [[ ! -s "${VCF_DIR}/sampled_chr${c}.norm.vcf.gz" ]]; then
            NEED_PREPARE=true
            break
        fi
    done
    if [[ "${NEED_PREPARE}" == "false" ]]; then
        echo "  [SKIP] All 22 autosomal VCFs already present in ${VCF_DIR}."
    fi
fi

if [[ "${NEED_PREPARE}" == "true" ]]; then
    echo "  Running utils/prepare_n2504_vcfs.py (--chromosomes ${CHROM_MODE})..."
    python3 "${PROJECT_ROOT}/utils/prepare_n2504_vcfs.py" \
        --project-root "${PROJECT_ROOT}" \
        --chromosomes "${CHROM_MODE}"
fi

echo ""
echo "================================================================="
echo " Large files verification complete!"
echo "================================================================="
