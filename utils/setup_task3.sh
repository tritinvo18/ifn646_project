#!/usr/bin/env bash
# ==============================================================================
# utils/setup_task3.sh
#
# Complete environment and data preparation script for Task 3:
# Population-scale variant analysis for CXCL11 (1000 Genomes Project, N=2,504).
#
# Steps performed:
#   1. Environment & tool check (Python, pysam, bcftools, samtools)
#   2. Directory structure verification (data, results)
#   3. Input metadata & BED file verification
#   4. Reference genome verification & preparation (GCF_000001405.26 / GRCh38.fa)
#   5. 1000 Genomes N=2,504 regional VCF preparation (sampled_chr4.norm.vcf.gz)
#
# Usage:
#   bash utils/setup_task3.sh [--all-autosomes] [--force-vcf]
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${PROJECT_ROOT}/data"
VCF_DIR="${DATA_DIR}/vcfs/n2504"
TABLES_DIR="${PROJECT_ROOT}/results/tables/n2504_cxcl11"
FIGURES_DIR="${PROJECT_ROOT}/results/figures/n2504_cxcl11"

# Options
CHROM_MODE="target"
FORCE_VCF=false

for arg in "$@"; do
    case "${arg}" in
        --all-autosomes)
            CHROM_MODE="autosomes"
            ;;
        --force-vcf)
            FORCE_VCF=true
            ;;
        -h|--help)
            echo "Usage: bash scripts/setup_task3.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --all-autosomes   Extract Chr 1-22 VCFs (default: target Chr 4 only)"
            echo "  --force-vcf       Force rebuild of regional VCF(s)"
            echo "  -h, --help        Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: ${arg}" >&2
            exit 1
            ;;
    esac
done

# Ensure conda env is on PATH if present
if [[ -d "/home/jimmy/.conda/envs/ifn646/bin" && ":${PATH}:" != *":/home/jimmy/.conda/envs/ifn646/bin:"* ]]; then
    export PATH="/home/jimmy/.conda/envs/ifn646/bin:${PATH}"
fi

echo "================================================================="
echo " Task 3 Preparation: CXCL11 N=2,504 Population Variant Analysis"
echo " Project root: ${PROJECT_ROOT}"
echo " Chromosome mode: ${CHROM_MODE}"
echo "================================================================="

# ------------------------------------------------------------------------------
# 1. Environment & Dependency Check
# ------------------------------------------------------------------------------
echo ""
echo "[Step 1/5] Checking CLI tools and Python dependencies..."

MISSING_TOOLS=()
for tool in python3 samtools bcftools; do
    if ! command -v "${tool}" &>/dev/null; then
        MISSING_TOOLS+=("${tool}")
    fi
done

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    echo "  ERROR: The following required CLI tools were not found in PATH:" >&2
    for t in "${MISSING_TOOLS[@]}"; do
        echo "    - ${t}" >&2
    done
    echo "  Please activate the conda environment (e.g., 'conda activate ifn646') or install them." >&2
    exit 1
fi
echo "  [OK] CLI tools available: python3, samtools, bcftools"

# Check Python packages
python3 -c "
import importlib.util
required = ['pysam', 'pandas', 'numpy', 'matplotlib', 'seaborn']
missing = [pkg for pkg in required if importlib.util.find_spec(pkg) is None]
if missing:
    print('  ERROR: Missing Python packages:', ', '.join(missing))
    exit(1)
print('  [OK] Python packages verified:', ', '.join(required))
"

# ------------------------------------------------------------------------------
# 2. Directory Structure Verification
# ------------------------------------------------------------------------------
echo ""
echo "[Step 2/5] Verifying directory structure..."
mkdir -p "${DATA_DIR}"
mkdir -p "${VCF_DIR}"
mkdir -p "${TABLES_DIR}"
mkdir -p "${FIGURES_DIR}"
echo "  [OK] Directories ready:"
echo "       - Data:    ${DATA_DIR}"
echo "       - VCFs:    ${VCF_DIR}"
echo "       - Tables:  ${TABLES_DIR}"
echo "       - Figures: ${FIGURES_DIR}"

# ------------------------------------------------------------------------------
# 3. Input Metadata & BED Verification
# ------------------------------------------------------------------------------
echo ""
echo "[Step 3/5] Verifying required metadata and region BED files..."

SAMPLE_META="${DATA_DIR}/n2504/sample_metadata.tsv"
SAMPLE_IDS="${DATA_DIR}/n2504/sample_ids.txt"
BED_FILE="${DATA_DIR}/task3_offtarget_regions.bed"

for req_file in "${SAMPLE_META}" "${SAMPLE_IDS}" "${BED_FILE}"; do
    if [[ ! -s "${req_file}" ]]; then
        echo "  ERROR: Required input file missing or empty: ${req_file}" >&2
        exit 1
    fi
done
echo "  [OK] 1000G Sample metadata: $(wc -l < "${SAMPLE_META}") records"
echo "  [OK] Sample IDs:            $(wc -l < "${SAMPLE_IDS}") samples"
echo "  [OK] Off-target BED:        $(wc -l < "${BED_FILE}") regions"

# ------------------------------------------------------------------------------
# 4. Reference Genome Check (NCBI GCF_000001405.26 / GRCh38.fa)
# ------------------------------------------------------------------------------
echo ""
echo "[Step 4/5] Verifying Reference Genome (GCF_000001405.26 / GRCh38.fa)..."
bash "${PROJECT_ROOT}/utils/download_large_files.sh" --${CHROM_MODE}

# ------------------------------------------------------------------------------
# 5. 1000 Genomes Regional VCFs Verification
# ------------------------------------------------------------------------------
echo ""
echo "[Step 5/5] Verifying N=2,504 Regional VCF files..."

CHR4_VCF="${VCF_DIR}/sampled_chr4.norm.vcf.gz"
NEED_REBUILD=false

if [[ "${FORCE_VCF}" == "true" ]]; then
    NEED_REBUILD=true
elif [[ "${CHROM_MODE}" == "target" ]]; then
    if [[ ! -s "${CHR4_VCF}" || ! -s "${CHR4_VCF}.n2504_regions_complete" ]]; then
        NEED_REBUILD=true
    fi
else
    for c in {1..22}; do
        if [[ ! -s "${VCF_DIR}/sampled_chr${c}.norm.vcf.gz" ]]; then
            NEED_REBUILD=true
            break
        fi
    done
fi

if [[ "${NEED_REBUILD}" == "true" ]]; then
    echo "  Rebuilding regional VCFs via utils/prepare_n2504_vcfs.py..."
    python3 "${PROJECT_ROOT}/utils/prepare_n2504_vcfs.py" \
        --project-root "${PROJECT_ROOT}" \
        --chromosomes "${CHROM_MODE}"
else
    echo "  [OK] Regional VCF(s) present and verified complete."
fi

echo ""
echo "================================================================="
echo " Task 3 setup is complete and ready!"
echo " You can now run the Task 3 notebook:"
echo "   notebooks/03_variant_analysis_cxcl11_n2504_genome_wide-KM.ipynb"
echo "================================================================="
