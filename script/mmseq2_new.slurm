#!/bin/bash -l
#SBATCH --job-name=IBK_MMseqs_fast
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=180G
#SBATCH --time=100:00:00
#SBATCH --partition=batch,guest
#SBATCH --output=/work/samodha/sachin/ShotgunM/IBK_over_1M/MEGAHIT_global_coassembly/Downstream_1000/LOGS/mmseqs_fast_%j.out
#SBATCH --error=/work/samodha/sachin/ShotgunM/IBK_over_1M/MEGAHIT_global_coassembly/Downstream_1000/LOGS/mmseqs_fast_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

module load mmseqs2/18.8cc5c


# ============================================================
# PATHS
# ============================================================

BASE="/work/samodha/sachin/ShotgunM/IBK_over_1M/MEGAHIT_global_coassembly"

CONTIGS="${BASE}/IBK_global_contigs_1000.fa"

DB="/work/samodha/sachin/databases/mmseqs2/UniRef90"

# IMPORTANT:
# Completely separate folder from the currently running job
OUTDIR="${BASE}/Downstream_1000/01_MMseqs2_taxonomy_s5.7"

TMPDIR="${OUTDIR}/tmp"

PREFIX="${OUTDIR}/IBK_global_1000_s5.7"

mkdir -p "${OUTDIR}"
mkdir -p "${TMPDIR}"
mkdir -p "${BASE}/Downstream_1000/LOGS"


# ============================================================
# INFORMATION
# ============================================================

echo "============================================================"
echo "IBK MMseqs2 FAST TAXONOMY RUN"
echo "============================================================"
echo "Date:        $(date)"
echo "Node:        $(hostname)"
echo "Contigs:     ${CONTIGS}"
echo "Database:    ${DB}"
echo "Output:      ${OUTDIR}"
echo "Temporary:   ${TMPDIR}"
echo "Threads:     ${SLURM_CPUS_PER_TASK}"
echo "Sensitivity: 5.7"
echo "============================================================"


# ============================================================
# CHECK INPUT
# ============================================================

if [[ ! -s "${CONTIGS}" ]]; then
    echo "ERROR: Contig file missing or empty:"
    echo "${CONTIGS}"
    exit 1
fi

echo
echo "Contig file:"
ls -lh "${CONTIGS}"


# ============================================================
# CHECK DATABASE
# ============================================================

if [[ ! -e "${DB}.dbtype" ]]; then
    echo "ERROR: MMseqs2 database not found:"
    echo "${DB}.dbtype"
    exit 1
fi

if [[ ! -e "${DB}_mapping" ]]; then
    echo "ERROR: taxonomy mapping file not found:"
    echo "${DB}_mapping"
    exit 1
fi

if [[ ! -e "${DB}_taxonomy" ]]; then
    echo "ERROR: taxonomy database not found:"
    echo "${DB}_taxonomy"
    exit 1
fi

echo
echo "MMseqs2 database found:"
mmseqs dbtype "${DB}"

echo
echo "Database files:"
ls -lh "${DB}.dbtype"
ls -lh "${DB}_mapping"
ls -lh "${DB}_taxonomy"


# ============================================================
# RUN MMSEQS2
# ============================================================

echo
echo "============================================================"
echo "Starting MMseqs2 easy-taxonomy with sensitivity 5.7"
echo "============================================================"
echo

mmseqs easy-taxonomy \
    "${CONTIGS}" \
    "${DB}" \
    "${PREFIX}" \
    "${TMPDIR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --tax-lineage 2 \
    -s 5.7


# ============================================================
# FINISHED
# ============================================================

echo
echo "============================================================"
echo "MMseqs2 FAST TAXONOMY FINISHED"
echo "============================================================"
echo "Finished: $(date)"
echo

echo "Output directory:"
echo "${OUTDIR}"

echo
echo "Output files:"
ls -lh "${OUTDIR}"

echo
echo "DONE"
