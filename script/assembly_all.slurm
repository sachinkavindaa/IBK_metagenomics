#!/bin/bash -l

#SBATCH --job-name=IBK_MEGAHIT_single
#SBATCH --array=1-209
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=30:00:00
#SBATCH --partition=guest,batch
#SBATCH --output=/work/samodha/sachin/ShotgunM/IBK_over_1M/logs/megahit_%A_%a.out
#SBATCH --error=/work/samodha/sachin/ShotgunM/IBK_over_1M/logs/megahit_%A_%a.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail


# ============================================================
# Load modules
# ============================================================

module load megahit/1.2
module load bbmap/39.06
module load quast/5.0


# ============================================================
# Base paths
# ============================================================

BASE="/work/samodha/sachin/ShotgunM/IBK_over_1M"

INPUT="${BASE}/04_Bovine_removal"

ASSEMBLY_DIR="${BASE}/05_Single_Assembly"

ALL_CONTIGS_DIR="${BASE}/06_All_Contigs"

CONTIGS_500_DIR="${BASE}/07_Contigs_500bp"

QUAST_BASE="${BASE}/08_MetaQUAST"

LOGDIR="${BASE}/logs"


mkdir -p "${ASSEMBLY_DIR}"
mkdir -p "${ALL_CONTIGS_DIR}"
mkdir -p "${CONTIGS_500_DIR}"
mkdir -p "${QUAST_BASE}"
mkdir -p "${LOGDIR}"


# ============================================================
# Select sample according to SLURM array task
#
# Example:
# ibk98_final_clean_paired.1.fastq
#
# becomes:
# ibk98
# ============================================================

SAMPLE=$(
    find "${INPUT}" \
        -maxdepth 1 \
        -type f \
        -name "*_final_clean_paired.1.fastq" \
        -printf "%f\n" \
    | sed 's/_final_clean_paired\.1\.fastq$//' \
    | sort -V \
    | sed -n "${SLURM_ARRAY_TASK_ID}p"
)


# ============================================================
# Make sure sample exists
# ============================================================

if [[ -z "${SAMPLE}" ]]; then

    echo "ERROR: No sample found for SLURM array task:"
    echo "${SLURM_ARRAY_TASK_ID}"

    exit 1
fi


echo "============================================================"
echo "SLURM Job ID      : ${SLURM_JOB_ID}"
echo "Array Task        : ${SLURM_ARRAY_TASK_ID}"
echo "Sample            : ${SAMPLE}"
echo "============================================================"


# ============================================================
# Input FASTQ files
# ============================================================

R1="${INPUT}/${SAMPLE}_final_clean_paired.1.fastq"

R2="${INPUT}/${SAMPLE}_final_clean_paired.2.fastq"

SGL="${INPUT}/${SAMPLE}_final_clean_single.fastq"


echo ""
echo "Input files:"
echo "R1  = ${R1}"
echo "R2  = ${R2}"
echo "SGL = ${SGL}"
echo ""


# ============================================================
# Check paired FASTQ files
# ============================================================

if [[ ! -s "${R1}" ]]; then

    echo "ERROR: Missing or empty R1:"
    echo "${R1}"

    exit 1
fi


if [[ ! -s "${R2}" ]]; then

    echo "ERROR: Missing or empty R2:"
    echo "${R2}"

    exit 1
fi


# ============================================================
# Sample-specific assembly directory
# ============================================================

SAMPLE_DIR="${ASSEMBLY_DIR}/${SAMPLE}"

ASM="${SAMPLE_DIR}/megahit_out"

mkdir -p "${SAMPLE_DIR}"


# ============================================================
# Final output files
# ============================================================

ALL_CONTIGS="${ALL_CONTIGS_DIR}/${SAMPLE}_all_contigs.fa"

CONTIGS_500="${CONTIGS_500_DIR}/${SAMPLE}_500_contigs.fa"

QUASTDIR="${QUAST_BASE}/${SAMPLE}"


# ============================================================
# Run MEGAHIT
# ============================================================

echo ""
echo "============================================================"
echo "Starting MEGAHIT"
echo "Sample: ${SAMPLE}"
echo "============================================================"
echo ""


if [[ -s "${ASM}/final.contigs.fa" ]]; then

    echo "Existing assembly found:"
    echo "${ASM}/final.contigs.fa"

    echo ""
    echo "Skipping MEGAHIT."

else

    # Remove incomplete directory from failed previous run
    if [[ -d "${ASM}" ]]; then

        echo "Removing incomplete MEGAHIT directory:"
        echo "${ASM}"

        rm -rf "${ASM}"

    fi


    # --------------------------------------------------------
    # Paired reads + singleton reads
    # --------------------------------------------------------

    if [[ -s "${SGL}" ]]; then

        echo "Singleton reads found."
        echo "Running paired + singleton assembly."

        megahit \
            -1 "${R1}" \
            -2 "${R2}" \
            -r "${SGL}" \
            -o "${ASM}" \
            -t "${SLURM_CPUS_PER_TASK}" \
            -m 0.9 \
            --min-contig-len 200

    else

        # ----------------------------------------------------
        # Paired reads only
        # ----------------------------------------------------

        echo "WARNING: Singleton file not found or empty."
        echo "Running paired reads only."

        megahit \
            -1 "${R1}" \
            -2 "${R2}" \
            -o "${ASM}" \
            -t "${SLURM_CPUS_PER_TASK}" \
            -m 0.9 \
            --min-contig-len 200

    fi

fi


# ============================================================
# Check MEGAHIT result
# ============================================================

if [[ ! -s "${ASM}/final.contigs.fa" ]]; then

    echo ""
    echo "ERROR: MEGAHIT failed."
    echo "final.contigs.fa was not produced."
    echo ""
    echo "Expected:"
    echo "${ASM}/final.contigs.fa"

    exit 1
fi


echo ""
echo "MEGAHIT assembly successfully found."


# ============================================================
# Save ALL MEGAHIT contigs
# ============================================================

echo ""
echo "============================================================"
echo "Saving ALL contigs"
echo "============================================================"
echo ""


if [[ ! -s "${ALL_CONTIGS}" ]]; then

    cp "${ASM}/final.contigs.fa" "${ALL_CONTIGS}"

else

    echo "All-contig file already exists:"
    echo "${ALL_CONTIGS}"

    echo "Skipping copy."

fi


# ============================================================
# Create >=500 bp contig file
# ============================================================

echo ""
echo "============================================================"
echo "Filtering contigs >=500 bp"
echo "============================================================"
echo ""


if [[ ! -s "${CONTIGS_500}" ]]; then

    reformat.sh \
        in="${ASM}/final.contigs.fa" \
        out="${CONTIGS_500}" \
        minlength=500

else

    echo "500-bp filtered file already exists:"
    echo "${CONTIGS_500}"

    echo "Skipping filtering."

fi


# ============================================================
# Check filtered contig file
# ============================================================

if [[ ! -s "${CONTIGS_500}" ]]; then

    echo ""
    echo "ERROR: No >=500 bp contigs were produced."
    echo "${CONTIGS_500}"

    exit 1
fi


# ============================================================
# Count contigs
# ============================================================

ALL_COUNT=$(grep -c "^>" "${ALL_CONTIGS}" || true)

COUNT_500=$(grep -c "^>" "${CONTIGS_500}" || true)


echo ""
echo "============================================================"
echo "Contig summary"
echo "============================================================"

echo ""
echo "Sample: ${SAMPLE}"
echo "All contigs:       ${ALL_COUNT}"
echo "Contigs >=500 bp:  ${COUNT_500}"


# ============================================================
# Run metaQUAST on >=500 bp contigs
# ============================================================

echo ""
echo "============================================================"
echo "Running metaQUAST"
echo "============================================================"
echo ""


if [[ -f "${QUASTDIR}/report.txt" ]]; then

    echo "metaQUAST report already exists:"
    echo "${QUASTDIR}/report.txt"

    echo "Skipping metaQUAST."

else

    if [[ -d "${QUASTDIR}" ]]; then
        rm -rf "${QUASTDIR}"
    fi


    metaquast.py \
        "${CONTIGS_500}" \
        -o "${QUASTDIR}" \
        --threads "${SLURM_CPUS_PER_TASK}"

fi


# ============================================================
# Finished
# ============================================================

echo ""
echo "============================================================"
echo "COMPLETED SUCCESSFULLY"
echo "============================================================"

echo ""
echo "Sample:"
echo "${SAMPLE}"

echo ""
echo "Original MEGAHIT assembly:"
echo "${ASM}/final.contigs.fa"

echo ""
echo "All contigs copied to:"
echo "${ALL_CONTIGS}"

echo ""
echo "Contigs >=500 bp:"
echo "${CONTIGS_500}"

echo ""
echo "metaQUAST:"
echo "${QUASTDIR}"

echo ""
echo "All contigs count:"
echo "${ALL_COUNT}"

echo ""
echo ">=500 bp contig count:"
echo "${COUNT_500}"

echo ""
echo "============================================================"
