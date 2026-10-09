#!/bin/bash -l
#SBATCH --job-name=IBK_Kraken2_C01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=batch,guest
#SBATCH --output=/work/samodha/sachin/ShotgunM/IBK_over_1M/MEGAHIT_global_coassembly/Downstream_500/LOGS/kraken2_conf01_%j.out
#SBATCH --error=/work/samodha/sachin/ShotgunM/IBK_over_1M/MEGAHIT_global_coassembly/Downstream_500/LOGS/kraken2_conf01_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ============================================================
# Load Kraken2
# ============================================================

module load kraken2

# ============================================================
# Paths
# ============================================================

BASE="/work/samodha/sachin/ShotgunM/IBK_over_1M/MEGAHIT_global_coassembly"

# Global coassembly contigs >=500 bp
CONTIGS="${BASE}/IBK_global_contigs_1000.fa"

# Separate directory for confidence = 0.1
OUT="${BASE}/Downstream_1000/01_Kraken2_conf05"

LOGDIR="${BASE}/Downstream_1000/LOGS"

# Use HCC PlusPF Kraken2 database
KRAKEN_DB="${KRAKEN2_PLUSPF_DB}"

# ============================================================
# Create output directories
# ============================================================

mkdir -p "$OUT" "$LOGDIR"

# ============================================================
# Check input contigs
# ============================================================

if [[ ! -f "$CONTIGS" ]]; then
    echo "ERROR: Input contig file not found:"
    echo "$CONTIGS"
    exit 1
fi

# ============================================================
# Check Kraken2 database
# ============================================================

if [[ ! -d "$KRAKEN_DB" ]]; then
    echo "ERROR: Kraken2 database not found:"
    echo "$KRAKEN_DB"
    exit 1
fi

# ============================================================
# Print run information
# ============================================================

echo "============================================================"
echo "IBK Global Coassembly - Kraken2"
echo "============================================================"
echo ""
echo "Start time : $(date)"
echo "Node       : $(hostname)"
echo "Job ID     : ${SLURM_JOB_ID:-NA}"
echo "Threads    : ${SLURM_CPUS_PER_TASK}"
echo ""
echo "Input:"
echo "$CONTIGS"
echo ""
echo "Database:"
echo "$KRAKEN_DB"
echo ""
echo "Confidence:"
echo "0.1"
echo ""
echo "Output directory:"
echo "$OUT"
echo ""
echo "============================================================"

# ============================================================
# Run Kraken2
# ============================================================

echo ""
echo "Starting Kraken2..."
echo ""

kraken2 \
    --db "$KRAKEN_DB" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --confidence 0.05 \
    --use-names \
    --report "$OUT/kraken2_report_conf01.txt" \
    --output "$OUT/kraken2_output_conf01.txt" \
    "$CONTIGS"

# ============================================================
# Finished
# ============================================================

echo ""
echo "============================================================"
echo "Kraken2 completed successfully."
echo "============================================================"
echo ""
echo "End time: $(date)"
echo ""
echo "Results:"
echo ""
echo "Report:"
echo "$OUT/kraken2_report_conf01.txt"
echo ""
echo "Per-contig classification:"
echo "$OUT/kraken2_output_conf01.txt"
echo ""
echo "============================================================"
