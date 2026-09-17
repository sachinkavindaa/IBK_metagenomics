#!/bin/bash -l
#SBATCH --job-name=quast_master
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --partition=batch,guest
#SBATCH --output=quast_master_%j.out
#SBATCH --error=quast_master_%j.err
#SBATCH --mail-user=dchuleta2@huskers.unl.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

BASE="/work/samodha/sachin/ShotgunM/IBK_over_1M/MEGAHIT_single/Assembly_singlepass"

OUT="${BASE}/QUAST_master_summary.tsv"

# Header
printf "Sample\tAll_Num_Contigs\tAll_Total_Length\tAll_Largest_Contig\tAll_N50\t500_Num_Contigs\t500_Total_Length\t500_Largest_Contig\t500_N50\t1000_Num_Contigs\t1000_Total_Length\t1000_Largest_Contig\t1000_N50\n" > "$OUT"

for SAMPLE_DIR in "$BASE"/ibk*; do

    [ -d "$SAMPLE_DIR" ] || continue

    SAMPLE=$(basename "$SAMPLE_DIR")

    REPORT="$SAMPLE_DIR/quast_comparison/report.tsv"

    if [ ! -f "$REPORT" ]; then
        echo "WARNING: report.tsv not found for $SAMPLE"
        continue
    fi

    # Extract rows from QUAST report
    NUM_CONTIGS=$(awk -F'\t' '$1=="# contigs" {print $2"\t"$3"\t"$4}' "$REPORT")
    TOTAL_LENGTH=$(awk -F'\t' '$1=="Total length" {print $2"\t"$3"\t"$4}' "$REPORT")
    LARGEST=$(awk -F'\t' '$1=="Largest contig" {print $2"\t"$3"\t"$4}' "$REPORT")
    N50=$(awk -F'\t' '$1=="N50" {print $2"\t"$3"\t"$4}' "$REPORT")

    # Split the three assembly columns
    read -r ALL_CONTIGS C500_CONTIGS C1000_CONTIGS <<< "$(echo "$NUM_CONTIGS" | tr '\t' ' ')"
    read -r ALL_LENGTH C500_LENGTH C1000_LENGTH <<< "$(echo "$TOTAL_LENGTH" | tr '\t' ' ')"
    read -r ALL_LARGEST C500_LARGEST C1000_LARGEST <<< "$(echo "$LARGEST" | tr '\t' ' ')"
    read -r ALL_N50 C500_N50 C1000_N50 <<< "$(echo "$N50" | tr '\t' ' ')"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$SAMPLE" \
        "$ALL_CONTIGS" \
        "$ALL_LENGTH" \
        "$ALL_LARGEST" \
        "$ALL_N50" \
        "$C500_CONTIGS" \
        "$C500_LENGTH" \
        "$C500_LARGEST" \
        "$C500_N50" \
        "$C1000_CONTIGS" \
        "$C1000_LENGTH" \
        "$C1000_LARGEST" \
        "$C1000_N50" >> "$OUT"

done

echo "========================================"
echo "Master QUAST summary created:"
echo "$OUT"
echo
echo "Number of samples:"
tail -n +2 "$OUT" | wc -l
echo "========================================"
