#!/bin/bash -l
#SBATCH --job-name=IBK_QC
#SBATCH --array=1-209
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=168:00:00
#SBATCH --partition=guest,batch
#SBATCH --output=run_QC_%A_%a.out
#SBATCH --error=run_QC_%A_%a.err
#SBATCH --mail-user=echandrasekara2@huskers.unl.edu
#SBATCH --mail-type=FAIL,END

set -euo pipefail

# ============================================================
# USER-CONFIGURABLE PATHS
# ============================================================

PROJECT="/work/samodha/sachin/ShotgunM/IBK_over_1M"

RAW_ROOT="${PROJECT}/raw_reads"

# Adapter FASTA
ADAPTERS="${PROJECT}/adaptor/adapters.fa"

# ------------------------------------------------------------
# Human Bowtie2 index
# GRCh38 - confirmed on HCC
# ------------------------------------------------------------
HUMAN_INDEX="/work/HCC/BCRF/indices/bowtie2/Homo_sapiens/NCBI/GRCh38/genome"

# ------------------------------------------------------------
# Bovine Bowtie2 index
# ARS-UCD1.2
# ------------------------------------------------------------
BOVINE_INDEX="/work/HCC/BCRF/indices/bowtie2/Bos_taurus/NCBI/ARS-UCD1.2/genome"

# Main QC output
OUT_ROOT="${PROJECT}/QC_cleaned"

# Final reads for downstream MEGAHIT assembly
MEGAHIT_INPUT="${PROJECT}/MEGAHIT_single"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# BBDuk quality settings
TRIMQ=20
MINLEN=75
MAQ=20


# ============================================================
# SOFTWARE
# ============================================================

module purge
module load bbmap/39.06
module load bowtie/2.5
module load fastqc/0.12
module load biodata


# ============================================================
# HELPER FUNCTIONS
# ============================================================

log() {
    echo "[$(date '+%F %T')] $*"
}


fail() {
    echo "ERROR: $*" >&2
    exit 1
}


count_fastq_reads() {

    local f="$1"

    if [[ ! -s "$f" ]]; then
        echo 0
        return
    fi

    if [[ "$f" == *.gz ]]; then
        gzip -cd "$f" | awk 'END {print int(NR/4)}'
    else
        awk 'END {print int(NR/4)}' "$f"
    fi
}


pct() {

    local numerator="$1"
    local denominator="$2"

    awk \
        -v n="$numerator" \
        -v d="$denominator" \
        'BEGIN {
            if (d == 0)
                printf "0.00";
            else
                printf "%.2f", 100*n/d
        }'
}


# ============================================================
# CHECK REQUIRED PATHS
# ============================================================

[[ -d "$RAW_ROOT" ]] || \
    fail "Raw-read directory not found: $RAW_ROOT"

[[ -s "$ADAPTERS" ]] || \
    fail "Adapter FASTA not found: $ADAPTERS"

# Check human index
if [[ ! -f "${HUMAN_INDEX}.1.bt2" && ! -f "${HUMAN_INDEX}.1.bt2l" ]]; then
    fail "Human Bowtie2 index not found: $HUMAN_INDEX"
fi

# Check bovine index
if [[ ! -f "${BOVINE_INDEX}.1.bt2" && ! -f "${BOVINE_INDEX}.1.bt2l" ]]; then
    fail "Bovine Bowtie2 index not found: $BOVINE_INDEX"
fi


mkdir -p \
    "$OUT_ROOT" \
    "$MEGAHIT_INPUT"


# ============================================================
# FIND RAW R1 FILES
# ============================================================
#
# Expected format:
#
# ibk_1_1.fq.gz
# ibk_1_2.fq.gz
#
# or
#
# ibk1_1.fq.gz
# ibk1_2.fq.gz
#
# ============================================================

mapfile -t R1_FILES < <(
    find "$RAW_ROOT" \
        -maxdepth 1 \
        -type f \
        -name "*_1.fq.gz" \
        | sort -V
)


NUM_SAMPLES=${#R1_FILES[@]}

(( NUM_SAMPLES > 0 )) || \
    fail "No *_1.fq.gz files found in $RAW_ROOT"


log "Number of R1 files found: $NUM_SAMPLES"


# ============================================================
# RESOLVE SAMPLE FOR THIS ARRAY TASK
# ============================================================

INDEX=$((SLURM_ARRAY_TASK_ID - 1))


if (( INDEX >= NUM_SAMPLES )); then

    log "Array index ${SLURM_ARRAY_TASK_ID} exceeds number of samples (${NUM_SAMPLES})."

    exit 0

fi


RAW_R1="${R1_FILES[$INDEX]}"

sample=$(basename "$RAW_R1" "_1.fq.gz")

RAW_R2="${RAW_ROOT}/${sample}_2.fq.gz"


[[ -s "$RAW_R1" ]] || \
    fail "R1 missing: $RAW_R1"

[[ -s "$RAW_R2" ]] || \
    fail "R2 missing: $RAW_R2"


# ============================================================
# SAMPLE OUTPUT DIRECTORIES
# ============================================================

OUTDIR="${OUT_ROOT}/${sample}"

TMPDIR_SAMPLE="${OUTDIR}/tmp"

FASTQC_RAW="${OUTDIR}/fastqc_raw"

FASTQC_CLEAN="${OUTDIR}/fastqc_clean"

mkdir -p \
    "$OUTDIR" \
    "$TMPDIR_SAMPLE" \
    "$FASTQC_RAW" \
    "$FASTQC_CLEAN"


log "============================================================"
log "Sample: $sample"
log "R1: $RAW_R1"
log "R2: $RAW_R2"
log "Output: $OUTDIR"
log "============================================================"


# ============================================================
# RAW READ COUNT
# ============================================================

RAW_PAIRS=$(count_fastq_reads "$RAW_R1")

RAW_READS=$(( RAW_PAIRS * 2 ))


log "Raw read pairs: $RAW_PAIRS"
log "Raw total reads: $RAW_READS"


# ============================================================
# FASTQC - RAW READS
# ============================================================

log "Running FastQC on raw reads"

fastqc \
    --threads "$THREADS" \
    --outdir "$FASTQC_RAW" \
    "$RAW_R1" \
    "$RAW_R2"


# ============================================================
# STEP 1
# ADAPTER + QUALITY TRIMMING
# ============================================================

TRIM_R1="${TMPDIR_SAMPLE}/${sample}_trim_R1.fastq.gz"

TRIM_R2="${TMPDIR_SAMPLE}/${sample}_trim_R2.fastq.gz"

TRIM_SINGLE="${TMPDIR_SAMPLE}/${sample}_trim_singletons.fastq.gz"


log "Running adapter and quality trimming with BBDuk"


bbduk.sh \
    in1="$RAW_R1" \
    in2="$RAW_R2" \
    out1="$TRIM_R1" \
    out2="$TRIM_R2" \
    outs="$TRIM_SINGLE" \
    ref="$ADAPTERS" \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    tpe \
    tbo \
    qtrim=rl \
    trimq="$TRIMQ" \
    minlen="$MINLEN" \
    maq="$MAQ" \
    threads="$THREADS" \
    stats="${OUTDIR}/${sample}_adapter_quality_stats.txt"


TRIM_PAIRS=$(count_fastq_reads "$TRIM_R1")

TRIM_SINGLE_READS=$(count_fastq_reads "$TRIM_SINGLE")

TRIM_READS=$(( TRIM_PAIRS * 2 + TRIM_SINGLE_READS ))


log "Reads after adapter/quality trimming: $TRIM_READS"


# ============================================================
# STEP 2
# PHIX + SEQUENCING ARTIFACT REMOVAL
# ============================================================

PHIX_R1="${TMPDIR_SAMPLE}/${sample}_phixclean_R1.fastq.gz"

PHIX_R2="${TMPDIR_SAMPLE}/${sample}_phixclean_R2.fastq.gz"

PHIX_SINGLE_FROM_PAIRS="${TMPDIR_SAMPLE}/${sample}_phix_pair_singletons.fastq.gz"

PHIX_SINGLE="${TMPDIR_SAMPLE}/${sample}_phix_singletons.fastq.gz"


# ------------------------------------------------------------
# Paired reads
# ------------------------------------------------------------

log "Removing PhiX and sequencing artifacts from paired reads"


bbduk.sh \
    in1="$TRIM_R1" \
    in2="$TRIM_R2" \
    out1="$PHIX_R1" \
    out2="$PHIX_R2" \
    outs="$PHIX_SINGLE_FROM_PAIRS" \
    ref=artifacts,phix \
    k=31 \
    hdist=1 \
    threads="$THREADS" \
    ordered \
    cardinality \
    stats="${OUTDIR}/${sample}_phix_artifact_paired_stats.txt"


# ------------------------------------------------------------
# Existing trimming singletons
# ------------------------------------------------------------

log "Removing PhiX and sequencing artifacts from singleton reads"


if [[ -s "$TRIM_SINGLE" ]]; then

    bbduk.sh \
        in="$TRIM_SINGLE" \
        out="$PHIX_SINGLE" \
        ref=artifacts,phix \
        k=31 \
        hdist=1 \
        threads="$THREADS" \
        ordered \
        cardinality \
        stats="${OUTDIR}/${sample}_phix_artifact_singleton_stats.txt"

else

    : > "$PHIX_SINGLE"

fi


# ------------------------------------------------------------
# Combine clean singleton streams
# ------------------------------------------------------------

ALL_SINGLE="${TMPDIR_SAMPLE}/${sample}_all_phixclean_singletons.fastq.gz"

cat \
    "$PHIX_SINGLE" \
    "$PHIX_SINGLE_FROM_PAIRS" \
    > "$ALL_SINGLE"


PHIX_PAIRS=$(count_fastq_reads "$PHIX_R1")

PHIX_SINGLE_READS=$(count_fastq_reads "$ALL_SINGLE")

PHIX_READS=$(( PHIX_PAIRS * 2 + PHIX_SINGLE_READS ))


log "Reads after PhiX/artifact removal: $PHIX_READS"


# ============================================================
# STEP 3
# HUMAN HOST REMOVAL - GRCh38
# ============================================================

HUMAN_R1="${TMPDIR_SAMPLE}/${sample}_nonhuman_R1.fastq.gz"

HUMAN_R2="${TMPDIR_SAMPLE}/${sample}_nonhuman_R2.fastq.gz"

HUMAN_SINGLE="${TMPDIR_SAMPLE}/${sample}_nonhuman_singletons.fastq.gz"

HUMAN_LOG_PAIRED="${OUTDIR}/${sample}_bowtie2_human_paired.log"

HUMAN_LOG_SINGLE="${OUTDIR}/${sample}_bowtie2_human_singletons.log"


log "Removing HUMAN reads using GRCh38"


# ------------------------------------------------------------
# Human removal - paired reads
# ------------------------------------------------------------

bowtie2 \
    -x "$HUMAN_INDEX" \
    -1 "$PHIX_R1" \
    -2 "$PHIX_R2" \
    -p "$THREADS" \
    --very-sensitive \
    --un-conc-gz "${TMPDIR_SAMPLE}/${sample}_nonhuman_R%.fastq.gz" \
    -S /dev/null \
    2> "$HUMAN_LOG_PAIRED"


[[ -e "$HUMAN_R1" ]] || \
    fail "Human removal did not create $HUMAN_R1"

[[ -e "$HUMAN_R2" ]] || \
    fail "Human removal did not create $HUMAN_R2"


# ------------------------------------------------------------
# Human removal - singleton reads
# ------------------------------------------------------------

if [[ -s "$ALL_SINGLE" ]] && \
   [[ $(count_fastq_reads "$ALL_SINGLE") -gt 0 ]]; then

    bowtie2 \
        -x "$HUMAN_INDEX" \
        -U "$ALL_SINGLE" \
        -p "$THREADS" \
        --very-sensitive \
        --un-gz "$HUMAN_SINGLE" \
        -S /dev/null \
        2> "$HUMAN_LOG_SINGLE"

else

    : > "$HUMAN_SINGLE"

    echo "No singleton reads to align." \
        > "$HUMAN_LOG_SINGLE"

fi


HUMAN_PAIRS=$(count_fastq_reads "$HUMAN_R1")

HUMAN_SINGLE_READS=$(count_fastq_reads "$HUMAN_SINGLE")

AFTER_HUMAN_READS=$(( HUMAN_PAIRS * 2 + HUMAN_SINGLE_READS ))


HUMAN_REMOVED=$(( PHIX_READS - AFTER_HUMAN_READS ))

HUMAN_REMOVED_PCT=$(pct "$HUMAN_REMOVED" "$PHIX_READS")


log "Reads after human removal: $AFTER_HUMAN_READS"

log "Human reads removed: $HUMAN_REMOVED (${HUMAN_REMOVED_PCT}%)"


# ============================================================
# STEP 4
# BOVINE HOST REMOVAL - ARS-UCD1.2
# ============================================================

FINAL_R1="${OUTDIR}/${sample}_final_clean_R1.fastq.gz"

FINAL_R2="${OUTDIR}/${sample}_final_clean_R2.fastq.gz"

FINAL_SINGLE="${OUTDIR}/${sample}_final_clean_singletons.fastq.gz"

BOVINE_LOG_PAIRED="${OUTDIR}/${sample}_bowtie2_bovine_paired.log"

BOVINE_LOG_SINGLE="${OUTDIR}/${sample}_bowtie2_bovine_singletons.log"


log "Removing BOVINE reads using ARS-UCD1.2"


# ------------------------------------------------------------
# Bovine removal - paired reads
# ------------------------------------------------------------

bowtie2 \
    -x "$BOVINE_INDEX" \
    -1 "$HUMAN_R1" \
    -2 "$HUMAN_R2" \
    -p "$THREADS" \
    --very-sensitive \
    --un-conc-gz "${OUTDIR}/${sample}_final_clean_R%.fastq.gz" \
    -S /dev/null \
    2> "$BOVINE_LOG_PAIRED"


[[ -e "$FINAL_R1" ]] || \
    fail "Bovine removal did not create $FINAL_R1"

[[ -e "$FINAL_R2" ]] || \
    fail "Bovine removal did not create $FINAL_R2"


# ------------------------------------------------------------
# Bovine removal - singleton reads
# ------------------------------------------------------------

if [[ -s "$HUMAN_SINGLE" ]] && \
   [[ $(count_fastq_reads "$HUMAN_SINGLE") -gt 0 ]]; then

    bowtie2 \
        -x "$BOVINE_INDEX" \
        -U "$HUMAN_SINGLE" \
        -p "$THREADS" \
        --very-sensitive \
        --un-gz "$FINAL_SINGLE" \
        -S /dev/null \
        2> "$BOVINE_LOG_SINGLE"

else

    : > "$FINAL_SINGLE"

    echo "No non-human singleton reads to align." \
        > "$BOVINE_LOG_SINGLE"

fi


FINAL_PAIRS=$(count_fastq_reads "$FINAL_R1")

FINAL_SINGLE_READS=$(count_fastq_reads "$FINAL_SINGLE")

FINAL_READS=$(( FINAL_PAIRS * 2 + FINAL_SINGLE_READS ))


BOVINE_REMOVED=$(( AFTER_HUMAN_READS - FINAL_READS ))

BOVINE_REMOVED_PCT=$(pct "$BOVINE_REMOVED" "$AFTER_HUMAN_READS")


log "Final reads after bovine removal: $FINAL_READS"

log "Bovine reads removed: $BOVINE_REMOVED (${BOVINE_REMOVED_PCT}%)"


# ============================================================
# STEP 5
# FINAL FASTQC
# ============================================================

log "Running FastQC on final non-host reads"


FASTQC_INPUTS=(
    "$FINAL_R1"
    "$FINAL_R2"
)


if [[ -s "$FINAL_SINGLE" ]] && \
   [[ $(count_fastq_reads "$FINAL_SINGLE") -gt 0 ]]; then

    FASTQC_INPUTS+=("$FINAL_SINGLE")

fi


fastqc \
    --threads "$THREADS" \
    --outdir "$FASTQC_CLEAN" \
    "${FASTQC_INPUTS[@]}"


# ============================================================
# STEP 6
# READ TRACKING TABLE
# ============================================================

TRACK="${OUTDIR}/${sample}_read_tracking.tsv"


pct_raw() {

    pct "$1" "$RAW_READS"

}


{
    printf "sample\tstage\tpaired_pairs\tsingleton_reads\ttotal_reads\tpercent_of_raw_reads\n"

    printf "%s\traw\t%d\t0\t%d\t100.00\n" \
        "$sample" \
        "$RAW_PAIRS" \
        "$RAW_READS"

    printf "%s\tafter_adapter_quality\t%d\t%d\t%d\t%s\n" \
        "$sample" \
        "$TRIM_PAIRS" \
        "$TRIM_SINGLE_READS" \
        "$TRIM_READS" \
        "$(pct_raw "$TRIM_READS")"

    printf "%s\tafter_phix_artifact\t%d\t%d\t%d\t%s\n" \
        "$sample" \
        "$PHIX_PAIRS" \
        "$PHIX_SINGLE_READS" \
        "$PHIX_READS" \
        "$(pct_raw "$PHIX_READS")"

    printf "%s\tafter_human_removal\t%d\t%d\t%d\t%s\n" \
        "$sample" \
        "$HUMAN_PAIRS" \
        "$HUMAN_SINGLE_READS" \
        "$AFTER_HUMAN_READS" \
        "$(pct_raw "$AFTER_HUMAN_READS")"

    printf "%s\tfinal_after_bovine_removal\t%d\t%d\t%d\t%s\n" \
        "$sample" \
        "$FINAL_PAIRS" \
        "$FINAL_SINGLE_READS" \
        "$FINAL_READS" \
        "$(pct_raw "$FINAL_READS")"

} > "$TRACK"


# ============================================================
# HOST REMOVAL SUMMARY
# ============================================================

HOST_SUMMARY="${OUTDIR}/${sample}_host_removal_summary.tsv"


{
    printf "sample\tphix_clean_reads\tafter_human_reads\thuman_removed\thuman_removed_percent\tafter_bovine_reads\tbovine_removed\tbovine_removed_percent\tfinal_retention_percent\n"

    printf "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n" \
        "$sample" \
        "$PHIX_READS" \
        "$AFTER_HUMAN_READS" \
        "$HUMAN_REMOVED" \
        "$HUMAN_REMOVED_PCT" \
        "$FINAL_READS" \
        "$BOVINE_REMOVED" \
        "$BOVINE_REMOVED_PCT" \
        "$(pct_raw "$FINAL_READS")"

} > "$HOST_SUMMARY"


# ============================================================
# STEP 7
# COPY FINAL READS TO MEGAHIT INPUT
# ============================================================

log "Copying final cleaned reads to MEGAHIT input directory"


cp -f \
    "$FINAL_R1" \
    "$MEGAHIT_INPUT/"


cp -f \
    "$FINAL_R2" \
    "$MEGAHIT_INPUT/"


if [[ -s "$FINAL_SINGLE" ]] && \
   [[ $(count_fastq_reads "$FINAL_SINGLE") -gt 0 ]]; then

    cp -f \
        "$FINAL_SINGLE" \
        "$MEGAHIT_INPUT/"

fi


# ============================================================
# STEP 8
# CLEAN TEMPORARY FILES
# ============================================================

rm -rf "$TMPDIR_SAMPLE"


# ============================================================
# FINAL SUMMARY
# ============================================================

log "============================================================"

log "QC COMPLETE: $sample"

log "Raw reads:                 $RAW_READS"

log "After adapter trimming:    $TRIM_READS"

log "After PhiX/artifacts:      $PHIX_READS"

log "After human removal:       $AFTER_HUMAN_READS"

log "Human removed:             $HUMAN_REMOVED (${HUMAN_REMOVED_PCT}%)"

log "After bovine removal:      $FINAL_READS"

log "Bovine removed:            $BOVINE_REMOVED (${BOVINE_REMOVED_PCT}%)"

log "Final retention:           $(pct_raw "$FINAL_READS")%"

log "Tracking table:            $TRACK"

log "Host summary:              $HOST_SUMMARY"

log "============================================================"
