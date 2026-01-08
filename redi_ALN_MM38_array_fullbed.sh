#!/bin/bash
#SBATCH --job-name=redi_ALN_MM38
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rotembc1@mskcc.org
#SBATCH --chdir=/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts
#SBATCH --output=/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts/logs/redi_ALN_MM38_%A_%a.out
#SBATCH --error=/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts/logs/redi_ALN_MM38_%A_%a.err

set -euo pipefail
umask 002
export LC_ALL=C
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

LOG_DIR="/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts/logs"
mkdir -p "$LOG_DIR"
RUNLOG="${LOG_DIR}/redi_ALN_MM38_runtime_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec > >(tee -a "$RUNLOG") 2>&1

echo "============================================================"
echo "[BOOT] $(date)"
echo "[JOB ] job_id=${SLURM_JOB_ID:-NA}"
echo "[TASK] array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
echo "[NODE] $(hostname)"
echo "[PWD ] $(pwd)"
echo "============================================================"
trap 'echo "[FATAL] line=$LINENO exit=$? cmd=$BASH_COMMAND"; exit 1' ERR

# -------------------- ENV (conda) --------------------
eval "$(conda shell.bash hook)"
conda activate reditools2

PY="$(command -v python)"
SAM="$(command -v samtools)"
BEDTOOLS="$(command -v bedtools)"

echo "[INFO] python   = $PY"
echo "[INFO] samtools = $SAM"
echo "[INFO] bedtools = $BEDTOOLS"

python - <<'PY'
import pysam
print("[SYSTEM] PYSAM VERSION", pysam.__version__)
print("[SYSTEM] PYSAM PATH", pysam.__path__)
PY

# -------------------- PATHS --------------------
PIPE_BASE="/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline"

IN_DIR="/data1/greenbab/users/suns3/MDanderson/bulkRNA_cellline/ALN/sorted_bams"
OUT_DIR="${PIPE_BASE}/output"
WORK="${OUT_DIR}/_work_ALN_TE"

# FINAL output folders (NO _v7)
OUT_SINE="${OUT_DIR}/SINE"
OUT_LINE="${OUT_DIR}/LINE"
OUT_LTR="${OUT_DIR}/LTR"
OUT_NONTE="${OUT_DIR}/NonTE"
OUT_WHOLE="${OUT_DIR}/Whole"

mkdir -p "$WORK" "$OUT_SINE" "$OUT_LINE" "$OUT_LTR" "$OUT_NONTE" "$OUT_WHOLE"

# mm10 reference + TE beds (your BAMs are labeled MM38 but you’re using mm10/mm38 mouse build refs here)
SEQ_REF="${PIPE_BASE}/seq_ref_mm10"
REF="${SEQ_REF}/mm10.fa"
SINE_BED="${SEQ_REF}/mm10_SINE.bed"
LINE_BED="${SEQ_REF}/mm10_LINE.bed"
LTR_BED="${SEQ_REF}/mm10_LTR.bed"

RT="/data1/greenbab/users/rotembc1/pipeline2.0/REDItools2/src/cineca/reditools.py"

echo "[INFO] REDItools = $RT"
echo "[INFO] IN_DIR    = $IN_DIR"
echo "[INFO] REF       = $REF"

[[ -s "$REF" ]] || { echo "[FATAL] missing REF: $REF"; exit 1; }
[[ -s "${REF}.fai" ]] || "$SAM" faidx "$REF"
[[ -s "$SINE_BED" ]] || { echo "[FATAL] missing SINE_BED: $SINE_BED"; exit 1; }
[[ -s "$LINE_BED" ]] || { echo "[FATAL] missing LINE_BED: $LINE_BED"; exit 1; }
[[ -s "$LTR_BED"  ]] || { echo "[FATAL] missing LTR_BED: $LTR_BED"; exit 1; }
[[ -f "$RT" ]] || { echo "[FATAL] missing REDItools script: $RT"; exit 1; }

# -------------------- pick ONLY MM38 BAMs --------------------
mapfile -t BAMS < <(ls -1 "${IN_DIR}"/*_MM38.sorted.bam 2>/dev/null | sort -V)
N=${#BAMS[@]}
echo "[INFO] N_BAMS = $N"
(( N > 0 )) || { echo "[FATAL] no *_MM38.sorted.bam found in $IN_DIR"; exit 1; }

IDX=$((SLURM_ARRAY_TASK_ID - 1))
if (( IDX < 0 || IDX >= N )); then
  echo "[FATAL] task id ${SLURM_ARRAY_TASK_ID} out of range (1..${N})"
  exit 1
fi

BAM="${BAMS[$IDX]}"
b="$(basename "$BAM")"
SAMPLE="${b%.sorted.bam}"

echo "------------------------------------------------------------"
echo "[INFO] SAMPLE=$SAMPLE"
echo "[INFO] BAM=$BAM"

"$SAM" quickcheck -v "$BAM" || { echo "[FATAL] BAM failed samtools quickcheck (corrupt?): $BAM"; exit 1; }

# -------------------- build REF-matching BAM (prevents KeyError contigs) --------------------
SAMPLE_WORK="${WORK}/${SAMPLE}"
mkdir -p "$SAMPLE_WORK"

REF_FAI="${REF}.fai"
REF_CONTIGS="${SAMPLE_WORK}/ref.contigs.txt"
BAM_CONTIGS="${SAMPLE_WORK}/bam.contigs.txt"
KEEP_CONTIGS="${SAMPLE_WORK}/${SAMPLE}.keep.refcontigs.txt"

cut -f1 "$REF_FAI" | sort -u > "$REF_CONTIGS"
"$SAM" idxstats "$BAM" | cut -f1 | grep -v '^\*$' | sort -u > "$BAM_CONTIGS"
comm -12 "$REF_CONTIGS" "$BAM_CONTIGS" > "$KEEP_CONTIGS"

NKEEP=$(wc -l < "$KEEP_CONTIGS" | awk '{print $1}')
echo "[INFO] REF∩BAM contigs = $NKEEP"
if (( NKEEP < 5 )); then
  echo "[FATAL] Too few contigs overlap between BAM and REF. Check build/species."
  echo "        head REF contigs:"; head -5 "$REF_CONTIGS" || true
  echo "        head BAM contigs:"; head -5 "$BAM_CONTIGS" || true
  exit 1
fi

REFMATCH_SORT="${SAMPLE_WORK}/${SAMPLE}.refmatch.sorted.bam"

echo "[INFO] Building REF-matching BAM -> $REFMATCH_SORT"
"$SAM" view -h "$BAM" \
| awk -v FS="\t" -v OFS="\t" -v keep="$KEEP_CONTIGS" '
BEGIN{ while((getline<keep)>0){ ok[$1]=1 } }
/^@/{
  if($1=="@SQ"){
    sn=""
    for(i=1;i<=NF;i++){ if($i ~ /^SN:/){ sn=substr($i,4); break } }
    if(sn=="" || ok[sn]) print
  } else print
  next
}
{
  r=$3; rn=$7
  if(!ok[r]) next
  if(rn!="=" && rn!="*" && !ok[rn]) next
  print
}' \
| "$SAM" view -b - \
| "$SAM" sort -@ "$OMP_NUM_THREADS" -o "$REFMATCH_SORT" -

"$SAM" index "$REFMATCH_SORT"

# -------------------- run REDItools (WRITE FULL TABLE AS .bed) --------------------
# This is the "i.bed" style output you want.
WHOLE_BED="${OUT_WHOLE}/${SAMPLE}_RNA_redi.bed"
REDI_STDERR="${SAMPLE_WORK}/${SAMPLE}.reditools.stderr.log"

echo "[INFO] Running REDItools -> $WHOLE_BED"
set +e
"$PY" "$RT" -f "$REFMATCH_SORT" -r "$REF" -o "$WHOLE_BED" -t "${OMP_NUM_THREADS}" 2> "$REDI_STDERR"
rc=$?
set -e
echo "[INFO] REDItools exit_code=$rc"

if [[ ! -s "$WHOLE_BED" ]]; then
  echo "[WARN] REDItools produced empty output for $SAMPLE. Writing empty class beds."
  : > "${OUT_SINE}/${SAMPLE}_SINE_RNA_redi.bed"
  : > "${OUT_LINE}/${SAMPLE}_LINE_RNA_redi.bed"
  : > "${OUT_LTR}/${SAMPLE}_LTR_RNA_redi.bed"
  : > "${OUT_NONTE}/${SAMPLE}_nonTE_RNA_redi.bed"
  exit 0
fi

# -------------------- Split FULL REDItools table by TE class (still FULL columns) --------------------
# Convert TE interval beds into (chr, position) keys; then filter the REDItools table by those keys.
# REDItools Position is 1-based; our TE beds are 0-based -> keypos = start+1

make_keys () {
  local bed="$1"
  local outkeys="$2"
  awk -v OFS="\t" '{print $1, $2+1}' "$bed" | sort -k1,1 -k2,2n -u > "$outkeys"
}

filter_reditools_by_keys () {
  local keys="$1"
  local inbed="$2"
  local outbed="$3"
  awk -v OFS="\t" '
    NR==FNR { key[$1 FS $2]=1; next }
    FNR==1 { print; next }
    { k=$1 FS $2; if (k in key) print }
  ' "$keys" "$inbed" > "$outbed"
}

LINE_KEYS="${SAMPLE_WORK}/${SAMPLE}.LINE.keys"
SINE_KEYS="${SAMPLE_WORK}/${SAMPLE}.SINE.keys"
LTR_KEYS="${SAMPLE_WORK}/${SAMPLE}.LTR.keys"

make_keys "$LINE_BED" "$LINE_KEYS"
make_keys "$SINE_BED" "$SINE_KEYS"
make_keys "$LTR_BED"  "$LTR_KEYS"

LINE_OUT="${OUT_LINE}/${SAMPLE}_LINE_RNA_redi.bed"
SINE_OUT="${OUT_SINE}/${SAMPLE}_SINE_RNA_redi.bed"
LTR_OUT="${OUT_LTR}/${SAMPLE}_LTR_RNA_redi.bed"
NONTE_OUT="${OUT_NONTE}/${SAMPLE}_nonTE_RNA_redi.bed"

echo "[INFO] Filtering REDItools table into TE classes (full columns)"
filter_reditools_by_keys "$LINE_KEYS" "$WHOLE_BED" "$LINE_OUT"
filter_reditools_by_keys "$SINE_KEYS" "$WHOLE_BED" "$SINE_OUT"
filter_reditools_by_keys "$LTR_KEYS"  "$WHOLE_BED" "$LTR_OUT"

# NonTE = Whole minus (SINE ∪ LINE ∪ LTR)
TE_KEYS="${SAMPLE_WORK}/${SAMPLE}.TE.keys"
cat "$SINE_KEYS" "$LINE_KEYS" "$LTR_KEYS" | sort -k1,1 -k2,2n -u > "$TE_KEYS"

awk -v OFS="\t" '
  NR==FNR { te[$1 FS $2]=1; next }
  FNR==1 { print; next }
  { k=$1 FS $2; if (!(k in te)) print }
' "$TE_KEYS" "$WHOLE_BED" > "$NONTE_OUT"

echo "[DONE] $(date) wrote:"
echo "  $WHOLE_BED"
echo "  $LINE_OUT"
echo "  $SINE_OUT"
echo "  $LTR_OUT"
echo "  $NONTE_OUT"

