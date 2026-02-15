#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_fast_tcr_pipeline.sh \
    --samples /abs/path/samples.tsv \
    --transcriptome /abs/path/refdata-gex \
    --vdj-reference /abs/path/refdata-vdj \
    --outdir /abs/path/results \
    [--total-cores 32] \
    [--total-mem-gb 384] \
    [--jobs 8]

Required TSV columns:
  sample_id  gex_fastqs  gex_fastq_id  vdj_fastqs  vdj_fastq_id  [chemistry]
EOF
}

SAMPLES=""
TRANSCRIPTOME=""
VDJ_REFERENCE=""
OUTDIR=""
TOTAL_CORES=32
TOTAL_MEM_GB=384
JOBS=8

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES="$2"; shift 2 ;;
    --transcriptome) TRANSCRIPTOME="$2"; shift 2 ;;
    --vdj-reference) VDJ_REFERENCE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --total-cores) TOTAL_CORES="$2"; shift 2 ;;
    --total-mem-gb) TOTAL_MEM_GB="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$SAMPLES" || -z "$TRANSCRIPTOME" || -z "$VDJ_REFERENCE" || -z "$OUTDIR" ]]; then
  usage
  exit 1
fi

if ! command -v cellranger >/dev/null 2>&1; then
  echo "ERROR: cellranger not found in PATH" >&2
  exit 1
fi

mkdir -p "$OUTDIR"/{cellranger_runs,multi_csvs,logs,postprocess}

CORES_PER_JOB=$(( TOTAL_CORES / JOBS ))
MEM_PER_JOB=$(( TOTAL_MEM_GB / JOBS ))
if (( CORES_PER_JOB < 1 )); then CORES_PER_JOB=1; fi
if (( MEM_PER_JOB < 8 )); then MEM_PER_JOB=8; fi

echo "total_cores=$TOTAL_CORES jobs=$JOBS cores_per_job=$CORES_PER_JOB mem_per_job_gb=$MEM_PER_JOB"

run_one() {
  local sample_id="$1"
  local gex_fastqs="$2"
  local gex_fastq_id="$3"
  local vdj_fastqs="$4"
  local vdj_fastq_id="$5"
  local chemistry="$6"

  local sample_dir="$OUTDIR/cellranger_runs/$sample_id"
  local csv="$OUTDIR/multi_csvs/${sample_id}.csv"
  local log="$OUTDIR/logs/${sample_id}.log"

  mkdir -p "$sample_dir"

  cat > "$csv" <<EOF
[gene-expression]
reference,$TRANSCRIPTOME
chemistry,$chemistry

[vdj]
reference,$VDJ_REFERENCE

[libraries]
fastq_id,fastqs,feature_types
$gex_fastq_id,$gex_fastqs,Gene Expression
$vdj_fastq_id,$vdj_fastqs,VDJ-T
EOF

  (
    cd "$OUTDIR/cellranger_runs"
    cellranger multi \
      --id "$sample_id" \
      --csv "$csv" \
      --localcores "$CORES_PER_JOB" \
      --localmem "$MEM_PER_JOB"
  ) >"$log" 2>&1
}

PIDS=()
ACTIVE=0

while IFS=$'\t' read -r sample_id gex_fastqs gex_fastq_id vdj_fastqs vdj_fastq_id chemistry; do
  if [[ "$sample_id" == "sample_id" ]]; then
    continue
  fi
  if [[ -z "$chemistry" ]]; then
    chemistry="auto"
  fi

  run_one "$sample_id" "$gex_fastqs" "$gex_fastq_id" "$vdj_fastqs" "$vdj_fastq_id" "$chemistry" &
  PIDS+=("$!")
  ACTIVE=$((ACTIVE + 1))

  if (( ACTIVE >= JOBS )); then
    wait -n
    ACTIVE=$((ACTIVE - 1))
  fi
done < "$SAMPLES"

for pid in "${PIDS[@]}"; do
  wait "$pid"
done

python3 scripts/postprocess_tcr_fast.py \
  --samples-tsv "$SAMPLES" \
  --cellranger-root "$OUTDIR/cellranger_runs" \
  --outdir "$OUTDIR/postprocess"

echo "Done. Outputs in: $OUTDIR/postprocess"
