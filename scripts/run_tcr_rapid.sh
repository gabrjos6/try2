#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_tcr_rapid.sh \
    --samples /abs/path/rapid_tcr_samples.tsv \
    --outdir /abs/path/rapid_tcr_out \
    [--tool auto|mixcr|trust4] \
    [--species hsa] \
    [--mixcr-preset generic-single-cell-gex] \
    [--total-cores 32] \
    [--jobs 16] \
    [--limit-input 0] \
    [--trust4-fasta /abs/path/human_IMGT+C.fa]

Required TSV columns:
  sample_id   r1   r2
EOF
}

SAMPLES=""
OUTDIR=""
TOOL="auto"
SPECIES="hsa"
MIXCR_PRESET="generic-single-cell-gex"
TOTAL_CORES=32
JOBS=16
LIMIT_INPUT=0
TRUST4_FASTA=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples) SAMPLES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --tool) TOOL="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --mixcr-preset) MIXCR_PRESET="$2"; shift 2 ;;
    --total-cores) TOTAL_CORES="$2"; shift 2 ;;
    --jobs) JOBS="$2"; shift 2 ;;
    --limit-input) LIMIT_INPUT="$2"; shift 2 ;;
    --trust4-fasta) TRUST4_FASTA="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$SAMPLES" || -z "$OUTDIR" ]]; then
  usage
  exit 1
fi

mkdir -p "$OUTDIR"/{raw,logs,postprocess}
CORES_PER_JOB=$(( TOTAL_CORES / JOBS ))
if (( CORES_PER_JOB < 1 )); then CORES_PER_JOB=1; fi

select_tool() {
  if [[ "$TOOL" == "mixcr" || "$TOOL" == "trust4" ]]; then
    echo "$TOOL"
    return
  fi
  if command -v mixcr >/dev/null 2>&1; then
    echo "mixcr"
  elif command -v run-trust4 >/dev/null 2>&1; then
    echo "trust4"
  else
    echo "none"
  fi
}

run_mixcr() {
  local sample="$1" r1="$2" r2="$3"
  local sdir="$OUTDIR/raw/$sample"
  local prefix="$sdir/$sample"
  local log="$OUTDIR/logs/$sample.log"
  mkdir -p "$sdir"

  local limit_arg=()
  if [[ "$LIMIT_INPUT" != "0" ]]; then
    limit_arg=(--limit-input "$LIMIT_INPUT")
  fi

  mixcr analyze "$MIXCR_PRESET" \
    --species "$SPECIES" \
    --threads "$CORES_PER_JOB" \
    --rna \
    --no-reports \
    --no-json-reports \
    --use-local-temp \
    "${limit_arg[@]}" \
    -f \
    "$r1" "$r2" "$prefix" >"$log" 2>&1

  local best=""
  if [[ -f "${prefix}.clna" ]]; then
    best="${prefix}.clna"
  elif [[ -f "${prefix}.clns" ]]; then
    best="${prefix}.clns"
  else
    best="$(ls "$sdir"/*.clna "$sdir"/*.clns 2>/dev/null | head -n1 || true)"
  fi

  if [[ -n "$best" && -f "$best" ]]; then
    mixcr exportAirr -f "$best" "$sdir/$sample.airr.tsv" >>"$log" 2>&1
  fi
}

run_trust4() {
  local sample="$1" r1="$2" r2="$3"
  local sdir="$OUTDIR/raw/$sample"
  local log="$OUTDIR/logs/$sample.log"
  mkdir -p "$sdir"
  if [[ -z "$TRUST4_FASTA" ]]; then
    echo "TRUST4 requires --trust4-fasta" >"$log"
    return 1
  fi

  run-trust4 \
    -1 "$r1" -2 "$r2" \
    -f "$TRUST4_FASTA" \
    -t "$CORES_PER_JOB" \
    --skipReadRealign \
    --od "$sdir" \
    -o "$sample" >"$log" 2>&1
}

ACTIVE=0
PIDS=()
TOOL_RESOLVED="$(select_tool)"
if [[ "$TOOL_RESOLVED" == "none" ]]; then
  echo "ERROR: neither mixcr nor run-trust4 available in PATH" >&2
  exit 1
fi

echo "tool=$TOOL_RESOLVED total_cores=$TOTAL_CORES jobs=$JOBS cores_per_job=$CORES_PER_JOB limit_input=$LIMIT_INPUT"

while IFS=$'\t' read -r sample_id r1 r2; do
  if [[ "$sample_id" == "sample_id" ]]; then
    continue
  fi
  if [[ -z "$sample_id" || -z "$r1" || -z "$r2" ]]; then
    continue
  fi

  if [[ "$TOOL_RESOLVED" == "mixcr" ]]; then
    run_mixcr "$sample_id" "$r1" "$r2" &
  else
    run_trust4 "$sample_id" "$r1" "$r2" &
  fi
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

python3 scripts/postprocess_tcr_rapid.py \
  --samples-tsv "$SAMPLES" \
  --result-root "$OUTDIR/raw" \
  --outdir "$OUTDIR/postprocess"

echo "Done. Outputs in: $OUTDIR/postprocess"
