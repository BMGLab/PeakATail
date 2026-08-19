#!/usr/bin/env bash
# ONE self-contained job: nextflow -resume for phase1 (B1_cohort_full switch),
# THEN nextflow -resume for phase3 (filter-effect 4 scenarios), THEN hub
# reindex, THEN a DONE marker. Launch once in tmux; runs to completion with
# no further ssh. Both resumes reuse every cached upstream task (peak-calling,
# clustering, GEX celltyping, SWITCH_COMBINE) unchanged -- only
# SWITCH_CELLTYPE's script body changed (mwu_percell added), so -resume
# re-executes ONLY that process, for every celltype in both sessions.
set -uo pipefail   # NOT -e: a non-zero nextflow exit on one phase must not
                    # skip cleanup/reporting of the rest

RR=/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed
FE=/mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08
MARKER_DIR=/mnt/ssd2/Laugney_Aligned/peakatail_experiments/SWITCH_RESUME_2026-08-13
mkdir -p "$MARKER_DIR"
LOG="$MARKER_DIR/run.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "=== run_switch_resume_all.sh START (pid $$) ==="

log "=== PHASE 1: nextflow -resume (B1_cohort_full switch, adds mwu_percell) ==="
bash "${RR}/pipeline/run_switch_resume.sh" >> "$LOG" 2>&1
p1rc=$?
log "PHASE 1 resume exit=${p1rc}"

log "=== PHASE 3: nextflow -resume (filter-effect 4 scenarios, adds mwu_percell) ==="
bash "${RR}/pipeline/run_phase3.sh" >> "$LOG" 2>&1
p3rc=$?
log "PHASE 3 resume exit=${p3rc}"

log "=== HUB REINDEX ==="
cd /mnt/ssd2/Laugney_Aligned/hub/peakatail-hub || log "WARN: hub compose dir not found, skipping reindex"
if [ -d /mnt/ssd2/Laugney_Aligned/hub/peakatail-hub ]; then
  DC="docker compose"; docker compose version >/dev/null 2>&1 || DC="docker-compose"
  $DC stop backend >> "$LOG" 2>&1
  if ! $DC run --rm backend hub index /runs >> "$LOG" 2>&1; then
    log "  'run --rm' indexing failed, falling back to exec-while-up"
    $DC up -d backend >> "$LOG" 2>&1
    for i in $(seq 1 30); do
      s=$($DC ps --format "{{.Service}} {{.Health}}" 2>/dev/null | awk '/backend/{print $2}')
      [ "$s" = "healthy" ] && break
      sleep 5
    done
    $DC exec -T backend hub index /runs >> "$LOG" 2>&1
  fi
  $DC up -d backend >> "$LOG" 2>&1
  for i in $(seq 1 30); do
    s=$($DC ps --format "{{.Service}} {{.Health}}" 2>/dev/null | awk '/backend/{print $2}')
    [ "$s" = "healthy" ] && break
    sleep 5
  done
  log "HUB REINDEX done (backend health: ${s:-unknown})"
fi

{
  echo "run_switch_resume_all.sh COMPLETE at $(date)"
  echo "phase1 (B1) resume exit : ${p1rc}"
  echo "phase3 (FE) resume exit : ${p3rc}"
  echo "B1 mwu_percell diff tsvs      : $(find "$RR/runs/B1_cohort_full/B3_switch/diff" -path '*mwu_percell/differential*' -name '*.tsv' 2>/dev/null | wc -l)"
  echo "B1 classic (rank-fixed) files : $(find "$RR/runs/B1_cohort_full/B3_switch/length" -path '*/classic/pdui_classic.tsv' 2>/dev/null | wc -l) / 24"
  echo "FE mwu_percell diff tsvs      : $(find "$FE/runs" -path '*mwu_percell/differential*' -name '*.tsv' 2>/dev/null | wc -l)"
  echo "FE classic (rank-fixed) files : $(find "$FE/runs" -path '*/classic/pdui_classic.tsv' 2>/dev/null | wc -l)"
  echo "nextflow logs: ${RR}/logs/nextflow_switch_resume.log , ${FE}/logs/nextflow.log"
} > "$MARKER_DIR/DONE_ALL.marker"
cat "$MARKER_DIR/DONE_ALL.marker" | tee -a "$LOG"
log "=== run_switch_resume_all.sh COMPLETE -- marker at $MARKER_DIR/DONE_ALL.marker ==="
