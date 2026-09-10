#!/usr/bin/env bash
# finish_all.sh — ONE self-contained, restart-safe job to run everything left
# before the BioLab server becomes unreachable (user leaving the lab kills
# SSH for everyone). Launch once in a dedicated tmux session; no further ssh
# needed. Runs SEQUENTIALLY:
#   1. B1_cohort_full: mwu_percell diff, 24 celltypes (new strategy alongside
#      fisher/nb_multi in the canonical layout).
#   2. B1_cohort_full: classic length recompute (rank fix), 24 celltypes.
#      Backs up the old canonical dir as .PRE_FIX_bak, then overwrites, then
#      regenerates the trend/ (3'UTR slope) that depends on classic PDUI.
#   3. FILTER_EFFECT_2026-08: 4-scenario diff(mwu_percell) + length(classic +
#      proportion) re-run, reusing the already-combined h5ads frozen out of
#      nextflow's work dir in step 0 (no reclustering/re-celltyping).
#   4. analyze_filter_effect() -> corrected comparison table (diff_strategy
#      mwu_percell).
#   5. Reindex the hub so corrected B1 data shows in the frontend.
#   6. DONE marker + summary log at a known ssd2 path.
#
# Safety: nothing written to $HOME; all scratch under /mnt/ssd2 (TMPDIR/XDG/
# MPL/NUMBA pinned); every ema invocation gets its own short-path TMPDIR
# (AF_UNIX 107-char limit); skip-if-exists on every unit of work so a restart
# (tmux respawn, resumed script) picks up where it left off; stays clear of
# TRIM_RERUN_2026-08 (separate TMPDIR root, separate output root, disjoint
# celltype-file work -- switch diff/length, not reannotate, so no
# get_n_jobs() oversubscription risk at 14-way x 4 threads = 56 cores).
set -uo pipefail   # NOT -e: one celltype's failure must not kill the whole run

ROOT=/mnt/ssd2/Laugney_Aligned
EXP=$ROOT/peakatail_experiments
B1=$EXP/RERUN_2026-08_fixed/runs/B1_cohort_full
FE=$EXP/FILTER_EFFECT_2026-08/runs
EMA=$ROOT/.peakatail_env/bin/ema
PY=$ROOT/.peakatail_env/bin/python
GTF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf
TMPBASE=$ROOT/rt
STAGE_ORDER='Normal,StageI,IVprimary,Met'
LOGROOT=$EXP/FINISH_ALL_2026-08-13
THREADS=4
CONCURRENCY=14

mkdir -p "$LOGROOT"/{logs,cache/mpl,cache/numba,config} "$TMPBASE"
export XDG_CACHE_HOME="$LOGROOT/cache"
export XDG_CONFIG_HOME="$LOGROOT/config"
export MPLCONFIGDIR="$LOGROOT/cache/mpl"
export NUMBA_CACHE_DIR="$LOGROOT/cache/numba"
export PEAKATAIL_NO_TIMESTAMP=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

MAINLOG="$LOGROOT/logs/finish_all.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MAINLOG"; }

has_tsv() { find "$1" -maxdepth 1 -name '*.tsv' -print -quit 2>/dev/null | grep -q .; }

log "=== finish_all.sh START (pid $$) ==="

########################################################################
# STEP 0: freeze the filter-effect combined/ snapshots out of nextflow's
# ephemeral work_v2fresh dir into a permanent canonical location (same
# layout as B1's on-disk B3_switch/combined/) -- survives any nextflow
# clean, and makes the rest of this script independent of nextflow state.
########################################################################
log "=== STEP 0: freeze FE combined/ snapshots ==="
declare -A FE_WORK=(
  [baseline]="$EXP/FILTER_EFFECT_2026-08/work_v2fresh/10/8a25143840295d1b2c08f87a9b18a6/combined"
  [atlas_filter]="$EXP/FILTER_EFFECT_2026-08/work_v2fresh/b3/dddac21bc4453a22513ee752f9c3f3/combined"
  [annot_filter_3utr]="$EXP/FILTER_EFFECT_2026-08/work_v2fresh/a5/e89477a83d1a402ff3c4fcc2fda153/combined"
  [ip_filter]="$EXP/FILTER_EFFECT_2026-08/work_v2fresh/1b/85b540fa54001f872a0d65f142b98a/combined"
)
for sc in baseline atlas_filter annot_filter_3utr ip_filter; do
  dest="$FE/$sc/B3_switch/combined"
  if [ -f "$dest/pasbed.bed" ]; then
    log "  $sc: already frozen ($(ls "$dest"/*.h5ad 2>/dev/null | wc -l) h5ads), skipping"
    continue
  fi
  mkdir -p "$dest"
  cp -n "${FE_WORK[$sc]}"/*.h5ad "$dest/" 2>>"$MAINLOG"
  cp -n "${FE_WORK[$sc]}/pasbed.bed" "$dest/" 2>>"$MAINLOG"
  log "  $sc: froze $(ls "$dest"/*.h5ad 2>/dev/null | wc -l) h5ads -> $dest"
done

########################################################################
# STEP 1: B1 mwu_percell diff, 24 celltypes.
########################################################################
log "=== STEP 1: B1 mwu_percell diff ==="
B1_H5DIR="$B1/B3_switch/combined"
B1_PB="$B1_H5DIR/pasbed.bed"

run_mwu_one() {
  local h5="$1" sl out tmp rc
  sl=$(basename "$h5" .h5ad)
  out="$B1/B3_switch/diff/${sl}/mwu_percell"
  if has_tsv "$out/differential" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP B1 mwu_percell $sl (already present)"
    return 0
  fi
  tmp="$TMPBASE/m$$_${RANDOM}"
  mkdir -p "$tmp"
  TMPDIR="$tmp" "$EMA" switch diff -i "$h5" --pasbed "$B1_PB" --gtf "$GTF" \
      --cluster-key stage --strategy mwu_percell --marker-top-n 0 --fdr 0.05 \
      --threads "$THREADS" --no-progress --no-plots -o "$out" \
      > "$LOGROOT/logs/b1_mwu_${sl}.log" 2>&1
  rc=$?
  echo "[$(date '+%H:%M:%S')] B1 mwu_percell $sl exit=$rc"
  rm -rf "$tmp"
  return $rc
}
export -f run_mwu_one has_tsv
export B1 B1_PB EMA GTF TMPBASE THREADS LOGROOT

ls "$B1_H5DIR"/*.h5ad | xargs -P "$CONCURRENCY" -I{} bash -c 'run_mwu_one "$@"' _ {} \
  >> "$MAINLOG" 2>&1
log "STEP 1 done ($(find "$B1/B3_switch/diff" -path '*mwu_percell/differential*' -name '*.tsv' | wc -l) tsv files across 24 celltypes)"

########################################################################
# STEP 2: B1 classic length recompute (rank fix) + trend regen, 24 celltypes.
########################################################################
log "=== STEP 2: B1 classic length recompute + trend regen ==="

run_classic_one() {
  local h5="$1" sl dir bak tmp rc pdui trendout
  sl=$(basename "$h5" .h5ad)
  dir="$B1/B3_switch/length/${sl}/classic"
  bak="$B1/B3_switch/length/${sl}/classic.PRE_FIX_bak"
  if [ -d "$dir" ] && [ ! -d "$bak" ]; then
    mv "$dir" "$bak"
    echo "[$(date '+%H:%M:%S')] B1 classic $sl: backed up old -> classic.PRE_FIX_bak"
  fi
  if [ -f "$dir/pdui_classic.tsv" ]; then
    echo "[$(date '+%H:%M:%S')] SKIP B1 classic $sl (already recomputed)"
  else
    tmp="$TMPBASE/c$$_${RANDOM}"
    mkdir -p "$tmp"
    TMPDIR="$tmp" "$EMA" switch length -i "$h5" --cluster-key stage --strategy classic \
        --pdui-pseudocount 1.0 --threads "$THREADS" --no-progress --no-plots -o "$dir" \
        > "$LOGROOT/logs/b1_classic_${sl}.log" 2>&1
    rc=$?
    echo "[$(date '+%H:%M:%S')] B1 classic $sl exit=$rc"
    rm -rf "$tmp"
    [ $rc -ne 0 ] && return $rc
  fi
  # trend regen off the fresh pdui_classic.tsv (depends on classic, must
  # follow it -- old trend/ would otherwise silently keep pre-fix numbers).
  pdui=$(find "$dir" -maxdepth 1 -name 'pdui_classic.tsv' | head -1)
  trendout="$B1/B3_switch/trend/${sl}"
  trendbak="$B1/B3_switch/trend/${sl}.PRE_FIX_bak"
  if [ -n "$pdui" ]; then
    if [ -d "$trendout" ] && [ ! -d "$trendbak" ]; then
      mv "$trendout" "$trendbak"
    fi
    if [ -f "$trendout/length_trend.json" ]; then
      echo "[$(date '+%H:%M:%S')] SKIP B1 trend $sl (already regenerated)"
    else
      "$EMA" switch trend --pdui "$pdui" --stage-order "$STAGE_ORDER" \
          --stage-col cluster --value-col pdui --no-progress -o "$trendout" \
          > "$LOGROOT/logs/b1_trend_${sl}.log" 2>&1
      echo "[$(date '+%H:%M:%S')] B1 trend $sl exit=$?"
    fi
  fi
}
export -f run_classic_one
export STAGE_ORDER

ls "$B1_H5DIR"/*.h5ad | xargs -P "$CONCURRENCY" -I{} bash -c 'run_classic_one "$@"' _ {} \
  >> "$MAINLOG" 2>&1
log "STEP 2 done ($(find "$B1/B3_switch/length" -path '*/classic/pdui_classic.tsv' -not -path '*PRE_FIX_bak*' | wc -l) classic files, $(find "$B1/B3_switch/trend" -name 'length_trend.json' -not -path '*PRE_FIX_bak*' | wc -l) trend files)"

########################################################################
# STEP 3: filter-effect 4-scenario re-run -- diff(mwu_percell) +
# length(classic + proportion) + trend regen off classic, reusing the
# frozen combined/ h5ads from step 0.
########################################################################
log "=== STEP 3: filter-effect 4-scenario diff+length+trend re-run ==="

run_fe_one() {
  local sc="$1" h5="$2" sl pb diffout tmp rc
  sl=$(basename "$h5" .h5ad)
  pb="$FE/$sc/B3_switch/combined/pasbed.bed"

  diffout="$FE/$sc/B3_switch/diff/${sl}/mwu_percell"
  if has_tsv "$diffout/differential" 2>/dev/null; then
    echo "[$(date '+%H:%M:%S')] SKIP FE $sc diff/mwu_percell $sl (already present)"
  else
    tmp="$TMPBASE/fd$$_${RANDOM}"
    mkdir -p "$tmp"
    TMPDIR="$tmp" "$EMA" switch diff -i "$h5" --pasbed "$pb" --gtf "$GTF" \
        --cluster-key stage --strategy mwu_percell --marker-top-n 0 --fdr 0.05 \
        --threads "$THREADS" --no-progress --no-plots -o "$diffout" \
        > "$LOGROOT/logs/fe_${sc}_diff_${sl}.log" 2>&1
    echo "[$(date '+%H:%M:%S')] FE $sc diff/mwu_percell $sl exit=$?"
    rm -rf "$tmp"
  fi

  local classic_dir="$FE/$sc/B3_switch/length/${sl}/classic"
  for S in classic proportion; do
    local lenout="$FE/$sc/B3_switch/length/${sl}/${S}"
    local lenbak="$FE/$sc/B3_switch/length/${sl}/${S}.PRE_FIX_bak"
    if [ -d "$lenout" ] && [ ! -d "$lenbak" ]; then
      mv "$lenout" "$lenbak"
      echo "[$(date '+%H:%M:%S')] FE $sc length/$S $sl: backed up old -> ${S}.PRE_FIX_bak"
    fi
    if [ -f "$lenout/pdui_${S}.tsv" ] || has_tsv "$lenout" 2>/dev/null; then
      echo "[$(date '+%H:%M:%S')] SKIP FE $sc length/$S $sl (already recomputed)"
      continue
    fi
    tmp="$TMPBASE/fl$$_${RANDOM}"
    mkdir -p "$tmp"
    TMPDIR="$tmp" "$EMA" switch length -i "$h5" --cluster-key stage --strategy "$S" \
        --pdui-pseudocount 1.0 --threads "$THREADS" --no-progress --no-plots -o "$lenout" \
        > "$LOGROOT/logs/fe_${sc}_length_${S}_${sl}.log" 2>&1
    echo "[$(date '+%H:%M:%S')] FE $sc length/$S $sl exit=$?"
    rm -rf "$tmp"
  done

  # trend regen off the fresh classic pdui.
  local pdui trendout trendbak
  pdui=$(find "$classic_dir" -maxdepth 1 -name 'pdui_classic.tsv' 2>/dev/null | head -1)
  trendout="$FE/$sc/B3_switch/trend/${sl}"
  trendbak="$FE/$sc/B3_switch/trend/${sl}.PRE_FIX_bak"
  if [ -n "$pdui" ]; then
    if [ -d "$trendout" ] && [ ! -d "$trendbak" ]; then
      mv "$trendout" "$trendbak"
    fi
    if [ -f "$trendout/length_trend.json" ]; then
      echo "[$(date '+%H:%M:%S')] SKIP FE $sc trend $sl (already regenerated)"
    else
      "$EMA" switch trend --pdui "$pdui" --stage-order "$STAGE_ORDER" \
          --stage-col cluster --value-col pdui --no-progress -o "$trendout" \
          > "$LOGROOT/logs/fe_${sc}_trend_${sl}.log" 2>&1
      echo "[$(date '+%H:%M:%S')] FE $sc trend $sl exit=$?"
    fi
  fi
}
export -f run_fe_one
export FE

JOBLIST="$LOGROOT/fe_joblist.txt"
: > "$JOBLIST"
for sc in baseline atlas_filter annot_filter_3utr ip_filter; do
  for h5 in "$FE/$sc/B3_switch/combined"/*.h5ad; do
    [ -e "$h5" ] || continue
    echo "${sc}|${h5}" >> "$JOBLIST"
  done
done
log "  FE job list: $(wc -l < "$JOBLIST") (scenario,celltype) pairs"

cat "$JOBLIST" | xargs -P "$CONCURRENCY" -I{} bash -c '
  IFS="|" read -r sc h5 <<< "{}"
  run_fe_one "$sc" "$h5"
' >> "$MAINLOG" 2>&1
log "STEP 3 done ($(find "$FE" -path '*mwu_percell/differential*' -name '*.tsv' | wc -l) FE mwu_percell tsvs, $(find "$FE" -name 'pdui_classic.tsv' -not -path '*PRE_FIX_bak*' | wc -l) FE classic files)"

########################################################################
# STEP 4: analyze_filter_effect() -> corrected comparison table
# (diff-strategy=mwu_percell so it reflects the fixed diff test, not fisher).
########################################################################
log "=== STEP 4: analyze_filter_effect() (diff-strategy=mwu_percell) ==="
ANALYSIS_OUT="$LOGROOT/analyze_filter_effect_corrected"
"$PY" -m ema.benchmark.sweep_analysis \
  --scenario baseline="$FE/baseline" \
  --scenario atlas_filter="$FE/atlas_filter" \
  --scenario annot_filter_3utr="$FE/annot_filter_3utr" \
  --scenario ip_filter="$FE/ip_filter" \
  --diff-strategy mwu_percell \
  --out "$ANALYSIS_OUT" \
  > "$LOGROOT/logs/analyze_filter_effect.log" 2>&1
log "STEP 4 done (exit=$?) -> $ANALYSIS_OUT/filter_effect_summary.md"

########################################################################
# STEP 5: reindex the hub (B1's corrected mwu_percell + classic now show up).
# DuckDB is single-writer -- stop the persistent service, index via a
# throwaway container, then bring the persistent service back up. Falls
# back to exec-while-up if `run --rm` doesn't work for any reason, so this
# step can't strand the hub down.
########################################################################
log "=== STEP 5: hub reindex ==="
cd /mnt/ssd2/Laugney_Aligned/hub/peakatail-hub || log "  WARN: hub compose dir not found, skipping reindex"
if [ -d /mnt/ssd2/Laugney_Aligned/hub/peakatail-hub ]; then
  DC="docker compose"; docker compose version >/dev/null 2>&1 || DC="docker-compose"
  $DC stop backend >> "$LOGROOT/logs/hub_reindex.log" 2>&1
  if ! $DC run --rm backend hub index /runs >> "$LOGROOT/logs/hub_reindex.log" 2>&1; then
    log "  'run --rm' indexing failed, falling back to exec-while-up"
    $DC up -d backend >> "$LOGROOT/logs/hub_reindex.log" 2>&1
    for i in $(seq 1 30); do
      s=$($DC ps --format "{{.Service}} {{.Health}}" 2>/dev/null | awk '/backend/{print $2}')
      [ "$s" = "healthy" ] && break
      sleep 5
    done
    $DC exec -T backend hub index /runs >> "$LOGROOT/logs/hub_reindex.log" 2>&1
  fi
  $DC up -d backend >> "$LOGROOT/logs/hub_reindex.log" 2>&1
  for i in $(seq 1 30); do
    s=$($DC ps --format "{{.Service}} {{.Health}}" 2>/dev/null | awk '/backend/{print $2}')
    [ "$s" = "healthy" ] && break
    sleep 5
  done
  log "STEP 5 done (backend health: ${s:-unknown})"
fi

########################################################################
# DONE marker
########################################################################
{
  echo "finish_all.sh COMPLETE at $(date)"
  echo "B1 mwu_percell diff tsvs      : $(find "$B1/B3_switch/diff" -path '*mwu_percell/differential*' -name '*.tsv' | wc -l)"
  echo "B1 classic recompute files    : $(find "$B1/B3_switch/length" -path '*/classic/pdui_classic.tsv' -not -path '*PRE_FIX_bak*' | wc -l) / 24"
  echo "B1 trend regenerated          : $(find "$B1/B3_switch/trend" -name 'length_trend.json' -not -path '*PRE_FIX_bak*' | wc -l) / 24"
  echo "FE mwu_percell diff tsvs      : $(find "$FE" -path '*mwu_percell/differential*' -name '*.tsv' | wc -l)"
  echo "FE classic recompute files    : $(find "$FE" -name 'pdui_classic.tsv' -not -path '*PRE_FIX_bak*' | wc -l)"
  echo "FE proportion recompute files : $(find "$FE" -path '*/proportion/*' -name '*.tsv' -not -path '*PRE_FIX_bak*' | wc -l)"
  echo "FE trend regenerated          : $(find "$FE" -name 'length_trend.json' -not -path '*PRE_FIX_bak*' | wc -l)"
  echo "analyze_filter_effect summary : $ANALYSIS_OUT/filter_effect_summary.md"
  echo "hub reindex log               : $LOGROOT/logs/hub_reindex.log"
} > "$LOGROOT/DONE_ALL.marker"
cat "$LOGROOT/DONE_ALL.marker" | tee -a "$MAINLOG"
log "=== finish_all.sh COMPLETE -- marker at $LOGROOT/DONE_ALL.marker ==="
