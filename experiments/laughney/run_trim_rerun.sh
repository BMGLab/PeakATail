#!/usr/bin/env bash
# Trim-axis (A2/A3) CLEAN re-run — off B1_cohort_full's saved peak-calls via
# `peakatail reannotate`, one invocation per experiments/laughney/trim_cluster_grid.tsv
# row, matching main.nf's REANNOTATE_BRANCH script body exactly. Bypasses
# nextflow deliberately: `--phase phase2` would ALSO re-launch GRID_RUN
# (peakcall method sweep -- expensive re-peak-calling), which isn't wanted
# here; this script runs ONLY the reannotate branches.
#
# Isolation (per team-lead's go-ahead to run concurrently once isolated):
#   - Each branch's --out is a distinct directory -> distinct gtf_cache/
#     (ema/reannotate.py's process_gtf_cached(output_dir=str(out)) already
#     writes there, same as FILTER_EFFECT_BRANCH's already-verified isolation).
#   - Each branch gets its own TMPDIR (short, PID-based).
#   - ema/outputs.py's atomic pas_gene.tsv/annotatedpas.bed write (temp+
#     os.replace) is deployed on the server venv -- even in the (now
#     structurally impossible, since --out never collides) worst case, no
#     writer can produce a torn file.
#   - Dedicated, previously-unused output root -- runs/reannotate/A2_trim_*
#     in RERUN_2026-08_fixed (the original, buggy-numbers run) is NEVER
#     touched; this is a fresh comparison target.
#
# Usage: bash run_trim_rerun.sh [max_concurrent] [threads_per_branch]
set -euo pipefail

RERUN_ROOT="/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed"
TRIM_ROOT="/mnt/ssd2/Laugney_Aligned/peakatail_experiments/TRIM_RERUN_2026-08"
BASE_RUN="${RERUN_ROOT}/runs/B1_cohort_full"
GTF="/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf"
EMA="/mnt/ssd2/Laugney_Aligned/.peakatail_env/bin/ema"
GRID_TSV="$(cd "$(dirname "$0")" && pwd)/trim_cluster_grid.tsv"

MAX_CONCURRENT="${1:-6}"
THREADS_PER_BRANCH="${2:-8}"

mkdir -p "${TRIM_ROOT}"/{runs,logs,cache/mpl,cache/numba,config} /mnt/ssd2/Laugney_Aligned/rt

export XDG_CACHE_HOME="${TRIM_ROOT}/cache"
export XDG_CONFIG_HOME="${TRIM_ROOT}/config"
export MPLCONFIGDIR="${TRIM_ROOT}/cache/mpl"
export NUMBA_CACHE_DIR="${TRIM_ROOT}/cache/numba"
export PEAKATAIL_NO_TIMESTAMP=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

test -f "${BASE_RUN}/posbed.bed" || { echo "FATAL: base run missing: ${BASE_RUN}/posbed.bed" >&2; exit 1; }

echo "== Trim-axis re-run -> ${TRIM_ROOT} =="
echo "   base run (READ-ONLY): ${BASE_RUN}"
echo "   max_concurrent=${MAX_CONCURRENT} threads_per_branch=${THREADS_PER_BRANCH}"
echo "   estimated core usage: $((MAX_CONCURRENT * THREADS_PER_BRANCH)) cores"

run_one() {
    local branch_name="$1" max_gene_distance="$2" utr_multiplier="$3" \
          include_extended="$4" cluster_method="$5" resolution="$6" n_neighbors="$7"
    local out="${TRIM_ROOT}/runs/${branch_name}"
    local tmp="/mnt/ssd2/Laugney_Aligned/rt/trim_$$_${branch_name}"
    mkdir -p "$tmp"
    local ext_flag=""
    [ "$include_extended" = "true" ] && ext_flag="--include-extended"
    local nn_flag=""
    [ -n "$n_neighbors" ] && [ "$n_neighbors" != "NA" ] && nn_flag="--n-neighbors ${n_neighbors}"

    echo "[$(date '+%H:%M:%S')] START ${branch_name}"
    TMPDIR="$tmp" "$EMA" reannotate \
        --base-run "$BASE_RUN" \
        --out "$out" \
        --gtf "$GTF" \
        --max-gene-distance "$max_gene_distance" \
        --utr-multiplier "$utr_multiplier" ${ext_flag} \
        --cluster-method "$cluster_method" \
        --resolution "$resolution" ${nn_flag} \
        --threads "$THREADS_PER_BRANCH" \
        > "${TRIM_ROOT}/logs/${branch_name}.log" 2>&1
    local rc=$?
    echo "[$(date '+%H:%M:%S')] DONE  ${branch_name} (exit=${rc})"
    return $rc
}
export -f run_one
export TRIM_ROOT BASE_RUN GTF EMA THREADS_PER_BRANCH TMPDIR

# Read the grid (skip header/comments), dispatch with a concurrency cap via xargs -P.
tail -n +2 "$GRID_TSV" | grep -v '^#' | \
    xargs -P "$MAX_CONCURRENT" -I{} -d '\n' bash -c '
        IFS=$'"'"'\t'"'"' read -r bn dist mult ext method res nn <<< "{}"
        run_one "$bn" "$dist" "$mult" "$ext" "$method" "$res" "$nn"
    '

echo "== TRIM RERUN DONE =="
