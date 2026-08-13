#!/usr/bin/env bash
# Re-run ONLY SWITCH_CELLTYPE (now with mwu_percell added + the engine-side
# classic-rank and proportion fixes already deployed) via `nextflow -resume`
# against the B1_cohort_full (phase1) session — reuses ALL cached upstream
# (peak-calling, clustering, GEX celltyping, SWITCH_COMBINE) since only
# SWITCH_CELLTYPE's script body changed. Deliberately skips run_sweep.sh's
# gen_configs.py regeneration step (not needed for a switch-only resume, and
# risks perturbing cached task hashes for unrelated processes).
#
#   bash run_switch_resume.sh
#
# Same ssd2-only pinning discipline as run_sweep.sh/run_phase3.sh — nothing
# to $HOME, TMPDIR short (AF_UNIX 107-char limit).
set -euo pipefail

SWEEP_ROOT="${SWEEP_ROOT:-/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed}"
MASTER="${MASTER:-/home/amiramiritabat/laughney_realign/peakatail_laughney.yaml}"
PIPE_DIR="${SWEEP_ROOT}/pipeline"   # deployed main.nf copy — same one run_phase3.sh reads

export NXF_HOME="${SWEEP_ROOT}/nxf_home"
export NXF_TEMP="${SWEEP_ROOT}/nxf_temp"
export NXF_WORK="${SWEEP_ROOT}/work"
export NXF_PLUGINS_DIR="${SWEEP_ROOT}/nxf_home/plugins"
export TMPDIR="/mnt/ssd2/Laugney_Aligned/rt"
export XDG_CACHE_HOME="${SWEEP_ROOT}/cache"
export XDG_CONFIG_HOME="${SWEEP_ROOT}/config"
export MPLCONFIGDIR="${SWEEP_ROOT}/cache/mpl"
export NUMBA_CACHE_DIR="${SWEEP_ROOT}/cache/numba"

mkdir -p "${SWEEP_ROOT}"/{runs,logs,work,nxf_home,nxf_temp,cache/mpl,cache/numba,config} /mnt/ssd2/Laugney_Aligned/rt

echo "== switch-only resume (phase1, B1_cohort_full) -> ${SWEEP_ROOT} =="
echo "   main.nf from: ${PIPE_DIR} (no gen_configs.py — configs already exist and are unchanged)"
nextflow -log "${SWEEP_ROOT}/logs/nextflow_switch_resume.log" \
    run "${PIPE_DIR}/main.nf" -profile biolab -resume \
    -work-dir "${NXF_WORK}" \
    --sweep_root "${SWEEP_ROOT}" \
    --master "${MASTER}" \
    --repo "${PIPE_DIR}" \
    --phase phase1 \
    "$@"

echo
echo "== DONE. B3_switch under: ${SWEEP_ROOT}/runs/B1_cohort_full/B3_switch =="
