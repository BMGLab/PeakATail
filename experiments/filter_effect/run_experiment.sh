#!/usr/bin/env bash
# Launch the FILTER-EFFECT experiment on BioLab — UNATTENDED-SAFE.
#
#   bash run_experiment.sh              # launch/resume all 4 scenarios
#
# HOME is nearly full and /tmp lives on it, so EVERY nextflow/cache/temp path
# is pinned to ssd2 below — nothing is written to $HOME. Dedicated,
# previously-unused output root (nothing under RERUN_2026-08_fixed touched):
#   /mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08
set -euo pipefail

EXP_ROOT="${EXP_ROOT:-/mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08}"
PIPE_DIR="$(cd "$(dirname "$0")" && pwd)"

export NXF_HOME="${EXP_ROOT}/nxf_home"
export NXF_TEMP="${EXP_ROOT}/nxf_temp"
export NXF_WORK="${EXP_ROOT}/work"
export NXF_PLUGINS_DIR="${EXP_ROOT}/nxf_home/plugins"
export TMPDIR="/mnt/ssd2/Laugney_Aligned/rt"   # SHORT — AF_UNIX socket path-length gotcha
export XDG_CACHE_HOME="${EXP_ROOT}/cache"
export XDG_CONFIG_HOME="${EXP_ROOT}/config"
export MPLCONFIGDIR="${EXP_ROOT}/cache/mpl"
export NUMBA_CACHE_DIR="${EXP_ROOT}/cache/numba"

mkdir -p "${EXP_ROOT}"/{runs,logs,work,nxf_home,nxf_temp,cache/mpl,cache/numba,config,refs} /mnt/ssd2/Laugney_Aligned/rt

echo "== launching FILTER_EFFECT sweep -> ${EXP_ROOT} =="
cd "${PIPE_DIR}"
nextflow -log "${EXP_ROOT}/logs/nextflow.log" \
    run "${PIPE_DIR}/main.nf" -profile biolab -resume \
    -work-dir "${NXF_WORK}" \
    --exp_root "${EXP_ROOT}" \
    "$@"

echo
echo "== DONE. Results under: ${EXP_ROOT} =="
echo "   per-scenario: ${EXP_ROOT}/runs/{baseline,atlas_filter,annot_filter_3utr,ip_filter}/{07_clustering,B2_gex_celltyping,B3_switch}"
echo "   provenance:   ${EXP_ROOT}/logs/{report.html,timeline.html,trace.txt}"
