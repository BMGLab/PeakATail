#!/usr/bin/env bash
# Launch Phase 3 (FILTER-EFFECT experiment) — ISOLATED from the RERUN sweep's
# own nextflow session/work-dir/logs/configs. Reads main.nf from the RERUN
# pipeline dir (where the corrected Phase 3 code is deployed) but never
# writes there: no gen_configs.py call (phase3 doesn't use params.configs),
# NXF_HOME/NXF_WORK/logs/trace/report/timeline all pinned under the
# dedicated FILTER_EFFECT_2026-08 root via --sweep_root, and
# --filter_effect_base_run is passed explicitly so it still resolves to the
# REAL (read-only) runs/grid/lg_annotate under RERUN_2026-08_fixed even
# though --sweep_root now points elsewhere.
set -euo pipefail

EXP_ROOT="${EXP_ROOT:-/mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08}"
RERUN_ROOT="${RERUN_ROOT:-/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed}"
PIPE_DIR="${RERUN_ROOT}/pipeline"           # read main.nf from here (deployed copy)
BASE_RUN="${RERUN_ROOT}/runs/grid/lg_annotate"   # READ-ONLY

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

test -f "${BASE_RUN}/posbed.bed" || { echo "FATAL: base run missing: ${BASE_RUN}/posbed.bed" >&2; exit 1; }

echo "== launching Phase 3 (FILTER-EFFECT) -> ${EXP_ROOT} =="
echo "   main.nf from: ${PIPE_DIR} (READ-ONLY: nothing under ${RERUN_ROOT} is written)"
echo "   base run:     ${BASE_RUN} (READ-ONLY)"
nextflow -log "${EXP_ROOT}/logs/nextflow.log" \
    run "${PIPE_DIR}/main.nf" -profile biolab -resume \
    -work-dir "${NXF_WORK}" \
    --sweep_root "${EXP_ROOT}" \
    --phase phase3 \
    --filter_effect_base_run "${BASE_RUN}" \
    "$@"

echo
echo "== DONE. Results under: ${EXP_ROOT} =="
echo "   per-scenario: ${EXP_ROOT}/runs/{baseline,atlas_filter,annot_filter_3utr,ip_filter}/{07_clustering,B2_gex_celltyping,B3_switch,branch_manifest.json}"
echo "   provenance:   ${EXP_ROOT}/logs/{report.html,timeline.html,trace.txt}"
