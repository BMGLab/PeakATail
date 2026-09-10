#!/usr/bin/env bash
# Copy real figures from peakatail_runs into docs/assets/figures/
# Run from the repository root: bash docs/assets/figures/copy_figures.sh

set -e
cd "$(dirname "$0")/../../.."

RUNDIR="peakatail_runs/full_v8_2026-05-11_152746"
DIFF_RUN="${RUNDIR}/switch_diff_2026-05-11_205015"
LEN_RUN="${RUNDIR}/switch_length_2026-05-11_205045"
PROP_RUN="${RUNDIR}/switch_length_2026-05-11_194149"
GV_RUN="${RUNDIR}/switch_geneview_2026-05-11_212341"
OUT="docs/assets/figures"

cp "${RUNDIR}/figures/umap_default.png"               "${OUT}/umap_default.png"
cp "${RUNDIR}/figures/clusters_default.png"           "${OUT}/clusters_default.png"
cp "${RUNDIR}/figures/peak_qc_default.png"            "${OUT}/peak_qc_default.png"
cp "${RUNDIR}/figures/resource_timeline.png"          "${OUT}/resource_timeline.png"
cp "${DIFF_RUN}/figures/volcano_0_vs_4.png"           "${OUT}/volcano_0_vs_4.png"
cp "${LEN_RUN}/figures/pdui_distribution.png"         "${OUT}/pdui_distribution.png"
cp "${PROP_RUN}/figures/pdui_distribution.png"        "${OUT}/proportion_distribution.png"
cp "${GV_RUN}/figures/gene_ENSG00000103275.png"       "${OUT}/gene_ENSG00000103275.png"
cp "${DIFF_RUN}/figures/figures_INDEX.md"             "${OUT}/switch_diff_figures_INDEX.md"

echo "Copied 9 figures to ${OUT}/"
ls -lh "${OUT}/"
