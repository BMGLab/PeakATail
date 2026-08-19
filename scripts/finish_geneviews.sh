#!/usr/bin/env bash
# finish_geneviews.sh — SEPARATE, self-contained job. Polls for finish_all.sh's
# DONE_ALL.marker (so it only starts once the corrected diff/length/trend +
# hub reindex are done), then renders real interactive geneviews (plotly
# .html + matplotlib .png) for the significant shorten/lengthening TREND
# genes per B1_cohort_full celltype, into a clean, professor-browsable tree:
#   FINISH_ALL_2026-08-13/geneviews/<celltype>/gene_<ENSG>.{html,png,svg}
#
# Does NOT touch the running finish_all tmux/session/files — read-only
# against B3_switch/{combined,trend}, writes only under its own geneviews/
# output dir. Safe to launch now; it just waits.
set -uo pipefail

ROOT=/mnt/ssd2/Laugney_Aligned
EXP=$ROOT/peakatail_experiments
B1=$EXP/RERUN_2026-08_fixed/runs/B1_cohort_full
EMA=$ROOT/.peakatail_env/bin/ema
GTF=/home/sharedFolder/humanSTARindex/Homo_sapiens.GRCh38.99.gtf
TMPBASE=$ROOT/rt
LOGROOT=$EXP/FINISH_ALL_2026-08-13
OUT=$LOGROOT/geneviews
TOP_N="${TOP_N:-5}"                 # top-N trend genes per celltype
SPEARMAN_THRESH="${SPEARMAN_THRESH:-0.8}"   # same cutoff sweep_analysis.py uses for "recurrent"
CONCURRENCY="${CONCURRENCY:-28}"    # load is low (~9/112) and geneview doesn't oversubscribe like reannotate

MARKER=$LOGROOT/DONE_ALL.marker
GV_MARKER=$LOGROOT/GENEVIEWS_DONE.marker
SELECTION_TSV=$LOGROOT/geneview_gene_selection.tsv
GENE_LISTS_DIR=$LOGROOT/geneview_genes

mkdir -p "$OUT" "$LOGROOT/logs" "$GENE_LISTS_DIR" "$TMPBASE"
export XDG_CACHE_HOME="$LOGROOT/cache"
export XDG_CONFIG_HOME="$LOGROOT/config"
export MPLCONFIGDIR="$LOGROOT/cache/mpl"
export NUMBA_CACHE_DIR="$LOGROOT/cache/numba"
export PEAKATAIL_NO_TIMESTAMP=1
mkdir -p "$XDG_CACHE_HOME" "$XDG_CONFIG_HOME" "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"

LOG="$LOGROOT/logs/finish_geneviews.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

log "=== finish_geneviews.sh START (pid $$) — waiting on $MARKER ==="
waited=0
while [ ! -f "$MARKER" ]; do
  sleep 30
  waited=$((waited + 30))
  if [ $((waited % 600)) -eq 0 ]; then
    log "  still waiting on finish_all DONE_ALL.marker (${waited}s elapsed)"
  fi
done
log "finish_all DONE_ALL.marker found (waited ${waited}s) -- proceeding"

########################################################################
# 1. Select top-N trend genes per celltype from the CORRECTED
#    length_trend_by_gene.tsv (built off the rank-fixed classic PDUI in
#    finish_all's step 2). |spearman| >= threshold, ranked by |slope|.
########################################################################
log "=== selecting top-${TOP_N} trend genes per celltype (|spearman|>=${SPEARMAN_THRESH}) ==="
python3 - "$B1/B3_switch/trend" "$TOP_N" "$SPEARMAN_THRESH" "$GENE_LISTS_DIR" > "$SELECTION_TSV" <<'PYEOF'
import csv, glob, os, sys

trend_root, top_n, sp_thresh, gene_lists_dir = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
print("celltype\tgene_id\tslope\tspearman\tdirection")
for ct_dir in sorted(glob.glob(os.path.join(trend_root, "*"))):
    ct = os.path.basename(ct_dir)
    if ct.endswith(".PRE_FIX_bak") or ct.endswith("_PRE_FIX_bak"):
        continue
    f = os.path.join(ct_dir, "length_trend_by_gene.tsv")
    if not os.path.isfile(f):
        continue
    rows = []
    with open(f) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                slope = float(row["slope"])
                sp = float(row["spearman"])
            except (KeyError, ValueError):
                continue
            if abs(sp) >= sp_thresh:
                rows.append((row["gene_id"], slope, sp, row.get("direction", "")))
    rows.sort(key=lambda x: -abs(x[1]))
    picked = rows[:top_n]
    if picked:
        with open(os.path.join(gene_lists_dir, f"{ct}.txt"), "w") as gf:
            for gene_id, slope, sp, direction in picked:
                gf.write(gene_id + "\n")
                print(f"{ct}\t{gene_id}\t{slope}\t{sp}\t{direction}")
PYEOF
n_pairs=$(($(wc -l < "$SELECTION_TSV") - 1))
log "gene selection written -> $SELECTION_TSV ($n_pairs celltype-gene pairs, $(ls "$GENE_LISTS_DIR" | wc -l) celltypes)"

########################################################################
# 2. Render geneviews, one `ema switch geneview` call per celltype (all its
#    selected genes in one call — the CLI renders a panel per gene), then
#    flatten figures/ up one level for a clean browsable tree.
########################################################################
render_celltype() {
  local ct="$1"
  local genelist="$GENE_LISTS_DIR/${ct}.txt"
  [ -f "$genelist" ] || return 0
  local h5="$B1/B3_switch/combined/${ct}.h5ad"
  local pb="$B1/B3_switch/combined/pasbed.bed"
  if [ ! -f "$h5" ]; then
    echo "[$(date '+%H:%M:%S')] SKIP $ct (no combined h5ad)"
    return 0
  fi
  local outdir="$OUT/${ct}"
  mkdir -p "$outdir"

  local need=0
  while IFS= read -r g; do
    [ -f "$outdir/gene_${g}.html" ] || need=1
  done < "$genelist"
  if [ "$need" -eq 0 ]; then
    echo "[$(date '+%H:%M:%S')] SKIP $ct (already rendered)"
    return 0
  fi

  local gene_flags=()
  while IFS= read -r g; do
    gene_flags+=(--gene-id "$g")
  done < "$genelist"

  local tmp="$TMPBASE/gv$$_${RANDOM}"
  mkdir -p "$tmp"
  TMPDIR="$tmp" "$EMA" switch geneview -i "$h5" --pasbed "$pb" --gtf "$GTF" \
      --cluster-key stage "${gene_flags[@]}" \
      --plot-engine both --subtitle "$ct" \
      --no-progress -o "$outdir" \
      > "$LOGROOT/logs/geneview_${ct}.log" 2>&1
  local rc=$?
  if [ -d "$outdir/figures" ]; then
    mv "$outdir"/figures/* "$outdir"/ 2>/dev/null
    rmdir "$outdir/figures" 2>/dev/null
  fi
  echo "[$(date '+%H:%M:%S')] geneview $ct ($(wc -l < "$genelist") genes) exit=$rc"
  rm -rf "$tmp"
}
export -f render_celltype
export B1 EMA GTF TMPBASE OUT LOGROOT GENE_LISTS_DIR

log "=== rendering geneviews (${CONCURRENCY}-way) ==="
ls "$GENE_LISTS_DIR" | sed 's/\.txt$//' | xargs -P "$CONCURRENCY" -I{} bash -c 'render_celltype "$@"' _ {} \
  >> "$LOG" 2>&1
log "=== geneview rendering done ==="

########################################################################
# 3. Manifest + DONE marker
########################################################################
python3 - "$OUT" "$SELECTION_TSV" > "$LOGROOT/geneviews_manifest.tsv" <<'PYEOF'
import csv, os, sys

out_root, sel_path = sys.argv[1], sys.argv[2]
print("celltype\tgene_id\tslope\tspearman\tdirection\thtml\tpng")
with open(sel_path) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        ct, g = row["celltype"], row["gene_id"]
        html = os.path.join(out_root, ct, f"gene_{g}.html")
        png = os.path.join(out_root, ct, f"gene_{g}.png")
        print(f"{ct}\t{g}\t{row['slope']}\t{row['spearman']}\t{row['direction']}\t{os.path.exists(html)}\t{os.path.exists(png)}")
PYEOF

{
  echo "finish_geneviews.sh COMPLETE at $(date)"
  echo "celltypes processed   : $(ls "$GENE_LISTS_DIR" | wc -l)"
  echo "gene-celltype pairs   : $(($(wc -l < "$SELECTION_TSV") - 1))"
  echo "html rendered         : $(find "$OUT" -name '*.html' | wc -l)"
  echo "png rendered          : $(find "$OUT" -name '*.png' | wc -l)"
  echo "manifest              : $LOGROOT/geneviews_manifest.tsv"
  echo "output tree           : $OUT"
} > "$GV_MARKER"
cat "$GV_MARKER" | tee -a "$LOG"
log "=== finish_geneviews.sh COMPLETE -- marker at $GV_MARKER ==="
