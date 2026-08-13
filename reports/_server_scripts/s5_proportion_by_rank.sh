#!/usr/bin/env bash
# Rank-stratified summary of the `proportion` length metric.
#
# WHY THIS EXISTS.  The mean of `proportion` is not a measurement.  Within each
# (gene, cell) the per-PAS proportions sum to 1, so the mean over all rows is
# exactly n_gene_cell_pairs / n_rows -- i.e. the reciprocal of the mean number of
# PAS per gene.  It is a property of the gene-to-PAS structure and carries no
# information about reads, stage or biology.  Measured on the corrected sweep it
# is constant to six decimal places across all four stages within every cell type
# (spread exactly 0.0), taking only a handful of distinct values genome-wide.
#
# Any comparison of "mean proportion" across conditions is therefore comparing a
# constant, will always report "no difference", and would be misread as evidence
# that APA does not change.
#
# The informative summary is proportion stratified by PAS rank: rank 1 is the
# proximal-most site, higher ranks are progressively distal.  A shift of read
# share from high to low ranks across stages is what 3'UTR shortening would look
# like in this metric, and unlike the mean it is free to vary.
#
# Ranks are capped into a 1..5+ bucket so the table stays small.  Everything is
# reported unconditionally and conditioned on the gene having reads in that cell,
# because section 12.1 showed the unconditional form tracks detection rate.
#
# Parallel across cell types; nothing written outside $OUT / $TMPDIR.
set -euo pipefail

SWEEP=/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed
LEN=$SWEEP/runs/B1_cohort_full/B3_switch/length
OUT=$SWEEP/analysis_extra
PARTS=$OUT/_prop_rank_parts
export TMPDIR=/mnt/ssd2/Laugney_Aligned/rt
mkdir -p "$OUT" "$TMPDIR" "$PARTS"

JOBS=${JOBS:-6}

one_celltype() {
  d=$1
  ct=$(basename "$(dirname "$(dirname "$d")")")
  out="$PARTS/$ct.tsv"
  [ -s "$out" ] && { echo "[skip] $ct" >&2; return 0; }
  echo "[rank] $ct" >&2
  # proportion.tsv columns:
  #   1 gene_id  4 rank  5 cell  6 proportion  8 total_reads_gene  9 cluster
  awk -F'\t' -v CT="$ct" '
    NR == 1 { next }
    {
      r = $4 + 0; if (r > 5) r = 5;          # bucket 5 = "5th PAS or more distal"
      k = $9 "\x1f" r;
      n[k]++; sum[k] += $6;
      if ($8 + 0 > 0) { nc[k]++; sumc[k] += $6 }
    }
    END {
      for (k in n) {
        split(k, a, "\x1f");
        printf "%s\t%s\t%s\t%d\t%.8f\t%d\t%.8f\n",
          CT, a[1], a[2], n[k], sum[k] / n[k],
          nc[k], (nc[k] ? sumc[k] / nc[k] : -1)
      }
    }
  ' "$d" > "$out.tmp" && mv "$out.tmp" "$out"
}
export -f one_celltype
export PARTS

find "$LEN" -mindepth 3 -maxdepth 3 -name proportion.tsv -print0 \
  | xargs -0 -P "$JOBS" -I{} bash -c 'one_celltype "$@"' _ {}

DEST=$OUT/proportion_by_rank.tsv
printf 'celltype\tstage\trank_bucket\tn_rows\tmean_prop_all\tn_cov\tmean_prop_cov\n' > "$DEST"
cat "$PARTS"/*.tsv >> "$DEST"
echo "[rank] wrote $DEST" >&2
wc -l "$DEST" >&2
