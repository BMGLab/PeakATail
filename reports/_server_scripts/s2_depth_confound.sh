#!/usr/bin/env bash
# Depth-confound cross-check for the "global 3'UTR shortening" claim.
#
# The B3 length-trend layer averages `pdui` over EVERY (gene, cell) row of
# classic/pdui_classic.tsv, including rows with total_reads == 0 (where pdui is
# reported as 0.0). If sequencing depth / detection rate falls across stages,
# mean PDUI falls too -- with no change in real 3'UTR usage. This script
# streams every per-cell PDUI table once (17 GB) and aggregates, per
# (celltype, stage), the quantities needed to separate the two explanations:
#
#   * mean PDUI over ALL rows            <- what the trend layer reports
#   * mean PDUI conditional on coverage  <- the depth-controlled estimate
#   * fraction of rows with zero coverage, and mean reads per row
#
# Writes one TSV. Nothing is written outside $OUT / $TMPDIR.
set -euo pipefail

SWEEP=/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed
LEN=$SWEEP/runs/B1_cohort_full/B3_switch/length
OUT=$SWEEP/analysis_extra
export TMPDIR=/mnt/ssd2/Laugney_Aligned/rt
mkdir -p "$OUT" "$TMPDIR"

DEST=$OUT/depth_confound_by_celltype_stage.tsv
printf 'celltype\tstage\tn_rows\tn_cells\tn_genes\tmean_total_reads\tfrac_zero_cov\tmean_pdui_all\tn_cov1\tmean_pdui_cov1\tn_cov5\tmean_pdui_cov5\tn_cov10\tmean_pdui_cov10\tn_cov20\tmean_pdui_cov20\n' > "$DEST"

for d in "$LEN"/*/classic/pdui_classic.tsv; do
  ct=$(basename "$(dirname "$(dirname "$d")")")
  echo "[depth] $ct" >&2
  # Columns: 1 gene_id .. 5 cell, 6 pdui, 7 prox_reads, 8 distal_reads,
  #          9 total_reads, 10 cluster(stage)
  awk -F'\t' -v CT="$ct" '
    NR == 1 { next }
    {
      s = $10; p = $6 + 0; tot = $9 + 0;
      n[s]++; sump[s] += p; sumr[s] += tot;
      if (tot == 0) zero[s]++;
      cells[s "\x1f" $5] = 1;
      genes[s "\x1f" $1] = 1;
      if (tot >= 1)  { n1[s]++;  p1[s]  += p }
      if (tot >= 5)  { n5[s]++;  p5[s]  += p }
      if (tot >= 10) { n10[s]++; p10[s] += p }
      if (tot >= 20) { n20[s]++; p20[s] += p }
    }
    END {
      for (k in cells) { split(k, a, "\x1f"); nc[a[1]]++ }
      for (k in genes) { split(k, a, "\x1f"); ng[a[1]]++ }
      for (s in n) {
        printf "%s\t%s\t%d\t%d\t%d\t%.6f\t%.6f\t%.8f\t%d\t%.8f\t%d\t%.8f\t%d\t%.8f\t%d\t%.8f\n",
          CT, s, n[s], nc[s], ng[s],
          sumr[s] / n[s], zero[s] / n[s], sump[s] / n[s],
          n1[s],  (n1[s]  ? p1[s]  / n1[s]  : -1),
          n5[s],  (n5[s]  ? p5[s]  / n5[s]  : -1),
          n10[s], (n10[s] ? p10[s] / n10[s] : -1),
          n20[s], (n20[s] ? p20[s] / n20[s] : -1)
      }
    }
  ' "$d" >> "$DEST"
done

echo "[depth] wrote $DEST" >&2
wc -l "$DEST" >&2
