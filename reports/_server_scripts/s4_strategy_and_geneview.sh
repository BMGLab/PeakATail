#!/usr/bin/env bash
# Regenerate the inputs for the three remaining pre-fix report figures, from the
# corrected sweep:
#
#   switch_strategy_real_data.png  -> differential-strategy comparison
#   length_strategy_real_data.png  -> length-strategy comparison
#   gene_walk_*.png                -> per-gene PAS tracks
#
# All three per-cell-type length tables are streamed once each (~91 GB total:
# classic 17 GB, proportion 57 GB, shannon 17 GB).  Nothing is loaded into
# memory and nothing is written outside $OUT / $TMPDIR.
#
# Every score is summarised BOTH unconditionally and conditioned on the gene
# having any reads in that cell.  Section 12.1 of the technical report showed the
# classic PDUI layer emits zero-coverage rows with a score of 0.0, which drove a
# spurious stage trend; the same question has to be asked of `proportion` and
# `shannon` before either is recommended, so the conditioning is carried through
# here rather than assumed away.
set -euo pipefail

SWEEP=/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed
B3=$SWEEP/runs/B1_cohort_full/B3_switch
OUT=$SWEEP/analysis_extra
export TMPDIR=/mnt/ssd2/Laugney_Aligned/rt
mkdir -p "$OUT" "$TMPDIR"

# Genes to render as PAS tracks. EZR is the most recurrent gene in the corrected
# recurrence table; DPYD and PHACTR1 are carried over from the superseded figures
# so the old and new versions can be compared like for like.
GENES="ENSG00000092820 ENSG00000188641 ENSG00000112137"

# ---------------------------------------------------------------------------
# 1. length-strategy summary  (score column differs per strategy)
# ---------------------------------------------------------------------------
LS=$OUT/length_strategy_summary.tsv
printf 'celltype\tstrategy\tstage\tn_rows\tn_units\tn_cells\tfrac_zero_cov\tmean_all\tsd_all\tn_cov\tmean_cov\tsd_cov\n' > "$LS"

# $1 gene_id, score col, reads col, cell col, stage col -- passed per strategy
summarise() {                      # <file> <celltype> <strategy> <score_col> <reads_col> <cell_col> <stage_col>
  awk -F'\t' -v CT="$2" -v ST="$3" -v SC="$4" -v RD="$5" -v CE="$6" -v CL="$7" '
    NR == 1 { next }
    {
      s = $CL; v = $SC + 0; r = $RD + 0;
      n[s]++; sum[s] += v; sq[s] += v * v;
      if (r == 0) zero[s]++;
      else { nc[s]++; sumc[s] += v; sqc[s] += v * v }
      units[s "\x1f" $1] = 1; cells[s "\x1f" $CE] = 1;
    }
    END {
      for (k in units) { split(k, a, "\x1f"); nu[a[1]]++ }
      for (k in cells) { split(k, a, "\x1f"); ncell[a[1]]++ }
      for (s in n) {
        m  = sum[s] / n[s];
        # clamp: E[x^2]-E[x]^2 can go marginally negative in floating point
        vv = (sq[s] / n[s]) - m * m; if (vv < 0) vv = 0;
        sd = (n[s] > 1) ? sqrt(vv) : 0;
        mc = (nc[s] ? sumc[s] / nc[s] : -1);
        vc = (nc[s] ? (sqc[s] / nc[s]) - mc * mc : 0); if (vc < 0) vc = 0;
        sc = (nc[s] > 1 ? sqrt(vc) : 0);
        printf "%s\t%s\t%s\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%d\t%.6f\t%.6f\n",
          CT, ST, s, n[s], nu[s], ncell[s], zero[s] / n[s], m, sd, nc[s], mc, sc
      }
    }
  ' "$1" >> "$LS"
}

for d in "$B3"/length/*/; do
  ct=$(basename "$d")
  echo "[length] $ct" >&2
  # classic:    6 pdui,        9 total_reads,      5 cell, 10 cluster
  [ -f "$d/classic/pdui_classic.tsv" ]      && summarise "$d/classic/pdui_classic.tsv"      "$ct" classic    6 9 5 10
  # proportion: 6 proportion,  8 total_reads_gene, 5 cell,  9 cluster
  [ -f "$d/proportion/proportion.tsv" ]     && summarise "$d/proportion/proportion.tsv"     "$ct" proportion 6 8 5 9
  # shannon:    6 normalized_entropy, 8 total_reads_gene, 4 cell, 9 cluster
  [ -f "$d/shannon/entropy_shannon.tsv" ]   && summarise "$d/shannon/entropy_shannon.tsv"   "$ct" shannon    6 8 4 9
done
echo "[length] wrote $LS" >&2

# ---------------------------------------------------------------------------
# 2. differential-strategy comparison: fisher vs nb_multi
# ---------------------------------------------------------------------------
DS=$OUT/diff_strategy_summary.tsv
printf 'celltype\tstrategy\tcontrast\tn_tested\tn_sig\tfrac_sig\n' > "$DS"
OV=$OUT/diff_strategy_overlap.tsv
printf 'celltype\tn_fisher_sig\tn_nbmulti_sig\tn_both\tjaccard\n' > "$OV"

for d in "$B3"/diff/*/; do
  ct=$(basename "$d")
  echo "[diff] $ct" >&2
  fsig=$TMPDIR/f_sig.$$; : > "$fsig"
  for f in "$d"/fisher/differential/fisher_*.tsv; do
    [ -f "$f" ] || continue
    contrast=$(basename "$f" .tsv | sed 's/^fisher_//')
    awk -F'\t' -v CT="$ct" -v C="$contrast" -v S="$fsig" '
      NR == 1 { next } { n++; if ($10 + 0 < 0.05) { sig++; print $1 >> S } }
      END { printf "%s\tfisher\t%s\t%d\t%d\t%.6f\n", CT, C, n, sig, (n ? sig / n : 0) }
    ' "$f" >> "$DS"
  done
  nb="$d/nb_multi/differential/nb_multi_omnibus.tsv"
  nsig=$TMPDIR/n_sig.$$; : > "$nsig"
  if [ -f "$nb" ]; then
    awk -F'\t' -v CT="$ct" -v S="$nsig" '
      NR == 1 { next } { n++; if ($3 + 0 < 0.05) { sig++; print $1 >> S } }
      END { printf "%s\tnb_multi\tomnibus\t%d\t%d\t%.6f\n", CT, n, sig, (n ? sig / n : 0) }
    ' "$nb" >> "$DS"
  fi
  # Jaccard over the union of PAS each strategy called significant anywhere
  sort -u "$fsig" > "$fsig.u" 2>/dev/null || : > "$fsig.u"
  sort -u "$nsig" > "$nsig.u" 2>/dev/null || : > "$nsig.u"
  nf=$(wc -l < "$fsig.u"); nn=$(wc -l < "$nsig.u")
  nb_both=$(comm -12 "$fsig.u" "$nsig.u" | wc -l)
  un=$(( nf + nn - nb_both ))
  jac=$(awk -v a="$nb_both" -v b="$un" 'BEGIN { printf "%.6f", (b ? a / b : 0) }')
  printf '%s\t%d\t%d\t%d\t%s\n' "$ct" "$nf" "$nn" "$nb_both" "$jac" >> "$OV"
  rm -f "$fsig" "$nsig" "$fsig.u" "$nsig.u"
done
echo "[diff] wrote $DS and $OV" >&2

# ---------------------------------------------------------------------------
# 3. gene-walk tracks: per-gene, per-stage read support at each PAS
# ---------------------------------------------------------------------------
GW=$OUT/gene_walk_tracks.tsv
printf 'gene_id\tcelltype\tstage\tpas_id\trank\tchrom\tstart\tend\tstrand\tn_cells\tn_cells_expr\treads\n' > "$GW"
for d in "$B3"/length/*/; do
  ct=$(basename "$d")
  f="$d/proportion/proportion.tsv"
  [ -f "$f" ] || continue
  awk -F'\t' -v CT="$ct" -v G="$GENES" '
    BEGIN { split(G, gg, " "); for (i in gg) want[gg[i]] = 1 }
    NR == 1 { next }
    ($1 in want) {
      k = $1 "\x1f" $9 "\x1f" $3;
      reads[k] += $7; cells[k]++; if ($7 + 0 > 0) expr[k]++;
      rank[k] = $4; chrom[k] = $10; st[k] = $11; en[k] = $12; sd[k] = $13;
    }
    END {
      for (k in reads) {
        split(k, a, "\x1f");
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%.1f\n",
          a[1], CT, a[2], a[3], rank[k], chrom[k], st[k], en[k], sd[k],
          cells[k], expr[k], reads[k]
      }
    }
  ' "$f" >> "$GW"
done
echo "[geneview] wrote $GW" >&2
wc -l "$LS" "$DS" "$OV" "$GW" >&2
