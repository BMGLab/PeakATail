# FILTER-EFFECT experiment

Compares how each PAS filter changes **clustering**, **differential APA**, and
**3'UTR shorten/lengthen** when the filter is applied to the downstream
matrix, on the Laughney 6-GSM method-comparison subset (spans all 4 stages:
Normal / StageI / IVprimary / Met), `lambda_gradient` peak-calling.

The pipeline is **Phase 3** of the existing corrected-rerun sweep
(`experiments/laughney/main.nf`) — there is no separate nextflow here. It is
opt-in: `bash run_sweep.sh --phase phase3` (never runs as part of the default
`--phase all`). See `main.nf`'s "PHASE 3" section for the processes
(`PREP_3UTR_BED`, `FILTER_EFFECT_BRANCH`, reusing `GEX_CELLTYPE` /
`SWITCH_COMBINE` / `SWITCH_CELLTYPE` verbatim).

## Mechanism: `peakatail reannotate`, NOT a fresh peak-call

Peak-calling is the expensive step and it's already done: Phase 3 branches
off **`runs/grid/lg_annotate`** — the 6-GSM subset, `lambda_gradient`, atlas
annotate + ip annotate run Phase 2's `GRID_RUN` already produced (row
`lg_annotate` in `peakcall_grid.tsv`, "keep-all atlas + fasta"). Each
scenario is an `peakatail reannotate --base-run runs/grid/lg_annotate` branch —
cheap (trim → annotate → preprocess → cluster only), no BAM streaming.

This is also what makes "PAS results keep everything, filter only applies
downstream" literally true here (unlike a fresh filtered `peakatail run` would be
— see the retracted first attempt at this experiment, corrected after
team-lead review): the PAS filters (new `--atlas-filter` / `--ip-filter` /
`--annot-filter` flags on `peakatail reannotate`, added in this PR) are applied to
a **copy** of `lg_annotate`'s `posbed.bed`/`negbed.bed` written into each
branch's own `--out` directory, immediately before `find_close()`.
`lg_annotate`'s own `posbed.bed`/`negbed.bed`/`annotatedpas.bed` are never
touched — it stays the complete, unfiltered PAS-results reference every
scenario's PAS count is compared against. Dropping a PAS from the filtered
copy only removes it from `find_close()`'s gene-assignment table, which is
what the clustering matrix is built from (`genes.index ∩ matrix_pas_ids`) —
so the drop is confined to the downstream matrix, exactly as asked.

## Scenarios (one-factor-at-a-time against `baseline`)

| scenario | `peakatail reannotate` flags (beyond `--base-run runs/grid/lg_annotate --gtf ...`) |
|---|---|
| `baseline` | none — reproduces `lg_annotate`'s downstream byte-for-byte |
| `atlas_filter` | `--atlas-filter` (drops PAS with no atlas match; reuses `lg_annotate`'s cached `unified/atlas_status.tsv`) |
| `annot_filter_3utr` | `--annot-filter --annotation-bed <3'UTR-only BED>` (keeps only PAS inside an annotated transcript 3'UTR) |
| `ip_filter` | `--ip-filter --genome-fasta <fa>` (drops PAS flagged as internal-priming artifacts) |

Then per scenario (unchanged from the original design, reusing Phase 1's
processes): `scripts/gex_celltyping.py` → `peakatail switch combine --group-key
stage --split-key celltype` → per-celltype `peakatail switch diff` (fisher,
nb_multi) + `peakatail switch length` (classic, shannon) + `peakatail switch trend
--stage-order Normal,StageI,IVprimary,Met`.

## 3'UTR-only annotation BED

`annot_filter_3utr` needs a 3'UTR-*only* BED, not `peakatail reannotate
--annot-filter`'s default fallback (no `--annotation-bed` given → the
GTF-derived gene BED, `ema/annotate/gtftobed.py`, `source_type="gene"` —
whole gene bodies, not 3'UTRs). Built once by `PREP_3UTR_BED` via
`scripts/build_3utr_bed.py`: scans the GTF for `three_prime_utr` feature rows
(col 3) directly, converts to BED, `bedtools merge -s` per strand. Output:
`<filter_effect_root>/refs/three_prime_utr.bed`.

## Output layout (dedicated root, `runs/grid/lg_annotate` is read-only)

```
/mnt/ssd2/Laugney_Aligned/peakatail_experiments/FILTER_EFFECT_2026-08/
  runs/<scenario>/{07_clustering,B2_gex_celltyping,B3_switch,pasbed.bed,posbed.bed,negbed.bed,branch_manifest.json,...}
  refs/three_prime_utr.bed
```

`branch_manifest.json`'s `"pas_filters"` key records `n_pas_before_filter` /
`n_pas_after_filter` / `n_pas_dropped` for that scenario, plus the raw
`apply_filters`/atlas-match stats.

## Analysis

`python -m ema.benchmark.sweep_analysis --scenario baseline=<root>/runs/baseline
--scenario atlas_filter=<root>/runs/atlas_filter --scenario
annot_filter_3utr=<root>/runs/annot_filter_3utr --scenario
ip_filter=<root>/runs/ip_filter --out <root>/logs/filter_effect_analysis`
(see `analyze_filter_effect()` in `ema/benchmark/sweep_analysis.py`) —
compares PAS-retained count/delta, clustering (n_clusters, GEX ARI/AMI),
Fisher differential hit counts, and 3'UTR trend across the 4 scenarios.
