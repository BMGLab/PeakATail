# FILTER-EFFECT experiment

Compares how each PAS filter changes **clustering**, **differential APA**, and
**3'UTR shorten/lengthen** when the filter is applied to the downstream
matrix, on the Laughney 6-GSM method-comparison subset (spans all 4 stages:
Normal / StageI / IVprimary / Met), `lambda_gradient` peak-calling.

## Scenarios (one-factor-at-a-time against `baseline`)

| scenario | atlas-mode | ip-filter-mode | annot-filter |
|---|---|---|---|
| `baseline` | annotate (keep+flag) | annotate (keep+flag) | off |
| `atlas_filter` | **filter** (drop non-atlas PAS) | annotate | off |
| `annot_filter_3utr` | annotate | annotate | **on**, `--annotation-bed` = GTF-derived 3'UTR-only BED |
| `ip_filter` | annotate | **filter** (drop internal-priming PAS) | off |

Exact `ema run` invocation per scenario (see `main.nf::RUN_SCENARIO`):

```
ema run -c <A2_base_lambda_gradient.yaml (reused from the corrected-rerun sweep)> \
    --atlas <polyasite_v3.bed> --atlas-distance 50 --atlas-mode <annotate|filter> \
    --ip-filter --genome-fasta <fa> --ip-filter-mode <annotate|filter> \
    [--annot-filter --annotation-bed <three_prime_utr.bed>] \
    --output <exp_root>/runs/<scenario> --threads 16 --no-progress
```

Then per scenario: `scripts/gex_celltyping.py` → `ema switch combine
--group-key stage --split-key celltype` → per-celltype `ema switch diff`
(fisher, nb_multi) + `ema switch length` (classic, shannon) + `ema switch
trend --stage-order Normal,StageI,IVprimary,Met`.

## IMPORTANT caveat — read before interpreting results

The user's framing was "PAS results always keep every PAS; the filter only
applies going into clustering." At the `ema run` code level this is **not**
what `--atlas-mode filter` / `--ip-filter-mode filter` / `--annot-filter` do:
`ema/main.py::_apply_pas_filters()` runs **before** `find_close()` builds the
gene-assignment table that `write_pas_gene_artifacts()` writes
`annotatedpas.bed` from. So a "filter" scenario's own `annotatedpas.bed` /
`pasbed.bed` already has fewer rows than `baseline`'s — the drop is not
confined to the clustering matrix. There is no existing `ema` flag that keeps
a PAS in the result list while excluding it only from clustering (that would
need new engine code — out of scope here).

This still answers the question asked: `baseline` is the one scenario whose
PAS list is provably complete (nothing dropped), and each filtered
scenario's PAS-count delta against `baseline` **is** the "removed for
downstream" number the comparison reports. Just don't describe a filtered
scenario's own `annotatedpas.bed` as "the full PAS results, filter applied
only downstream" — it isn't, by construction of the current filter wiring.

## 3'UTR-only annotation BED

`annot_filter_3utr` needs a 3'UTR-*only* BED, not the default `--annot-filter`
fallback (`ema/annotate/gtftobed.py`'s gene BED, `source_type="gene"` —
whole gene bodies, not 3'UTRs). Built once by `BUILD_3UTR_BED` via
`scripts/build_3utr_bed.py`: scans the GTF for `three_prime_utr` feature rows
(col 3) directly, converts to BED, `bedtools merge -s` per strand. Output:
`<exp_root>/refs/three_prime_utr.bed`.

## Output layout

```
<exp_root>/runs/<scenario>/{07_clustering,B2_gex_celltyping,B3_switch,annotatedpas.bed,pasbed.bed,...}
<exp_root>/logs/{nextflow.log,report.html,timeline.html,trace.txt}
<exp_root>/refs/three_prime_utr.bed
```

## Analysis

`python -m ema.benchmark.sweep_analysis --filter-effect <exp_root>/runs
--scenario baseline=baseline --scenario atlas_filter=atlas_filter
--scenario annot_filter_3utr=annot_filter_3utr --scenario ip_filter=ip_filter
--out <exp_root>/logs/filter_effect_analysis` (see `analyze_filter_effect()`).
