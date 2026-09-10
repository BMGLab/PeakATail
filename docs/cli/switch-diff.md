# `peakatail switch diff`

`peakatail switch diff` tests for differential alternative polyadenylation (APA)
between every pair of Leiden clusters in one or more `clusters.h5ad` files
produced by `peakatail run`. For each cluster pair (c1, c2) and each PAS that passes
the cell-count filter, the command runs the selected statistical strategy
(default: Fisher exact test) and writes one TSV per pair under a
`switch_diff_<timestamp>/differential/` subdirectory.

When the input `--h5ad` files come from a recognisable
`peakatail_runs/<run>/` path and `--output` is left at its default, the output
is automatically routed **inside** the originating run directory as
`peakatail_runs/<run>/switch_diff_<timestamp>/`. This keeps all results for a
run self-contained. Source: `ema/cli/common.py::resolve_subcommand_output_dir`
and `detect_run_dir`.

!!! note "When to use it"
    - You have finished `peakatail run` and want to identify PAS that are
      differentially used between cell types or conditions.
    - You want to rank genes by how strongly their 3' isoform choice differs
      between two cluster populations.
    - You are feeding results into `peakatail switch geneview` to visualise
      per-cluster PAS distributions for the top hits.

!!! warning "When NOT to use it"
    - You have not yet run `peakatail run` — you need `clusters.h5ad` first.
    - You want to quantify the absolute level of 3' UTR shortening per cluster,
      not the pairwise difference. Use `peakatail switch length` for that.
    - You have more than ~20 clusters and want to test all pairwise combinations
      with a strategy other than `fisher`. NB-based strategies scale as O(n_pairs);
      use `--cluster-pairs` to limit to biologically meaningful contrasts.

## Quick example

```bash
uv run peakatail switch diff \
  --h5ad peakatail_runs/emaout_2026-05-11_120000/per_dataset/sample1/clusters.h5ad \
  --strategy fisher \
  --fdr 0.05 \
  --min-cells-per-group 10
```

What lands on disk after this command (inside the originating run dir):

- `peakatail_runs/emaout_.../switch_diff_<ts>/differential/fisher_0_vs_1.tsv` — per-pair result TSV with augmented schema (see Output files).
- `peakatail_runs/emaout_.../switch_diff_<ts>/markers.tsv` — top marker PAS per cluster used for pre-filtering (only when `--marker-top-n > 0`; the default 0 writes no markers file).
- `peakatail_runs/emaout_.../switch_diff_<ts>/peakatail_<ts>.log` — run log.
- Volcano plot figures in `switch_diff_<ts>/figures/` (when plotting is enabled).

## Full `--help` output

```text
Usage: peakatail switch diff [OPTIONS]

  Differential APA test (Fisher / NB regression) across cluster pairs.

Options:
  --list-strategies               Print available diff strategies and exit.
  --threads INTEGER               Max parallel workers (auto-detected if not
                                  set). Respected by ResourceManager as an
                                  absolute ceiling.
  -v, --verbose                   Increase verbosity. -v = DEBUG for ema.*;
                                  -vv = DEBUG everywhere.
  -q, --quiet                     WARNING and up only. Overrides --verbose.
  --log-level TEXT                Explicit logger level (DEBUG/INFO/WARNING
                                  /ERROR) or `logger.name=LEVEL` (repeatable:
                                  comma-separated).
  --no-log-file                   Don't write peakatail_<ts>.log next to the
                                  outputs.
  --no-progress                   Suppress Rich progress bars.
  -c, --config PATH               YAML config; CLI flags override individual
                                  keys.
  -o, --output PATH               Output directory (timestamp suffix added
                                  automatically).
  --plot-engine TEXT              Engines: 'matplotlib' (default), 'plotly',
                                  'both', 'none', or comma list.
  --plot-format TEXT              Restrict output formats. Default 'all' =
                                  png+svg+html as appropriate.
  --no-plots                      Disable all plotting (alias for --plot-
                                  engine none).
  -i, --h5ad PATH                 Per-dataset clusters.h5ad. Repeat for
                                  multi-dataset.  [required]
  --pasbed PATH                   Optional PAS BED for context.
  --gtf PATH
  --cluster-pairs TEXT            `c1,c2;c3,c4` — limit to specific pairs.
  --cluster-key TEXT              [default: leiden]
  --marker-top-n INTEGER          Pre-filter the tested PAS to the union of the
                                  top-N marker PAS per cluster. 0 (default)
                                  disables pre-selection and is the only FDR-
                                  controlled setting (issue #94): markers are
                                  ranked with the SAME cluster labels the
                                  differential test then contrasts, so any non-
                                  zero value double-dips on the labels, making
                                  every strategy anti-conservative. The
                                  restriction no longer changes the within-gene
                                  Fisher denominator (that is computed from the
                                  full matrix), but it still selects what is
                                  tested. Speed-only; not a statistical filter.
                                  For speed WITHOUT the double-dip use
                                  --prefilter-min-cells instead.  [default: 0]
  --prefilter-min-cells INTEGER   LABEL-INDEPENDENT speed pre-filter (issue
                                  #94): test only the PAS detected (count > 0)
                                  in at least N cells, counted over ALL cells
                                  POOLED. 0 (default) disables it, leaving
                                  behaviour unchanged. This is the safe
                                  alternative to --marker-top-n: the criterion
                                  never looks at --cluster-key, so the PAS kept
                                  are identical under any permutation of the
                                  group labels and the null stays calibrated.
                                  Like --marker-top-n it gates only WHICH PAS
                                  are tested -- the within-gene Fisher
                                  denominator still comes from the full matrix.
                                  [default: 0]
  --marker-method TEXT            [default: wilcoxon]
  -s, --strategy TEXT             Differential APA strategy (run --list-
                                  strategies to see).  [default: fisher]
  --fdr FLOAT                     [default: 0.05]
  --per-worker-mb INTEGER         [default: 300]
  --min-cells-per-group INTEGER   Minimum cells per group for a PAS to enter
                                  differential testing.  [default: 10]
  --log2fc-thresh FLOAT           log2 fold-change threshold drawn on the
                                  volcano plot. Default 1.0.  [default: 1.0]
  --help                          Show this message and exit.
```

## Flags

### Inputs

| Flag | Type | Default | Description |
|---|---|---|---|
| `--h5ad` / `-i` | PATH (repeatable) | — | One or more `clusters.h5ad` files from `peakatail run`. Repeating this flag accumulates all h5ads into one run; all-pairs testing is performed within each h5ad independently, then results are merged. Required. |
| `--pasbed` | PATH | — | Optional PAS BED file. When present it is used to auto-discover coordinate columns for the output TSV (chrom, start, end, strand). If not given, the runner looks for `pasbed.bed` next to each `--h5ad` file (walking up to 4 parent directories). |
| `--gtf` | PATH | — | Optional GTF annotation file. Currently passed through to the runner but not used by the `fisher` strategy. Accepted for forward compatibility. |
| `--cluster-key` | TEXT | `leiden` | The `adata.obs` column containing cluster labels. Change this when using `--external-clusters` in `peakatail run` or a custom labelling scheme. |
| `--cluster-pairs` | TEXT | — | Restrict testing to specific cluster pairs. Format: `c1,c2;c3,c4` (semicolon-separated pairs, comma-separated within each pair). When omitted, all pairwise combinations are tested. |

### Strategy options

| Flag | Type | Default | Description |
|---|---|---|---|
| `--strategy` / `-s` | TEXT | `fisher` | Differential APA strategy. Run `peakatail switch diff --list-strategies` to see registered names. `fisher` applies a within-gene Fisher exact test (see Within-gene Fisher framing below). |
| `--marker-top-n` | INT | `0` (disabled) | **Speed shortcut, not a statistical filter — leave it at 0.** Pre-filters the PAS matrix to the union of the top-N marker PAS per cluster before differential testing. Any non-zero value ranks those markers with the **same cluster labels** the differential test then contrasts (a label double-dip), so the reported q-values are **not FDR-calibrated** — see [Why `--marker-top-n` defaults to 0](#why---marker-top-n-defaults-to-0). When set, the markers TSV is saved to `markers.tsv` for inspection. |
| `--prefilter-min-cells` | INT | `0` (disabled) | **Label-independent speed knob — the safe alternative to `--marker-top-n`.** Tests only the PAS detected (count > 0) in at least N cells, counted over **all cells pooled**. The criterion never looks at `--cluster-key`, so the PAS kept are identical under any permutation of the group labels and the null stays calibrated — see [Cutting the tested PAS set without a double-dip](#cutting-the-tested-pas-set-without-a-double-dip). Like `--marker-top-n` it gates only *which* PAS are tested; the within-gene Fisher denominator still comes from the full matrix. Source: `ema/switch_test/prefilter.py`. |
| `--marker-method` | TEXT | `wilcoxon` | Marker ranking method passed to `scanpy.tl.rank_genes_groups`. Options include `wilcoxon`, `t-test`, `logreg`. |
| `--min-cells-per-group` | INT | 10 | Minimum number of cells (with non-zero counts for NB strategies) in each cluster group for a PAS to be included in differential testing. PAS failing this filter in either cluster of a pair are dropped. Source: `ema/cli/config_schema.py`, `ema/switch_test/runner.py::run_diff`. |

### Filtering

| Flag | Type | Default | Description |
|---|---|---|---|
| `--fdr` | FLOAT | 0.05 | Benjamini–Hochberg FDR threshold. Rows with `qvalue < fdr` are considered significant. Used both to count significant hits in the log and to shade the volcano plot. |

### Output

| Flag | Type | Default | Description |
|---|---|---|---|
| `--output` / `-o` | PATH | `switch_out` | Output directory base name. Auto-routed inside the originating run dir when `--output` is left at default and `--h5ad` files come from a `peakatail_runs/` path. |
| `--log2fc-thresh` | FLOAT | 1.0 | log2 fold-change threshold drawn as vertical lines on the volcano plot. Does not filter the TSV output. |

### Diagnostics and performance

| Flag | Type | Default | Description |
|---|---|---|---|
| `--per-worker-mb` | INT | 300 | Estimated peak RAM per parallel worker in MB. Used by `ResourceManager` to cap outer parallelism: `n_outer = available_RAM / per_worker_mb`. Lower this to run more workers on memory-constrained machines; raise it if workers are crashing with OOM errors. |
| `--threads` | INT | auto | Absolute thread ceiling. See [Common flags](index.md#common-flags). |

## Why `--marker-top-n` defaults to 0

`--marker-top-n` used to default to `200`: before testing, the PAS matrix was
restricted to the union of the top-200 marker PAS per cluster, ranked by
`scanpy.tl.rank_genes_groups` on `--cluster-key`. That is a **label
double-dip** — the markers are chosen with the *same* labels the differential
test then contrasts, so the PAS that enter the test are exactly the ones that
already look cluster-associated by chance. Restricting the p-value *set* this
way alone put 17.4 % of null p-values below 0.05 (nominal 5 %).

Restricting the matrix used to also shrink the within-gene Fisher denominator,
because the "rest of the gene" background became the same label-selected
subset — the same PAS scored an `n_reads_gene` of 1,746 restricted vs 5,289
unrestricted, and only 629 of 6,453 p-values agreed between a marker-on and a
marker-off run. That half of issue #94 is **fixed**: `fisher` is now handed the
unrestricted matrix for the denominator, so `--marker-top-n N` changes only
*which* PAS are tested and reported, and each reported p-value is bit-identical
to the one the unrestricted run produces.

> **This holds for every strategy, but it did not always.** `fisher` takes the
> unrestricted matrix for its within-gene denominator. `nb_pairwise` and
> `nb_multi` additionally derive a per-cell library-size offset
> (`offset = log(total counts in that cell)`); before 0.3.0 they computed it
> from whatever matrix they were handed, so a pre-selection silently changed
> every cell's offset and therefore every fitted coefficient, dispersion and
> p-value — by up to two orders of magnitude on a small fixture. The offset is
> now always taken from the full matrix, because a cell's sequencing depth
> cannot depend on which hypotheses you chose to test.

The same holds under
`--isoform-agg within_utr` / `between_utr`: each group's background is built
from **all** of its PAS, and the selection only decides which rows are
reported. The numbers in the table below were
measured before that fix; the label double-dip they are driven by is unchanged.

Measured on a correctly-keyed matrix under a 20-run **label-permutation null**
(cluster labels shuffled, so there is nothing true to find; issue #94):

| Configuration | Null p < 0.05 | Runs with a q < 0.05 "hit" |
|---|---|---|
| `fisher --count-mode reads --marker-top-n 200` | 20.3 % | 20 / 20 |
| `fisher --count-mode cells --marker-top-n 200` | 13.0 % | 19 / 20 |
| `nb_pairwise --marker-top-n 200` | 24.7 % | 19 / 20 |
| `fisher --count-mode cells --marker-top-n 0` | **3.0 %** | **0 / 20** |

Only `--marker-top-n 0` controls the FDR, so it is now the default: a flagless
`peakatail switch diff` is calibrated. Any non-zero value still works but logs a loud
warning — use it as a **speed shortcut / ranking screen** on large datasets
(NB strategies scale badly in the number of PAS), never as evidence of
significance. If you need both speed and calibration, cut the search space
with something independent of the labels instead — `--prefilter-min-cells N`
(see [Cutting the tested PAS set without a double-dip](#cutting-the-tested-pas-set-without-a-double-dip)),
`--cluster-pairs`, `--min-cells-per-group`, or a PAS list from a separate
dataset.

## Cutting the tested PAS set without a double-dip

`--marker-top-n 0` is calibrated but tests everything, which is slow on a large
matrix (the NB strategies scale badly in the number of PAS). `--prefilter-min-cells N`
is the speed knob to reach for instead:

```bash
peakatail switch diff -i clusters.h5ad --cluster-key celltype --prefilter-min-cells 25
```

It keeps only the PAS detected (count > 0) in at least `N` cells, counted over
**all cells pooled** — the labels are never grouped, split or read. That is
structural, not a promise: the criterion is computed by
`ema.switch_test.prefilter.select_expressed_pas(count_matrix, min_cells)`, whose
signature has no `cluster_key`, no label vector and no AnnData parameter, so
there is nothing for it to double-dip on. `run_diff` calls it before the label
vector is even built.

Consequences, measured on the same label-permutation null the marker flag was
measured on (2000 PAS, 120 cells, `fisher --count-mode cells`, 5 permutations):

| Configuration | PAS tested | Null p < 0.05 |
|---|---|---|
| unfiltered (`--marker-top-n 0`) | 2000 | 3.6 % |
| `--prefilter-min-cells 80` | 469 | **3.5 %** |
| `--marker-top-n 200` | 400 | 16.6 % |

Two further properties hold exactly, not approximately, and are pinned by
`tests/test_label_independent_prefilter_i94.py`:

* the pre-filtered PAS set is **identical under every permutation of the
  labels** (the marker set is different every time), and
* every surviving p-value is **bit-identical** to the one the unfiltered run
  reports for that PAS — the pre-filter removes hypotheses, it never changes a
  test. The q-values are then plain Benjamini–Hochberg over that smaller,
  label-blind set of hypotheses.

  For `nb_pairwise` / `nb_multi` this holds only because their library-size
  offset is taken from the full matrix (see the note under `--marker-top-n`
  above); before 0.3.0 it did not hold for them.

**Choosing N.** The criterion counts *cells*, not reads, because the per-cell
contingency table (`--count-mode cells`) is built from exactly that number: a
PAS detected in fewer than `--min-cells-per-group` cells cannot populate a
usable table in either group anyway. A read-total threshold would mostly rank
sequencing depth, and a variance/abundance threshold starts to correlate with
the between-group difference being tested even without reading the labels. A
safe starting point is roughly `2 × --min-cells-per-group`; raise it until the
run is fast enough, and note that (as with any independent filter) testing
fewer hypotheses makes the BH threshold less stringent — that is the intended
multiplicity saving, not selection bias.

## `nb_pairwise` and the dispersion floor

`nb_pairwise` clips its per-PAS Negative-Binomial dispersion to `[1e-4, 10]`.
A PAS that lands on the **lower** clip is fitted as a Poisson GLM, whose Wald
standard error is a lower bound — so its p-value is anti-conservative. In the
same label-permutation null as above, 8.5 % of null tests hit that floor and
those tests produced **67 % of `nb_pairwise`'s false `q < 0.05` hits**.

Such rows are now marked `dispersion_floored = True` in the per-pair TSV and
their `qvalue` column is left **empty (`NaN`)**, so a `qvalue < fdr` filter can
never call them significant. The row and its raw `pvalue` are still written —
read that p-value as a *lower bound*, not as an error rate.

!!! warning "nb_pairwise q-values need permutation calibration"

    Withholding the floored rows removes the dominant source of
    anti-conservatism, but it does **not** make the remaining `nb_pairwise`
    q-values calibrated (`nb_pairwise --marker-top-n 0` still put 5.1 % of null
    p-values below 0.05 and produced a hit in 20/20 permutation runs, the
    dispersion-floor tail being the bulk of it). If a `nb_pairwise` q-value has
    to carry an error-rate claim, calibrate it against a label-permutation null
    of your own data. `fisher --count-mode cells --marker-top-n 0` is the
    configuration measured to control the FDR out of the box.

## Within-gene Fisher framing

As of commit `f5ed80d`, the `fisher` strategy uses a **within-gene** framing
when `adata.var["gene_id"]` is present. For each gene, it groups all PAS
belonging to that gene and tests whether the read distribution across PAS
differs between cluster 1 and cluster 2 using a Fisher exact test on the
contingency table:

```
           | PAS_1  | PAS_2  | ... | PAS_N  |
cluster_1  |  r_11  |  r_12  | ... |  r_1N  |
cluster_2  |  r_21  |  r_22  | ... |  r_2N  |
```

This within-gene design tests for differential usage of a gene's own PAS
repertoire, rather than comparing a PAS against all other PAS genome-wide. It
is more APA-appropriate because it removes confounding from differential gene
expression.

When `adata.var` has no `gene_id` column (older h5ad files), the strategy logs
a warning and falls back to the global (cross-gene) path.

## Output files

Output is written to `<out_dir>/differential/` (created automatically).

**`differential/<strategy>_<c1>_vs_<c2>.tsv`**

One TSV per cluster pair. Columns (in order) for the default `--strategy fisher`:

| Column | Type | Description |
|---|---|---|
| `pas_id` | str | PAS identifier matching `adata.var_names`. Under `--isoform-agg between_utr` the row unit is a 3'UTR isoform and this reads `GENE::TRANSCRIPT`. |
| `gene_id` | str | Gene annotation from `adata.var["gene_id"]` (empty string if unavailable). Under `--isoform-agg between_utr` it is the group's gene. |
| `chrom` | str | Chromosome from `pasbed.bed` (empty if pasbed not found). **Omitted under `--isoform-agg between_utr`.** |
| `start` | str | Genomic start position (0-based). **Omitted under `--isoform-agg between_utr`.** |
| `end` | str | Genomic end position. **Omitted under `--isoform-agg between_utr`.** |
| `strand` | str | `+` or `-`. **Omitted under `--isoform-agg between_utr`.** |
| `cluster1` | str | First cluster label of this pair. |
| `cluster2` | str | Second cluster label of this pair. |
| `pvalue` | float | Raw two-sided Fisher exact p-value for this PAS. |
| `qvalue` | float | Benjamini–Hochberg adjusted p-value (FDR) across all PAS tested in this pair. |
| `n_cells` | int | Cells in the pair (`n_cells_cluster1 + n_cells_cluster2`). |
| `n_cells_cluster1` | int | Cells carrying the `cluster1` label. |
| `n_cells_cluster2` | int | Cells carrying the `cluster2` label. |
| `n_cells_expr_cluster1` | int | Cells of cluster 1 with ≥1 read at this PAS. |
| `n_cells_expr_cluster2` | int | Cells of cluster 2 with ≥1 read at this PAS. |
| `n_reads_pas_cluster1` | int | Reads at this PAS in cluster 1. |
| `n_reads_pas_cluster2` | int | Reads at this PAS in cluster 2. |
| `n_reads_gene_cluster1` | int | Total reads for this gene in cluster 1 (fisher within-gene framing). Summed over **all** PAS of the gene, including any excluded by `--marker-top-n`. |
| `n_reads_gene_cluster2` | int | Total reads for this gene in cluster 2, on the same basis. |
| `odds_ratio` | float | Odds ratio of the 2×2 table (this PAS vs the gene's other PAS, cluster 1 vs cluster 2), in the unit chosen by `--count-mode` — cells by default, reads under `--count-mode reads`. |
| `delta_proportion` | float | `prop(cluster1) - prop(cluster2)` of the within-gene usage proportion. **Positive ⇒ the PAS is used more in `cluster1`.** |
| `log2fc` | float | `log2(prop(cluster2) / prop(cluster1))` of those same proportions (plus a small `eps` so an empty cluster stays finite). **Positive ⇒ the PAS is used more in `cluster2`.** |

!!! warning "`delta_proportion` and `log2fc` use opposite sign conventions"
    `delta_proportion` is `prop(cluster1) - prop(cluster2)` while `log2fc` is
    `log2(prop(cluster2) / prop(cluster1))`, so on the same row the two
    normally carry **opposite signs**: `delta_proportion > 0` means the PAS is
    used more in `cluster1`, whereas `log2fc > 0` means it is used more in
    `cluster2`. Filter on one of them, never on both with the same inequality
    — and note the volcano plot's x-axis is `log2fc`, i.e. cluster2-positive.

The statistical columns are strategy-specific. `--strategy nb_pairwise` writes
`pvalue`, `qvalue`, `log2fc`, `dispersion`, `dispersion_floored`,
`n_cells`, `test_stat`; the
`nb_multi` omnibus (written to `<strategy>_omnibus.tsv`, not to a per-pair
file) writes `pvalue`, `qvalue`, `test_stat`, `df`, `dispersion`, `n_cells`.

The augmented column order (pas_id, gene_id, chrom, start, end, strand, cluster1, cluster2, then statistical columns) is produced by the `_augment_diff_df` helper in `ema/switch_test/runner.py`.

The UTR-scoped values of `--isoform-agg` add a `diff_group_id` column naming the
group each row was tested within, and `between_utr` changes the row unit:

| `--isoform-agg` | Row unit | `pas_id` | `gene_id` | `diff_group_id` | `chrom`/`start`/`end`/`strand` |
|---|---|---|---|---|---|
| `per_gene` (default) | PAS | PAS id | gene of the PAS | *(column absent)* | present |
| `within_utr` (alias `per_isoform`) | PAS | PAS id | gene of the PAS | `GENE::TRANSCRIPT` of the tested 3'UTR (`GENE::_gene_` for the `--utr-unmatched gene` fallback bucket) | present |
| `between_utr` | 3'UTR isoform | `GENE::TRANSCRIPT` | the gene (same value as `diff_group_id`) | the gene | **absent** |

Under `between_utr` a row is a whole 3'UTR, not a single cleavage site, so there
is no one position to report: the four coordinate columns are **omitted from the
TSV entirely** rather than written as empty strings ([issue #110](https://github.com/BMGLab/PeakATail/issues/110)).
A join keyed on them therefore fails with a missing-column error instead of
silently matching nothing. `gene_id` *is* populated for these rows, so joining
`between_utr` output to `per_gene` output on `gene_id` works as expected.

**`markers.tsv`**

Written only when `--marker-top-n > 0` (not at the default 0). Two-column TSV: `cluster` and `pas_id`. Lists the top-N marker PAS per cluster used as pre-filter for differential testing.

**`figures/volcano_<c1>_vs_<c2>.*`**

Volcano plot (log2FC vs -log10 qvalue) per cluster pair. Written by
`ema.viz.pipeline_hooks::render_switch_diff_outputs`. Format depends on
`--plot-engine` and `--plot-format`.

## How it relates to other commands

- **[`peakatail run`](run.md)** — produces the `clusters.h5ad` and `pasbed.bed` inputs.
- **[`peakatail switch geneview`](switch-geneview.md)** — consumes the `differential/*.tsv` files via `--diff-tsv` to auto-rank genes for per-cluster track plots.
- **[`peakatail switch length`](switch-length.md)** — complementary quantification; results can be overlaid in `geneview`.

## See also

- Strategy details in [`../strategies/`](../strategies/) — how Fisher and NB regression are implemented.
- Tutorial in [`../tutorials/`](../tutorials/) — end-to-end differential APA walkthrough.
