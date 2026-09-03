# `ema switch diff`

`ema switch diff` tests for differential alternative polyadenylation (APA)
between every pair of Leiden clusters in one or more `clusters.h5ad` files
produced by `ema run`. For each cluster pair (c1, c2) and each PAS that passes
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
    - You have finished `ema run` and want to identify PAS that are
      differentially used between cell types or conditions.
    - You want to rank genes by how strongly their 3' isoform choice differs
      between two cluster populations.
    - You are feeding results into `ema switch geneview` to visualise
      per-cluster PAS distributions for the top hits.

!!! warning "When NOT to use it"
    - You have not yet run `ema run` — you need `clusters.h5ad` first.
    - You want to quantify the absolute level of 3' UTR shortening per cluster,
      not the pairwise difference. Use `ema switch length` for that.
    - You have more than ~20 clusters and want to test all pairwise combinations
      with a strategy other than `fisher`. NB-based strategies scale as O(n_pairs);
      use `--cluster-pairs` to limit to biologically meaningful contrasts.

## Quick example

```bash
uv run ema switch diff \
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
Usage: ema switch diff [OPTIONS]

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
| `--h5ad` / `-i` | PATH (repeatable) | — | One or more `clusters.h5ad` files from `ema run`. Repeating this flag accumulates all h5ads into one run; all-pairs testing is performed within each h5ad independently, then results are merged. Required. |
| `--pasbed` | PATH | — | Optional PAS BED file. When present it is used to auto-discover coordinate columns for the output TSV (chrom, start, end, strand). If not given, the runner looks for `pasbed.bed` next to each `--h5ad` file (walking up to 4 parent directories). |
| `--gtf` | PATH | — | Optional GTF annotation file. Currently passed through to the runner but not used by the `fisher` strategy. Accepted for forward compatibility. |
| `--cluster-key` | TEXT | `leiden` | The `adata.obs` column containing cluster labels. Change this when using `--external-clusters` in `ema run` or a custom labelling scheme. |
| `--cluster-pairs` | TEXT | — | Restrict testing to specific cluster pairs. Format: `c1,c2;c3,c4` (semicolon-separated pairs, comma-separated within each pair). When omitted, all pairwise combinations are tested. |

### Strategy options

| Flag | Type | Default | Description |
|---|---|---|---|
| `--strategy` / `-s` | TEXT | `fisher` | Differential APA strategy. Run `ema switch diff --list-strategies` to see registered names. `fisher` applies a within-gene Fisher exact test (see Within-gene Fisher framing below). |
| `--marker-top-n` | INT | `0` (disabled) | **Speed shortcut, not a statistical filter — leave it at 0.** Pre-filters the PAS matrix to the union of the top-N marker PAS per cluster before differential testing. Any non-zero value ranks those markers with the **same cluster labels** the differential test then contrasts (a label double-dip), so the reported q-values are **not FDR-calibrated** — see [Why `--marker-top-n` defaults to 0](#why---marker-top-n-defaults-to-0). When set, the markers TSV is saved to `markers.tsv` for inspection. |
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
to the one the unrestricted run produces. The numbers in the table below were
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
`ema switch diff` is calibrated. Any non-zero value still works but logs a loud
warning — use it as a **speed shortcut / ranking screen** on large datasets
(NB strategies scale badly in the number of PAS), never as evidence of
significance. If you need both speed and calibration, cut the search space
with something independent of the labels instead (e.g. `--cluster-pairs`,
`--min-cells-per-group`, or a PAS list from a separate dataset).

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

One TSV per cluster pair. Columns (in order):

| Column | Type | Description |
|---|---|---|
| `pas_id` | str | PAS identifier matching `adata.var_names`. |
| `gene_id` | str | Gene annotation from `adata.var["gene_id"]` (empty string if unavailable). |
| `chrom` | str | Chromosome from `pasbed.bed` (empty if pasbed not found). |
| `start` | str | Genomic start position (0-based). |
| `end` | str | Genomic end position. |
| `strand` | str | `+` or `-`. |
| `cluster1` | str | First cluster label of this pair. |
| `cluster2` | str | Second cluster label of this pair. |
| `n_reads_gene_cluster1` | int | Total reads for this gene in cluster 1 (fisher within-gene framing). Summed over **all** PAS of the gene, including any excluded by `--marker-top-n`. |
| `n_reads_gene_cluster2` | int | Total reads for this gene in cluster 2, on the same basis. |
| `statistic` | float | Test statistic (odds ratio for Fisher). |
| `pvalue` | float | Raw p-value. |
| `qvalue` | float | Benjamini–Hochberg adjusted p-value (FDR). |

The augmented column order (pas_id, gene_id, chrom, start, end, strand, cluster1, cluster2, then statistical columns) is produced by the `_augment_diff_df` helper in `ema/switch_test/runner.py`. Source: lines 316–384.

**`markers.tsv`**

Written only when `--marker-top-n > 0` (not at the default 0). Two-column TSV: `cluster` and `pas_id`. Lists the top-N marker PAS per cluster used as pre-filter for differential testing.

**`figures/volcano_<c1>_vs_<c2>.*`**

Volcano plot (log2FC vs -log10 qvalue) per cluster pair. Written by
`ema.viz.pipeline_hooks::render_switch_diff_outputs`. Format depends on
`--plot-engine` and `--plot-format`.

## How it relates to other commands

- **[`ema run`](run.md)** — produces the `clusters.h5ad` and `pasbed.bed` inputs.
- **[`ema switch geneview`](switch-geneview.md)** — consumes the `differential/*.tsv` files via `--diff-tsv` to auto-rank genes for per-cluster track plots.
- **[`ema switch length`](switch-length.md)** — complementary quantification; results can be overlaid in `geneview`.

## See also

- Strategy details in [`../strategies/`](../strategies/) — how Fisher and NB regression are implemented.
- Tutorial in [`../tutorials/`](../tutorials/) — end-to-end differential APA walkthrough.
