# Output files

Every file PeakATail writes is listed here with its canonical path, format, and the pipeline stage that produces it. The paths use `<run>/` as shorthand for the timestamped output directory (e.g. `peakatail_runs/emaout_2026-05-11_152746/`).

---

## Run root

### `run_config.json`

**Path:** `<run>/run_config.json`  
**Produced by:** `OutputManager.save_run_config()` at the start of every pipeline run.  
**Format:** JSON object.

Records the full resolved configuration (YAML values + CLI overrides) plus a timestamp. Useful for reproducing a run or understanding what parameters were active.

### `peakatail_<timestamp>.log`

**Path:** `<run>/peakatail_<timestamp>.log`  
**Produced by:** the logging configuration in `ema/logging_config.py`.  
**Format:** plain text, structured log lines.

Full DEBUG-level log of the run. Includes per-stage timings, cell/PAS counts at each filter step, and any warnings about missing data. Rotate or delete after inspection; large BAMs can produce multi-MB logs.

### `resources.jsonl`

**Path:** `<run>/resources.jsonl`  
**Produced by:** `_ResourceSampler` thread in `ema/main.py`, started when `--plot-engines` is set.  
**Format:** newline-delimited JSON, one record per 5-second sample.

Each record: `{"elapsed_s": float, "rss_gb": float, "cpu_pct": float}`. Read by the `resource_timeline` visualisation strategy.

---

## Stage stats JSONs

One stats file per pipeline stage, under numbered subdirectories.

| File | Stage | Key fields |
|------|-------|------------|
| `01_peak_calling/peak_calling_stats.json` | Peak calling | `strategy`, `lambda_fold_change`, `lambda_window` |
| `02_cb_filter/cb_filter_stats.json` | Cell barcode filter | `min_read` |
| `03_gtf_annotation/gtf_annotation_stats.json` | GTF annotation | `gtf_path`, `utr_count` |
| `04_pas_gene_assignment/pas_gene_stats.json` | PAS-gene assignment | `max_gene_distance`, `assigned_pas_count` |
| `05_annotated_matrix/annotated_stats.json` | Annotated matrix | `input_pas_count`, `annotated_pas_count`, `cell_count` |
| `06_preprocessing/preprocessing_stats.json` | Preprocessing | `cells_after_filter`, `pas_after_filter` |
| `07_clustering/clustering_stats.json` | Clustering | `final_cells`, `final_pas` |

---

## Per-dataset directory: `per_dataset/<id>/`

For each dataset ID listed in `datasets:` (or `"only"` for a single-BAM run), PeakATail writes these files.

### `raw/pos.bed`, `raw/neg.bed`, `raw/pas.bed`

**Produced by:** `write_raw_peak_outputs()` in `ema/outputs.py`.  
**Format:** BED6.

Raw peak-calling output before any filtering. `pos.bed` = positive-strand PAS; `neg.bed` = negative-strand PAS; `pas.bed` = union of both strands. These are the unfiltered call set; downstream tools that need the broadest possible PAS set should read from here.

### `raw/pos.mtx`, `raw/neg.mtx`

**Produced by:** `write_raw_peak_outputs()`.  
**Format:** MatrixMarket sparse integer matrix.

Raw count matrices before cell-barcode or PAS filtering. Rows are PAS (matching the BED row order), columns are cell barcodes (matching `raw/cb.tsv`).

### `raw/cb.tsv`

**Produced by:** `write_raw_peak_outputs()`.  
**Format:** plain text, one barcode per line.

Cell barcodes in the order corresponding to the columns of `raw/pos.mtx` and `raw/neg.mtx`.

### `posbed.bed`, `negbed.bed`, `pasbed.bed`

**Produced by:** `write_per_dataset_beds()` in `ema/outputs.py`.  
**Format:** BED6. The 4th column is the integer PAS ID; the 5th column is always 0 (unused score field); the 6th column is the strand.

The canonical post-merge BED files for this dataset. `pasbed.bed` is the union of positive and negative strand PAS and is the primary coordinate reference for all downstream analysis. `peakatail switch diff` reads it to annotate result rows with `chrom`, `start`, `end`, and `strand`.

### `filtered_cb.tsv`

**Produced by:** `write_filtered_cb()` in `ema/outputs.py`.  
**Format:** TSV, header `barcode\tmin_read=<n>`.

Cell barcodes that passed the `min_read` threshold. One barcode per line after the header. Pass this to other tools (e.g. `samtools view -D CB:`) to subset a BAM to the same cells.

### `annotatedpas.bed`

**Produced by:** `write_pas_gene_artifacts()` in `ema/outputs.py`.  
**Format:** BED6 + gene_id column (7-column tab-delimited).

`pasbed.bed` extended with a trailing gene_id column. Useful for intersecting PAS coordinates with gene models in external tools.

### `pas_gene.tsv`

**Produced by:** `write_pas_gene_artifacts()` in `ema/outputs.py`.  
**Format:** TSV, two columns: `pas_id` and `gene_id`.

Explicit mapping from integer PAS ID to Ensembl gene ID. Produced by `find_close()`, which searches for the nearest gene 3' end within `max_gene_distance` bp (default 5000 bp).

### `annotated_matrix.mtx`

**Produced by:** `write_annotated_matrix()` in `ema/outputs.py`.  
**Format:** MatrixMarket sparse integer matrix.

PAS-by-cell count matrix after gene annotation but before any cell or PAS filtering. Rows correspond to `annotated_pas_ids.tsv`; columns correspond to `annotated_cells.tsv`.

### `annotated_pas_ids.tsv`

**Produced by:** `write_annotated_matrix()`.  
**Format:** single-column TSV, header `pas_id`.

Row index for `annotated_matrix.mtx`. Integer PAS IDs that have a gene assignment.

### `annotated_cells.tsv`

**Produced by:** `write_annotated_matrix()`.  
**Format:** single-column TSV, header `barcode`.

Column index for `annotated_matrix.mtx`. Cell barcodes in matrix column order (before any cell-level filtering).

### `preprocessed.h5ad`

**Produced by:** `write_preprocessed_h5ad()` in `ema/outputs.py`.  
**Format:** HDF5/AnnData.

The count matrix after applying `min_cells` (PAS filter) and `min_pas_per_cell` (cell filter), but before clustering. Contains `adata.var['gene_id']` for each PAS. Useful for inspecting the filter step in isolation or loading into Scanpy for alternative analyses.

### `clusters.h5ad`

**Produced by:** `clustering()` in `ema/clustering/clustering.py`.  
**Format:** HDF5/AnnData.

The primary output of `peakatail run`. Contains:

- `adata.obs['leiden']` — cluster label for each cell (string).
- `adata.var['gene_id']` — gene assignment for each PAS.
- `adata.obsm['X_pca']` or `adata.obsm['X_lsi']` — dimensionality-reduced embedding.
- `adata.obsm['X_umap']` — UMAP coordinates.

All `peakatail switch` subcommands take this file as their primary input via `--h5ad`.

---

## Peakcalling directory: `peakcalling/`

One pair of files per BAM per direction:

| File | Contents |
|------|----------|
| `<id>_<n>.pos.bed` | Positive-strand PAS calls for dataset `<id>`, BAM index `<n>` |
| `<id>_<n>.neg.bed` | Negative-strand PAS calls |
| `<id>_<n>.pos.mtx` | Positive-strand count matrix |
| `<id>_<n>.neg.mtx` | Negative-strand count matrix |
| `<id>_<n>.cb.tsv` | Unified cell barcode list for this BAM |
| `<id>_<n>.peak_jobs.json` | Per-(contig, strand) peak-calling worker timings, PAS/cell counts and peak RSS. Only written when peak calling ran in parallel (see `--peak-workers`). |

These are intermediate files. The canonical per-dataset view is under `per_dataset/<id>/`.

---

## Unified directory: `unified/` (multi-sample only)

| File | Contents |
|------|----------|
| `unified.bed` | Merged PAS coordinates (no atlas) |
| `atlas_snapped.bed` | Atlas-snapped PAS coordinates |
| `atlas_mapping.tsv` | Called PAS → atlas PAS mapping (3 columns: `called_pas_id`, `atlas_id`, `distance_bp`) |
| `atlas_pas_id_lookup.tsv` | Atlas ID → unified integer PAS ID |
| `concatenated.mtx` | Full joint count matrix across all datasets |
| `concatenated_cbs.tsv` | Joint cell barcode list (barcodes prefixed with `<dataset_id>_`) |

---

## Cross-dataset directory: `cross_dataset/` (multi-sample only)

| File | Contents |
|------|----------|
| `canonical_cluster_map.tsv` | Maps each (dataset_id, cluster_label) pair to a canonical cluster label |

---

## GTF cache: `gtf_cache/`

| File | Contents |
|------|----------|
| `gene_end.bed` | 3' gene boundaries extracted from GTF |
| `raw_feature.tsv` | All features from the GTF used for annotation |
| `utr_lengths.tsv` | Per-gene 3'UTR length table |
| `utr_lengths.json` | Same data as JSON (used by Python code) |
| `manifest.json` | Cache metadata: GTF path, mtime, PeakATail version |

The cache is invalidated and rebuilt when the GTF path or its modification time changes.

---

## Figures directory: `figures/` (run root)

Produced by `ema/viz/pipeline_hooks.py::render_run_outputs()` after the pipeline completes. Only generated when `--plot-engines` is set (default: matplotlib only).

| File | Strategy | Contents |
|------|----------|---------|
| `umap_<dataset>.png/.svg/.html` | `umap_matplotlib` / `umap_plotly` | UMAP coloured by Leiden cluster |
| `clusters_<dataset>.png/.svg/.html` | `cluster_sizes_matplotlib` | Bar chart of cells per cluster |
| `peak_qc_<dataset>.png/.svg/.html` | `peak_qc_matplotlib` | Four-panel QC: peaks per chrom, peak widths, PAS per cell, reads per cell |
| `resource_timeline.png/.svg/.html` | `resource_timeline_matplotlib` | RAM and CPU over pipeline wall time |
| `tile_timing.png/.svg/.html` | `tile_timing_matplotlib` | Per-tile timing when `--tiles` was used |
| `run_report.html` | `run_report` | Standalone HTML with all figures embedded |

Every figure has a `.meta.json` sidecar with metadata:

```json
{
  "figure_name": "umap_default",
  "viz_strategy": "umap_matplotlib",
  "generated_at": "2026-05-11T15:32:00+00:00",
  "peakatail_version": "0.3.1",
  "n_cells": 1051,
  "n_clusters": 12
}
```

---

## Switch diff outputs: `switch_diff_<timestamp>/`

Produced by `peakatail switch diff`. The subdirectory is created inside the run directory of the input h5ad, keeping all analysis of a run in one place.

| File | Contents |
|------|----------|
| `markers.tsv` | **Only written when `--marker-top-n N` is passed with `N > 0`** (the default is `0`, so a flagless run produces no `markers.tsv`). Top-N marker PAS per cluster, used to subset testing. Columns: `cluster`, `pas_id`, `score`. Marker pre-selection ranks PAS with the same cluster labels the test then contrasts, so it is a speed shortcut, not calibrated inference — see [issue #94](https://github.com/BMGLab/PeakATail/issues/94) |
| `differential/<strategy>_<c1>_vs_<c2>.tsv` | Per-pair differential results. Columns: `pas_id`, `gene_id`, `chrom`, `start`, `end`, `strand`, `cluster1`, `cluster2`, then the strategy's statistical columns. Under `--isoform-agg within_utr`/`between_utr` a `diff_group_id` column names the group each row was tested within. See [`ema switch diff` — Output files](../cli/switch-diff.md#output-files) for the full schema and [Tutorial 04 — switch analysis](../tutorials/04-switch-analysis.md) for a walkthrough |
| `figures/volcano_<c1>_vs_<c2>.png/.svg/.html` | Volcano plot for one cluster pair |
| `figures/figures_INDEX.md` | Human-readable index of all figures with per-pair statistics |
| `figures/figures_INDEX.json` | Machine-readable version of the index |

Under `--isoform-agg between_utr` a row is a whole 3'UTR isoform rather than a
single PAS: `pas_id` reads `GENE::TRANSCRIPT`, `gene_id` holds that UTR's gene
(the same value as `diff_group_id`), and the four PAS-level coordinate columns
`chrom`, `start`, `end` and `strand` are **not written at all** — a whole UTR has
no single cleavage position, so omitting them makes a join keyed on them fail
loudly instead of matching nothing ([issue #110](https://github.com/BMGLab/PeakATail/issues/110)).
`within_utr` keeps PAS as the row unit and so keeps the full column set.

---

## Switch length outputs: `switch_length_<timestamp>/`

Produced by `peakatail switch length`. The output filename depends on the strategy class's `output_filename` attribute:

| Strategy | Output file | Contents |
|----------|-------------|----------|
| `classic` | `pdui_classic.tsv` | Per-cell PDUI (proximal/distal pair). Columns: `gene_id`, `transcript_id`, `cell`, `pdui`, `cluster` |
| `proportion` | `pdui_proportion.tsv` | Per-cell per-PAS proportion within gene. Columns: `gene_id`, `transcript_id`, `pas_id`, `rank`, `cell`, `proportion`, `cluster` |
| `shannon` | `entropy_shannon.tsv` | Shannon entropy of PAS usage distribution per cell per gene. Columns: `gene_id`, `cell`, `entropy`, `cluster` |

All strategy outputs include a `cluster` column (from `adata.obs['leiden']`) and coordinate columns (`chrom`, `start`, `end`, `strand`) looked up from `pasbed.bed`.

Figures:

| File | Contents |
|------|----------|
| `figures/pdui_distribution.png/.svg/.html` | Violin plot of the strategy metric per cluster |
| `figures/length_shifts.png/.svg/.html` | Bar chart of mean PDUI shift per gene across a selected cluster pair |

---

## Switch geneview outputs: `switch_geneview_<timestamp>/`

Produced by `peakatail switch geneview`.

| File | Contents |
|------|----------|
| `figures/gene_<id>.png/.svg` | Gene track figure for one gene |
| `figures/gene_<id>.meta.json` | Rendering metadata: gene coordinates, n_pas, n_clusters_rendered, top proportion per cluster |
| `figures/figures_INDEX.md` | Index of all rendered genes |
| `figures/figures_INDEX.json` | Machine-readable index |
