# Strategy reference

PeakATail is built around pluggable strategy registries. Every algorithm
choice — how cells are clustered, how APA usage is quantified, how
differential tests are run, how clusters are matched across datasets, and
which figures are produced — is a named strategy that can be swapped without
touching the rest of the pipeline. This page lists every registered strategy
and links to its full reference.

## Clustering

Clustering groups cells by their 3' UTR usage so that downstream APA tests
have biologically meaningful contrasts. See [clustering strategies](clustering.md).

| Strategy | CLI flag | When to use |
|---|---|---|
| `leiden_tfidf` | `--cluster-method leiden_tfidf` | Sparse PAS counts (10x Chromium, Drop-seq). Recommended default. |
| `leiden_libsize` | `--cluster-method leiden_libsize` | Dense gene-count-like data; when TF-IDF produces too many small clusters. |
| `external` | `--cluster-method external` | Import Seurat / Scanpy clusters produced from gene expression. |

## Quantification (PDUI)

Quantification converts the PAS count matrix into a per-cell APA summary
score. See [quantification strategies](quantification.md).

| Strategy | CLI flag | Output file | When to use |
|---|---|---|---|
| `classic` | `--pdui-method classic` | `pdui_classic.tsv` | Binary proximal/distal question on a per-gene or per-isoform basis. |
| `proportion` | `--pdui-method proportion` | `proportion.tsv` | Full per-PAS usage vector per cell; feeds `proportion_heatmap` viz. |
| `shannon` | `--pdui-method shannon` | `entropy_shannon.tsv` | Heterogeneity of PAS usage; genes with many PAS, or pseudotime analyses. |

## Differential APA

Differential testing identifies PAS that shift between two or more clusters.
See [diff strategies](diff.md).

| Strategy | CLI flag | Multi-condition | When to use |
|---|---|---|---|
| `fisher` | `--diff-method fisher` | No | Fast exploratory screen; within-gene framing (DEXSeq-style). |
| `nb_pairwise` | `--diff-method nb_pairwise` | No | Cell-level NB GLM; corrects for library size; no pseudo-replication. |
| `nb_multi` | `--diff-method nb_multi` | Yes | Omnibus test across all clusters at once. |

## Cross-dataset matching

Matching assigns canonical cluster IDs when the same experiment is processed
in multiple batches or from multiple samples. See [match strategies](match.md).

| Strategy | CLI flag | When to use |
|---|---|---|
| `marker_overlap` | `--match-method marker_overlap` | Default. Datasets share PAS features; reasonably similar protocols. |
| `mnn` | `--match-method mnn` | Low marker overlap; different capture technologies or library sizes. |
| `jaccard` | `--match-method jaccard` | Same physical cells in both runs (regression testing or split-BAM experiments). |

## Visualization

Every plot type is a registered strategy with optional engine variants
(matplotlib, plotly, scanpy). See [viz strategies](viz.md).

| Plot type | Registered names | What it shows |
|---|---|---|
| `umap` | `umap_matplotlib`, `umap_plotly`, `umap_scanpy` | 2-D UMAP colored by cluster. |
| `cluster_sizes` | `cluster_sizes_matplotlib`, `cluster_sizes_plotly` | Bar chart of cell counts per cluster. |
| `peak_qc` | `peak_qc_matplotlib`, `peak_qc_plotly` | 4-panel peak and cell QC distributions. |
| `volcano` | `volcano_matplotlib`, `volcano_plotly` | Differential APA volcano: log2fc vs -log10(q). |
| `pdui_distribution` | `pdui_distribution_matplotlib`, `pdui_distribution_plotly`, `pdui_distribution_scanpy` | Per-cluster PDUI violin plots. |
| `proportion_heatmap` | `proportion_heatmap_matplotlib`, `proportion_heatmap_plotly` | PAS x cluster mean proportion heatmap. |
| `entropy_distribution` | `entropy_distribution_matplotlib`, `entropy_distribution_plotly` | Per-cluster Shannon entropy violin plots. |
| `diff_agreement` | `diff_agreement_matplotlib`, `diff_agreement_plotly` | Jaccard overlap between significant PAS sets from different test strategies. |
| `length_shifts` | `length_shifts_matplotlib`, `length_shifts_plotly` | Gene x cluster-pair PDUI delta heatmap. |
| `pas_overlap` | `pas_overlap_matplotlib`, `pas_overlap_plotly` | UpSet / Venn of PAS detected across datasets. |
| `atlas_snap_diag` | `atlas_snap_diag_matplotlib`, `atlas_snap_diag_plotly` | Atlas-snap result: snapped vs unsnapped counts + distance histogram. |
| `cluster_match_sankey` | `cluster_match_sankey_matplotlib`, `cluster_match_sankey_plotly` | Alluvial chart of original → canonical cluster correspondence. |
| `match_confidence` | `match_confidence_matplotlib`, `match_confidence_plotly` | Heatmap of match confidence scores per (dataset, cluster) pair. |
| `tile_timing` | `tile_timing_matplotlib`, `tile_timing_plotly` | Bar chart of per-tile peak-calling wall time, outliers highlighted. |
| `resource_timeline` | `resource_timeline_matplotlib`, `resource_timeline_plotly` | Dual-axis RAM (GB) and CPU% over elapsed seconds. |
| `gene_track` | `gene_track_matplotlib`, `gene_track_plotly` | Per-gene: isoform structure + per-cluster PAS coverage bars. |
