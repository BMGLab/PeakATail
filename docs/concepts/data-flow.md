# Data flow

This page shows how data moves through PeakATail from raw BAM files to differential APA results. The diagram covers both the single-sample path (one BAM in `datasets:`) and the multi-sample path (two or more datasets), and the switch subcommands that run after the main pipeline.

## Pipeline diagram

```mermaid
flowchart TD
    BAM["BAM file(s)\n(10x Chromium, sorted+indexed)"]
    MERGE["samtools merge\n(merge_strategy: before)"]
    PEAK["Peak calling\nema/countmatrix/peackcalling.py\n→ pos.bed, neg.bed, pos.mtx, neg.mtx, cb.tsv\n(per BAM, under peakcalling/)"]
    RAWSNAP["Per-dataset raw snapshot\nema/outputs.py::write_raw_peak_outputs\n→ per_dataset/<id>/raw/\n  pos.bed, neg.bed, pas.bed,\n  pos.mtx, neg.mtx, cb.tsv"]
    PERDATABED["Per-dataset BEDs\nema/outputs.py::write_per_dataset_beds\n→ per_dataset/<id>/\n  posbed.bed, negbed.bed, pasbed.bed"]
    UNIFY{"Single sample\nor multi-sample?"}
    ATLAS["Atlas snap\nema/datasets/atlas_snap.py\n→ unified/atlas_snapped.bed\n   unified/atlas_mapping.tsv"]
    MERGE2["PAS coordinate merge\nema/datasets/pas_merge.py\n→ unified/unified.bed\n   unified/atlas_mapping.tsv"]
    CONCAT["Matrix concatenation\nema/datasets/pas_merge.py::concat_matrices\n→ unified/concatenated.mtx\n   unified/concatenated_cbs.tsv"]
    CBFILTER["Cell barcode filter\nema/matrixfilter.py::filter_cb\nmin_read threshold\n→ filtered_cb.tsv"]
    GTF["GTF annotation\nema/annotate/gtf_cache.py\n→ gtf_cache/gene_end.bed\n   gtf_cache/utr_lengths.json\n(parallel with peak calling)"]
    FINDCLOSE["PAS-gene assignment\nema/annotate/find_close.py\n→ annotatedpas.bed\n   pas_gene.tsv"]
    ANNOTMATRIX["Annotated count matrix\nema/outputs.py::write_annotated_matrix\n→ per_dataset/<id>/\n  annotated_matrix.mtx\n  annotated_pas_ids.tsv\n  annotated_cells.tsv"]
    PREPROCESS["Preprocessing\nema/matrixfilter.py::preprocessing\nmin_cells, min_pas_per_cell filters\n→ per_dataset/<id>/preprocessed.h5ad"]
    CLUSTER["Clustering\nema/clustering/clustering.py\nTF-IDF + LSI + Leiden\n→ per_dataset/<id>/clusters.h5ad"]
    MATCH["Cross-dataset cluster matching\nema/clustering/cross_dataset/\n→ cross_dataset/canonical_cluster_map.tsv\n(only when >1 dataset)"]
    VIZ["Post-pipeline visualisation\nema/viz/pipeline_hooks.py\n→ figures/umap_*.png\n   figures/clusters_*.png\n   figures/peak_qc_*.png\n   figures/resource_timeline.png\n   figures/run_report.html"]

    DIFF["ema switch diff\nFisher exact / NB regression\n→ switch_diff_<ts>/\n  differential/fisher_A_vs_B.tsv\n  markers.tsv (only if --marker-top-n > 0)\n  figures/volcano_*.png\n  figures/figures_INDEX.md"]
    LENGTH["ema switch length\nPDUI / proportion / entropy\n→ switch_length_<ts>/\n  pdui_classic.tsv\n  pdui_proportion.tsv\n  entropy_shannon.tsv\n  figures/pdui_distribution.png"]
    GENEVIEW["ema switch geneview\nPer-gene coverage tracks\n→ switch_geneview_<ts>/\n  figures/gene_<id>.png\n  figures/figures_INDEX.md"]

    BAM -->|one per BAM| MERGE
    BAM -->|merge_strategy: none| PEAK
    MERGE --> PEAK
    PEAK --> RAWSNAP
    PEAK --> PERDATABED
    PERDATABED --> UNIFY
    UNIFY -->|"atlas: set"| ATLAS
    UNIFY -->|no atlas| MERGE2
    UNIFY -->|single sample| CBFILTER
    ATLAS --> CONCAT
    MERGE2 --> CONCAT
    CONCAT --> CBFILTER
    CBFILTER --> GTF
    GTF -->|"utr_lengths"| FINDCLOSE
    PERDATABED -->|"pasbed.bed"| FINDCLOSE
    FINDCLOSE --> ANNOTMATRIX
    ANNOTMATRIX --> PREPROCESS
    PREPROCESS --> CLUSTER
    CLUSTER --> MATCH
    CLUSTER --> VIZ
    CLUSTER -->|"clusters.h5ad"| DIFF
    CLUSTER -->|"clusters.h5ad"| LENGTH
    DIFF -->|"fisher_*.tsv"| GENEVIEW
    LENGTH -.->|"optional overlay"| GENEVIEW
```

## Stage annotations

### BAM input

PeakATail reads sorted, indexed BAM files produced by Cell Ranger or STARsolo. Each BAM must have cell barcodes in the tag specified by `barcode_tag` (default `CB`). UMI tags are not used; PeakATail counts raw read 3' ends, not UMI-deduplicated molecules.

### Peak calling

`ema/countmatrix/peackcalling.py` iterates reads strand by strand. For each position where a read's 3' end falls, it increments a per-position coverage array. Peaks are detected from this coverage profile using the strategy class set by `--peak-strategy` (default `original`). Each peak emits one or more PAS coordinates (up to `max_pas` per peak, separated by at least `pas_gap` bp). Output is one BED file and one MatrixMarket sparse count matrix per strand per BAM.

### Atlas snap vs coordinate merge

In the multi-sample path, all per-BAM BEDs must be unified to a common coordinate space before the matrices can be added together.

- **Atlas snap** (`atlas:` key set): each called coordinate is replaced by the nearest entry in the reference BED within `atlas_distance` bp. Coordinates that do not match any atlas entry are discarded. This is the recommended approach when a validated poly(A) site database is available (e.g. PolyASite 2.0 or GENCODE).
- **Coordinate merge** (no atlas): nearby PAS from different datasets are merged if they fall within `pas_gap` bp of each other. A representative coordinate (median position) is selected. This is less precise but works without an external reference.

### Clustering

`ema/clustering/clustering.py` dispatches to one of two built-in strategies:

- `leiden_tfidf` (default): TF-IDF normalisation (Signac Method 1) → SVD/LSI dimensionality reduction with depth-correlation filtering (ArchR-style) → k-NN graph → Leiden community detection.
- `leiden_libsize`: library-size normalisation → PCA on highly variable PAS → k-NN graph → Leiden.

The resolution parameter controls cluster granularity; higher values produce more clusters.

### Switch subcommands

`ema switch diff`, `ema switch length`, and `ema switch geneview` are independent commands that read `clusters.h5ad` as input. They write their results inside a timestamped subdirectory under the original run directory, keeping every analysis self-contained.
