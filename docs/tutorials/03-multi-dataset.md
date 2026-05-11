# Multi-dataset workflow

PeakATail natively handles multiple BAM files or multiple sample groups in a single run. This tutorial covers:

1. Configuring two datasets in one YAML
2. Optional: snapping called PAS to a reference atlas
3. Aligning clusters across datasets with `ema switch match`

---

## Step 1 — Configure two datasets

```yaml
# multi_run.yaml
datasets:
  - id: control
    merge_strategy: before        # merge two replicate BAMs, then call peaks
    bams:
      - /data/control_rep1.bam
      - /data/control_rep2.bam

  - id: treated
    merge_strategy: none          # single BAM, no merging
    bams:
      - /data/treated.bam

gtf: /refs/Homo_sapiens.GRCh38.99.gtf
seqlen: 150
cb_len: 16
barcode_tag: CB
min_read: 1500
min_cells: 3
min_pas_per_cell: 50
pas_gap: 100
```

### `merge_strategy` — which one to pick

PeakATail supports three values for `merge_strategy` per dataset. Pick by what
your BAMs represent biologically:

| `merge_strategy` | Use when… | What it does |
|---|---|---|
| `before` | Multiple BAMs are **technical replicates / lanes of the same library** — same cells, same barcodes, different sequencing runs. | `samtools merge` all BAMs into one before peak calling. More reads per PAS → better detection in lowly-covered regions. |
| `after` | Multiple BAMs are **separate libraries with potentially different cells** but you still want to call them as one logical dataset (rare). | Peak-call each BAM independently, then merge the PAS coordinate sets (union with gap-collapse). PAS counts are NOT summed across BAMs. |
| `none` | **Single BAM per dataset**. The typical case. | Skip the merge step entirely. |

!!! warning "Don't merge libraries with different cell-barcode spaces"
    `merge_strategy: before` does NOT correct barcodes. If you merge BAMs from
    two distinct 10x runs, identical-looking barcodes from different cells
    will collide. Use multi-dataset mode (one `datasets:` entry per library)
    + `ema switch match` for cross-library analysis instead.

---

## Step 2 — Run with atlas snap (optional but recommended for multi-sample)

When you have two datasets, the same PAS may be called at slightly different coordinates in each BAM. The `atlas` key resolves this by snapping every called coordinate to the nearest entry in a reference poly(A) site database.

```yaml
# Add to multi_run.yaml
atlas: /refs/polyasite_2.0_GRCh38_sorted.bed
atlas_distance: 50   # maximum snap distance in bp (default 50)
```

Download the PolyASite 2.0 human atlas from [polyasite.unibas.ch](https://polyasite.unibas.ch/atlas) and sort it with `bedtools sort`. With the atlas enabled, the pipeline:

1. Calls PAS per dataset independently.
2. For each called PAS, finds the nearest atlas entry within `atlas_distance` bp.
3. Replaces the called coordinate with the atlas coordinate, producing a unified coordinate space.
4. Writes `unified/atlas_mapping.tsv` so you can trace which called PAS became which atlas PAS.

Without the atlas, coordinates are merged by proximity (gap controlled by `pas_gap`), but subtle inter-sample differences in peak summits remain.

---

## Step 3 — Run the pipeline

```bash
uv run ema run --config multi_run.yaml --threads 8
```

The multi-sample path runs each dataset through its own downstream pipeline (annotation, matrix filtering, preprocessing, clustering) in parallel worker processes. Progress is reported per dataset.

Output layout for two datasets:

```
peakatail_runs/emaout_<timestamp>/
  peakcalling/
    control_0.pos.bed
    control_0.neg.bed
    treated_0.pos.bed
    treated_0.neg.bed
    ...
  unified/
    atlas_snapped.bed          # unified coordinate set (when --atlas is set)
    atlas_mapping.tsv          # called PAS -> atlas PAS mapping
    concatenated.mtx           # joint count matrix
    concatenated_cbs.tsv       # joint cell barcode list (prefixed with dataset id)
  per_dataset/
    control/
      pasbed.bed
      clusters.h5ad
      ...
    treated/
      pasbed.bed
      clusters.h5ad
      ...
  cross_dataset/
    canonical_cluster_map.tsv  # cluster correspondence across datasets
```

---

## Step 4 — Run cross-dataset cluster matching

After clustering, cluster labels are independent integers within each dataset. Cluster 0 in `control` may correspond to cluster 3 in `treated`. `ema switch match` aligns them using marker PAS overlap:

```bash
uv run ema switch match \
  --h5ad peakatail_runs/emaout_<timestamp>/per_dataset/control/clusters.h5ad \
  --h5ad peakatail_runs/emaout_<timestamp>/per_dataset/treated/clusters.h5ad \
  --strategy marker_overlap \
  --n-top-markers 50
```

The `marker_overlap` strategy:

1. For each dataset independently, ranks PAS by their Wilcoxon score within each cluster (top 50 per cluster).
2. For each pair of clusters across datasets, computes Jaccard similarity of their marker PAS sets.
3. Assigns a canonical cluster label to the cluster pair with the highest Jaccard score.

Output: `cross_dataset/canonical_cluster_map.tsv` with columns `dataset_id`, `cluster_label`, and `canonical_cluster`. Use this table to group clusters for differential testing:

```bash
# Test control cluster 0 (= canonical cluster A) vs treated cluster 3 (= canonical cluster A)
uv run ema switch diff \
  --h5ad .../per_dataset/control/clusters.h5ad \
  --h5ad .../per_dataset/treated/clusters.h5ad \
  --cluster-pairs "0,3" \
  --strategy fisher
```

!!! note "Atlas snap is strongly recommended for multi-sample differential testing"
    Without atlas coordinates, a PAS at chr1:22639574 in control and chr1:22639571 in treated is treated as two separate PAS. With atlas snap both resolve to the same database entry and are directly comparable.
