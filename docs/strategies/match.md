# Cross-dataset matching strategies

Cross-dataset matching assigns canonical cluster IDs when the same experiment
is processed in multiple batches, or when you want to align cell populations
across independent samples. Each strategy reads a list of `.h5ad` files
(one per dataset), each of which must have `obs["leiden"]` populated by a
prior `peakatail run` or clustering step. The output is a DataFrame mapping every
`(dataset_id, original_cluster)` pair to a `canonical_cluster` integer.

For CLI usage see [../cli/switch-match.md](../cli/switch-match.md).

All three strategies share the same output schema:

| Column | Type | Description |
|---|---|---|
| `dataset_id` | str | Human-readable dataset identifier |
| `original_cluster` | str | Leiden cluster label from that dataset |
| `canonical_cluster` | int | Canonical cluster ID shared across datasets |
| `match_confidence` | float | Mean Jaccard similarity for all edges within the component (0.0 for singletons) |
| `matched_to` | str | JSON list of `[[dataset_id, original_cluster], ...]` matched partners |

---

## `marker_overlap` {#marker_overlap}

**Question:** Which clusters in dataset A correspond to which clusters in
dataset B based on the overlap of their top-ranked differentially expressed
PAS?

### Algorithm

1. For each dataset, load the `.h5ad` file and compute per-cluster marker PAS
   using Wilcoxon rank-sum (`scanpy.tl.rank_genes_groups`, method="wilcoxon").
   The top-`n_top_markers` PAS by score form the cluster "fingerprint". If
   the data looks raw (integers, max value > 20), log-normalization is applied
   first.
2. For each pair of datasets (A, B), compute the Jaccard similarity matrix
   between every cluster in A and every cluster in B:

    ```
    J(A_i, B_j) = |markers_Ai ∩ markers_Bj| / |markers_Ai ∪ markers_Bj|
    ```

3. Solve the maximum-weight 1-to-1 assignment using the Hungarian algorithm
   (`scipy.optimize.linear_sum_assignment`) on the negated similarity matrix.
   Only pairs with Jaccard > 0 are accepted.
4. Build a graph where nodes are `(dataset_id, original_cluster)` tuples and
   edges connect matched pairs. Extract connected components using union-find
   (transitive closure): if A.0 ~ B.2 and B.2 ~ C.5, all three share one
   canonical ID.
5. Assign canonical integer IDs (1, 2, ...) to each component. Unmatched
   singletons get unique IDs. `match_confidence` is the mean Jaccard of all
   edges within the component.

Marker computation is parallelized across datasets via `joblib.Parallel`
(loky backend, per_worker_mb=200).

### Inputs and outputs

**Inputs:**

- `h5ad_paths`: list of paths to `.h5ad` files
- `dataset_ids`: list of human-readable identifiers (same order)
- `n_top_markers`: number of top marker PAS per cluster

**Outputs:** canonical mapping DataFrame (schema above).

### Tunable hyperparameters

| Parameter | Default | CLI flag | YAML key | Description |
|---|---|---|---|---|
| `n_top_markers` | 50 | `--n-top-markers` | `n_top_markers` | Top-N marker PAS per cluster used as the cluster fingerprint. |

### Interpretation

`match_confidence` near 1.0 indicates near-identical marker sets between
the matched clusters (strong biological correspondence). Values < 0.1 indicate
weak overlap; inspect whether the datasets share the same PAS feature space.
Unmatched clusters (confidence = 0.0) are singletons; they may represent
populations unique to that dataset.

A sanity check: run with `n_top_markers=10` and `n_top_markers=100`; stable
matches (same canonical ID assignments) indicate robust correspondence.

### Limitations

- Requires datasets to share a meaningful fraction of PAS features. Datasets
  processed with different reference BED files or different capture
  technologies may have few shared markers even for the same cell type;
  use `mnn` in that case.
- Wilcoxon rank-sum marker computation inside scanpy requires the h5ad to
  have at least 2 clusters. Single-cluster datasets fall back to returning
  all features as markers.
- The Hungarian assignment enforces 1-to-1 matching per dataset pair. When
  one dataset has more clusters than another, the extra clusters become
  singletons.

---

## `mnn` {#mnn}

**Question:** Which clusters in dataset A correspond to which clusters in
dataset B based on how their cells link in a shared embedding space, even
when marker gene overlap is low?

### Algorithm

1. Load all `.h5ad` files and concatenate their cell × feature matrices along
   the obs axis, intersecting to the common set of PAS features
   (`var_names`). A `batch` column records the dataset origin. Raises
   `ValueError` if no common features exist.
2. Log-normalize if the data looks raw, then run Truncated SVD (LSI) with
   `n_components` components on the concatenated matrix. L2-normalize the
   embedding rows.
3. For each pair of datasets (A, B), find Mutual Nearest Neighbors in the
   shared embedding:
   - For each cell in A, find the `k_neighbors` nearest cells in B
     (cosine distance).
   - For each cell in B, find the `k_neighbors` nearest cells in A.
   - A pair (i, j) is mutual if i is in B's neighbors of j AND j is in A's
     neighbors of i.
4. For each cluster A_i, count MNN votes to each cluster B_j. Compute:

    ```
    sim(A_i, B_j) = (MNN votes from A_i to B_j) / |cells in A_i|
    ```

5. Hungarian assignment + transitive closure, identical to `marker_overlap`.

### Inputs and outputs

**Inputs:**

- `h5ad_paths`, `dataset_ids`
- `n_components`: LSI components for shared embedding
- `k_neighbors`: MNN search neighborhood size

**Outputs:** canonical mapping DataFrame (schema above).

### Tunable hyperparameters

| Parameter | Default | CLI flag | YAML key | Description |
|---|---|---|---|---|
| `n_components` | 30 | `--mnn-components` | `mnn_components` | LSI components for the shared embedding. More components capture finer structure at higher cost. |
| `k_neighbors` | 10 | `--mnn-k-neighbors` | `mnn_k_neighbors` | Nearest neighbors for MNN search. Larger values produce denser MNN graphs. |

### Interpretation

MNN is most useful when datasets come from different protocols (e.g., Smart-seq2
and 10x Chromium) or have very different library sizes that reduce the marker
gene overlap. MNN `match_confidence` is the fraction of cluster A_i cells
whose dominant MNN partner is in cluster B_j; values above 0.3 indicate
reliable matches.

Check the shared embedding by inspecting whether cells from different datasets
intermingle in a UMAP of the concatenated matrix. Strong batch separation in
the embedding will cause false negatives (no MNN links across the boundary).

### Limitations

- Re-embedding ALL datasets together is memory-intensive: memory scales with
  total cell count x n_features. For > 100,000 cells, consider subsampling.
- No batch correction (no Harmony, no Scanorama). Strong batch effects may
  dominate the embedding, causing biologically similar clusters to fail to
  link.
- All datasets must share the same PAS feature space (same unified PAS BED
  file). Datasets processed with different BED files will have no common
  features and the strategy raises `ValueError`.
- MNN search cost is O(n_cells_A * n_cells_B) per pair; for large datasets
  this is slower than `marker_overlap`.

---

## `jaccard` {#jaccard}

**Question:** Do the same physical cells appear in the same cluster across two
runs of the same BAM (for regression testing or split-BAM experiments)?

### Algorithm

1. Load each `.h5ad` and build a mapping from cluster label to the set of
   stripped cell barcodes in that cluster.
2. Strip dataset prefixes from barcodes: `sample#ACGT-1` → `ACGT-1`, and
   `sampleA_ACGTGCATAGCT` → `ACGTGCATAGCT` (only when the suffix looks like
   a real DNA barcode; `cellA_0` is kept as-is to avoid misidentifying
   synthetic test names).
3. For each pair of datasets (A, B), compute the Jaccard similarity matrix
   between cluster barcode sets:

    ```
    J(A_i, B_j) = |CB(A_i) ∩ CB(B_j)| / |CB(A_i) ∪ CB(B_j)|
    ```

4. Hungarian assignment + transitive closure, identical to `marker_overlap`.

### Inputs and outputs

**Inputs:**

- `h5ad_paths`, `dataset_ids`
- `n_top_markers`, `n_jobs`: accepted for interface compatibility but unused

**Outputs:** canonical mapping DataFrame (schema above).

### Tunable hyperparameters

None. This strategy performs pure set operations on barcodes; there are no
algorithm parameters.

### Interpretation

When the same physical cells are clustered twice (deterministic run with the
same seed), `match_confidence` should approach 1.0. Values below 0.9 for
identical runs indicate non-determinism in clustering (e.g., different igraph
versions) or barcode format mismatches. Use `match_confidence` as a
reproducibility metric.

### Limitations

- Biologically distinct samples will share no barcodes, producing no matches.
  Do not use `jaccard` for cross-sample comparison; use `marker_overlap` or
  `mnn` instead.
- The DNA-barcode heuristic (`[ACGTN]{6,}(-\d+)?`) may incorrectly strip
  prefixes from synthetic test barcodes that happen to contain long DNA-like
  strings. Verify barcode stripping on your dataset with a small test before
  relying on it.
