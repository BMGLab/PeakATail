# `peakatail switch match`

`peakatail switch match` assigns canonical cluster identities across two or more
`clusters.h5ad` files produced by `peakatail run`. When the same experiment is
processed in multiple datasets (e.g. different time points, conditions, or
replicates), cluster labels from Leiden are arbitrary integers that do not
correspond between datasets. This command finds which cluster in dataset A
best corresponds to which cluster in dataset B and writes a mapping table
with a confidence score.

Three strategies are available:

- **`marker_overlap`** — Jaccard similarity of top-N marker PAS per cluster (default, fastest).
- **`mnn`** — Mutual Nearest Neighbours in a shared LSI embedding.
- **`jaccard`** — Direct Jaccard overlap of cell barcodes (only meaningful when the same physical cells appear in multiple datasets).

When `--output` is left at its default and the `--h5ad` files come from a
recognisable `peakatail_runs/<run>/` directory, output is routed **inside** the
originating run directory as `peakatail_runs/<run>/switch_match_<timestamp>/`.
Source: `ema/cli/common.py::resolve_subcommand_output_dir`.

!!! note "When to use it"
    - You have multiple datasets from the same experiment and want to align
      cluster labels for cross-condition comparison.
    - You want to verify that clusters are reproducible across biological
      replicates before pooling differential APA results.
    - You are running `peakatail switch diff` on multiple datasets and need canonical
      cluster IDs to interpret shared vs. unique cluster populations.

!!! warning "When NOT to use it"
    - You have only one dataset. The command produces an identity mapping but
      it is not informative.
    - You are using `jaccard` for datasets from different biological samples
      that share no cell barcodes. The result will contain no matches.
    - You are using `mnn` for datasets with no shared PAS features (different
      genomic regions captured). The runner raises a `ValueError`.

## Quick example

```bash
uv run peakatail switch match \
  --h5ad peakatail_runs/emaout_.../per_dataset/sampleA/clusters.h5ad \
  --h5ad peakatail_runs/emaout_.../per_dataset/sampleB/clusters.h5ad \
  --strategy marker_overlap \
  --n-top-markers 50
```

What lands on disk:

- `peakatail_runs/emaout_.../switch_match_<ts>/cluster_match.tsv` — canonical mapping table.
- `peakatail_runs/emaout_.../switch_match_<ts>/peakatail_<ts>.log` — run log.
- Cluster similarity heatmap figures (when plotting is enabled).

## Full `--help` output

```text
Usage: peakatail switch match [OPTIONS]

  Cross-dataset cluster matching (marker_overlap / mnn / jaccard).

Options:
  --list-strategies               Print available match strategies and exit.
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
  -i, --h5ad PATH                 Per-dataset clusters.h5ad files.  [required]
  -s, --strategy TEXT             [default: marker_overlap]
  --n-top-markers INTEGER         [default: 50]
  --mnn-components INTEGER        SVD components for MNN shared embedding
                                  (--strategy mnn). Default 30.  [default: 30]
  --mnn-k-neighbors INTEGER       Nearest neighbours for MNN search (--
                                  strategy mnn). Default 10.  [default: 10]
  --help                          Show this message and exit.
```

## Flags

### Inputs

| Flag | Type | Default | Description |
|---|---|---|---|
| `--h5ad` / `-i` | PATH (repeatable) | — | Two or more `clusters.h5ad` files. The order determines the synthetic dataset IDs (`"0"`, `"1"`, `"2"`, ...) used in the output table. Required. |

### Strategy options

| Flag | Type | Default | Description |
|---|---|---|---|
| `--strategy` / `-s` | TEXT | `marker_overlap` | Matching algorithm. Run `peakatail switch match --list-strategies` to see all registered names. |
| `--n-top-markers` | INT | 50 | Number of top marker PAS per cluster to use for `marker_overlap`. Higher values include more PAS in the Jaccard fingerprint; lower values are faster. Has no effect on `mnn` or `jaccard` (accepted for API compatibility). |
| `--mnn-components` | INT | 30 | Number of TruncatedSVD (LSI) components for the shared embedding used by `--strategy mnn`. More components capture finer structure but increase cost. Source: `ema/clustering/cross_dataset/mnn.py::MNNStrategy`. |
| `--mnn-k-neighbors` | INT | 10 | Number of nearest neighbours searched per cell in the MNN step. Larger k produces denser MNN graphs and slower kNN search but may improve matching for sparse clusters. Source: `ema/clustering/cross_dataset/mnn.py::MNNStrategy`. |

### Output

| Flag | Type | Default | Description |
|---|---|---|---|
| `--output` / `-o` | PATH | `switch_out` | Output directory base name. Auto-routed inside the originating run dir when inputs are from a `peakatail_runs/` path. |

## Strategies

### `marker_overlap` (default)

Computes per-cluster marker PAS using Wilcoxon rank-sum test via
`scanpy.tl.rank_genes_groups`. The top-N markers form the cluster's
"fingerprint". Jaccard similarity between fingerprints is computed for all
cluster pairs across datasets. Hungarian assignment solves the maximum-weight
1-to-1 match; transitive closure propagates canonical IDs across three or more
datasets.

**When to use:** Datasets from different samples with the same PAS feature
space. Fails gracefully for datasets with very different marker gene sets.

**Parallelism:** Marker computation is parallelised across datasets using
`joblib` with the `loky` backend. Worker count respects `--threads`.

Source: `ema/clustering/cross_dataset/marker_overlap.py::MarkerOverlapStrategy`.

### `mnn`

Concatenates cells from all datasets into a joint matrix, computes a shared
TruncatedSVD (LSI) embedding without explicit batch correction, then finds
Mutual Nearest Neighbour (MNN) pairs across datasets. For each cluster pair
(A_i, B_j), the similarity is the fraction of A_i cells whose dominant MNN
partner is in B_j. Hungarian assignment + transitive closure follows.

**When to use:** Datasets with low marker PAS overlap due to different capture
technologies or library size distributions. More sensitive to batch effects than
`marker_overlap`; use when marker overlap is insufficient.

**Limitations:** Memory scales with total cell count (all datasets concatenated).
All datasets must share the same feature space (same PAS BED). For >100K
total cells, consider subsampling. Source: `ema/clustering/cross_dataset/mnn.py::MNNStrategy`.

### `jaccard`

Computes Jaccard similarity between the sets of cell barcodes in each cluster,
with dataset prefixes stripped (`sampleA#ACGT` → `ACGT`; `sampleA_ACGTACGT` →
`ACGTACGT` for DNA-like suffixes). Useful only when the same physical cells
appear in multiple datasets.

**When to use:** Regression testing (same BAM processed twice should give
confidence ≈ 1.0) or two analysis passes of the same experiment. Not
appropriate for distinct biological samples. Source:
`ema/clustering/cross_dataset/jaccard.py::JaccardCellSetStrategy`.

## Output files

**`cluster_match.tsv`**

Written by all three strategies. Canonical mapping table.

| Column | Type | Description |
|---|---|---|
| `dataset_id` | str | Synthetic dataset identifier (`"0"`, `"1"`, ...) matching the order of `--h5ad` arguments. |
| `original_cluster` | str | Leiden cluster label as it appears in `adata.obs["leiden"]`. |
| `canonical_cluster` | int | Integer ID assigned to the connected component this cluster belongs to across datasets. Clusters that match across datasets share the same `canonical_cluster`. |
| `match_confidence` | float | Mean Jaccard similarity across all edges in this component. 1.0 = perfect match; 0.0 = isolated (no match found in other datasets). |
| `matched_to` | str | JSON-encoded list of `[dataset_id, original_cluster]` pairs that were matched to this entry within the same canonical component. Empty list `[]` for singletons. |

## How it relates to other commands

- **[`peakatail run`](run.md)** — produces the `clusters.h5ad` inputs. The `leiden` column in `obs` is required.
- **[`peakatail switch diff`](switch-diff.md)** — use `cluster_match.tsv` to interpret which clusters are comparable across datasets before running pairwise differential APA.
- **[`peakatail switch geneview`](switch-geneview.md)** — visualise PAS usage for canonically matched clusters.

## See also

- Strategy details in [`../strategies/`](../strategies/) — algorithm descriptions for cross-dataset matching.
- Tutorial in [`../tutorials/`](../tutorials/) — multi-dataset matching walkthrough.
