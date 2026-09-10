# Differential APA strategies

Differential APA testing identifies poly(A) sites that are used at different
rates between two or more cell populations. All three strategies operate on a
per-cell count matrix of shape `(n_cells, n_pas)` and a `cluster_labels` Series
mapping cell barcodes to cluster names. Each strategy applies Benjamini-Hochberg
FDR correction across tested PAS and returns a DataFrame indexed by `pas_id`.

Choose a strategy based on the trade-off between speed, statistical rigor, and
whether you need to compare two clusters or all clusters at once.

For CLI usage see [../cli/switch-diff.md](../cli/switch-diff.md).

---

## `fisher` {#fisher}

**Question:** For each PAS within a gene, does the within-gene proportion of
reads at that PAS differ significantly between cluster 1 and cluster 2?

### Algorithm

1. Aggregate reads per PAS across all cells in cluster 1 and cluster 2
   separately (`agg1`, `agg2`).
2. Group PAS by gene using `pas_gene_map`. When `pas_gene_map` is not
   supplied, a warning is logged and all PAS are treated as a single global
   group (legacy behavior, loses the APA-specific framing).
3. For each gene G with at least 2 PAS, and for each PAS p in G, construct a
   2x2 contingency table:

   ```
                         cluster 1        cluster 2
   reads at p           reads_p_c1       reads_p_c2
   reads at other PAS   other_c1         other_c2
   ```

   where `other_cX = total gene reads in cX - reads_p_cX`.

4. Run two-sided Fisher exact test via `scipy.stats.fisher_exact`.
5. Compute within-gene proportions `prop1 = reads_p_c1 / gene_reads_c1`
   and `prop2 = reads_p_c2 / gene_reads_c2`.
6. Compute `delta_proportion = prop1 - prop2` and
   `log2fc = log2((prop2 + eps) / (prop1 + eps))` where
   `eps = 1 / (max(gene_reads_c1, gene_reads_c2) + 1)` prevents log(0).
7. Apply Benjamini-Hochberg FDR to all p-values using
   `scipy.stats.false_discovery_control`.

This framing is the canonical APA differential-usage design matching the
DEXSeq-for-splicing paradigm: the test asks "for this gene, does this PAS
change its within-gene share?" rather than "does this PAS change relative to
all other PAS genome-wide?".

### Inputs and outputs

**Inputs:**

- `count_matrix`: `pd.DataFrame` shape `(n_cells, n_pas)`
- `cluster_labels`: `pd.Series` indexed by cell barcode
- `cluster1`, `cluster2`: cluster label strings (both required)
- `pas_gene_map`: `dict[str, str]` mapping PAS ID string to gene ID (recommended)
- `min_cells_per_group`: minimum cells per cluster (total count, not non-zero)

**Output columns (indexed by `pas_id`):**

| Column | Type | Description |
|---|---|---|
| `pvalue` | float | Two-sided Fisher p-value |
| `qvalue` | float | BH FDR-corrected q-value |
| `n_cells` | int | Total cells in both clusters |
| `n_cells_cluster1` | int | Cells in cluster 1 |
| `n_cells_cluster2` | int | Cells in cluster 2 |
| `n_cells_expr_cluster1` | int | Cells with non-zero count at this PAS in cluster 1 |
| `n_cells_expr_cluster2` | int | Cells with non-zero count at this PAS in cluster 2 |
| `n_reads_pas_cluster1` | int | Total reads at this PAS in cluster 1 |
| `n_reads_pas_cluster2` | int | Total reads at this PAS in cluster 2 |
| `n_reads_gene_cluster1` | int | Total reads at all PAS of this gene in cluster 1 |
| `n_reads_gene_cluster2` | int | Total reads at all PAS of this gene in cluster 2 |
| `odds_ratio` | float | Fisher odds ratio |
| `delta_proportion` | float | prop1 - prop2 (within-gene proportions) |
| `log2fc` | float | log2(prop2 / prop1), positive = higher in cluster 2 |

Results are sorted by `qvalue` ascending.

### Tunable hyperparameters

| Parameter | Default | CLI flag | YAML key | Description |
|---|---|---|---|---|
| `min_cells_per_group` | 10 | `--min-cells-per-group` | `min_cells_per_group` | Minimum cells in each cluster for a PAS to be testable. |
| `fdr` | 0.05 | `--fdr` | `fdr` | BH FDR threshold applied by the caller when filtering results. |

### Interpretation

A significant result (qvalue < FDR threshold) with `delta_proportion > 0`
means PAS p captures a larger fraction of its gene's reads in cluster 1 than
in cluster 2. The `volcano` figure plots `log2fc` vs `-log10(qvalue)`;
significant PAS appear in the upper-left (down in cluster 2) or upper-right
(up in cluster 2) quadrants beyond the dashed threshold lines.

Audit a result row: if `n_reads_gene_cluster1 = 5` the test was run on very
few reads and the odds ratio is unreliable. Use `n_reads_pas_cluster1` and
`n_reads_gene_cluster1` to filter low-coverage results before reporting.

### Limitations

- Pseudo-replication: Fisher operates on aggregated (pseudo-bulk) read counts.
  Treating all cells in a cluster as a single observation inflates statistical
  power. For rigorous cell-level testing use `nb_pairwise`.
- Genes with fewer than 2 PAS are silently excluded (no within-gene comparison
  is possible).
- Fisher is anti-conservative for overdispersed count data. A q-value of 0.01
  from Fisher should be treated as directional evidence, not as a precise
  error rate.

---

## `nb_pairwise` {#nb_pairwise}

**Question:** For each PAS individually, is the per-cell read count in cluster
1 significantly different from cluster 2 after correcting for library size?

### Algorithm

1. Restrict the count matrix to cells in cluster 1 and cluster 2.
2. For each PAS, filter cells with fewer than `min_cells_per_group` non-zero
   counts in either group; skip that PAS.
3. Estimate per-PAS dispersion using a Negative Binomial MLE
   (`statsmodels.discrete.discrete_model.NegativeBinomial`, method="bfgs",
   maxiter=200). If the solver fails to converge, fall back to a
   method-of-moments estimate:

    ```
    alpha_hat = (sample_variance - mean) / mean^2
    ```

   Dispersion is clipped to `[1e-4, 10.0]` for numerical stability. PAS that
   land on the **lower** clip are flagged `dispersion_floored` and get no
   q-value — see [The dispersion floor](#dispersion-floor).

4. Fit a NB GLM with fixed dispersion:

   `count ~ cluster_indicator + offset(log(library_size))`

   using IRLS (method="irls", maxiter=200). The cluster indicator is 0 for
   cluster 1 and 1 for cluster 2.

5. Extract the Wald z-statistic and p-value for the cluster coefficient.
   `log2fc = coefficient / log(2)`.

6. Apply BH FDR across all tested PAS. Dispersion-floored PAS stay in the BH
   input (so the family size *m* is unchanged) but their reported `qvalue` is
   set to `NaN`.

Per-PAS GLM fits are embarrassingly parallel. Dense numpy slices (cells x 1
per PAS) are pre-extracted once before dispatch so workers never pickle large
sparse matrices. Parallelism is capped by `ResourceManager` with
`per_worker_mb=300`.

### Inputs and outputs

**Inputs:**

- `count_matrix`: `pd.DataFrame` shape `(n_cells, n_pas)`
- `cluster_labels`: `pd.Series`
- `cluster1`, `cluster2`: cluster label strings (both required)
- `min_cells_per_group`: minimum non-zero cells per group per PAS

**Output columns (indexed by `pas_id`):**

| Column | Type | Description |
|---|---|---|
| `pvalue` | float | Wald test p-value |
| `qvalue` | float | BH FDR-corrected q-value — **`NaN` when `dispersion_floored`** |
| `log2fc` | float | log2 fold change (cluster2 / cluster1) |
| `dispersion` | float | Estimated NB dispersion alpha |
| `dispersion_floored` | bool | `True` when alpha hit the `1e-4` lower clip |
| `n_cells` | int | Total cells in both clusters |
| `test_stat` | float | Wald z-statistic |

Results are sorted by `qvalue` ascending.

### Tunable hyperparameters

| Parameter | Default | CLI flag | YAML key | Description |
|---|---|---|---|---|
| `min_cells_per_group` | 10 | `--min-cells-per-group` | `min_cells_per_group` | Minimum cells with non-zero counts in each cluster per PAS. Raise to reduce noisy low-count tests. |
| `fdr` | 0.05 | `--fdr` | `fdr` | BH FDR threshold applied when filtering results. |

### Interpretation

`nb_pairwise` operates at the cell level and correctly accounts for library
size differences between clusters via the `offset(log(library_size))` term.
A positive `log2fc` means the PAS is more abundant in cluster 2. The
`dispersion` column shows the per-PAS overdispersion.

The `fallback_used` flag (dropped from the output DataFrame but visible in
logs) is set when the MLE solver failed and the method-of-moments path was
used. Filter by `dispersion < 9` to exclude rows where the dispersion hit
the upper clip and the GLM may be unreliable.

### The dispersion floor {#dispersion-floor}

When the per-PAS dispersion estimate collapses onto the **lower** clip
(`alpha = 1e-4`) the NB GLM is effectively a **Poisson** fit: the Wald standard
error is then a lower bound and the p-value is anti-conservative. Worse, alpha
is estimated under the *full* model (`count ~ cluster_indicator + offset`), so
sampling noise that happens to mimic a group difference is absorbed by the
fitted coefficient and pushes alpha *down* — exactly the tests that then look
most significant. Under the 20-run label-permutation null of
[issue #94](https://github.com/BMGLab/PeakATail/issues/94), 8.5 % of null tests
hit the floor and those tests accounted for **67 % of all false `q < 0.05`
calls**.

Those rows are therefore reported with `dispersion_floored = True` and
`qvalue = NaN`, so no `qvalue < fdr` filter can call them significant. The row
and its raw `pvalue` are kept so you can still inspect the evidence — treat
that p-value as a *lower bound* on the true one, never as a calibrated error
rate. `NbPairwiseStrategy.test` logs a warning naming how many PAS were
affected.

!!! warning "nb_pairwise q-values are not permutation-calibrated"

    Withholding the floored rows removes the dominant source of
    anti-conservatism, but it does not make the remaining q-values calibrated.
    If a `nb_pairwise` q-value has to carry an error-rate claim (e.g. in a
    manuscript), calibrate it against a label-permutation null of your own
    data. A permutation-calibrated q — or dispersion shrinkage toward a
    common trend, as edgeR/DESeq2 do — is still open on issue #94.

### Limitations

- Significantly slower than Fisher (10-100x depending on cell count and
  parallelism). On 50,000 cells x 20,000 PAS with 8 workers, expect 30-90
  minutes.
- Does not implement a within-gene normalization. A PAS that simply has
  more reads because its gene is more expressed will appear significant.
  Combine with `proportion` output to distinguish expression changes from
  APA shifts.
- Not available for multi-condition comparisons; use `nb_multi` when
  comparing three or more clusters simultaneously.

---

## `nb_multi` {#nb_multi}

**Question:** Does the read count at each PAS differ significantly across all
clusters at once?

### Algorithm

1. Use all clusters in `cluster_labels` (not just a pair). Build a one-hot
   dummy matrix (K-1 columns, drop-first) plus an intercept column as the
   full model design matrix.
2. Per-PAS, require `min_cells_per_group` cells with non-zero counts total,
   and signal in at least 2 clusters; otherwise skip.
3. Estimate dispersion using a full-model NB MLE (same procedure as
   `nb_pairwise`, method="bfgs", maxiter=200, MoM fallback).
4. Fit full model: `count ~ C(cluster) + offset(log(library_size))`.
5. Fit null model: `count ~ 1 + offset(log(library_size))`.
6. Likelihood ratio test statistic: `LRT = 2 * (LL_full - LL_null)`, which
   follows chi-squared with `K - 1` degrees of freedom.
7. Apply BH FDR across all tested PAS.

Parallelized over PAS batches identically to `nb_pairwise`.

### Inputs and outputs

**Inputs:**

- `count_matrix`: `pd.DataFrame` shape `(n_cells, n_pas)`
- `cluster_labels`: `pd.Series` — all unique values define the groups
- `cluster1`, `cluster2`: accepted for interface compatibility but ignored
- `min_cells_per_group`: minimum total non-zero cells across all clusters per PAS

**Output columns (indexed by `pas_id`):**

| Column | Type | Description |
|---|---|---|
| `pvalue` | float | LRT chi-squared p-value |
| `qvalue` | float | BH FDR-corrected q-value |
| `test_stat` | float | LRT statistic (2 * delta log-likelihood) |
| `df` | int | Degrees of freedom (K - 1) |
| `dispersion` | float | NB dispersion alpha from full model |
| `n_cells` | int | Total cells across all clusters |

Results are sorted by `qvalue` ascending. There is no `log2fc` column because
the test is omnibus (K > 2 levels).

### Tunable hyperparameters

| Parameter | Default | CLI flag | YAML key | Description |
|---|---|---|---|---|
| `min_cells_per_group` | 10 | `--min-cells-per-group` | `min_cells_per_group` | Minimum total non-zero cells for a PAS to be tested. Requires signal in at least 2 clusters. |
| `fdr` | 0.05 | `--fdr` | `fdr` | BH FDR threshold applied when filtering results. |

### Interpretation

`nb_multi` tells you which PAS change across any cluster contrast without
specifying which pair. It is the right first pass when you have many clusters
and want a ranked list of the most dynamically regulated PAS. Follow up
significant results from `nb_multi` with `nb_pairwise` on the cluster pairs
of interest to get directional `log2fc` values.

A significant result with `df = K - 1` degrees of freedom: the `test_stat`
on a chi-squared distribution with K-1 df gives the p-value. Check that
`test_stat` >> 0 for significant PAS; test_stat near 0 despite low p-value
is a numerical artifact.

### Limitations

- No `log2fc` column. You must re-run `nb_pairwise` or examine cluster-level
  coefficients from the full model to identify direction and magnitude.
- Requires at least 2 cluster levels; raises `ValueError` if fewer are found.
- Runtime scales with K (number of clusters) and n_PAS. With K=10 clusters
  and 20,000 PAS, expect 2-5x the runtime of `nb_pairwise` for a two-cluster
  test.
