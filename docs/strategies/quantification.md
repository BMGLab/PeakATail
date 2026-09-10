# Quantification strategies

Quantification converts the PAS count matrix into a per-cell summary of
alternative polyadenylation (APA) usage. The three strategies answer different
biological questions. Choose based on what you want to measure.

Each strategy implements a `compute(count_matrix, pas_isoform_map, ...)` method
and writes one TSV to the output directory. The count matrix has shape
`(n_pas, n_cells)` with integer PAS IDs as the row index. The isoform map
encodes `{pas_id: [(gene_id, transcript_id, transcript_pos, rank, total_pas_in_transcript), ...]}`.

All three strategies support two aggregation modes:

- `per_gene` (genomic-coordinate ordering): proximal = first PAS by genomic
  position on the gene, distal = last PAS. Transcript ID column is set to
  the sentinel `"_gene_"`.
- `per_isoform` (transcript-coordinate ordering): proximal/distal are
  determined by the rank field in the isoform map for each transcript
  individually.

Strategies parallelize using `joblib.Parallel` when `n_genes > 5000`.

---

## `classic` {#classic}

**Question:** For each gene and cell, what fraction of reads go to the distal
poly(A) site (3' UTR lengthening)?

**Output file:** `pdui_classic.tsv`

### Algorithm

1. For each gene, identify the proximal PAS (rank 0) and the distal PAS
   (highest rank) from the isoform map. For genes with only one PAS, skip.
2. Extract raw counts for proximal and distal PAS across all cells.
3. Compute per-cell PDUI:

    ```
    PDUI = distal / (proximal + distal + pseudocount)
    ```

4. When `proximal + distal + pseudocount = 0`, PDUI is set to `NaN`.
5. Collect results in long format; one row per (gene, transcript, cell).

With `aggregation="per_gene"`, proximal/distal are selected using min/max
rank across all transcripts of the gene (genomic-coordinate proxy). With
`aggregation="per_isoform"`, selection is done independently per transcript
so a transcript with only 2 PAS always uses those 2 directly.

### Inputs and outputs

**Inputs:**

- `count_matrix`: `pd.DataFrame` shape `(n_pas, n_cells)`
- `pas_isoform_map`: `dict[int, list[tuple[str, str, int, int, int]]]`

**Output columns in `pdui_classic.tsv`:**

| Column | Type | Description |
|---|---|---|
| `gene_id` | str | Gene identifier |
| `transcript_id` | str | Transcript ID, or `"_gene_"` when `aggregation="per_gene"` |
| `proximal_pas_id` | int | PAS ID selected as proximal |
| `distal_pas_id` | int | PAS ID selected as distal |
| `cell` | str | Cell barcode |
| `pdui` | float64 | Distal usage index in [0, 1] or NaN |
| `proximal_reads` | float | Raw read count at proximal PAS |
| `distal_reads` | float | Raw read count at distal PAS |
| `total_reads` | float | proximal_reads + distal_reads (excludes pseudocount) |

### Tunable hyperparameters

| Parameter | Default | CLI flag | YAML key | Description |
|---|---|---|---|---|
| `pseudocount` | 0.0 | `--pdui-pseudocount` | `pdui_pseudocount` | Added to denominator. Set to 1.0 to avoid NaN for zero-count cells. |
| `aggregation` | `"per_isoform"` | `--isoform-agg` | `isoform_agg` | Level at which proximal/distal PAS are selected. |

### Interpretation

PDUI = 0 means all reads at the proximal site (shorter 3' UTR). PDUI = 1
means all reads at the distal site (longer 3' UTR). A gene with consistently
higher PDUI in one cluster compared to another is a candidate for 3' UTR
lengthening in that cluster.

Sanity checks: the `proximal_reads + distal_reads` should be non-zero for
most expressed genes. A high fraction of NaN rows usually means `pseudocount=0`
and many zero-count cells; set `pseudocount=1.0` if you want to retain those
rows.

### Limitations

- Uses only the two extreme PAS (proximal and distal). Genes with intermediate
  PAS that shift between clusters are invisible to this metric.
- For genes with many transcript isoforms, the `per_gene` mode collapses to a
  single pair which may not reflect any real transcript.
- Pseudo-replication: PDUI is computed from aggregated read counts per PAS,
  not from per-cell distributions. For cell-level statistical tests, use the
  `proportion` output as input to a per-cell model.

---

## `proportion` {#proportion}

**Question:** For each gene and cell, what fraction of reads falls at each
individual poly(A) site?

**Output file:** `proportion.tsv`

### Algorithm

1. For each gene (or transcript when `per_isoform`), extract the count
   sub-matrix: rows = PAS belonging to this gene, columns = all cells.
2. For each cell column, divide each PAS count by the gene-level total
   (+ pseudocount):

    ```
    p[i, cell] = (x[i, cell] + pseudocount)
                 / sum_over_j(x[j, cell] + pseudocount)
    ```

3. Proportions sum to 1.0 per (gene, cell). Cells where the gene total is
   zero receive NaN.
4. Output one row per (gene, transcript, PAS, cell).

The implementation avoids full-matrix densification: it operates on
column-wise numpy slices per gene and returns column-oriented numpy arrays
that are concatenated into the final DataFrame. This bounds peak memory to
the final DataFrame size regardless of gene count.

### Inputs and outputs

**Inputs:**

- `count_matrix`: `pd.DataFrame` shape `(n_pas, n_cells)`
- `pas_isoform_map`: `dict[int, list[tuple[str, str, int, int, int]]]`

**Output columns in `proportion.tsv`:**

When `aggregation="per_gene"`:

| Column | Type | Description |
|---|---|---|
| `gene_id` | str | Gene identifier |
| `transcript_id` | str | `"_gene_"` sentinel |
| `pas_id` | int | PAS identifier |
| `rank` | int | Proximal-to-distal rank (0 = most proximal) |
| `cell` | str | Cell barcode |
| `proportion` | float64 | Fraction of gene reads at this PAS for this cell |
| `reads_at_pas` | float | Raw count at this PAS (before pseudocount) |
| `total_reads_gene` | float | Sum across all gene PAS for this cell |

When `aggregation="per_isoform"`, the last column is named
`total_reads_transcript` and `transcript_id` holds the actual transcript ID.

### Tunable hyperparameters

| Parameter | Default | CLI flag | YAML key | Description |
|---|---|---|---|---|
| `pseudocount` | 0.0 | `--pdui-pseudocount` | `pdui_pseudocount` | Added to each count before per-gene normalization. |
| `aggregation` | `"per_isoform"` | `--isoform-agg` | `isoform_agg` | `"per_gene"` or `"per_isoform"`. |

### Interpretation

This output feeds the `proportion_heatmap` visualization, which selects the
top-N most variable PAS across clusters and renders a PAS x cluster heatmap.
A PAS with high proportion in cluster A and low proportion in cluster B is a
candidate for differential testing with `fisher` or `nb_pairwise`.

The `reads_at_pas` and `total_reads_gene` columns let you audit a proportion
value: `proportion = 0.8` with `reads_at_pas = 2, total_reads_gene = 2.5`
indicates very low coverage.

### Limitations

- Proportions are not corrected for library size differences between clusters.
  A gene with 100x more reads in cluster A than cluster B will have the same
  proportions if the PAS distribution is identical, which is the intended
  behavior.
- The full long-format output can be very large for datasets with many cells
  and many PAS. At 10,000 cells x 20,000 PAS x 3 PAS per gene on average,
  the TSV can reach several GB. Consider using `per_gene` aggregation and
  filtering to expressed genes before downstream use.

---

## `shannon` {#shannon}

**Question:** How heterogeneous (diverse) is the poly(A) site usage for each
gene in each cell?

**Output file:** `entropy_shannon.tsv`

### Algorithm

1. For each gene (or transcript), extract counts and compute per-PAS
   proportions within each cell.
2. Compute Shannon entropy per cell:

    ```
    H = -sum_i( p_i * log2(p_i) )    (in bits)
    ```

   where `p_i = x_i / sum_j(x_j)` and the convention `0 * log2(0) = 0`
   is applied.

3. Normalize to [0, 1]:

    ```
    H_norm = H / log2(N)
    ```

   where `N` is the number of PAS in the gene or transcript. When `N = 1`,
   `log2(1) = 0` and `H_norm` is set to NaN.

4. Boundary cases from the code:
   - All reads at one PAS: `H = 0`, `H_norm = 0`
   - Uniform usage across N PAS: `H = log2(N)`, `H_norm = 1`
   - Zero-count cell: `H = NaN`, `H_norm = NaN`

### Inputs and outputs

**Inputs:**

- `count_matrix`: `pd.DataFrame` shape `(n_pas, n_cells)`
- `pas_isoform_map`: `dict[int, list[tuple[str, str, int, int, int]]]`

**Output columns in `entropy_shannon.tsv`:**

| Column | Type | Description |
|---|---|---|
| `gene_id` | str | Gene identifier |
| `transcript_id` | str | Transcript ID or `"_gene_"` |
| `pas_ids` | str | Semicolon-joined PAS IDs included in entropy computation |
| `cell` | str | Cell barcode |
| `entropy` | float64 | H in bits |
| `normalized_entropy` | float64 | H / log2(N) in [0, 1] or NaN |
| `n_pas` | int | Number of PAS in the gene/transcript |
| `total_reads_gene` | float | Sum of reads across all gene PAS (per_gene mode) |

When `aggregation="per_isoform"`, the last column is `total_reads_transcript`.

### Tunable hyperparameters

| Parameter | Default | CLI flag | YAML key | Description |
|---|---|---|---|---|
| `pseudocount` | 0.0 | `--pdui-pseudocount` | `pdui_pseudocount` | Added to each count before computing proportions. Non-zero values shift entropy toward uniformity. |
| `aggregation` | `"per_isoform"` | `--isoform-agg` | `isoform_agg` | `"per_gene"` or `"per_isoform"`. |

### When entropy is informative

Entropy is most informative for genes with three or more PAS. For single-PAS
genes, `H = 0` always (there is no choice) and `H_norm = NaN`. For two-PAS
genes, entropy reduces to a single number that is monotonically related to the
classic PDUI, so `classic` is the simpler choice.

Entropy becomes interesting in pseudotime analysis or when comparing
progenitor vs. differentiated cells: progenitors often show higher entropy
(broader PAS usage) while terminally differentiated cells show lower entropy
(focused usage at one site). Compare `entropy_distribution` violin plots
across clusters as a first screen.

### Limitations

- Adding a pseudocount shifts entropy toward uniformity. A pseudocount of 1.0
  on a gene with 2 reads (1 at each of 2 PAS) will produce the same entropy
  as a gene with 10,000 reads (5,000 at each). Filter low-coverage cells
  before relying on entropy values.
- Entropy does not indicate which PAS is used more, only how evenly usage
  is spread. Combine with `proportion` output to identify the dominant PAS in
  low-entropy cells.
