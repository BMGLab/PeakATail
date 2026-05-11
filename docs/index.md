# PeakATail

PeakATail detects poly(A) sites (PAS) at single-cell resolution from 10x Chromium scRNA-seq BAM files,
clusters cells by their 3' UTR usage patterns, and tests for differential alternative polyadenylation (APA)
between cell types or conditions. The tool is packaged as the `ema` CLI, installable via `pip` or `uv`.

## Why PeakATail

- **Discover APA isoforms per cell type.** Peak-calling runs per-strand on sorted BAM files to identify PAS
  coordinates; each PAS is associated with its nearest gene via GTF-based UTR lookup.
- **Cluster cells by 3' UTR usage.** PAS-by-cell count matrices are processed through TF-IDF + LSI
  dimensionality reduction and Leiden community detection so clusters reflect 3' isoform choice, not
  total expression level.
- **Identify cell-type-specific PAS switching.** `ema switch diff` runs Fisher exact tests or negative
  binomial regression across every cluster pair and reports differentially used PAS with FDR control.

## Install

PeakATail requires **samtools** and **bedtools** on your system path, and Python 3.11 or later.

=== "pip (PyPI)"

    ```bash
    pip install peakatail
    ```

=== "uv (recommended)"

    ```bash
    uv pip install peakatail
    ```

=== "From source"

    ```bash
    git clone https://github.com/BMGLab/PeakATail.git
    cd PeakATail
    uv sync          # installs all deps from uv.lock
    # or: pip install -e .
    ```

!!! note "System dependencies"
    ```bash
    sudo apt-get install samtools bedtools   # Debian / Ubuntu
    ```

## Quick start

**Step 1 — Run the full pipeline** (peak-calling + clustering):

```bash
ema run --config example.yaml
```

The `example.yaml` at the repo root shows the full schema. At minimum, provide a `datasets` block
pointing to your BAM file(s), a `gtf` path, and read-geometry parameters:

```yaml
datasets:
  - id: sample1
    merge_strategy: none
    bams:
      - /path/to/cellranger_output/possorted_genome_bam.bam
gtf: /path/to/gencode.v44.annotation.gtf
seqlen: 150
cb_len: 16
barcode_tag: CB
```

**Step 2 — Inspect per-dataset clusters** in the timestamped output directory (default: `emaout/`).
Each dataset gets a `clusters.h5ad` inside `per_dataset/<id>/`.

**Step 3 — Test for differential APA** between cluster pairs:

```bash
ema switch diff \
  --h5ad emaout/per_dataset/sample1/clusters.h5ad \
  --strategy fisher \
  --fdr 0.05
```

Results land in `switch_diff_<timestamp>/` inside the same run directory.

## Documentation map

<div class="grid cards" markdown>

-   **CLI reference**

    ---

    Every `ema` subcommand, flag, and option documented with types,
    defaults, and examples.

    [:octicons-arrow-right-24: CLI reference](cli/)

-   **Strategies**

    ---

    The algorithms behind peak calling, clustering, quantification,
    differential testing, and visualisation.

    [:octicons-arrow-right-24: Strategies](strategies/)

-   **Tutorials**

    ---

    Step-by-step guides: single-sample run, multi-sample atlas mode,
    cluster-pair differential APA, and more.

    [:octicons-arrow-right-24: Tutorials](tutorials/)

-   **Concepts**

    ---

    Data flow through the pipeline, output file layout, and the
    YAML configuration schema.

    [:octicons-arrow-right-24: Concepts](concepts/)

</div>

## Status and citation

PeakATail is developed at the BMG Lab. Source code and issue tracker:
[github.com/BMGLab/PeakATail](https://github.com/BMGLab/PeakATail).

If you use PeakATail in your research, please cite the repository until a
peer-reviewed publication is available.
