# PeakATail

<p align="center" markdown>
  <img src="assets/logos/main_logo.png" alt="PeakATail — coiled snake wordmark with poly(A) tail" width="480" />
</p>

PeakATail detects poly(A) sites (PAS) at single-cell resolution from any scRNA-seq BAM
that carries a cell-barcode (`CB`) tag — STARsolo, CellRanger, Alevin-fry,
or any aligner that emits the standard 10x-style tag schema. UMI (`UB`) tags are
**not** required or used: PeakATail counts raw read 3′ ends, not UMI-deduplicated
molecules. It clusters cells by their
3' UTR usage patterns and tests for differential alternative polyadenylation (APA) between
cell types or conditions. The tool is packaged as the `ema` CLI, installable via `pip` or
`uv`.

!!! tip "Input requirements"
    PeakATail does **not** correct cell barcodes — your aligner must already have applied a
    barcode whitelist. The BAM must carry a `CB:Z` (corrected barcode) tag. UMI (`UB:Z`)
    tags are not read — PeakATail counts raw read 3′ ends rather than deduplicating
    molecules. See [Preparing your BAM](#preparing-your-bam) below for the recommended
    STAR/STARsolo command.

## Why PeakATail

- **Discover APA isoforms per cell type.** Peak-calling runs per-strand on sorted BAM files to identify PAS
  coordinates; each PAS is associated with its nearest gene via GTF-based UTR lookup.
- **Cluster cells by 3' UTR usage.** PAS-by-cell count matrices are processed through TF-IDF + LSI
  dimensionality reduction and Leiden community detection so clusters reflect 3' isoform choice, not
  total expression level.
- **Identify cell-type-specific PAS switching.** `peakatail switch diff` runs Fisher exact tests or negative
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

## Preparing your BAM

PeakATail reads only what's already in your BAM — it does **not** correct
barcodes, demultiplex reads, or align FASTQs. You need a coord-sorted BAM with
the standard 10x-style `CB:Z` (corrected cell barcode) tag. UMI (`UB:Z`) tags
are not read (PeakATail counts raw read 3′ ends, not UMI-deduplicated
molecules), but the STARsolo command below still emits them for compatibility
with other tools. The recommended STARsolo invocation (matches the protocol the
reference benchmarks were run on):

```bash
STAR \
    --runThreadN 16 \
    --genomeDir <STAR_INDEX> \
    --readFilesIn R2.fastq.gz R1.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix sample_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
    --soloType CB_UMI_Simple \
    --soloCBstart 1  --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 10 \
    --soloCBwhitelist /path/to/737K-august-2016.txt
```

!!! warning "Barcode whitelist is required"
    `--soloCBwhitelist` is mandatory. PeakATail trusts the `CB:Z` tag — if you
    skip the whitelist, sequencing errors will appear as thousands of
    spurious "cells". For 10x v2/v3 the whitelists ship with CellRanger
    (`737K-august-2016.txt`, `3M-february-2018.txt`). For Drop-seq /
    inDrops / smart-seq3 / etc., use the protocol-specific whitelist that
    your aligner accepts (see your aligner's docs).

!!! note "What about CellRanger / Alevin-fry / kallisto|bustools?"
    Any aligner that emits a standard `CB` cell-barcode tag works (a `UB` UMI tag
    may also be present but is ignored). CellRanger
    BAMs work out-of-the-box. Alevin-fry emits the same tags via its
    `--sketch` / `--rad`-then-`convert` flow. For salmon/kallisto-bustools you
    need to convert the busfile back to a tagged BAM before passing it to
    PeakATail.

## Quick start

**Step 1 — Run the full pipeline** (peak-calling + clustering):

```bash
peakatail run --config example.yaml
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
peakatail switch diff \
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
