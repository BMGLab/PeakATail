# PeakATail

<p align="center">
    <img src="https://raw.githubusercontent.com/BMGLab/PeakATail/develop/docs/assets/logos/main_logo.png" alt="PeakATail — coiled snake with poly(A) tail" width="640">
</p>

PeakATail is a Python tool for single-cell poly(A) site (PAS) detection and alternative polyadenylation (APA) analysis. It works with any scRNA-seq BAM that carries a `CB:Z` (corrected cell barcode) tag — STARsolo, CellRanger, Alevin-fry, or any aligner that emits the standard 10x-style tag schema. UMI (`UB:Z`) tags are NOT required or used: PeakATail counts raw read 3'ends, not UMI-deduplicated molecules. PeakATail does NOT correct barcodes; your aligner must apply a barcode whitelist (e.g. STARsolo's `--soloCBwhitelist`). From the input BAM it calls polyadenylation sites at the read level, builds a per-cell PAS count matrix, clusters cells by their APA profiles, and provides downstream analyses including differential APA testing between clusters, 3'UTR length quantification, and cross-dataset cluster matching. Developed at BMGLab.

[![Docs](https://img.shields.io/badge/docs-bmglab.github.io%2FPeakATail-blue)](https://bmglab.github.io/PeakATail/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](pyproject.toml)

<!-- TODO[verify]: A CI workflow badge is not included because the repository currently contains only a docs-deploy workflow (.github/workflows/docs.yml). Add a badge once a CI workflow (e.g. .github/workflows/ci.yml) is committed. -->

---

## Quick install

### System prerequisites

PeakATail calls `samtools` and `bedtools` as subprocesses. Both must be on your `PATH` before running the pipeline.

```bash
# Debian/Ubuntu
sudo apt-get install samtools bedtools

# macOS (Homebrew)
brew install samtools bedtools
```

The code has been tested against samtools >= 1.10 (enforced at runtime in `ema/datasets/manager.py`).

### Python version

Python **3.11 or later** is required (declared in `pyproject.toml` as `requires-python = ">=3.11"`).

### Install from source

```bash
git clone https://github.com/BMGLab/PeakATail.git
cd PeakATail

# Create and activate a virtual environment
python3.11 -m venv .venv
source .venv/bin/activate

# Install with uv (recommended — matches the development environment)
pip install uv
uv pip install -e .

# Or with plain pip
pip install -e .
```

All Python dependencies are declared in `pyproject.toml` and installed automatically.

---

## Quick start

```bash
# Step 1: copy the example config and edit paths
cp example.yaml my_run.yaml
# Edit my_run.yaml — set datasets[].bams, gtf, seqlen, cb_len, barcode_tag

# Step 2: run the full pipeline (peak calling → annotation → clustering)
uv run ema run --config my_run.yaml --threads 4

# Step 3: differential APA between every cluster pair
uv run ema switch diff \
    -i peakatail_runs/<run>/per_dataset/<ds>/clusters.h5ad \
    --pasbed peakatail_runs/<run>/per_dataset/<ds>/pasbed.bed \
    --strategy fisher

# Step 4: 3'UTR length quantification across clusters
uv run ema switch length \
    -i peakatail_runs/<run>/per_dataset/<ds>/clusters.h5ad \
    --strategy classic
```

The pipeline writes all output under `peakatail_runs/<name>_<timestamp>/`. A typical run directory looks like:

```
peakatail_runs/emaout_20240501_143022/
    run_config.json                  # full resolved parameters
    peakatail_<ts>.log               # structured run log
    per_dataset/
        <dataset_id>/
            raw/                     # pre-filter BEDs and matrices
            posbed.bed               # filtered positive-strand PAS BED
            negbed.bed               # filtered negative-strand PAS BED
            pasbed.bed               # combined filtered PAS BED (both strands)
            filtered_cb.tsv          # barcodes that passed min_read filter
            pas_gene.tsv             # PAS-to-gene mapping table
            annotatedpas.bed         # pasbed.bed extended with gene_id column
            annotated_matrix.mtx     # post-annotation count matrix (MatrixMarket)
            annotated_pas_ids.tsv    # row index for annotated_matrix.mtx
            annotated_cells.tsv      # column index for annotated_matrix.mtx
            preprocessed.h5ad        # filtered AnnData before clustering
            clusters.h5ad            # AnnData with leiden cluster labels
    switch_diff_<ts>/                # output of `ema switch diff` (auto-routed)
    switch_length_<ts>/              # output of `ema switch length` (auto-routed)
    switch_match_<ts>/               # output of `ema switch match` (auto-routed)
    switch_geneview_<ts>/            # output of `ema switch geneview` (auto-routed)
```

---

## What PeakATail produces

### Peak calling stage

| File | Location | Description |
|---|---|---|
| `raw/pos.bed` | `per_dataset/<ds>/raw/` | Unfiltered positive-strand PAS calls (concatenated across BAM replicates) |
| `raw/neg.bed` | `per_dataset/<ds>/raw/` | Unfiltered negative-strand PAS calls |
| `raw/pas.bed` | `per_dataset/<ds>/raw/` | Union of raw pos + neg BED (pre-filter) |
| `raw/pos.mtx` | `per_dataset/<ds>/raw/` | Raw count matrix, positive strand (MatrixMarket) |
| `raw/neg.mtx` | `per_dataset/<ds>/raw/` | Raw count matrix, negative strand (MatrixMarket) |
| `raw/cb.tsv` | `per_dataset/<ds>/raw/` | Raw cell-barcode index aligned to raw MTX columns |
| `posbed.bed` | `per_dataset/<ds>/` | Filtered positive-strand PAS BED |
| `negbed.bed` | `per_dataset/<ds>/` | Filtered negative-strand PAS BED |
| `pasbed.bed` | `per_dataset/<ds>/` | Combined filtered PAS BED (both strands); used as input for downstream commands |
| `filtered_cb.tsv` | `per_dataset/<ds>/` | Barcodes that passed the `min_read` filter, with filter threshold in header |

### Annotation stage

| File | Location | Description |
|---|---|---|
| `pas_gene.tsv` | `per_dataset/<ds>/` | Two-column table: `pas_id`, `gene_id` |
| `annotatedpas.bed` | `per_dataset/<ds>/` | `pasbed.bed` extended with a trailing `gene_id` column |
| `annotated_matrix.mtx` | `per_dataset/<ds>/` | MatrixMarket sparse count matrix (rows = annotated PAS, cols = cells) |
| `annotated_pas_ids.tsv` | `per_dataset/<ds>/` | Row index for `annotated_matrix.mtx` |
| `annotated_cells.tsv` | `per_dataset/<ds>/` | Column index for `annotated_matrix.mtx` |

### Clustering stage

| File | Location | Description |
|---|---|---|
| `preprocessed.h5ad` | `per_dataset/<ds>/` | AnnData after cell/PAS filtering, before cluster labels are assigned |
| `clusters.h5ad` | `per_dataset/<ds>/` | AnnData with `leiden` cluster labels in `.obs`; primary input for all `switch` subcommands |

### Differential APA stage (`ema switch diff`)

| File | Location | Description |
|---|---|---|
| `diff_<c1>_vs_<c2>.tsv` | `switch_diff_<ts>/` | Per-cluster-pair differential APA results table |
| `cluster_match.tsv` | `switch_match_<ts>/` | Cross-dataset cluster correspondence scores (`ema switch match`) |
| `pdui_classic.tsv` | `switch_length_<ts>/` | Per-cell PDUI scores when strategy is `classic` |
| `proportion.tsv` | `switch_length_<ts>/` | Per-cell per-PAS proportion scores when strategy is `proportion` |
| `entropy_shannon.tsv` | `switch_length_<ts>/` | Per-cell Shannon entropy scores when strategy is `shannon` |

---

## Available CLI commands

All commands are accessed through the `ema` entry point installed by `pip install -e .`.

| Command | Description | Docs |
|---|---|---|
| `ema run` | Run the full pipeline: peak calling, annotation, clustering | [cli/run](https://bmglab.github.io/PeakATail/cli/run/) |
| `ema reannotate` | Branch a finished run into a new trim/filter/clustering variant without re-peak-calling | [cli/reannotate](https://bmglab.github.io/PeakATail/cli/reannotate/) |
| `ema switch diff` | Differential APA test across cluster pairs (Fisher / NB regression) | [cli/switch-diff](https://bmglab.github.io/PeakATail/cli/switch-diff/) |
| `ema switch length` | 3'UTR shortening/lengthening quantification (PDUI variants) | [cli/switch-length](https://bmglab.github.io/PeakATail/cli/switch-length/) |
| `ema switch trend` | Ordered-covariate (e.g. stage-progression) APA-length trend: slope + Spearman + direction | [cli](https://bmglab.github.io/PeakATail/cli/) |
| `ema switch combine` | Stitch stage/celltype-labelled `clusters.h5ad` files into one grouped h5ad for cross-group testing | [cli](https://bmglab.github.io/PeakATail/cli/) |
| `ema switch match` | Cross-dataset cluster matching | [cli/switch-match](https://bmglab.github.io/PeakATail/cli/switch-match/) |
| `ema switch geneview` | Gene-track visualisation: per-cluster PAS coverage and proportions | [cli/switch-geneview](https://bmglab.github.io/PeakATail/cli/switch-geneview/) |
| `ema collapse` | Pool `samtools merge` RG-suffixed run tags back into per-library cells | [cli](https://bmglab.github.io/PeakATail/cli/) |
| `ema merge` | Merge multiple BAM files into one sorted and indexed BAM | [cli/merge](https://bmglab.github.io/PeakATail/cli/merge/) |
| `ema parse-gtf` | Pre-warm the GTF cache so subsequent runs start immediately | [cli/parse-gtf](https://bmglab.github.io/PeakATail/cli/parse-gtf/) |
| `ema wizard` | Interactive setup wizard (also invoked by bare `ema`) | [cli/wizard](https://bmglab.github.io/PeakATail/cli/wizard/) |

Use `--list-strategies` on `ema run`, `ema switch diff`, `ema switch length`, and `ema switch match` to see the strategies registered in the current installation.

Every command accepts `--help` for full flag documentation.

---

## Available strategies

### Peak calling (`ema run --peak-strategy`)

| Name | Description |
|---|---|
| `original` | Wraps the original `pasfind()` logic unchanged; used as a baseline for comparison |
| `lambda_poisson` | MACS2-style: Poisson p-value against a local lambda estimated from floor positions |
| `sierra_iterative` | Sierra-style iterative peak subtraction; finds multiple PAS per UTR |
| `lambda_gradient` | Local-lambda estimation combined with gradient-based peak delineation |

### Clustering (`ema run --cluster-method`)

| Name | Description |
|---|---|
| `leiden_tfidf` | TF-IDF normalisation followed by LSI dimensionality reduction and Leiden community detection |
| `leiden_libsize` | Library-size normalisation followed by PCA and Leiden community detection |
| `external` | Load pre-computed cluster labels from a file instead of running clustering |

### Cross-dataset cluster matching (`ema switch match --strategy` / `ema run --match-method`)

| Name | Description |
|---|---|
| `marker_overlap` | Match clusters by overlap of top marker PAS sets |
| `jaccard` | Match clusters by Jaccard similarity of cell sets |
| `mnn` | Mutual nearest neighbours in a shared LSI embedding |

### 3'UTR length quantification (`ema switch length --strategy`)

| Name | Output file | Description |
|---|---|---|
| `classic` | `pdui_classic.tsv` | Classic 2-PAS PDUI: distal_count / (proximal + distal); proximal/distal selected by genomic or transcript-coordinate rank |
| `proportion` | `proportion.tsv` | Per-PAS proportion vector: reads per PAS as fraction of gene total per cell |
| `shannon` | `entropy_shannon.tsv` | Shannon entropy of the per-PAS proportion distribution per (gene, cell) |

### Differential APA testing (`ema switch diff --strategy`)

| Name | Description |
|---|---|
| `fisher` | Fisher's exact test on per-PAS read counts between two cluster groups |
| `nb_pairwise` | Negative binomial regression, pairwise cluster comparison |
| `nb_multi` | Negative binomial regression, multi-condition |

---

## Documentation

Full docs: **https://bmglab.github.io/PeakATail/**

The MkDocs documentation site is built automatically from the `develop` and `main` branches via the `.github/workflows/docs.yml` GitHub Actions workflow and published to GitHub Pages.

---

## Repository layout

```
PeakATail/
    ema/                    # core Python package (entry point: ema.cli:main)
        cli/                # Click subcommands and config schema
        strategies/         # peak-calling strategy registry
        clustering/         # clustering strategies and cross-dataset matching
        quantification/     # PDUI / proportion / entropy strategies
        switch_test/        # differential APA testing strategies
        outputs.py          # all file I/O for the pipeline
        main.py             # pipeline orchestration
        downstream_runner.py# per-dataset worker (safe for multiprocessing)
    data/                   # example data files
    test/                   # test suite (pytest)
    other_repos/            # reference implementations (Sierra, SCAPE, scTail, etc.)
    example.yaml            # canonical YAML config template
    combined_polya_scrna_methods.csv  # comparison of poly(A) scRNA methods
    pyproject.toml          # package metadata, dependencies, entry points
    ROADMAP.md              # development and publication roadmap
    LICENSE                 # MIT License
    Dockerfile              # container build
```

---

## Citation

If you use PeakATail, please cite the software entry below for now. The `author`
field mirrors the package author declared in `pyproject.toml`; the full
manuscript author list and a machine-readable `CITATION.cff` will be added with
the paper release (tracked in the release-engineering issue).

<!-- TODO[verify]: No CITATION.cff file exists yet. Replace this software entry
with the published journal citation — and expand `author` to the full manuscript
author list — once the paper is out. -->

```bibtex
@software{peakatail,
  author  = {Amiri Tabat, Amir},
  title   = {{PeakATail}: single-cell poly(A) site detection and APA analysis},
  url     = {https://github.com/BMGLab/PeakATail},
  version = {0.3.0},
  note    = {Preprint in preparation; replace with the journal citation when available}
}
```

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for the full text. `pyproject.toml` declares the matching SPDX identifier (`license = "MIT"`).
