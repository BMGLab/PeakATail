# Installation

## System prerequisites

PeakATail calls `samtools` and `bedtools` at runtime via subprocess. Both must be on your `PATH` before you run `ema run`.

=== "Debian / Ubuntu"

    ```bash
    sudo apt-get install samtools bedtools
    ```

=== "macOS (Homebrew)"

    ```bash
    brew install samtools bedtools
    ```

=== "Conda"

    ```bash
    conda install -c bioconda samtools bedtools
    ```

Verify that the binaries are visible:

```bash
samtools --version   # expect 1.17 or later
bedtools --version   # expect 2.30 or later
```

## Python version

PeakATail requires **Python 3.11 or later** (`requires-python = ">=3.11"` in `pyproject.toml`). Check your active interpreter:

```bash
python --version   # must print 3.11.x or higher
```

Using [uv](https://github.com/astral-sh/uv) (recommended) pins the interpreter automatically. If you are on an older system Python, install 3.11 first:

```bash
uv python install 3.11
```

## Install PeakATail

=== "From PyPI (production)"

    ```bash
    pip install peakatail
    # or with uv:
    uv pip install peakatail
    ```

=== "From source (development)"

    ```bash
    git clone https://github.com/BMGLab/PeakATail.git
    cd PeakATail
    uv pip install -e .
    ```

    The `-e` flag installs in editable mode so source changes take effect immediately without reinstalling.

=== "Sync all extras (full dev env)"

    ```bash
    cd PeakATail
    uv sync
    ```

    `uv sync` reads `uv.lock` and installs the exact pinned set of packages, including the `dev` extra (`pytest`, `ruff`, `mypy`).

## Verify the installation

```bash
uv run ema --version
```

Expected output (version number matches `pyproject.toml`):

```
ema, version 0.3.0
```

Running `ema --help` lists every subcommand:

```
Usage: ema [OPTIONS] COMMAND [ARGS]...

  PeakATail — single-cell poly(A) site detection and APA analysis.

Options:
  --version  Show the version and exit.
  -h, --help  Show this message and exit.

Commands:
  merge      Merge BAM files before peak calling.
  parse-gtf  Pre-warm the GTF annotation cache.
  run        Full pipeline: peak-calling → clustering.
  switch     Per-cluster differential APA analyses.
  wizard     Interactive configuration wizard.
```

## Optional: pre-warm the GTF cache

The first `ema run` invocation parses your GTF to build a gene-end BED and a UTR-lengths table. For large GTFs (Ensembl GRCh38 is ~1 GB) this takes a few minutes. You can do it up front so it does not delay your first pipeline run:

```bash
uv run ema parse-gtf --gtf /path/to/Homo_sapiens.GRCh38.99.gtf
```

The cache is stored next to the GTF file in a `gtf_cache/` directory and is reused automatically on subsequent runs.
