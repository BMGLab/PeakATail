# Tutorials

These tutorials walk through real PeakATail workflows from start to finish. Every command shown here was verified against the `feature/product-cli` branch and the `full_v8_2026-05-11_152746` run directory; copy-paste should work on any machine that has the prerequisites installed.

## Where to start

| Your situation | Start here |
|---|---|
| First time using PeakATail | [Installation](01-installation.md), then [Quickstart](02-quickstart.md) |
| Have multiple BAM files or samples | [Multi-dataset workflow](03-multi-dataset.md) |
| Want to understand a differential result biologically | [Switch analysis](04-switch-analysis.md) |
| Need a single-gene visual summary | [Gene track deep dive](05-gene-track-deep-dive.md) |

## Tutorial map

**[01 — Installation](01-installation.md)**
Prerequisites, Python version, install command, and install verification.

**[02 — Quickstart](02-quickstart.md)**
End-to-end single-sample run in five steps: configure, run, inspect outputs, run differential, read the volcano.

**[03 — Multi-dataset workflow](03-multi-dataset.md)**
Two BAMs in one YAML, atlas snap mode, cross-dataset cluster matching with `ema switch match`.

**[04 — Switch analysis](04-switch-analysis.md)**
Full biological story: locate a significant differential row, decode every column, run `ema switch length`, interpret PDUI values per cluster, and conclude whether 3'UTR shortening or lengthening is happening.

**[05 — Gene track deep dive](05-gene-track-deep-dive.md)**
Use `ema switch geneview` to generate a per-gene figure showing PAS positions across the gene body and per-cluster read proportions. Understand what each panel in the figure means.
