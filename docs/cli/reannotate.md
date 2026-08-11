# `ema reannotate`

`ema reannotate` branches a completed `ema run` into a **new trim / filter /
clustering variant without re-running peak calling**. Peak calling (streaming
the BAMs) is the expensive stage; the trim (`find_close`) and everything
downstream of it (annotate → preprocess → cluster) is cheap and depends only on
artifacts a base run already wrote to disk. So a parameter sweep over
gene-distance, cell/PAS filters, or clustering settings can reuse **one** set of
raw peak calls and branch it many times.

It reuses the exact same internals as `ema run`'s downstream section
(`find_close` + `run_one_dataset_downstream` + provenance reconcile + manifest
writing), so a branch is behaviourally identical to having run `ema run` with
those parameters — it just skips peak calling.

!!! note "When to use it"
    - You have a finished `ema run` and want to test different
      `--max-gene-distance`, `--min-cells`, `--min-pas-per-cell`,
      `--resolution`, `--n-neighbors`, or `--cluster-method` values **without
      paying for peak calling again**.
    - You are running a sensitivity sweep (OFAT) over filter/clustering
      parameters from a single base peak-call.

!!! warning "When NOT to use it"
    - You want to change a **peak-calling** parameter (strategy, lambda,
      prominence, merge). Those change the raw peaks, so you must re-run
      `ema run`.
    - The base run has no `unified/concatenated.mtx` + `posbed.bed`/`negbed.bed`
      (i.e. it wasn't produced by a recent `ema run`).

## The output is a complete, chainable run

`--out` is written as a full run directory — indistinguishable from a base
run's downstream output — so the **new matrix continues straight into the next
steps**:

- per-dataset `07_clustering/<ds>/clusters.h5ad` (the new matrix + clustering)
- E3 provenance: `provenance/by_dataset/{pas,cell}_ledger.tsv` + `reconcile_summary.json`
- E2 `run_manifest.json` (hub-indexable)
- `03_gtf_annotation/<ds>/annotatedpas.bed`, root `pasbed.bed`, `branch_manifest.json`

So `ema switch diff`, `ema switch length`, `ema switch trend`, and the hub all
consume a reannotated run exactly as they would a base run.

## Quick example

```bash
# Branch a base run into a 3000 bp gene-distance + resolution-0.5 variant
uv run ema reannotate \
  --base-run peakatail_runs/emaout_2026-05-11_120000 \
  --out      peakatail_runs/branch_gd3000_res0p5 \
  --gtf      Homo_sapiens.GRCh38.99.gtf \
  --max-gene-distance 3000 \
  --resolution 0.5 \
  --threads 8
```

## Flags

Every trim / filter / clustering knob `run_one_dataset_downstream` accepts is a
CLI flag — nothing is pinned internally.

### Inputs
| Flag | Type | Default | Meaning |
|---|---|---|---|
| `--base-run` | DIR | — (required) | A completed `ema run` output dir to branch from. |
| `--out` | DIR | — (required) | Where to write the new variant run. |
| `--gtf` | FILE | — (required) | Same GTF as the base run. |

### Trim (`find_close`)
| Flag | Type | Default | Meaning |
|---|---|---|---|
| `--max-gene-distance` | INT | 5000 | Max distance (bp) for gene-end annotation. |
| `--utr-multiplier` | FLOAT | 2.0 | 3'UTR extension multiplier. |
| `--include-extended` | flag | off | Keep PAS in the extended (TIER_3) region. |

### Cell / PAS filters
| Flag | Type | Default | Meaning |
|---|---|---|---|
| `--min-read` | INT | 1500 | Minimum reads per cell barcode. |
| `--min-cells` | INT | 3 | Minimum cells expressing a PAS. |
| `--min-pas-per-cell` | INT | 50 | Minimum PAS per cell. |

### Clustering
| Flag | Type | Default | Meaning |
|---|---|---|---|
| `--cluster-method` | TEXT | leiden_tfidf | `leiden_tfidf`, `leiden_libsize`, or `external`. |
| `--resolution` | FLOAT | 1.0 | Leiden resolution. |
| `--n-neighbors` | INT | strategy default | kNN graph neighbours (30 tfidf / 10 libsize). |
| `--n-pcs` | INT | 40 | Principal components. |
| `--n-svd-components` | INT | 50 | SVD components (LSI). |
| `--n-top-hvg` | INT | 2000 | Top HVGs (leiden_libsize only). |
| `--tfidf-scale-factor` | FLOAT | 10000 | TF-IDF scale factor (leiden_tfidf only). |
| `--depth-corr-threshold` | FLOAT | 0.75 | Depth-correlation drop threshold (leiden_tfidf only). |
| `--external-clusters` | PATH | — | Cluster labels TSV (`--cluster-method external`). |
| `--random-seed` | INT | 42 | Clustering RNG seed. |

Plus the shared `--threads`, `--config`, and logging flags documented on the
[CLI index](index.md).
