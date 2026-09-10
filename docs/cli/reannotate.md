# `peakatail reannotate`

`peakatail reannotate` branches a completed `peakatail run` into a **new trim / filter /
clustering variant without re-running peak calling**. Peak calling (streaming
the BAMs) is the expensive stage; the trim (`find_close`) and everything
downstream of it (annotate → preprocess → cluster) is cheap and depends only on
artifacts a base run already wrote to disk. So a parameter sweep over
gene-distance, cell/PAS filters, or clustering settings can reuse **one** set of
raw peak calls and branch it many times.

It reuses the exact same internals as `peakatail run`'s downstream section
(`find_close` + `run_one_dataset_downstream` + provenance reconcile + manifest
writing), so a branch is behaviourally identical to having run `peakatail run` with
those parameters — it just skips peak calling.

!!! note "When to use it"
    - You have a finished `peakatail run` and want to test different
      `--max-gene-distance`, `--min-cells`, `--min-pas-per-cell`,
      `--resolution`, `--n-neighbors`, or `--cluster-method` values **without
      paying for peak calling again**.
    - You are running a sensitivity sweep (OFAT) over filter/clustering
      parameters from a single base peak-call.

!!! warning "When NOT to use it"
    - You want to change a **peak-calling** parameter (strategy, lambda,
      prominence, merge). Those change the raw peaks, so you must re-run
      `peakatail run`.
    - The base run has no `unified/concatenated.mtx` + `posbed.bed`/`negbed.bed`
      (i.e. it wasn't produced by a recent `peakatail run`).

## The output is a complete, chainable run

`--out` is written as a full run directory — indistinguishable from a base
run's downstream output — so the **new matrix continues straight into the next
steps**:

- per-dataset `07_clustering/<ds>/clusters.h5ad` (the new matrix + clustering)
- E3 provenance: `provenance/by_dataset/{pas,cell}_ledger.tsv` + `reconcile_summary.json`
- E2 `run_manifest.json` (hub-indexable)
- `03_gtf_annotation/<ds>/annotatedpas.bed`, root `pasbed.bed`, `branch_manifest.json`

So `peakatail switch diff`, `peakatail switch length`, `peakatail switch trend`, and the hub all
consume a reannotated run exactly as they would a base run.

## One branch, one `--out` (exit code 1 if not)

Every branch **must** have its own `--out` directory. Two `peakatail reannotate`
processes pointed at the same `--out` write the same
`04_pas_gene_assignment/<ds>/pas_gene.tsv`, `05_annotated_matrix/<ds>/*` and
`07_clustering/<ds>/clusters.h5ad` paths, so the artifacts that survive are an
arbitrary interleaving of two different parameter sets. This is not
hypothetical: five branches of a sweep grid that shared a `branch_name` (hence
an `--out`) ran concurrently under Nextflow and produced three different
`pas_gene.tsv` row counts for **identical** declared parameters, grouped by
write time.

So `peakatail reannotate` now **refuses to start** when another live `peakatail reannotate`
already holds the same `--out`. It takes an exclusive `flock` on
`<out>/.ema_reannotate.lock` before any work begins; if that claim fails it
exits **immediately with exit code 1**, writing nothing, and prints who holds
the directory:

```console
$ peakatail reannotate --base-run runs/base --out runs/branch_gd3000 --gtf ref.gtf
Error: --out is already in use by another running `peakatail reannotate`:
/path/runs/branch_gd3000 (held by pid=48211 host=node07 started=2026-09-03T23:34:14)
— two branches writing one output dir interleave their artifacts and silently
corrupt both. Give every branch its OWN --out; a duplicate branch_name in a
sweep grid is the usual cause. (Lock file: /path/runs/branch_gd3000/.ema_reannotate.lock)
$ echo $?
1
```

!!! tip "What to do when you hit this"
    - **Sweep grids are the usual cause.** Two rows resolving to the same
      `branch_name` — and therefore the same `${out_root}/reannotate/<branch_name>`
      — is a duplicate `--out`, not two branches. De-duplicate `branch_name`
      (or fold the differing parameters into it) so every row gets a distinct
      directory, then re-run.
    - **Genuinely want both variants?** Give the second one a different
      `--out`. Branches are cheap; they share the base run's peak calls.
    - **Want to overwrite an earlier branch?** Wait for the holder named in the
      message to finish (or stop it), then re-run into the same `--out`.
    - **Not a stale lock.** The claim is held by the running process's open
      file description, so the kernel releases it the moment that process exits
      — including a crash or a `kill`. A branch that died never blocks the
      re-run, and you never need to delete `.ema_reannotate.lock` by hand. If
      the refusal fires, some `peakatail reannotate` really is alive on that
      directory.

    The lock is advisory and only guards against a second `peakatail reannotate`; it
    does not stop an unrelated program from writing into `--out`.

`.ema_reannotate.lock` is a zero-value bookkeeping file (dot-prefixed, never
listed in `run_manifest.json`/`branch_manifest.json`) and is safe to leave in
place.

## Artifacts are written atomically

Text artifacts — including `pas_gene.tsv`, `annotatedpas.bed`, `pasbed.bed`,
the annotated-matrix `pas_ids`/`barcodes` sidecars, `run_config.json`,
`run_manifest.json` and the per-stage `*_stats.json` — are written to a temp
file in the **same directory** and then `os.replace()`d onto the final path.
A reader therefore never observes a partially written or truncated artifact:
at any instant the final path holds either the complete previous content or
the complete new content, and a failed write leaves the previous file
untouched. Nothing about the file contents changes — only when they become
visible.

## Quick example

```bash
# Branch a base run into a 3000 bp gene-distance + resolution-0.5 variant
uv run peakatail reannotate \
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
| `--base-run` | DIR | — (required) | A completed `peakatail run` output dir to branch from. |
| `--out` | DIR | — (required) | Where to write the new variant run. Must be unique per branch — a second `peakatail reannotate` on the same `--out` is refused with exit code 1, see [One branch, one `--out`](#one-branch-one-out-exit-code-1-if-not). |
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
