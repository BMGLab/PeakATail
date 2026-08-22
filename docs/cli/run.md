# `ema run`

`ema run` executes the complete PeakATail pipeline: it reads one or more BAM
files, calls poly(A) sites (PAS) per strand, builds cell-by-PAS count matrices,
annotates PAS with gene identities from a GTF file, filters low-quality
barcodes, and clusters cells using TF-IDF + LSI + Leiden. One `clusters.h5ad`
file is produced per dataset and placed inside a timestamped output directory.
All switch subcommands (`diff`, `length`, `match`, `geneview`) consume these
files as their primary input.

!!! note "When to use it"
    - You have one or more 10x Chromium scRNA-seq BAM files and want to discover
      PAS and cluster cells by 3' UTR usage.
    - You want a single command that produces all intermediate files (BEDs, count
      matrices, AnnData) for downstream analysis.
    - You are running a new sample or reprocessing with changed peak-calling or
      clustering parameters.

!!! warning "When NOT to use it"
    - You already have `clusters.h5ad` files from a previous run and only want
      to test differential APA. Use `ema switch diff` directly.
    - You want to change only the clustering resolution without re-running peak
      calling. Re-running `ema run` repeats the entire pipeline; there is
      currently no checkpoint resume.

## Quick example

```bash
uv run ema run \
  --config example.yaml \
  --threads 8 \
  --peak-strategy lambda_gradient \
  --cluster-method leiden_tfidf \
  --resolution 0.8 \
  -vv
```

After this command completes, the following files are on disk (example with
`--output emaout` and dataset id `sample1`):

- `peakatail_runs/emaout_2026-05-11_120000/per_dataset/sample1/clusters.h5ad` — AnnData with Leiden cluster labels in `obs["leiden"]` and gene annotations in `var["gene_id"]`.
- `peakatail_runs/emaout_2026-05-11_120000/per_dataset/sample1/pasbed.bed` — Filtered BED6 PAS coordinates consumed by `switch diff` and `switch geneview`.
- `peakatail_runs/emaout_2026-05-11_120000/per_dataset/sample1/pas_gene.tsv` — Two-column TSV mapping every `pas_id` to a `gene_id`.
- `peakatail_runs/emaout_2026-05-11_120000/run_config.json` — Resolved parameters for reproducibility.
- `peakatail_runs/emaout_2026-05-11_120000/resources.jsonl` — Per-second RSS + CPU samples (requires `psutil`).

## Full `--help` output

```text
Usage: ema run [OPTIONS]

  Run the full PeakATail pipeline.

Options:
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
                                  png+svg+html as appropriate. Examples: 'svg'
                                  / 'png,svg' / 'html'.
  --no-plots                      Disable all plotting (alias for --plot-
                                  engine none).
  --bam-dir PATH                  Single-BAM convenience.
  --bam-files TEXT                Comma-separated multi-BAM list.
  --gtf PATH                      GTF annotation file.
  --atlas PATH                    Reference PAS atlas BED.
  --atlas-distance INTEGER        Atlas snap distance (bp).  [default: 50]
  --seq-len INTEGER               Sequencing read length.
  --cb-len INTEGER                Cell-barcode length (bp).
  --barcode-tag TEXT              BAM tag holding the cell barcode (default
                                  CB).
  --bam-threads INTEGER           pysam decompression threads per BAM.
                                  [default: 4]
  --pipeline                      Pipeline mode (currently only used by
                                  --tiles).
  --batch-size INTEGER            Worker batch size for streaming reads.
                                  [default: 10000]
  --tiles                         Enable tile-based peak calling (parallel
                                  pool).
  --tile-size INTEGER             Tile size in bp (auto if unset).
  --tile-overlap INTEGER          Tile overlap in bp.  [default: 10000]
  --peak-strategy TEXT            Peak-calling strategy (run --list-strategies
                                  to see).  [default: lambda_gradient]
  --lambda-window INTEGER         Background lambda estimation window (bp).
                                  [default: 5000]
  --lambda-method TEXT            Lambda estimator (median / mean / ...).
                                  [default: median]
  --lambda-fold-change FLOAT      Lambda fold-change cutoff for peak
                                  detection.  [default: 2.0]
  --max-pas INTEGER               Maximum PAS sites kept per peak.  [default:
                                  5]
  --smoothing-window INTEGER      Coverage smoothing window (bp).  [default:
                                  50]
  --min-prominence FLOAT          scipy.signal.find_peaks prominence
                                  threshold.  [default: 5.0]
  --dynamic-threshold             Use a per-window dynamic peak threshold.
  --floor-threshold INTEGER       Minimum peak height (clamps dynamic
                                  threshold).  [default: 3]
  --pas-gap INTEGER               Minimum gap between PAS within a peak (bp).
                                  [default: 100]
  --min-pas-spacing INTEGER       Tier-1 (distance) of the post-detection PAS
                                  merger.  Adjacent PAS within one peak whose
                                  gap is below this value are merged uncon-
                                  ditionally.  -1 (default) auto-detects the
                                  median read length per BAM; 0 disables the
                                  distance tier.  [default: -1]
  --min-pas-prominence FLOAT      Tier-2 (valley) static fallback for the
                                  post-detection PAS merger.  Lambda strategies
                                  (lambda_poisson, lambda_gradient) ignore
                                  this and use their own compute_lambda(...).
                                  Non-lambda strategies (original,
                                  sierra_iterative) use this value as the
                                  valley-depth threshold.  Negative disables
                                  Tier 2.  [default: 5.0]
  --ip-filter                     Enable internal-priming filter. Defaults to
                                  --ip-filter-mode annotate (keeps every PAS,
                                  flags A-stretch ones); requires
                                  --genome-fasta.
  --ip-filter-mode [annotate|filter]
                                  annotate KEEPS every PAS and stamps the
                                  internal_priming flag (on annotatedpas.bed);
                                  filter DROPS flagged PAS.  [default:
                                  annotate]
  --genome-fasta PATH             Genome FASTA (.fai indexed) for --ip-filter.
  --annot-filter                  Enable annotation-region filter (drops PAS
                                  not overlapping a gene region; needs --gtf or
                                  --annotation-bed).
  --ip-a-stretch INTEGER          Min consecutive genomic A's downstream of a
                                  PAS to flag it as internal priming.
                                  [default: 6]
  --min-pas-per-cell INTEGER      Minimum PAS per cell (also bridges to
                                  filter_config.min_genes).  [default: 50]
  --min-read INTEGER              Minimum reads per cell barcode.  [default:
                                  1500]
  --min-cells INTEGER             Minimum cells expressing a PAS.  [default:
                                  3]
  --max-gene-distance INTEGER     Max distance for gene-end annotation (bp).
                                  [default: 5000]
  --utr-multiplier FLOAT          3'UTR length multiplier for extended-3'
                                  annotation.  [default: 2.0]
  --include-extended              Include extended-3' annotations.
  --cluster-method TEXT           Clustering strategy (leiden_tfidf /
                                  leiden_libsize / external).  [default:
                                  leiden_tfidf]
  --resolution FLOAT              Leiden resolution.  [default: 1.0]
  --n-pcs INTEGER                 Number of principal components.  [default:
                                  40]
  --external-clusters PATH        Path to external cluster labels TSV.
  --random-seed INTEGER           RNG seed for clustering reproducibility.
                                  [default: 42]
  --match-method TEXT             Cross-dataset cluster match strategy.
                                  [default: marker_overlap]
  --n-top-markers INTEGER         Number of top marker PAS per cluster.
                                  [default: 50]
  --n-neighbors INTEGER           Number of nearest neighbours for the kNN
                                  graph used by Leiden (leiden_tfidf default:
                                  30; leiden_libsize default: 10 — set
                                  --n-neighbors explicitly to override).
                                  [default: 30]
  --tfidf-scale-factor FLOAT      Scale factor for Signac Method 1 TF-IDF
                                  (leiden_tfidf strategy). Default 10000.
                                  [default: 10000.0]
  --depth-corr-threshold FLOAT    Pearson |r| threshold for removing LSI
                                  components correlated with sequencing depth
                                  (leiden_tfidf, ArchR-style). Default 0.75.
                                  [default: 0.75]
  --n-svd-components INTEGER      Number of SVD/PCA components computed before
                                  filtering/neighbor graph (leiden_tfidf and
                                  leiden_libsize strategies). Default 50.
                                  [default: 50]
  --n-top-hvg INTEGER             Number of highly variable genes/PAS selected
                                  before PCA (leiden_libsize strategy only).
                                  Default 2000.  [default: 2000]
  --benchmark                     No-op (use scripts/validate_strategies.py).
  --validate-db PATH              No-op (use scripts/validate_strategies.py).
  --log-level TEXT                Logger level or `name=LEVEL` (comma-
                                  separated).
  --no-log-file                   Don't write peakatail_<ts>.log next to
                                  outputs.
  --no-progress                   Suppress Rich progress bars.
  --list-strategies               Print available strategies and exit.
  --help                          Show this message and exit.
```

## Flags

### Inputs

| Flag | Type | Default | Description |
|---|---|---|---|
| `--config` / `-c` | PATH | — | YAML config file. CLI flags override individual YAML keys. Requires `datasets:` block when used as the primary input method. |
| `--bam-dir` | PATH | — | Convenience flag: treats a single BAM file as a one-dataset run with `id="default"`. Mutually exclusive with `--bam-files` when no `--config` is given. |
| `--bam-files` | TEXT | — | Comma-separated list of BAM paths for a one-dataset run. Same effect as `--bam-dir` but accepts multiple files. |
| `--gtf` | PATH | — | Ensembl or GENCODE GTF annotation file. Used for gene-end and UTR annotation of PAS sites. |
| `--atlas` | PATH | — | Reference PAS atlas BED (e.g. PolyASite v3.0). **Off unless provided.** When provided, PAS are *annotated* against it by default (see `--atlas-mode`), never dropped. |
| `--atlas-distance` | INT | 50 | Match distance in bp: a PAS whose 3' summit is within this distance of an atlas entry is `atlas_match=True`. |
| `--atlas-mode` | TEXT | `annotate` | `annotate` (default) keeps ALL PAS and adds `atlas_match`/`atlas_distance_bp` columns — the unified PAS set is built exactly like a no-atlas run (a pure overlay; novel PAS are never fragmented or dropped). `filter` restores the legacy snap-and-drop (drops PAS with no atlas hit). Use `annotate` when hunting alternative/novel polyA. |
| `--seq-len` | INT | 150 | Sequencing read length. Warns and defaults to 150 if not set. |
| `--cb-len` | INT | 16 | Cell-barcode length in bp. Warns and defaults to 16 if not set. |
| `--barcode-tag` | TEXT | CB | BAM tag carrying the cell barcode. Defaults to `CB` (Cell Ranger convention). |
| `--read-geometry` | `fixed`/`keep`/`true` | `fixed` | How a read's genomic interval is derived. See below. |
| `--read-exclude-flags` | INT | 0 | SAM flag mask vetoed on the coverage/count channel, like `samtools view -F`. `0` = no filtering. See below. |
| `--pas-features` | `off`/`on` | `on` | Append per-site scoring features to `pas_support.tsv`. Adds, drops and moves no PAS. See [`--pas-features`](#pas-features--per-site-scoring-covariates). |
| `--pas-score` | `none`/`calibrated`/`select` | `none` | Calibrated per-site PAS probability. `calibrated` appends a `pas_score` column and moves no call; `select` additionally uses it in place of the molecule-count threshold, inside tier-1 and inside the internal-priming veto. See [`--pas-score`](#pas-score--the-calibrated-per-site-score). |
| `--pas-score-model` | TEXT | `prime1` | Which model `--pas-score` evaluates: a name shipped with the package, or a path to a model JSON. |
| `--pas-score-min` | FLOAT | `-1` | Threshold for `--pas-score select`. Negative means "the threshold the model was shipped with". |

#### `--read-geometry` — how a read becomes an interval

Historically the caller normalised **every** accepted read to exactly
`--seq-len` bp of *reference* span: a read whose span was longer was
**discarded**, and a shorter one had its end rewritten to `start + seq_len`.
A spliced alignment's reference span includes its introns, so the discard fell
almost entirely on spliced reads — before the poly(A) clip detector **and**
before the count matrix.

Measured on the PBMC 10k v3 chr19+21 development slice (50,898,456 records,
`--seq-len 91`):

| | reads | share of valid-CB reads |
|---|---:|---:|
| valid-CB reads reaching the rule | 48,647,964 | 100 % |
| **discarded** because reference span > `--seq-len` | 11,768,752 | **24.19 %** |
| …of those, spliced | 11,542,490 | 98.08 % of the discard |
| **end rewritten** because span < `--seq-len` | 4,752,308 | 9.77 % |
| qualifying poly(A) clip reads lost with the discard | 12,427 | +4.65 % on top of the 267,520 kept |

The rewrite moves a read's 3' end a mean **13.53 bp downstream** of its real
last aligned base (up to 50 bp). Genome-wide the discard is smaller than on
this slice — 13.74 % of valid-CB reads, 96 % of them spliced, **88.8 M reads**
on the full PBMC BAM — because chr19+21 are spliced-richer than average.

| value | what it does |
|---|---|
| `fixed` | **Default.** The historical rule, unchanged: discard span > `--seq-len`, pad shorter reads to `start + seq_len`. Reproduces pre-existing output byte-for-byte. |
| `keep` | Stop discarding, change nothing else — the read still becomes a `seq_len`-long interval. The ablation arm: it isolates "stop throwing reads away" from "stop fabricating the 3' end". |
| `true` | The read's real aligned reference footprint. |

**Why `fixed` is still the default.** `true` was measured against `fixed` on
two development slices. On the default precision arm it is worth **+31.2 %**
(PBMC) and **+23.5 %** (mouse) of raw count-matrix mass — the reads it
restores really were missing from your quantification — but on detection it
costs **2.9 P@100 points on the mouse slice** for +0.6 points of recall, and
gains only +0.3 points of recall on PBMC. The calls it adds do agree with the
curated atlas far better than chance (0.37–0.52 vs a genic null of
0.014–0.031), just less well than the calls already being made, so it slides
along the precision/recall curve rather than lifting it.

**Turn it on when quantification is what you care about** — differential
usage, PDUI, anything reading the count matrix — and leave it off when you
are optimising the called PAS set. It changes coverage, so it changes peak
boundaries, the coverage-only tier and both matrices; expect a full re-run,
about 1.3× the wall time and 1.5× the peak RSS.

What `true` does with each CIGAR operation, exactly:

* **soft clip (`S`)** — excluded at *both* ends. Soft-clipped bases are not
  aligned to the reference, and the terminal poly(A) clip in particular is the
  clip detector's evidence, not coverage: counting it would push a clipped
  read's 3' end *past* its own cleavage site.
* **hard clip (`H`)** — never part of the reference span; nothing to do.
* **deletion (`D`)** — advances the reference and stays **inside** the span.
* **skipped region / intron (`N`)** — **removed** from the span. A spliced read
  does not cover its intron. On 7.82 M valid-CB reads of the slice the introns
  carry 10.02 Gb, about **14×** the real read mass over the same 105 Mb, so
  admitting them would turn the coverage state machine into a gene-body
  detector.

Acceptance is on that **de-introned footprint**, not on the query length. The
distinction matters: a short alignment carrying a long terminal soft clip has
reference span ≤ `--seq-len` but query length > `--seq-len`, so a query-length
rule would *drop* a read the old rule kept — and that shape is precisely a
poly(A) clip read. The footprint rule accepts **every** read `fixed` accepts,
plus the spliced ones whose aligned length fits; a read whose real aligned
footprint genuinely exceeds `--seq-len` (a very long deletion) is still
rejected, exactly as before.

#### `--read-exclude-flags` — multimappers on the coverage channel

The coverage/count channel applies **no** SAM-flag filter by default, so a read
aligned to *N* places contributes *N* reads of coverage and *N* matrix counts.
On the PBMC chr19+21 slice **10.42 %** of valid-CB reads are secondary
alignments (`0x100`); on GSE104556 (STARsolo) the rate is higher still.
`--read-exclude-flags 256` drops secondary alignments; `3844` is samtools'
`unmapped+secondary+qcfail+duplicate+supplementary`. Measured on the default
precision arm, `256` buys **+1.2 P@100 points on PBMC and +2.2 on mouse** for
**−0.1 / −0.3 points of recall** — a clean precision-for-recall trade, which
is a move *along* the curve rather than a lift, so it is **off by default** — and note that PCR duplicates
cannot inflate a *molecule* count in the first place, since duplicate reads
share their `(CB, UMI)` key. The independent clip-evidence channel has its own
switch, `--polya-clip-filter`.

### Output

| Flag | Type | Default | Description |
|---|---|---|---|
| `--output` / `-o` | PATH | `emaout` | Base name for the output directory. A timestamp is appended: `peakatail_runs/<name>_<timestamp>/`. When `--config` provides `output_dir`, the YAML value wins unless `--output` is explicitly set. |

### Concurrency

| Flag | Type | Default | Description |
|---|---|---|---|
| `--threads` | INT | auto | Absolute worker ceiling passed to `ResourceManager`. It bounds **peak calling** (one worker per contig and strand, see `--peak-workers`), the `--tiles` pool and the per-dataset downstream pool. Omitting it lets `ResourceManager` detect available cores. It is *not* a BLAS/OpenMP setting: the clustering libraries read `OMP_NUM_THREADS` from the environment. |
| `--peak-workers` | INT | auto | Worker processes for peak calling. Each job is one (contig, strand) pair called through a region fetch and merged deterministically, so the output is byte-identical to the single-process caller. Defaults to `ResourceManager.get_n_jobs(per_worker_mb=2500)`, i.e. the `--threads` ceiling capped by free RAM. `1` runs the legacy single-process two-pass caller. Requires a BAM index (`.bai`); without one the legacy caller runs and a warning says so. |
| `--bam-threads` | INT | 4 | pysam BGZF decompression threads per BAM reader. Budgeted per worker so `peak_workers × bam_threads` stays within the `--threads` ceiling (with 16 workers and `--threads 16` each worker gets 1). |
| `--batch-size` | INT | 10000 | Worker batch size when streaming reads in tile mode. Lower this if workers are hitting memory limits on large chromosomes. |
| `--tiles` | FLAG | off | Enable tile-based parallel peak calling. Splits chromosomes into overlapping tiles processed by a multiprocessing pool. Recommended for very large BAMs (>5 GB). |
| `--tile-size` | INT | auto | Tile size in bp. When unset, the pipeline picks a value based on chromosome lengths. |
| `--tile-overlap` | INT | 10000 | Overlap between adjacent tiles in bp. Overlap ensures PAS near tile boundaries are not missed. |
| `--pipeline` | FLAG | off | Pipeline mode; currently used only when `--tiles` is active. |

### Peak calling

| Flag | Type | Default | Description |
|---|---|---|---|
| `--peak-strategy` | TEXT | `lambda_gradient` | Algorithm used to call PAS peaks. Run `ema run --list-strategies` to see registered names. `lambda_gradient` is the recommended production strategy (highest precision in benchmark runs); pass `--peak-strategy original` for the unfiltered baseline. |
| `--lambda-window` | INT | 5000 | Window size in bp used to estimate local background signal lambda. Increase for sparse data where the default window may include too few reads. |
| `--lambda-method` | TEXT | `median` | Estimator for lambda within the window. `median` is robust to outliers; `mean` may be inflated by nearby peaks. |
| `--lambda-fold-change` | FLOAT | 2.0 | A region must exceed `lambda * fold_change` to be called as a peak. Raise to 3.0–4.0 to reduce false positives in noisy data. |
| `--max-pas` | INT | 5 | Maximum PAS sites retained per peak region. Rarely needs changing; increase for loci with complex alternative polyadenylation. |
| `--smoothing-window` | INT | 50 | Gaussian smoothing window in bp applied to per-strand coverage before peak finding. Larger values suppress noise at the cost of resolution. |
| `--min-prominence` | FLOAT | 5.0 | `scipy.signal.find_peaks` prominence threshold. Lower to 2.0 to rescue low-coverage PAS; raise to 10.0 to keep only prominent peaks. |
| `--dynamic-threshold` | FLAG | off | Use a per-window dynamic peak height threshold rather than a fixed cutoff. Useful for samples with highly variable library depth across chromosomes. |
| `--floor-threshold` | INT | 3 | When `--dynamic-threshold` is on, this is the minimum peak height. Prevents the dynamic threshold from falling so low that noise is called. |
| `--pas-gap` | INT | 100 | Minimum bp gap between two PAS within the same peak. Increase to merge closely-spaced PAS that likely represent the same site. |
| `--min-pas-spacing` | INT | `-1` | Tier-1 (distance) of the post-detection PAS merger. Adjacent PAS within one peak whose gap < this value are merged unconditionally. `-1` auto-detects the median read length per BAM (e.g. ~98 bp for 10x v2, ~150 bp for v3). `0` disables Tier 1. See [Post-Detection PAS Merger](../strategies/pas-merger.md). |
| `--min-pas-prominence` | FLOAT | `5.0` | Tier-2 (valley depth) of the post-detection PAS merger. Lambda strategies (`lambda_poisson`, `lambda_gradient`) **ignore** this value and use their own `compute_lambda(heights)` instead — fully dynamic. Non-lambda strategies (`original`, `sierra_iterative`) treat this as a static coverage-depth threshold. Negative disables Tier 2. |
| `--cleavage-offset` | INT | `0` | **3' cleavage-site offset correction** (issue #72). Called peak 3' ends stop ~90–105 nt short of the true cleavage site because 10x R2 coverage runs out before the poly(A) junction. When `> 0`, the reported PAS 3' end is shifted **downstream** (strand-aware) by this many bp after peak calling, so tight-cutoff benchmarks and atlas annotation score the inferred cleavage position rather than the coverage edge. A sane data-driven constant is ~90–100 (try `95`). `0` (default) preserves legacy behaviour (no shift). See [3' cleavage offset](#3-cleavage-site-offset-issue-72) below. |

### 3' cleavage-site offset (issue #72)

Atlas-independent motif analysis of PeakATail's calls (Laughney cohort,
22,629 PAS) showed the reported peak 3' end systematically stops **~90–105 nt
short** of the true cleavage site: AATAAA positional density peaks at +75 nt
downstream of the peak end (canonical AATAAA→cleavage spacing 15–30 nt), and
genomic A-fraction crests at +98 nt then cliffs to background — exactly where
10x R2 coverage runs out. Under tight-cutoff benchmarks this offset is scored
as a miss, punishing the *offset* rather than the calls.

`--cleavage-offset N` shifts each reported PAS 3' end downstream by `N` bp
(in the direction of transcription: increasing coordinate on `+`, decreasing,
clamped at 0, on `-`). The correction is applied in place to the per-dataset
strand BEDs immediately after peak calling, so every downstream artifact —
`pasbed.bed`, `annotatedpas.bed`, gene assignment, atlas matching, and the
benchmark harness — uses the inferred cleavage coordinate consistently. The
5' end of each peak (where R2 coverage is real) is preserved.

The offset is chemistry-dependent (R2 read length), so a per-run data-driven
estimate (AATAAA-mode + canonical spacing, or the A-fraction cliff) is
preferred over a constant; that estimator is stubbed in
`ema/countmatrix/cleavage_offset.py::estimate_cleavage_offset` with the
constant `DEFAULT_CLEAVAGE_OFFSET = 95` as the current fallback (see the
`TODO(issue #72)` there). Leave the flag at `0` for legacy behaviour.

### Filters

| Flag | Type | Default | Description |
|---|---|---|---|
| `--min-pas-per-cell` | INT | 50 | Minimum number of distinct PAS detected per cell barcode. Cells below this threshold are excluded. Also bridges to `filter_config.min_genes` in the legacy interface. |
| `--min-read` | INT | 1500 | Minimum total read count per cell barcode. Cells below this are discarded before count matrix construction. Reduce to 500 for low-depth protocols. |
| `--min-cells` | INT | 3 | Minimum number of cells a PAS must be expressed in to survive preprocessing. |

### Poly(A) read evidence

PeakATail reads the non-templated poly(A) tail off the reads themselves: a
read sequenced through the cleavage site carries the tail as a terminal soft
clip (A on `+`, T on `-`). This is **on by default and annotate-only** — each
PAS's clip-read support is written into **BED column 5** of `pasbed.bed`
(historically a hardcoded `0`), and `annotatedpas.bed` inherits it. No
coordinate or count changes unless you also change `--polya-mode` or pick the
[`clip_seeded`](../strategies/peak-calling.md#clip_seeded) strategy, which
*seeds* PAS candidates from clip clusters instead of coverage summits.

At startup the caller reports the observed clip rate and **warns loudly below
0.3% of CB reads** — a pipeline that trims poly(A) before alignment destroys
this evidence, and the run should not be read as clip-supported when it fires.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--polya-evidence` | `on`/`off` | `on` | Collect read-level poly(A) clip evidence during peak calling. `off` restores byte-identical pre-feature output. |
| `--polya-mode` | TEXT | `annotate` | `annotate` keeps every PAS and records support in BED column 5. `filter` drops PAS with zero support (the coverage-only tier) at the same seam as `--ip-filter`. `require` does the same and **fails the run** if nothing clip-supported survives. |
| `--polya-min-clip` | INT | 6 | Minimum terminal soft-clip length, and minimum A/T run flush against the alignment boundary. The adjacency requirement is what buys the measured 92x wrong-end specificity. |
| `--polya-min-purity` | FLOAT | 0.8 | Minimum A (`+`) / T (`-`) fraction across the clipped bases. |
| `--polya-window` | INT | 100 | Radius in bp around a PAS's strand-aware 3' base within which clip reads count as its support. |
| `--polya-seed-window` | INT | 25 | Single-linkage gap for clustering clip sites into candidates (`clip_seeded` only). |
| `--polya-min-umis` | INT | 1 | Minimum distinct `(cell barcode, UMI)` molecules for a clip cluster to be called, and the unit of BED column 5 (`clip_seeded` only). A read with no `UB` tag counts as one molecule; PCR duplicates of one molecule count once. **`--polya-min-reads` is a deprecated alias** — the gate always counted molecules; only the flag name and the score column said "reads". |
| `--polya-clip-filter` | `none`/`f3844` | `none` | Alignment filter on the clip-evidence channel. `f3844` counts only reads passing samtools `-F 3844` (drops secondary / supplementary / duplicate / qcfail / unmapped), which also **drops clusters whose evidence is entirely such alignments** — a call-set change (measured on PBMC: −8.3% of chr19 (+), −27.0% of chr21 (+) tier-1 clusters), so it is opt-in. `none` still de-duplicates by `(barcode, UMI)`, which is what makes PCR duplicates uncountable. Both counts are always in `pas_support.tsv`. |
| `--polya-count-window` | TEXT | `auto,25` | `UP,DOWN` bp, transcript orientation, around a tier-1 cluster's cleavage site. Read ends in `[site-UP, site+DOWN]` that belong to no coverage candidate are counted on the tier-1 row (a cluster inside a coverage peak also takes that peak's counts). `auto` == `--seq-len`, because R2 3' ends pile up just upstream of cleavage (`clip_seeded` only). |

#### `pas_support.tsv` — where the raw counts went

`pasbed.bed` stays plain **BED6**, so column 5 can carry exactly one number
and it carries the one the gate uses (molecules). Everything else is written
to a sidecar next to each caller BED (`<bed>.support.tsv`) and, for
single-BAM runs, merged into `<run>/pas_support.tsv`:

| column | meaning |
|---|---|
| `pas_id` | BED column 4 — join key back to `pasbed.bed` / `annotatedpas.bed` |
| `clip_reads` | raw poly(A) clip reads supporting this PAS |
| `clip_umis` | distinct `(barcode, UMI)` molecules — **equals BED column 5** |
| `clip_reads_f3844` | clip reads passing samtools `-F 3844` |
| `clip_umis_f3844` | distinct molecules among those reads |
| `window_reads` | reads counted into this PAS's count-matrix row |
| `tier` | `1` = clip-seeded cluster, `2` = coverage candidate |

With `--pas-features on` (the default) another 24 columns are **appended**
after these seven — see [`--pas-features`](#pas-features--per-site-scoring-covariates).
Columns are only ever added at the end: the seven above keep their names,
their order and their values, so a reader that indexes them positionally is
unaffected.

The sidecar is a superset of the run-root `pasbed.bed`, which is rewritten
after the cell/count filters; join on `pas_id`.

Comparing `clip_reads` with `clip_umis` is the honest way to see PCR
duplication at a site, and `clip_umis_f3844` shows what the stricter
`--polya-clip-filter f3844` gate would keep — without re-running.

#### `--pas-features` — per-site scoring covariates

`on` (default) appends 24 columns to `pas_support.tsv`. It **changes nothing
else**: no PAS is added, dropped or moved, `pasbed.bed` stays BED6, the count
matrix is untouched, and every pre-existing sidecar column keeps its position
and its bytes. `off` restores the seven-column sidecar exactly.

They exist so a per-site score can be **fitted offline** and applied as a
re-ranker inside the caller's existing gates. Nothing in the pipeline reads
them; they are an annotation.

**Cost.** Two of the columns are computed by the caller from numbers it
already has. The other 22 are computed at the internal-priming stage, inside
the pass that already walks the PAS BED with the genome open — the genome is
opened exactly as many times with the features on as with them off. The
sequence columns need `--genome-fasta`; without one they are written `NA` and
a warning says so. If `--genome-fasta` is supplied but `--ip-filter` is not,
the feature scan *is* that single pass (it drops nothing and writes no BED).

**Written by the caller, per emitted PAS**

| column | meaning |
|---|---|
| `clip_positions` | distinct poly(A) clip **positions** backing this PAS. Tier 1: the members of the single-linkage cluster. Tier 2: the clip positions inside the same ±`--polya-window` its four clip counts come from. |
| `clip_span` | bp between the first and the last of them (`0` when there is at most one). Bounded above by `--polya-seed-window × (clip_positions − 1)` for tier 1. |

**Written from the genome, in transcript orientation**

`r` is the offset from the cleavage base `c` (BED `end − 1` on `+`, BED
`start` on `−`). `r > 0` is **downstream in transcript orientation**, which on
`−` runs toward *lower* genomic coordinates; the window is
reverse-complemented there. The fetched window is `r ∈ [−40, +30]`.

| column | meaning |
|---|---|
| `seq_ok` | `1` if the whole `[−40, +30]` window was readable. `0` at a contig edge or a contig missing from the FASTA (every other sequence column is then `NA`). `NA` — not `0` — when no `--genome-fasta` was supplied at all, so "no genome" stays distinguishable from "edge". |
| `ip_tool_flag` | the caller's **own** internal-priming call for this site — the same boolean `--ip-filter` vetoes on. Emitted **in addition to** the veto, never as a replacement for it. |
| `ip_tool_afrac` | A fraction over the caller's internal-priming window (`--ip-window-left`/`--ip-window-right`, default `r ∈ [−9, +30]`), computed from the very string the veto tested. |
| `ip_tool_arun` | longest A run over that same window. |
| `a_count_d18` | A count in `r +1..+18`. |
| `a_frac_d18` | `a_count_d18 / 18`. |
| `a_run_d18` | longest A run in `r +1..+18`. |
| `a_frac_d30` | A fraction in `r +1..+30`. |
| `a_run_d30` | longest A run in `r +1..+30`. |
| `kin_ip_flag` | `1` when `a_count_d18 ≥ 12` — the "≥12 of 18 downstream A" rule the long-read internal-priming decoy set uses. |
| `hex_strong` | `AATAAA` or `ATTAAA` present in `r −40..−5`. |
| `hex_any12` | any of the 12 canonical hexamers (`AATAAA ATTAAA TATAAA AGTAAA AATACA CATAAA GATAAA AATATA AATAGA AAAAAG ACTAAA AAGAAA`) present in `r −40..−5`. |
| `hex_n_types` | how many distinct canonical hexamers matched. |
| `hex_best_off` | `r` of the **last base** of the 3′-most hexamer hit (always negative; `0` = no hit). |
| `hex_strong_off` | the same, restricted to `AATAAA`/`ATTAAA`. |

`ip_tool_afrac` is the one to look at first: on separating genuine long-read
3′ termini from internal-priming decoys, this single number — inverted — is
a stronger discriminator than the whole 49-feature model that motivated the
column set.

**Written from the PAS BED — same contig, same strand**

The neighbourhood is the candidate set present when the features are
collected: both tiers, before the internal-priming veto drops anything.
(`--polya-mode filter` removes the coverage-only tier *before* this point, so
its context columns describe the smaller set; `run_config.json` records which
you ran.)

| column | meaning |
|---|---|
| `d_prev_cand` | bp to the previous same-strand candidate (`1000000` = none). |
| `d_next_cand` | bp to the next same-strand candidate (`1000000` = none). |
| `n_cand_100` | **other** same-strand candidates within ±100 bp. |
| `n_cand_500` | **other** same-strand candidates within ±500 bp. |
| `mol_500_sum` | sum of BED column 5 over that ±500 bp neighbourhood, this candidate included. |
| `is_local_mol_max` | `1` when no candidate in the neighbourhood has more molecules (ties count as max). |
| `mol_frac_local` | this candidate's share of `mol_500_sum` (`0.0` when the neighbourhood has no molecules at all). |

**What is deliberately *not* emitted.** A second BAM pass for per-cell clip
statistics, end counts or pileup sharpness, and any molecule-end pileup
feature: both were measured and are worth nothing here (0.000–0.002 held-out
AUC for the first; after matching on local read depth, a molecule-end pileup
at a true missed site is as likely as at a random position of the same depth).
Nothing that needs a wider sequence window than `r ∈ [−40, +30]` is emitted
either.

**Multi-BAM runs** re-key their PAS ids when the datasets are merged, so there
is no run-root `pas_support.tsv` to extend; those runs get a standalone
`<run>/pas_features.tsv` with `pas_id` plus the 22 seam columns, in the merged
id space. The per-caller `<bed>.support.tsv` files still carry the two
call-time columns.

#### `--pas-score` — the calibrated per-site score

**Default `none`, which is v2.** Nothing is loaded, nothing is computed and no
column is written unless you ask.

`calibrated` evaluates a model that was fitted **offline** and ships as a table
of constants (`ema/countmatrix/models/pas_score_model_prime1.json`), and
appends one column to `pas_support.tsv`:

| column | meaning |
|---|---|
| `pas_score` | probability that this candidate is a real polyadenylation site, in `[0, 1]`; `NA` when the site could not be scored |

It **changes nothing else**: no PAS is added, dropped or moved, and every other
column keeps its bytes. `select` additionally uses the score as the selection
rule — see below.

**What the score replaces, and what it does not touch.** PeakATail's
precision-first arm is `tier 1 ∩ internal-priming-pass ∩ ≥ 2 clip molecules`.
The score replaces **only the `≥ 2 molecules` threshold**. Tier-1 membership
and the internal-priming veto stay hard gates in front of it, so
`--pas-score select` can only ever *remove* a tier-1 candidate: it cannot
promote a coverage-only tier-2 candidate and it cannot rescue a site the veto
dropped. That is not a stylistic choice — every configuration in which a score
was allowed to override the veto looked spectacular on the curated atlas and no
better against long-read truth.

**What it is made of.** Twenty-one covariates, every one of them a column of
`pas_support.tsv` verbatim: the clip channel (`clip_reads`, `clip_umis`, their
`-F 3844` twins, `tier`, `clip_positions`, `clip_span`), `window_reads`, the
canonical-hexamer block, and the local candidate context. The downstream
A-content columns are deliberately **excluded** — the internal-priming veto
already uses that evidence as a hard gate, and a model that recites the same
rule behind it is double-counting. Because every feature is a sidecar column,
any `pas_score` can be recomputed from the row it sits on.

**Requirements.** `--pas-features on` (the default) and `--genome-fasta`; both
are checked before the run starts, not warned about afterwards. A candidate
whose sequence window could not be read (`seq_ok` `0` — a contig missing from
the FASTA, or a site too close to a contig edge) gets `pas_score` `NA` and is
**exempt** from `select`: "we could not score it" is not spelled the same way
as "we scored it and it lost".

**Cost.** The score rides the pass the internal-priming filter already makes —
no extra pass over the FASTA, over the BAM or over the BED. Evaluating the
shipped 161-tree ensemble is a few seconds per 100 k candidates, in the parent
process, once per run.

**No scikit-learn at run time.** Models are fitted offline by
`scripts/prime/taskD_fit_model.py` and shipped as node arrays evaluated with
numpy; the exporter refuses to write a model whose numpy evaluation differs
from scikit-learn's by more than 1e-9 on any training row (the shipped one
agrees to 2.2e-16). A test asserts, in a subprocess, that loading a model and
scoring with it imports no scikit-learn.

**Where the threshold came from.** `prime1` was fitted on GSE104556 testis
mouse 1 and its threshold is the calibrated decision boundary `p ≥ 0.50` —
chosen without reading any precision, recall or call count, and never on a
dataset it is reported against.

### Internal-priming annotation (D9)

Like the atlas, the internal-priming filter is **off unless enabled** and
**annotates rather than drops** by default — false poly(A) sites caused by the
oligo-dT primer mis-binding genomic A-stretches get an `internal_priming` flag
column (on the PAS ledger + `annotatedpas.bed`) so you can filter downstream.

**`annotate` mode does not change `pasbed.bed`.** It rewrites the internal
pos/neg BEDs unchanged (keep-all, zero rows dropped) and writes the
`internal_priming` flag only to `annotatedpas.bed` — never to `pasbed.bed`.
So a run with `--ip-filter --ip-filter-mode annotate` produces a `pasbed.bed`
(and any benchmark computed from it) **byte-identical** to a run with no
`--ip-filter` at all. The only IP setting that changes `pasbed.bed` is
`--ip-filter-mode filter`, which drops the flagged PAS. Treat the IP axis as a
two-way contrast — keep-all (`annotate` ≡ off) vs `filter` — not three-way; an
explicit "off" arm alongside an `annotate` arm is a duplicate (see issue #69).

| Flag | Type | Default | Description |
|---|---|---|---|
| `--ip-filter` | FLAG | off | Enable the internal-priming check. Requires `--genome-fasta`. |
| `--ip-filter-mode` | TEXT | `annotate` | `annotate` flags A-stretch PAS but keeps them; `filter` drops them. |
| `--genome-fasta` | PATH | — | Genome FASTA (`.fai` indexed) — required with `--ip-filter`; used to read the sequence downstream of each PAS. |
| `--ip-a-stretch` | INT | 6 | Minimum consecutive genomic A's downstream of a PAS (in transcript orientation) to flag it as internal priming. |
| `--ip-a-fraction` | FLOAT | 0.7 | Alternative trigger: flag when the A-fraction of the window reaches this value. |
| `--ip-window-left` / `--ip-window-right` | INT | 10 / 30 | Window (bp) **upstream / downstream of the cleavage site in transcript orientation** examined for the A-stretch, on both strands. |
| `--annot-filter` | FLAG | off | Enable the annotation-region filter (drops PAS that do not overlap a gene region). Requires `--gtf` or `--annotation-bed`. Distinct from `--ip-filter`; it always drops non-overlapping peaks. |

**Strand handling.** Internal priming comes from a genome-encoded A-stretch
*downstream* of the cleavage site in the direction of transcription, so the
window is defined relative to the transcript and mirrored in genomic
coordinates on the `-` strand (`pos` = BED `end` on `+`, BED `start` on `-`):

```text
+ strand:  genomic [pos-left,  pos+right)   ...UUUUU|cleavage>AAAAAA...   scan for A-run / A-fraction
- strand:  genomic [pos-right, pos+left )   ...TTTTTT<cleavage|UUUUU...   reverse-complemented, same scan
```

!!! warning "Versions up to the 4efeb12 line tested the wrong side on `-`"
    Earlier builds applied the `+` genomic window to both strands, so on `-`
    the check covered 30 bp *upstream* / 10 bp downstream of the cleavage
    site in transcript orientation. Runs made with those builds carry an
    `internal_priming` flag (and, in `--ip-filter-mode filter`, a drop set)
    that is wrong for roughly 3–4 % of `-`-strand sites in each direction
    (sites missed and sites wrongly flagged). Re-run the filter
    (`ema reannotate --genome-fasta`) if you use the per-site flag.

### Annotation

| Flag | Type | Default | Description |
|---|---|---|---|
| `--max-gene-distance` | INT | 5000 | Maximum distance in bp from a PAS to a gene 3' end for annotation assignment. PAS farther than this from any gene are left unannotated. |
| `--utr-multiplier` | FLOAT | 2.0 | The 3' UTR region is extended by `utr_multiplier * annotated_UTR_length` for the "extended-3'" annotation category. |
| `--include-extended` | FLAG | off | Include PAS falling in the extended 3' UTR region in the output. By default these are excluded. |

### Clustering

| Flag | Type | Default | Description |
|---|---|---|---|
| `--cluster-method` | TEXT | `leiden_tfidf` | Clustering algorithm. `leiden_tfidf` uses Signac Method 1 TF-IDF + ArchR-style LSI + Leiden. `leiden_libsize` uses library-size normalisation + HVG selection + PCA + Leiden. `external` reads labels from `--external-clusters`. |
| `--resolution` | FLOAT | 1.0 | Leiden resolution parameter. Higher values produce more, smaller clusters. Start at 0.5 for coarse cell types, increase to 1.5–2.0 to split sub-populations. |
| `--n-pcs` | INT | 40 | Number of PCA components used to build the kNN graph for Leiden. Reduce to 20 for small datasets (<1000 cells). |
| `--n-neighbors` | INT | 30 | Number of nearest neighbours in the kNN graph. `leiden_tfidf` defaults to 30; `leiden_libsize` defaults to 10 internally but uses this value when explicitly set. |
| `--n-svd-components` | INT | 50 | Number of SVD/TruncatedSVD components computed before depth-correlation filtering and kNN graph construction. |
| `--tfidf-scale-factor` | FLOAT | 10000 | Scale factor for Signac Method 1 TF-IDF normalisation (`leiden_tfidf` only). Rarely needs changing unless count distributions are unusually skewed. |
| `--depth-corr-threshold` | FLOAT | 0.75 | Pearson \|r\| threshold for removing LSI components correlated with sequencing depth (ArchR-style, `leiden_tfidf` only). Set to 1.0 to disable depth-correlation filtering entirely. |
| `--n-top-hvg` | INT | 2000 | Number of highly variable PAS selected before PCA (`leiden_libsize` only). |
| `--external-clusters` | PATH | — | Path to a TSV with barcode→cluster label assignments. Required when `--cluster-method=external`. |
| `--random-seed` | INT | 42 | RNG seed for Leiden and TruncatedSVD reproducibility. |

### Cross-dataset matching

| Flag | Type | Default | Description |
|---|---|---|---|
| `--match-method` | TEXT | `marker_overlap` | Strategy for assigning canonical cluster IDs across datasets. See [`ema switch match`](switch-match.md). |
| `--n-top-markers` | INT | 50 | Number of top marker PAS per cluster used by `marker_overlap` for cross-dataset Jaccard comparison. |

### Logging and observability

| Flag | Type | Default | Description |
|---|---|---|---|
| `-v` / `--verbose` | COUNT | 0 | `-v` sets DEBUG for the `ema.*` logger namespace; `-vv` sets DEBUG everywhere including third-party libraries. |
| `-q` / `--quiet` | FLAG | off | Suppress INFO messages; show WARNING and above only. Overrides `-v`. |
| `--log-level` | TEXT | INFO | Global level string (`DEBUG`, `INFO`, `WARNING`, `ERROR`) or per-logger override syntax `ema.clustering=DEBUG,ema.annotate=WARNING` (comma-separated). |
| `--no-log-file` | FLAG | off | Do not write `peakatail_<ts>.log` inside the output directory. Useful for CI pipelines where logs are captured via stdout. |
| `--no-progress` | FLAG | off | Suppress Rich live-progress bars. Useful when running in non-interactive shells or capturing output. |
| `--list-strategies` | FLAG | off | Print all registered strategy names across the peak-calling, clustering, and match registries, then exit. |

### Plotting

| Flag | Type | Default | Description |
|---|---|---|---|
| `--plot-engine` | TEXT | `matplotlib` | Comma-separated list of rendering engines. `matplotlib` produces PNG/SVG; `plotly` produces interactive HTML; `both` produces all three; `none` disables plotting. |
| `--plot-format` | TEXT | `all` | Restrict output formats: `png`, `svg`, `html`, or comma list. `all` writes whatever each engine supports by default. |
| `--no-plots` | FLAG | off | Alias for `--plot-engine none`. Disables all figure output and skips the viz hooks entirely. |

### Validation flags (no-ops)

| Flag | Type | Default | Description |
|---|---|---|---|
| `--benchmark` | FLAG | off | No-op. Use `scripts/validate_strategies.py` instead. Emits a warning if set. |
| `--validate-db` | PATH | — | No-op. Use `scripts/validate_strategies.py` instead. Emits a warning if set. |

## Output files

All paths below are relative to the run root
`peakatail_runs/<name>_<timestamp>/`.

**`per_dataset/<id>/raw/pos.bed`** and **`raw/neg.bed`**
: Raw peak-calling output for the positive and negative strands respectively, concatenated across BAM replicates. Unfiltered.

**`per_dataset/<id>/raw/pas.bed`**
: Union of `pos.bed` and `neg.bed`. BED6 format: chrom, start, end, pas_id, score, strand.

**`per_dataset/<id>/pasbed.bed`**
: Filtered and (optionally atlas-snapped) PAS coordinates. BED6. This is the canonical PAS BED consumed by `ema switch diff` and `ema switch geneview`. Source: `ema/outputs.py::write_per_dataset_beds`.

**`per_dataset/<id>/filtered_cb.tsv`**
: Single-column TSV of barcodes passing the `--min-read` filter. Header: `barcode\tmin_read=<n>`. Source: `ema/outputs.py::write_filtered_cb`.

**`per_dataset/<id>/annotated_matrix.mtx`**, **`annotated_pas_ids.tsv`**, **`annotated_cells.tsv`**
: MatrixMarket sparse count matrix (PAS rows, cell columns) after PAS-gene annotation. Row and column indices in the parallel TSV files. Source: `ema/outputs.py::write_annotated_matrix`.

**`per_dataset/<id>/preprocessed.h5ad`**
: Filtered AnnData before clustering. Inspect this to confirm the cell and PAS counts after quality filtering. Source: `ema/outputs.py::write_preprocessed_h5ad`.

**`per_dataset/<id>/clusters.h5ad`**
: Final AnnData with Leiden cluster labels in `obs["leiden"]` and gene annotations in `var["gene_id"]`. Primary input for all `ema switch` subcommands.

**`per_dataset/<id>/pas_gene.tsv`**
: Two-column TSV (`pas_id`, `gene_id`) mapping every PAS to its annotated gene. Source: `ema/outputs.py::write_pas_gene_artifacts`.

**`per_dataset/<id>/annotatedpas.bed`**
: Extends `pasbed.bed` with a trailing `gene_id` column. Useful for IGV inspection.

**`run_config.json`**
: Resolved run parameters with a timestamp. Used for reproducibility and by the pipeline viz hooks.

**`resources.jsonl`**
: Newline-delimited JSON records with schema `{"elapsed_s": float, "rss_gb": float, "cpu_pct": float}`. Written by a background `_ResourceSampler` thread at 5-second intervals. Requires `psutil`; silently absent if not installed.

**`tile_timings.json`**
: Per-tile peak-calling wall times. Only written when `--tiles` is active.

## Performance and resources

Peak calling runs one spawned worker per (contig, strand) and the merge is
deterministic, so **`--threads` changes the wall time and nothing else** —
every output file is byte-identical to a `--peak-workers 1` run (pinned by
`tests/test_chrom_parallel_identity.py` and re-verified end-to-end on the
runs below).

Measured on one machine (dual-socket, `/usr/bin/time -v` peak RSS = the
largest single process, wall = whole `ema run`, `--peak-strategy
clip_seeded`, plots off):

| Dataset | BAM | Cells × PAS | Flags | Wall | Peak RSS |
|---|---|---|---|---|---|
| PBMC 10k v3, chr19+21 slice | 3 GB | 7,121 × 17,968 | `--threads 16` | 6 min 11 s | 1.2 GB |
| Mouse testis (GSE104556, STARsolo) | 15.7 GB | 10,339 × 65,075 | `--threads 12 --ip-filter` | 9 min 03 s | 3.7 GB |
| PBMC 10k v3, full CellRanger BAM | 44 GB | 23,303 × 390,493 | `--threads 16` | 27 min 43 s | 12.4 GB |

Where the time and the memory go on that last (largest) run: peak calling
11 min in 252 parallel jobs plus a 3 min single-threaded merge, cell-barcode
filter 8 min, annotation 1 min, clustering 4 min; the 12.4 GB peak is the
clustering stage, and no peak-calling worker exceeded 3.0 GB.

Rules of thumb:

- **Peak calling** costs 1–3 GB per worker (the deepest contig sets the
  ceiling) and its wall time is bounded by the largest contig — 10 min for
  human chr1 — so more workers than contigs buys nothing. Budget
  `--threads × 2.5 GB`.
- **Everything after it is single-process** and scales with the number of
  non-zeros, not with cells × PAS: the cell-barcode filter holds ~12 bytes
  per non-zero (2.5 GB for the 200 M non-zeros of the full PBMC run) and
  clustering ~3 × nnz × 8 bytes for the annotated matrix.
- **Disk** is the bigger constraint at scale: the raw and filtered
  MatrixMarket files of the full PBMC run are ~7 GB together.

If a run is memory-bound rather than CPU-bound, lower `--peak-workers`
(peak calling is the only stage that scales with it); the outputs do not
change.

## How it relates to other commands

After `ema run` completes, use the per-dataset `clusters.h5ad` and `pasbed.bed`
files as inputs to the switch subcommands:

- **[`ema switch diff`](switch-diff.md)** — differential APA between Leiden cluster pairs.
- **[`ema switch length`](switch-length.md)** — per-cluster PDUI / proportion / entropy quantification.
- **[`ema switch match`](switch-match.md)** — align cluster identities across multiple datasets.
- **[`ema switch geneview`](switch-geneview.md)** — visualise per-cluster PAS usage for specific genes.

## See also

- Strategy overview in [`../strategies/`](../strategies/) — detailed algorithm descriptions for peak calling and clustering.
- Tutorial in [`../tutorials/`](../tutorials/) — step-by-step single-sample and multi-sample walkthroughs.
