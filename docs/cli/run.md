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

### Output

| Flag | Type | Default | Description |
|---|---|---|---|
| `--output` / `-o` | PATH | `emaout` | Base name for the output directory. A timestamp is appended: `peakatail_runs/<name>_<timestamp>/`. When `--config` provides `output_dir`, the YAML value wins unless `--output` is explicitly set. |

### Concurrency

| Flag | Type | Default | Description |
|---|---|---|---|
| `--threads` | INT | auto | Absolute ceiling passed to `ResourceManager`. Wired directly into the singleton before pipeline dispatch; omitting this flag lets `ResourceManager` detect available cores automatically. |
| `--bam-threads` | INT | 4 | pysam decompression threads per BAM file. Increase to 8–16 on machines with fast storage to reduce BAM-read I/O time. |
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

The sidecar is a superset of the run-root `pasbed.bed`, which is rewritten
after the cell/count filters; join on `pas_id`.

Comparing `clip_reads` with `clip_umis` is the honest way to see PCR
duplication at a site, and `clip_umis_f3844` shows what the stricter
`--polya-clip-filter f3844` gate would keep — without re-running.

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
| `--ip-a-fraction` | FLOAT | 0.7 | Alternative: flag if the A-fraction in the window exceeds this. |
| `--annot-filter` | FLAG | off | Enable the annotation-region filter (drops PAS that do not overlap a gene region). Requires `--gtf` or `--annotation-bed`. Distinct from `--ip-filter`; it always drops non-overlapping peaks. |

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
