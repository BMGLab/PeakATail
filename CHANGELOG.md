# Changelog

## Unreleased — caller memory and CPU

Peak RSS and wall time only: **every output file is byte-identical** at the
default settings. Verified end-to-end by re-running two full benchmark arms
on the same inputs and comparing bytes against the runs they are supposed to
reproduce: the chr19+21 PBMC slice (21 files), GSE104556 mouse1 testis with
`--ip-filter` (18 files) and the full PBMC 10k v3 CellRanger BAM (18 files)
— `pasbed.bed`, `pas_support.tsv`, the raw and filtered matrices, the cell
list, `annotatedpas.bed`, `annotated_matrix.mtx`, `pas_gene.tsv` — plus
`preprocessed.h5ad` and `clusters.h5ad` (`X`, `X_lsi`, Leiden labels) on all
three.

| Run | Wall before | Wall after | Peak RSS before | Peak RSS after |
|---|---|---|---|---|
| PBMC 10k v3, chr19+21 slice (`--threads 16`) | 14 min 41 s | 6 min 11 s | 5.56 GB | 1.16 GB |
| GSE104556 mouse1 testis (`--threads 12 --ip-filter`) | 1 h 01 min | 9 min 03 s | 23.12 GB | 3.68 GB |
| PBMC 10k v3, full BAM (`--threads 16`) | 3 h 45 min 53 s | 27 min 43 s | 293.74 GB | 12.45 GB |

Not all of it is parallelism: the slice re-run with `--peak-workers 1` (the
legacy single-process caller) takes 11 min 27 s / 1.71 GB, with peak calling
at 613 s instead of 733 s and the cell-barcode filter at 25 s instead of 64 s.

### Fixed

- **TF-IDF no longer densifies the count matrix
  (`ema/clustering/strategies/leiden_tfidf.py`).** `_tfidf_signac_method1`
  called `X.toarray()` and then held four dense float64 copies at once
  (`astype`, `tf`, `tf * idf`, `log1p`), i.e. `4 × n_cells × n_PAS × 8`
  bytes. That single function was **>98 % of peak RSS on every dataset
  measured** — the formula predicts 104.3 GB / 140.7 GB / 236.8 GB /
  291.2 GB / 21.5 GB against measured peaks of 105.6 / 140.2 / 239.5 /
  293.7 / 23.1 GB — and most of the 31-minute clustering stage of the PBMC
  run. The same scalar sequence applied to the stored entries needs
  `~3 × nnz × 8` bytes and produces the **same bits**: structural zeros map
  to `log1p(0) == 0` and were dropped by the trailing `csr_matrix()` anyway,
  and the row sums / per-PAS cell counts are integer sums, exact in float64
  in any order. Checked against the previous implementation on the real
  `preprocessed.h5ad` of two runs (CSR structure, data, `X_lsi`, Leiden
  labels all identical): 0.14 s / 0.56 GB vs 3.1 s / 4.5 GB on the slice,
  0.76 s / 1.55 GB vs 17.4 s / 22.2 GB on mouse1.
- **The cell-barcode filter streams integers instead of a per-row string
  frame (`ema/matrixfilter.py`).** The vectorised reader kept an object
  `cb_str` column per non-zero (~95 B/row): 2.71 GB for the slice's 14.9 M
  rows, an estimated 20–27 GB for the ~200 M rows of the full PBMC run —
  the next peak once the TF-IDF was fixed. Files are now parsed once in
  4 M-row chunks into three integer arrays (int32 where the values fit,
  12 B/row) and both passes run on those. Row semantics are unchanged; the
  Pass-1 grouping moves from the raw token to its integer, which is the same
  partition while every token is canonical, so a non-canonical token
  (`"007"`, `"+7"`, `"7.0"`), a >3-token row, a non-integer token or an
  absurd column index now defers the whole call to `_filter_cb_legacy` —
  the byte-for-byte copy of the reference algorithm. On such pathological
  input the result is therefore the reference's, which the previous fast
  path only approximated (it truncated float-valued tokens instead of
  skipping the row).

### Added

- **Peak calling runs one worker per (contig, strand)
  (`ema/countmatrix/chrom_parallel.py`, `--peak-workers`).** It was a single
  pure-Python thread streaming the whole BAM twice — 3 h 01 min of the
  3 h 46 min PBMC run at ~143 % CPU — and `--threads` never reached it
  (it fed only the `--tiles` pool and the downstream pool). Each job calls
  `peak_calling(region=(contig, 0, length))` in a spawned worker with its own
  `BarcodeIndex` and `PeakCallingState`; the merge then walks the jobs in the
  order the sequential run emitted them (`+` strand first, contigs in BAM
  header order, then `-`), renumbers `pas_id` 1..N across both strands, and
  rebuilds the shared barcode index by appending each job's local `cb.tsv` in
  local first-write order — which reproduces the sequential singleton's
  column assignment exactly. Matrices and support sidecars are re-keyed
  through those two maps.
  - Default when the BAM has an index; without one (or with
    `--peak-workers 1`) the legacy single-process caller runs unchanged and
    logs why. Worker count comes from `ResourceManager` (2.5 GB per worker;
    measured 1.5 GB on PBMC chr19 (+)), so `--threads` now governs peak
    calling, and BGZF threads are budgeted per worker so
    `workers × bam_threads` stays inside the ceiling.
  - Wall time is bounded by the largest contig, not by the BAM: 253 s for
    130 jobs on mouse1 (12 workers), and 674 s for the 252 jobs of the full PBMC BAM (16 workers, largest
    worker 3.0 GB) against 3 h 01 min of sequential streaming.
  - Per-job timings, PAS/cell counts and peak RSS land in
    `peakcalling/<id>_<n>.peak_jobs.json`.
- **`--peak-workers`** (YAML `peak_workers`).

### Changed

- `read_check` tests the strand — one flag bit — before the `get_tag(CB)`
  lookup, the `reference_end` CIGAR walk and the composite build; each
  strand pass rejects roughly half the BAM on exactly that test. Every
  rejection returns the same sentinel, so this reorder cannot change a
  result (pinned against a verbatim copy of the previous implementation over
  the full grid of read shapes).
- The `"<RG>_<CB>"` composite is interned per distinct pair, so every
  structure holding it shares one string object per cell — tracemalloc found
  102 MB of duplicate composites alive at a single chromosome flush of the
  slice.
- Region-mode `peak_calling` honours `bam_threads` instead of forcing
  `threads=1`.
- `ema/main.py` passes `default_threshold` / `merge_len` from
  `variable_config` explicitly on the single-process path. They were
  previously left to `peak_calling`'s default arguments, which froze the
  values `variable_config` held when `peackcalling.py` was imported — before
  the CLI/YAML bridge ran — so `--default-threshold` / `--merge-len` were
  silently ignored there (the `--tiles` path always read them live). At the
  defaults (5 / 100) nothing changes; a non-default value is now honoured on
  every path.

### Notes

- `--tiles` remains unusable for a real run and is untouched by this work:
  its workers crash under spawn because `barcode_tag` / `cb_len` / `seqlen`
  are never propagated, it builds jobs for every header contig, and
  `merge_tiles` sorts contigs as strings and unions barcodes per strand.
  The new path deliberately passes every config value in the job spec.
- Still single-process after peak calling: `make_dataframe`'s
  `mmread` → CSC → CSR → CSC chain and `preprocessing`'s transpose are the
  remaining multi-GB copies on a whole-genome run.

## Unreleased — read-level poly(A) evidence

### Added
- **`ema/countmatrix/polya.py`: PeakATail now reads poly(A) evidence off the
  reads.** `clip_site()` is a strand-aware terminal soft-clip A/T detector
  returning the inferred cleavage coordinate. Measured on the pbmc_10k_v3
  BAM: 1.152% of CB-bearing reads carry a qualifying clip, the wrong-end
  control fires on 0.0125% (92x specificity), and 73.9% of clip sites fall
  within 100 bp of a PolyASite 2.0 site (13.9% for arbitrary read 3' ends).
- **BED column 5 of `pasbed.bed` now carries per-PAS clip-read support.** It
  was hardcoded to `0`; every downstream reader already names it `score` and
  ignores it, and `annotatedpas.bed` inherits it. `--polya-evidence off`
  restores the previous byte-identical output.
- **New `clip_seeded` peak strategy — the accuracy fix.** Clip-site clusters
  (single-linkage, `--polya-seed-window`, read-weighted modal position) are
  PRIMARY PAS candidates; coverage peaks that overlap no cluster are emitted
  as an explicit SECOND TIER with `score == 0`. Seeding rather than filtering
  is required because ~60% of the clip evidence lies outside every coverage
  peak window.

  Measured end-to-end on chr19+chr21 of pbmc_10k_v3 and scored with
  `scripts/benchmark_tools/score_tool.py` (detected-gene-restricted
  PolyASite 2.0 atlas, both arms taken through the identical gene-assignment
  stage the shipped benchmark reports at):

  | arm | n | P@100 | R@100 | F1@100 |
  |---|---:|---:|---:|---:|
  | shipped caller (`lambda_gradient`) | 12,928 | 0.1610 | 0.1956 | **0.1766** |
  | `clip_seeded`, both tiers | 19,354 | 0.2371 | 0.3782 | **0.2914** |
  | `clip_seeded`, clip-supported tier only | 11,318 | 0.3591 | 0.3379 | **0.3482** |
  | `clip_seeded`, coverage-only tier | 8,036 | 0.0652 | 0.0404 | 0.0499 |

  Precision and recall both roughly double. The clip-supported tier's recall
  (0.3379) sits at the measured ceiling for this slice — 33.56% of its
  detected-gene atlas sites carry any clip read within 100 bp — so that tier
  extracts essentially all the recall this evidence type can provide, and the
  coverage-only tier is what carries the rest. **Always report the two tiers
  separately, with their n.**
- New flags alongside the `ip_*` block: `--polya-evidence`, `--polya-mode`
  (`annotate` default / `filter` / `require`), `--polya-min-clip`,
  `--polya-min-purity`, `--polya-window`, `--polya-seed-window`,
  `--polya-min-umis` (formerly `--polya-min-reads`, still accepted),
  `--polya-clip-filter`. `filter` and `require` are enforced at the existing
  `_apply_pas_filters()` seam and drop exactly the coverage-only tier.
- **`pas_support.tsv` — per-PAS poly(A) support sidecar.** Written next to
  every caller BED (`<bed>.support.tsv`) and, for single-BAM runs, merged to
  the run root as `pas_support.tsv`, keyed by the BED's own `pas_id`:
  `clip_reads` (raw), `clip_umis` (== BED column 5), `clip_reads_f3844`,
  `clip_umis_f3844`, `window_reads` (the reads counted into that PAS's
  matrix row) and `tier`. `pasbed.bed` stays plain BED6 — no parser sees a
  new column and the name field is untouched. All three peak-calling paths
  write it (the tile merge re-keys it with the same pasnumber remapping the
  BED gets).
- **Startup warning when the observed clip rate is below 0.3% of CB reads** —
  a pipeline that trims poly(A) before alignment destroys this evidence
  channel, and the caller now says so instead of silently emitting an
  unsupported call set.

### Fixed
- **BED column 5 is now DISTINCT MOLECULES, the unit the flag always gated
  on.** `--polya-min-reads` documented "distinct molecules" and
  `ClipSeeder.flush()` did gate on them, but the score column wrote the raw
  clip-READ count — PCR duplicates and secondary alignments included. On the
  Stage-2 PBMC run that made 18,864 of the 79,751 "≥2-read" tier-1 sites
  (23.7%) a single molecule counted twice, and 47.7% of the exactly-2-read
  sites. The gate is unchanged, so **no PAS moves and none is added or
  dropped at the default `--polya-min-umis 1`**: only column 5 changes value
  (proven on the chr19+21 PBMC slice — columns 1–4 and 6 byte-identical to
  the previous commit, all four region BEDs). Raw clip reads remain
  available in `pas_support.tsv`.
  - `--polya-min-reads` → **`--polya-min-umis`**; the old spelling still
    works (CLI and YAML) and logs a deprecation warning.
  - A molecule is `(cell barcode, UMI)`. `read_check` keys cells by
    `RG + "_" + barcode`, so molecules deliberately drop the read-group
    prefix (`molecule_cb`): a BAM with one read group per lane would
    otherwise count one molecule once per lane (+4.7% on the PBMC chr21 (+)
    slice: 13,100 vs 12,483 molecules over identical clusters). Matrix cell
    identity is unchanged.
  - Reads with no `UB` tag still count as one molecule each.
- **Alignment filtering on the clip-evidence channel (`--polya-clip-filter`,
  `clip_read_ok` == samtools `-F 3844`).** `read_check` applies no
  secondary / supplementary / duplicate / qcfail / unmapped filter (it is
  the coverage path's contract; widening it would move every peak), so every
  clip read is now classified and BOTH counts are reported per PAS.
  `--polya-clip-filter f3844` makes the filtered reads stop being evidence
  (score, gate and tier tag all use the `-F 3844` molecule count).
  **The default stays `none`, deliberately**: UMI de-duplication already
  makes PCR duplicates uncountable (duplicates share their `(CB, UMI)` key),
  while dropping flagged reads outright also drops every cluster whose
  evidence is entirely such alignments — measured on PBMC as −8.3% of
  chr19 (+) and −27.0% of chr21 (+) tier-1 clusters (STARsolo emits
  multimappers as secondary; CellRanger flags ~50% duplicates). That is a
  call-set change and must be scored, not slipped into a default; the
  `*_f3844` sidecar columns let it be measured post hoc without a re-run.
- **The tier-1 clip fallback no longer counts a read the neighbouring
  cluster already counted.** A cluster whose candidate partition and
  midpoint-clipped cleavage window are both empty falls back to its own clip
  reads; those reads' `end1` could lie past the midpoint, inside the
  neighbour's window, and were counted twice (verified at `13:43135248+`,
  pas#23403 on mouse1: 37 clip reads with `end1` 43,135,267–43,135,294 past
  its clipped window's 43,135,262 upper bound; the row summed to 371 where
  the shipped candidate + remainder is 334). The fallback is now
  midpoint-clipped and candidate-excluded exactly like the window, so
  "no read end is ever counted twice" is true for every row; a fallback row
  can legitimately end up empty (its reads belong to the neighbour) and the
  new `empty_tier1_rows` stat counts those. On mouse1 this touched 5,250
  rows and ~8k of 140.7M counts (~0.006%).
- **`clip_seeded` tier-1 PAS are now counted from the reads that pile up at
  the cleavage site, not from their poly(A)-clipped reads.** The first
  full-data run of the strategy (GSE104556 mouse1 testis) surfaced a ~10x
  quantification regression: every tier-1 row of the count matrix held only
  the cluster's clip reads (~1% of the reads), so the annotated matrix carried
  4.1M counts against 41.9M from the shipped caller and the `min_read` cell
  filter lost 23% of the real STARsolo cells. The BED coordinates and tier
  tags were — and remain — byte-identical; only the matrix changes:
  - a cluster that suppresses a coverage candidate takes that candidate's
    accumulated `cb_positions` counts, restricted to its partition of the
    candidate (several clusters split a candidate at the midpoints between
    their anchors — the same `partition_peak_region` rule `find_pas` applies
    to multi-PAS peaks — each share counted with `reconstruct_cb_dict`
    semantics on a compact `_CandidateSlice` of the candidate's own
    `cb_positions`, which is what every positional coverage strategy's
    `get_cb_dict_for_pas` delegates to; a non-positional strategy such as
    `original` is detected by the shares not summing to the candidate and
    the candidate then goes whole to the nearest cluster), so the shipped
    mass is conserved exactly;
  - every cluster additionally counts the accepted read ends in its cleavage
    window (`--polya-count-window`, default `auto,25` ==
    `[site - seq_len, site + 25]` in transcript orientation) that belong to
    no coverage candidate, clipped at the midpoint to neighbouring clusters,
    so nothing is counted twice. For a cluster outside every coverage peak
    this is its whole count; inside a peak it recovers the reads the
    coverage caller counted into `cb_positions` but never put in a candidate
    (the streaming loop drops the ends still in `data_array` when a peak
    closes, so every candidate interval stops ~`seq_len` bp short of the
    pile's 3' end — on `+` those are exactly the reads at the cleavage site);
  - the clip-read count stays in BED column 5 as the support annotation.
  Tier-2 rows are untouched. Implemented once in `ClipSeeder.flush()` and
  fed by all three peak-calling paths (`ClipStream` records every accepted
  read end as two int32 arrays per chromosome — 8 bytes/read, released at
  each chromosome flush); `tests/test_polya_three_path_agreement.py` now pins
  per-PAS counts as well as coordinates, and
  `tests/test_polya_tier1_counts.py` pins the accounting on a synthetic BAM
  (a cluster inside a 200-read peak carries the peak, not ~2; an isolated
  cluster with 30 in-window reads + 3 clip reads counts 33; tier-2 rows are
  byte-for-byte the coverage strategy's). Each strand pass logs a mass
  accounting (`clip_seeded counting: {...}`) so a run can show where every
  count came from.

  Measured on GSE104556 mouse1 (`--seq-len 98`, 12 threads, 55 min wall):
  every caller BED (`posbed.bed`, `negbed.bed`, `01_peak_calling/*`,
  `annotatedpas.bed`) is byte-identical to the Stage-2 run. Raw matrices:
  140.7M counts vs 132.9M shipped (1.059x) — the suppressed candidates'
  partitions plus the tier-2 rows reproduce the shipped per-strand totals
  to the read (`+` 58,221,832 + 5,131,499 = 63,353,331; `-` 64,246,477 +
  5,330,874 = 69,577,351), the +7.8M is the cleavage-window remainder.
  `min_read` keeps 10,339 cells including 1,294/1,294 STARsolo cells (the
  clip-only run kept 995). Gene-assigned mass on the shipped run's cells:
  51.8M vs 49.2M shipped when the shipped `filterdmatrix.mtx` is keyed
  correctly (1.053x; per-gene Spearman 0.991, median per-gene ratio 1.039).
  The shipped `annotated_matrix.mtx` itself is mis-keyed (bug 0a: 226 of
  45,921 rows aligned) and sums to 31.8M counts — the often-quoted
  41,944,353 is that figure plus its MatrixMarket size line. Post-filter
  `pasbed.bed` grows from 58,848 to 85,993 PAS because 27k tier-1 rows now
  carry enough counts to pass `min_cells`; 120 one-to-two-read clusters
  drop out because their counts moved to cells below `min_read`.

  **Measured on GSE104556 mouse1** (`--seq-len 98`, 12 threads, 54 min wall,
  26.9 GB peak RSS — unchanged from the previous run's 55 min / 26.9 GB;
  `results/benchmark_tools/gse104556/peakatail_clipseeded_v3/mouse1`):
  - every caller BED (`posbed.bed`, `negbed.bed`, `01_peak_calling/*`,
    `peakcalling/*`, `raw/pas.bed`) is identical in columns 1–4 and 6 to the
    previous run; only column 5 changes (`annotatedpas.bed` is byte-identical
    — it does not carry the score column);
  - clip reads 2,181,334 → 2,091,279 molecules (0.959×) over 148,715 tier-1
    rows; 13.3% (+) / 13.2% (−) of tier-1 rows change value; **3,269 of the
    62,159 "≥2 clip read" sites (5.3%) are a single molecule** (the same
    statistic is 23.7% on PBMC, where the duplication rate is far higher);
  - counting is otherwise unchanged: `reads_from_candidates`,
    `reads_from_window`, `reads_tier2`, `tier1`, `tier2` and
    `suppressed_candidates` are identical to the read; the raw caller
    matrices lose exactly the double-counted fallback reads
    (140,744,000 → 140,743,255, −745 counts, −0.0005%) and the annotated
    matrix 52,703,320 → 52,703,105 (−215, −0.0004%, 85,816 vs 85,993 PAS —
    the emptied fallback rows fall below `min_cells`);
    562 of the 5,250 fallback rows (98 on `+`, all 464 on `−`) were entirely
    reads a neighbouring cluster had already counted;
  - `min_read` keeps the same 10,339 cells, including 1,294/1,294 STARsolo
    cells;
  - scoring is unmoved: tier-1 P@100 0.4992, R_det 0.2991, F1 0.3741 —
    identical to 4 dp before and after (both tiers 0.4093/0.3231/0.3611).
  - **post-hoc sensitivity, labelled as such**: gating tier-1 at ≥2
    MOLECULES gives n 30,346, P@100 0.6923, R_det 0.2164, F1 0.3298. On the
    chr19+21 PBMC slice the same comparison, at equal nominal thresholds,
    separates the two units: ≥2 reads n 4,749 / P 0.5519 / R 0.2442 /
    F1 0.3386 versus ≥2 molecules n 3,664 / P 0.6395 / R 0.2245 / F1 0.3323
    (≥1 is unchanged at n 11,318 / 0.3591 / 0.3379 / 0.3482). A read
    threshold and a molecule threshold are not the same operating point.

### Notes
- The evidence is computed in all three peak-calling paths (monolithic,
  `--pipeline`, `--tiles`) from the `AlignedSegment` at the call site;
  `read_check`'s 5-tuple return is deliberately unchanged.
  `tests/test_polya_three_path_agreement.py` pins that the three paths
  produce identical BED output and identical per-PAS counts.
## Unreleased

### Added
- **`--cleavage-offset` — data-driven 3' cleavage-site offset correction
  (issue #72).** Called peak 3' ends stop ~90–105 nt short of the true
  cleavage site because 10x R2 coverage runs out before the poly(A) junction
  (AATAAA density peaks +75 nt downstream of the peak end; genomic A-fraction
  crests at +98 nt). Tight-cutoff benchmarks then punish the offset, not the
  calls. `--cleavage-offset N` (default `0` = legacy/off) shifts each reported
  PAS 3' end downstream by `N` bp (strand-aware, clamped at 0) after peak
  calling, so `pasbed.bed`, `annotatedpas.bed`, gene assignment, atlas
  matching, and the benchmark harness all use the inferred cleavage position.
  A sane fallback constant is `95`. The per-run data-driven estimator
  (`ema/countmatrix/cleavage_offset.py::estimate_cleavage_offset`) is stubbed
  with a clear `TODO(issue #72)`; it currently returns the constant.

## 0.2.0 (2026-05-10) — Product CLI

### Breaking changes
- The `ema_switch`, `ema_merge`, `ema_parse_gtf` console scripts are removed.
  Use `ema switch {diff,length,match}`, `ema merge`, `ema parse-gtf` instead.
- `ema --pdui-*` and `ema --diff-method` flags are dropped (moved to `ema switch length` / `ema switch diff` in 0.1.b).
- Output directories now have a timestamp suffix (e.g. `emaout_2026-05-10_181523/`).
- **BREAKING (bug fix): YAML `min_read` / `min_cells` / `min_genes` /
  `min_pas_per_cell` now actually take effect.** Previously the legacy
  argparse shim only bridged `min_pas_per_cell` from YAML — the other three
  filter thresholds silently kept their dataclass defaults (2000, 3, 50)
  regardless of YAML values. The atlas baseline in
  `reports/baseline_apa_completeness/` was regenerated to reflect this
  fix: per-dataset clusters now have shape (393, 5755) instead of
  (256, 5467) when `min_read: 1500` is set in
  `test_run/atlas_full.yaml`. Peak counts, BAM column counts, and
  atlas-snapped unified peaks remain byte-identical to the prior
  baseline — only downstream cluster sizes change because more cells
  survive the looser min_read threshold.
- **BREAKING (bug fix): `run_merge` now actually writes to `--output`.**
  Previously the function called `merge(..., threads=...)` without
  forwarding the `output` argument, so pysam wrote to the default
  `Aligned.sortedByCoord.merged.out.bam` and the user's `--output` flag
  was silently ignored.
- **BREAKING (bug fix): `--cluster-method` / `--resolution` / `--n-pcs` /
  `--random-seed` / `--external-clusters` now take effect.** Previously
  the pipeline called `clustering(adata=adata)` with no kwargs in both
  the single-sample and multi-sample paths, so all five flags were
  silently dropped. Default behavior is unchanged (defaults match the
  old hard-coded values), but explicitly setting any of these now
  affects the output.
- **BREAKING (bug fix): strategy hyperparameters `--lambda-method` /
  `--max-pas` / `--smoothing-window` / `--min-prominence` now take
  effect.** Previously these were declared in argparse but never
  forwarded into `get_strategy(name, **kwargs)`. Now passed through
  using `inspect.signature` filtering so each strategy receives only
  the kwargs it accepts.

### Loud no-op warnings
- `--ip-filter` / `--genome-fasta` / `--annot-filter` / `--ip-a-stretch` /
  `--benchmark` / `--validate-db` are accepted by `ema run` but emit a
  WARNING when set, telling the user the flag is currently a no-op
  (the underlying `apply_filters` and `run_benchmark` integrations are
  not wired into the new pipeline body). No silent fallbacks.

### Added
- Click-based CLI with subcommand structure
- Interactive wizard launched by bare `ema`
- Rich logging + per-run timestamped log file
- Nested progress bars, spawn-worker safe
- pyproject.toml + uv.lock for both pip and uv users
- `ema run --list-strategies` discovery
- New `--pas-gap` CLI flag (was YAML-only)
- New common options on every subcommand: `-v/-vv/-q`, `--log-level`,
  `--no-log-file`, `--no-progress`

### Changed
- Package PyPI name: `peakatail` (was `ema`); import path stays `ema`
- Build backend: setuptools via pyproject.toml (was setup.py)
