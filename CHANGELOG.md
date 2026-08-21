# Changelog

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
- **Internal-priming filter now tests the correct side on the `-` strand.**
  `ema/experimental/internal_priming.py` applied the `+` genomic window
  `[pos-10, pos+30)` to both strands, so on `-` it scanned 30 bp *upstream* /
  10 bp downstream of the cleavage site in transcript orientation — mostly
  the wrong side for an oligo-dT priming artifact. The window is now defined
  in transcript orientation (`--ip-window-left` upstream, `--ip-window-right`
  downstream, both strands) and mirrored in genomic coordinates on `-`
  (`[pos-30, pos+10)`, reverse-complemented before the A-run / A-fraction
  scan). `+`-strand calls are byte-identical; flags, defaults and the
  `--ip-*` CLI are unchanged. New pure helpers `ip_window()`,
  `call_internal_priming()` and `check_internal_priming()` are unit-tested
  in `tests/test_internal_priming_strand.py` (implanted runs on a C/G-only
  genome, both strands, fraction rule, contig edges, end-to-end through
  `filter_internal_priming` and `apply_filters`).

  Measured post hoc on the Stage-2 final run (4efeb12) by re-applying both
  rules to the emitted `pasbed.bed` without re-running the caller: `+`
  strand, 0 changes on every dataset. On `-`, the corrected rule flags
  5,599 / 160,464 (3.5 %) of the PBMC sites that had passed the old filter
  (testis mouse1 1,505 / 35,487 = 4.2 %, mouse2 1,302 / 35,946 = 3.6 %),
  and of the 39,244 `-` sites the old rule removed from the PBMC candidate
  set, 11,217 would be kept by the corrected rule. On the pre-registered
  precision default (PBMC, tier-1, ≥2 molecules; PolyASite 2.0 within
  100 bp) removing the 386 newly-flagged sites moves P@100 0.7167 → 0.7184;
  the full corrected filter (newly flagged removed AND wrongly removed sites
  restored, 46,524 vs 44,394 sites) gives P@100 0.7062, R@100 0.1754 (was
  0.1707). Every pre-registered gate still passes; benchmark conclusions
  are unchanged. Any use of the per-site `internal_priming` flag from a
  pre-fix run should re-run the filter (`ema reannotate --genome-fasta`).
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
