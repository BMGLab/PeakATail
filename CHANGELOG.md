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
  `--polya-min-reads`. `filter` and `require` are enforced at the existing
  `_apply_pas_filters()` seam and drop exactly the coverage-only tier.
- **Startup warning when the observed clip rate is below 0.3% of CB reads** —
  a pipeline that trims poly(A) before alignment destroys this evidence
  channel, and the caller now says so instead of silently emitting an
  unsupported call set.

### Fixed
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
    to multi-PAS peaks — through the coverage strategy's own
    `get_cb_dict_for_pas`), so the shipped mass is conserved exactly;
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
