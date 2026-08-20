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

### Notes
- The evidence is computed in all three peak-calling paths (monolithic,
  `--pipeline`, `--tiles`) from the `AlignedSegment` at the call site;
  `read_check`'s 5-tuple return is deliberately unchanged.
  `tests/test_polya_three_path_agreement.py` pins that the three paths
  produce identical BED output.
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
