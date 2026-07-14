# Engine branch handoff — `engine/bugfixes-e-series`

Branch off `d8db1e9`. **33 atomic commits**, all with unit tests. Worktree:
`trees/engine-fixes`. Venv: `/home/user/PeakATail/.venv`.

**Test status:** `802 passed, 4 skipped, 5 failed` (pytest, full suite).
The 5 failures are **pre-existing at base `d8db1e9`** (verified against a temp
base worktree) and unrelated to this branch:

| Failing test | Why (pre-existing) |
|---|---|
| `test_cli_switch::test_isoform_agg_per_gene_dispatches_to_per_gene_branch` | BLOCKER-1 strategy-dispatch, base bug |
| `test_downstream_parallel::test_raises_on_empty_sub_indices` | test passes a `per_dataset_dir` kwarg the fn never accepted |
| `test_downstream_parallel::test_full_pipeline_calls_with_mocks` | same kwarg drift |
| `test_pyproject_install::test_entry_point_is_installed` | needs `pip install -e .` |
| `test_pyproject_install::test_entry_point_runs` | needs `pip install -e .` |

---

## What landed (all with synthetic unit tests, contract-verified where applicable)

### Bug fixes
- **B0** serialize resolved RunConfig (not argparse defaults) → `run_config.json`.
- **B1** strand-safe merge key `(dataset_id, strand, pasnumber)`; mapping TSV gains a `strand` column.
- **B2** `annotate()` resolves output paths at call time (kills `emaout/` shadow files).
- **B3** always write `filtered_cb.tsv` (dropped the `.exists()` gate).
- **B4 (data-controller)** `switch diff/length` resolve `pasbed` via `--pasbed` + run layout.
- **B5** per-dataset stage funnel (`stage_stats.json`) in multi-sample runs.
- **B6** round-trip `canonical_cluster` into `clusters.h5ad` obs.
- **B7** repair NaN `var['gene_id']` from the pas_gene map on `var_names` at concat.
- **D1** atlas snap-rate counts both strands + drops the `max(0,·)` clamp.
- **D2** atlas distance measured to the PAS 3' summit, not the peak interval.
- **D3** flag atlas snap_rate as coverage (circularity caveat), not precision.

### Write-side / contract
- **E1** content-addressed `pas_uid = chrom:pos:strand` sidecar (`unified/pas_uid.tsv`).
- **E2** `OutputManager → run_manifest.json`, conforming to frozen `peakatail_contract.RunManifest` (validated in the contract venv).
- **E3** append-only pas/cell drop ledgers (`ema/provenance.py`, `record_drop`), **wired at atlas-snap AND the per-dataset cb_filter/pas_gene/preprocess drop sites**, with a reusable invariant checker `check_survivor_invariant()` + `surviving_count_from_tsv()`, self-checked per worker. **Sidecar-then-reconcile:** each multiprocessing worker writes its own `provenance/` sidecar (no shared state); the parent (`main.py`, after the pool joins) calls `reconcile_dataset_ledgers()` — per-dataset `surviving PAS == n_vars(clusters.h5ad)` self-check + cohort concatenation into `provenance/by_dataset/` + `provenance/reconcile_summary.json`. Two-level model: the run-level atlas ledger and per-dataset ledgers answer different questions and are kept separate (never force-merged). Conforms to `PasLedgerRow`/`CellLedgerRow` column order.
- **E5** `findings_long` (FindingRow) + `length_long` (LengthRow) producers; parquet with `.tsv` fallback (`pyarrow>=15` added to deps, fallback kept).
- **E2 per-artifact content hash** — `run_manifest.json` artifact entries now carry `content_hash` (`"sha256:<hex>"`, streamed) so a consumer can detect a stale index when a referenced artifact changes but the manifest bytes don't. Rides as an extra key (contract ignores it on validate); **consuming it needs a contract `Artifact.content_hash` field bump — hub-team**.

### Statistics
- **D4** fisher `count_mode` — **opt-in `"cells"`** (de-pseudoreplicated per-cell contingency among gene-expressing cells). **DEFAULT is `"reads"` (legacy, unchanged)**; `test_nb_regression` restored verbatim to base. See flag #3.
- **D5** nb_multi `sample_split` — **opt-in, default `False`** (select PAS on half A, infer on disjoint half B). Default unchanged.
- **D8** shorten/lengthen polarity lives in `length_long` (LengthRow), computed by `structural_length_direction()` — deterministic one-vs-rest Δdistal-usage sign per (gene, cluster). `classic` uses PDUI (distal fraction); `proportion` uses the max-rank (distal) PAS; `shannon` → `undetermined`. Never NA-silent; always a contract `Direction` enum. Carries `direction_basis="structural"`. `findings_long` intentionally stays differential-usage (flat/undetermined is fine).

### Net-new subcommands
- **A4** `ema collapse` — samtools-merge RG-suffix pooling (pure `canonical_cb` + count-preserving sparse collapse).
- **A3** `ema switch trend` — ordered-stage APA length trend (slope + Spearman + direction, overall & per-gene), generalized beyond hardcoded stage sets.
- **A2** `ema switch combine` — stitch stage-labelled `clusters.h5ad` into per-group (optional split-by-celltype) h5ads for `--cluster-key stage`.
- **A1 (core)** `ema/celltype/scoring.py` — marker M1(cluster-argmax)/M2(cell-majority) assignment + agreement + ARI/AMI concordance.

---

## FLAGGED — needs real-run data or is net-new (documented, NOT blocking)

1. **A1 GEX-I/O + `ema celltype` CLI.** The STARsolo `Solo.out/Gene` pooling,
   scanpy clustering, and `scanpy.tl.score_genes` that PRODUCE the score matrix,
   plus the `ema celltype` command wiring them. The prototype
   (`scripts/gex_celltyping.py`) is explicitly "prove the biology first, then
   port"; the port needs real Laughney STARsolo matrices to validate before it
   lands. The deterministic math half (assignment + concordance) is landed +
   tested. **Next step:** point the porter at a real `--star-root` + master YAML;
   wrap `ema/celltype/scoring.py` behind a loader + CLI; validate ARI/AMI on the
   real cohort.

2. **E3 live cross-stage invariant on a REAL cohort run.** The
   sidecar-then-reconcile machinery is fully wired and tested on a synthetic
   2-dataset mini-run: per-dataset sidecars, `reconcile_dataset_ledgers()`, the
   `surviving PAS == n_vars(clusters.h5ad)` self-check, cohort concatenation,
   and `reconcile_summary.json`. What remains is **running it on a real
   multi-sample cohort** to confirm the invariant holds on live data (and to
   surface any drop site the accounting still misses). D10 `stratum_to_label`
   population is tied to A1 celltyping. **Next step:** execute a real
   multi-sample run; inspect `provenance/reconcile_summary.json` (`all_ok`) and
   any per-dataset warning in the run log.

3. **D4/D5 number-effect + default flip.** The per-cell fisher and split-select
   nb_multi are textbook-correct and landed as opt-in with synthetic tests, but
   their effect on **reported** significance must be characterized on the
   **no-atlas re-run** before flipping the defaults. **Next step:** after the
   re-run, compare `count_mode="cells"` vs `"reads"` (and `sample_split=True`
   vs `False`) FDR calls; if calibrated as expected, flip defaults + update
   `test_nb_regression`.

4. **D9** (`differs_across_stages` / `omnibus_q`) lives in the harvest scripts,
   not `ema/` — nothing to wire or drop in the engine.

---

## Coordination notes for hub-team

- `run_manifest.json`, `findings_long`, `length_long` all validate against the
  frozen contract in `/home/user/D/peakatail-hub/packages/contract/.venv`.
- **`length_long` direction (geneview overlays):** `direction` is a per-(gene,
  canonical_cluster) one-vs-rest structural call across the clusters present in
  the length output; `direction_basis="structural"` is an extra informational
  column the contract ignores on validate. Geneview length overlays can read
  `direction` directly for shorten/lengthen coloring.
- Ledger `dropped_at` stage ids in use: `atlas_snap`, `cb_filter`, `pas_gene`,
  `preprocess`. Survivor rows have `dropped_at == ""`.
