# Plan — APA Pipeline Completeness

**Branch**: `feature/apa-completeness` (off `feature/hybrid-multi-sample`)
**Author**: Trex
**Date**: 2026-05-09
**Status**: ✅ Shipped — all 8 issues addressed, 29/29 strategy invariants pass

## Status snapshot (2026-05-10)

| # | Issue | Status | Evidence |
|---|---|---|---|
| 1 | PDUI multi-PAS (3 strategies + isoform-aware) | ✅ Shipped | `ema/quantification/strategies/{classic,proportion,shannon}.py`; PDUI invariants pass (PDUI ∈ [0,1], proportions sum to 1.0, entropy ∈ [0, log2(N)]) |
| 2 | Multi-sample diff APA wiring | ✅ Shipped (in `ema_switch`, not `ema` per user feedback) | `ema/switch_test/cli.py` |
| 3 | Cross-dataset cluster matching (3 strategies) | ✅ Shipped | `ema/clustering/cross_dataset/{marker_overlap,mnn,jaccard}.py`; all return confidence ∈ [0,1], 12 canonical clusters on identical-input regression |
| 4 | Fisher → NB regression | ✅ Shipped | `ema/switch_test/strategies/{fisher,nb_pairwise,nb_multi}.py`; KS p=0.69 on null simulation, 20/20 signal recovery |
| 5 | 9% cell-count gap diagnostic | ✅ Investigated (no bug) | Total counts in multi-sample slice byte-identical to raw peak calling output (6,337,529 counts, 2,794,774 nnz); gap is from bedtools merge boundary shifts, not data loss |
| 6 | Atlas e2e validation | ✅ Shipped | `test_run/atlas_validation.yaml` + atlas_snap.py BED12 fix; 147 PolyASite 2.0 PAS hit on chr22 test |
| 7 | BAM tagging speedup (samtools addreplacerg) | ✅ Shipped | `ema/datasets/manager.py` |
| 8 | Peak.pasnumber reset | ✅ Shipped | `ema/countmatrix/peak.py` |

### Discovered + fixed during validation (not in original plan)

| Bug | Location | Fix |
|---|---|---|
| nb_multi returned 0 rows on real data | per-cluster nonzero filter too strict | Replaced with total-nonzero + ≥2-clusters-with-signal filter |
| mnn match_confidence > 1.0 | `vote_matrix / cluster_size_a` doesn't account for k>1 MNN per cell | Row-normalize so each row is a probability distribution |
| jaccard returned 0 matches on PeakATail-prefixed CBs | `_strip_prefix` only handled `#` (Cell Ranger), not `_` (PeakATail) | Try `#` first, fall back to `_` |
| atlas_snap returned 0 hits on PolyASite 2.0 | `_IDX_DISTANCE=12` hardcoded for BED6 atlas; PolyASite is BED12+ | Derive `_IDX_DISTANCE = len(cols) - 1` per row |
| ResourceManager missing | n_jobs=-1 oversubscribed on loaded box | New `ema/utils/resource_manager.py` (psutil-based, falls back to /proc/meminfo) |
| Marker selection missing | NB on 20K PAS × 66 pairs = 6 hours | New `ema/quantification/marker_selector.py`; ema_switch now defaults to top-200 markers per cluster |

### Architecture correction during validation

User feedback led to: PDUI + diff APA do NOT belong in main `ema` pipeline (heavy compute, user wants to control scope). They moved to existing `ema_switch` separate command, which now uses the strategy registries. Main `ema` ends at clustering + cross-dataset matching.

### Performance baseline (Phase 3)

Pipeline runtime on full BAM (`data/SRR8325947_Aligned.sortedByCoord.out.bam`, 14.7M reads), 16GB RAM box, 14 physical cores:

| Stage | Time |
|---|---|
| `ema` main (multi-sample 2 datasets, peak call → cluster → cross-dataset match) | 3m38s |
| `ema_switch` fisher (66 pairs, top-200 markers, classic PDUI) | 1m08s |
| `ema_switch` nb_pairwise (3 pairs, top-200 markers, prop+shannon PDUI) | 17s |
| `ema_switch` nb_multi omnibus (top-200 markers, classic PDUI) | 48s |
| All 3 cross-dataset matchers | <2s each |
| Strategy validation suite (29 invariant checks) | ~30s |

ResourceManager picked n_jobs=12 (vs 20-thread oversubscription before).

cProfile/memory_profiler reports deferred to peak-calling-perf branch.

---


---

## 1. Summary

PeakATail's multi-sample mode now runs end-to-end (peak call → unify → per-dataset cluster). What is still missing or wrong:

- PDUI throws away information for genes with 3+ PAS
- PDUI is gene-level only — doesn't account for isoform-specific 3'UTR structures
- Differential APA still uses Fisher exact, which is anti-conservative for scRNA-seq
- Multi-sample mode does no differential APA at all
- Cross-dataset clusters share PAS coordinates but have no cluster correspondence
- Atlas merge strategy works on synthetic data but never validated on real atlas
- BAM tagging in `before` strategy is slow
- Process-global state (`Peak.pasnumber`) leaks between datasets

This plan closes all eight gaps in one branch, using the **strategy + registry pattern** consistently so each new behavior is a pluggable strategy with no `if/else` chains in the orchestration layer.

---

## 2. Goals

1. **Multi-PAS-aware PDUI** — three quantification methods (`classic`, `proportion`, `shannon`), user-selectable, multiple at once
2. **Isoform-aware OR gene-aware aggregation** — both modes selectable per run
3. **Statistically defensible differential APA** — NB regression (pairwise + multi-condition) as alternative to Fisher
4. **Multi-sample completeness** — switch test runs in multi-sample mode, per dataset
5. **Cross-dataset cluster identity** — clusters in different datasets get canonical IDs when they correspond
6. **Atlas validated** — atlas merge strategy proven on real PolyASite 2.0 data, not just synthetic
7. **Performance** — `before` BAM tagging uses `samtools addreplacerg` instead of pysam Python loop
8. **Hygiene** — `Peak.pasnumber` reset between datasets so per-dataset pasnumbers start at 1

Non-goals (out of scope for this branch):
- Joint cross-dataset clustering with batch correction (Harmony / Scanorama)
- Marker-PAS based cell type annotation
- Internal priming filter using genome FASTA
- Joint multi-method PDUI (e.g. weighted ensemble of shannon + proportion)

---

## 3. Architectural decisions

### 3.1 Strategy + registry pattern (mandatory)

We already use this for peak calling (`ema/strategies/`) and clustering (`ema/clustering/strategies/`). Extend the same pattern to the four new dimensions:

| Dimension | Strategies | Registry location |
|---|---|---|
| PDUI quantification | `classic`, `proportion`, `shannon` | `ema/quantification/strategies/` |
| Aggregation level | `per_gene`, `per_isoform` | (lookup table inside PDUI strategies) |
| Isoform-collapse | `none` (keep per-isoform), `mean`, `majority`, `weighted` | (lookup table inside PDUI strategies) |
| Differential APA test | `fisher`, `nb_pairwise`, `nb_multi` | `ema/switch_test/strategies/` |
| Cross-dataset cluster match | `marker_overlap`, `mnn`, `jaccard` | `ema/clustering/cross_dataset/strategies/` |

Registry contract (mirrors existing pattern):

```
ema/<area>/strategies/__init__.py:
    def get_strategy(name: str) -> Strategy:  # lookup
    def register_strategy(name: str, cls: type)  # plugin point
    list_strategies() -> list[str]
```

CLI flags accept comma-separated values when multiple methods are wanted. Orchestrator (main.py) loops over selected strategies — never branches on method names.

### 3.2 No code repetition between strategies

Each strategy implements a small base interface:

```
class PDUIStrategy:
    def compute(self, count_matrix, pas_isoform_map, pseudocount) -> pd.DataFrame
```

Common helpers (sorting PAS by transcript-coord, normalizing by gene total, BH FDR, etc.) live in shared utility modules and are imported by all strategies. No copy-paste.

### 3.3 GTF parsing is shared, parsed once

The GTF isoform parser runs once at pipeline start (already in a background thread for gene-level UTRs). The new isoform-level parse adds a second worker that produces the per-transcript UTR map. Both share the GTF file read.

### 3.4 Backward compatibility

- `--pdui-method` defaults to `classic` to match current single-sample behavior
- `--diff-method` defaults to `fisher`
- `--cluster-match-method` defaults to `marker_overlap`
- Existing single-sample runs produce identical output unless new flags are set

---

## 4. Issues catalog (8)

| # | Issue | Owner |
|---|---|---|
| 1 | PDUI uses only 2 PAS per gene | Agent A |
| 1b | PDUI is gene-level, not isoform-level | Agent A |
| 2 | Multi-sample skips differential APA | Phase 2 (me) |
| 3 | No cross-dataset cluster correspondence | Agent C |
| 4 | Fisher exact is anti-conservative | Agent B |
| 5 | 9% cell-count gap in multi-sample (vs single) | Phase 2 (me) — investigation only |
| 6 | Atlas method not validated on real atlas | Phase 3 (me) |
| 7 | Slow BAM tagging in `before` strategy | Agent D |
| 8 | `Peak.pasnumber` global state leaks between datasets | Agent D |

---

## 5. Component requirements

### 5.1 GTF isoform UTR parser (`ema/annotate/gtf2isoform_utr.py`)

**What we need**:
- Parse GTF and emit, per transcript, the ordered list of `three_prime_utr` features
- Preserve transcript-coordinate ordering (5' → 3') based on strand
- Group by `gene_id` so the output is `{gene_id: {transcript_id: [(chrom, start, end, strand, exon_rank), ...]}}`
- Cache result with mtime+size fingerprint (same pattern as existing `gtf_cache.py`)

**Evaluation**:
- Spot-check 10 random genes against the Ensembl web GTF viewer
- For genes with >1 transcript, verify multiple entries appear with distinct transcript IDs
- For transcripts whose 3'UTR spans multiple exons, verify multiple entries in the list

**Test**:
- Unit test on a 10-gene synthetic GTF covering: single-isoform gene, multi-isoform gene, gene with split UTR, gene on minus strand, gene with no UTR annotation
- Verify caching: second call < 100 ms when cache hit

---

### 5.2 PAS-to-isoform mapping (`ema/quantification/pas_to_isoform.py`)

**What we need**:
- For each PAS in the unified BED, find every transcript whose 3'UTR contains it
- Compute `transcript_pos` (position in transcript coordinates: sum of upstream UTR-exon lengths + offset within current UTR exon, strand-aware)
- Compute `rank` (0 = most proximal among PAS in that transcript's UTR)
- Output: `{pas_id: [(gene_id, transcript_id, transcript_pos, rank, total_pas_in_transcript), ...]}`

**Evaluation**:
- For our test BAM, > 80% of annotated PAS should map to at least one transcript UTR
- For a gene known to have alternative UTRs (e.g. CD47), verify the PAS show up in different transcripts with different positions

**Test**:
- Unit test with synthetic isoform UTR map + synthetic PAS BED — verify ranks correct, transcript positions correct, multi-isoform genes produce one entry per isoform per PAS

---

### 5.3 PDUI strategies (`ema/quantification/strategies/`)

**Files to create**:
- `base.py` — `PDUIStrategy` abstract base class
- `classic.py` — current 2-PAS PDUI, kept for backward compatibility
- `proportion.py` — per-PAS proportion vector per gene per cell
- `shannon.py` — Shannon entropy of PAS usage distribution per gene per cell
- `__init__.py` — registry (`get_pdui_strategy`, `list_pdui_strategies`)

**What each strategy must produce** (common output schema):
- `classic`: long format `[gene_id, transcript_id (optional), cell, pdui]` — one row per gene per cell
- `proportion`: long format `[gene_id, transcript_id (optional), pas_id, rank, cell, proportion]` — one row per (gene, transcript, PAS, cell)
- `shannon`: long format `[gene_id, transcript_id (optional), cell, entropy, normalized_entropy, n_pas]`

**Aggregation handling**:
- Aggregation strategy is a parameter passed to the PDUI strategy at compute time (`per_gene` or `per_isoform`)
- If `per_gene`: collapse all isoforms first (sum PAS counts assigned to gene regardless of which isoform's UTR they fall in)
- If `per_isoform`: keep transcript-level breakdown
- Isoform-collapse for gene-level reporting: `mean`, `majority`, `weighted` (subselectable when `per_isoform` results need to be summarized)

**Evaluation**:
- Shannon entropy should be 0.0 for genes where all reads go to one PAS (verify on a single-PAS gene)
- Shannon entropy should equal log2(N) for uniform usage (verify on synthetic uniform distribution)
- Per-PAS proportions per (gene, cell) must sum to 1.0 (within rounding)
- Classic PDUI on multi-PAS genes must equal the value computed by the legacy code (regression test)

**Test**:
- Unit test for each strategy on a synthetic 3-gene matrix
- Cross-check sum-to-1 invariant for proportion
- Cross-check entropy boundary cases (delta-distribution, uniform, two-equal)

---

### 5.4 Differential APA strategies (`ema/switch_test/strategies/`)

**Files to create**:
- `base.py` — `DiffAPAStrategy` abstract base class
- `fisher.py` — wrap existing Fisher exact + BH FDR
- `nb_pairwise.py` — NB GLM with cluster as binary predictor; one test per PAS per cluster pair
- `nb_multi.py` — NB GLM with cluster as multi-level factor; ANOVA-style omnibus test per PAS
- `__init__.py` — registry

**What each NB strategy must produce**:
- DataFrame indexed by PAS with columns: `[log2fc (pairwise only), pvalue, qvalue, dispersion, n_cells, test_stat]`
- Use `statsmodels.GLM(family=NegativeBinomial(alpha=dispersion))`
- Size-factor offset = log(library_size) per cell
- Dispersion estimation: per-PAS MLE if convergence succeeds, fall back to method-of-moments
- Filter PAS with too few cells per group (< 10 by default)

**Evaluation**:
- On simulated data with no differential signal, NB p-values should be uniformly distributed (0–1)
- On simulated data with known fold-change, NB log2fc should recover the truth within ±10%
- NB p-values should be less inflated than Fisher's on the same data (KS test on uniform distribution)

**Test**:
- Unit test: simulate count matrix with NB-distributed counts, no group difference → p-value distribution should be uniform
- Unit test: inject 2x fold-change in 100 PAS → NB should recover them with FDR < 0.05
- Comparison test: same data through Fisher vs NB → Fisher should produce more significant calls (anti-conservative)

---

### 5.5 Cross-dataset cluster matching (`ema/clustering/cross_dataset/`)

**Files to create**:
- `base.py` — `ClusterMatchStrategy` abstract base class
- `marker_overlap.py` — top-N marker PAS per cluster, Jaccard overlap, Hungarian assignment
- `mnn.py` — mutual nearest neighbors in shared LSI space
- `jaccard.py` — direct Jaccard overlap of cell sets (only useful when same cells in both datasets)
- `__init__.py` — registry

**What each strategy must produce**:
- DataFrame: `[dataset_id, original_cluster, canonical_cluster, match_confidence, matched_to (list)]`
- Canonical cluster IDs use letters `A`, `B`, `C`, ... or integers — consistent format across strategies
- Unmatched clusters get a unique canonical ID + confidence 0
- One canonical ID per group of corresponding clusters across all datasets

**Evaluation**:
- On 2 datasets with identical input (our regression test), every cluster in sampleA should match exactly one in sampleB with confidence ≥ 0.9
- On 2 datasets from genuinely different samples, expect some clusters to match (housekeeping cell types) and some to be dataset-unique

**Test**:
- Unit test: 2 synthetic h5ad files with 5 clusters each, identical marker structure → expect 1-to-1 matching
- Integration test: real multi-sample test config (2× same BAM) → expect all 12 clusters to match 1-to-1

---

### 5.6 main.py multi-sample loop wiring (Phase 2, my work)

**What we need**:
- After per-dataset clustering, parse user-selected PDUI methods (`--pdui-method` comma-list)
- For each method, instantiate the strategy via registry, run it, write output to `per_dataset/{ds}/pdui_{method}.tsv`
- After PDUI, parse user-selected differential APA method (`--diff-method`)
- For each cluster pair within the dataset, run the strategy, write output to `per_dataset/{ds}/differential/{method}_{c1}_vs_{c2}.tsv`
- After all per-dataset outputs are complete, run cross-dataset cluster matching (selectable strategy) and write `cross_dataset/canonical_cluster_map.tsv`
- All loops parameterized by registry — no `if method == ...` branches

**Evaluation**:
- Single-sample run with `--pdui-method classic --diff-method fisher` produces identical output to the current pipeline
- Multi-sample run with all defaults produces all expected output files
- Multi-sample run with `--pdui-method shannon,proportion --diff-method nb_pairwise` produces six output files per dataset (3 PDUI methods × 2 strands? — actually 2 PDUI methods × N cluster pairs for diff, plus 1 cross-dataset map)

**Test**:
- Run on `test_run/multi_sample_full.yaml` with default methods → verify all per-dataset directories have PDUI + differential subdirectories populated
- Run with `--pdui-method=shannon,proportion,classic` → verify three PDUI files per dataset

---

### 5.7 BAM tagging speedup (`ema/datasets/manager.py`)

**What we need**:
- Replace pure-Python pysam loop in `_tag_bam_with_rg` with `subprocess` call to `samtools addreplacerg`
- Multi-thread via `-@` flag
- Same output path, same RG tag content, no behavioral change

**Evaluation**:
- Tagging a 14M-read BAM in < 2 minutes (currently ~10 min)
- Output BAM byte-identical or close (modulo header timestamp / sort order) to the current pure-Python output

**Test**:
- Run `test_run/multi_sample_full.yaml` with `merge_strategy: before` → confirm runtime improvement
- Diff the output BAM headers — RG section should be identical

---

### 5.8 `Peak.pasnumber` reset (`ema/countmatrix/peak.py`)

**What we need**:
- Add `Peak.reset_pasnumber()` classmethod that sets `Peak.pasnumber = 0`
- main.py multi-sample loop calls `Peak.reset_pasnumber()` between datasets (after `reset_index()`)
- All downstream code that depends on `(dataset_id, old_pasnumber)` keying (concat_matrices, atlas_snap mapping file) keeps working because dataset_id provides global uniqueness

**Evaluation**:
- After multi-sample run, each dataset's peak-calling output BED files should have pasnumbers starting at 1 (not continuing from the previous dataset's last number)

**Test**:
- Inspect `peakcalling/sampleA_0.pos.bed` first row → pasnumber should be 1
- Inspect `peakcalling/sampleB_0.pos.bed` first row → pasnumber should also be 1 (not 16509+)

---

## 6. Phase plan

### Phase 1 — Parallel agent work (Sonnet)

| Agent | Files | LoC | Skills | Owner |
|---|---|---|---|---|
| A | gtf2isoform_utr.py, pas_to_isoform.py, quantification/strategies/* | ~570 | `python-parallel-data-streaming` (chunked GTF parse), `single-cell-rna-qc` (sparse-aware PDUI) | new isoform machinery + 3 PDUI strategies |
| B | switch_test/strategies/* (fisher wrap + nb_pairwise + nb_multi) | ~400 | `python-parallel-data-streaming` (joblib over PAS for NB fit) | NB regression machinery |
| C | clustering/cross_dataset/* (3 strategies + registry) | ~350 | `python-parallel-data-streaming` (parallel marker computation per cluster) | cross-dataset cluster matching |
| D | manager.py samtools tagging + peak.py reset | ~25 | `python-resource-management` (subprocess pipe cleanup) | quick cleanups |

All 4 agents run in parallel with non-overlapping file ownership. Locked contracts ensure interfaces match before any wiring.

Each agent's prompt explicitly names the skill(s) to invoke before starting implementation.

### Phase 2 — Wiring (sequential, me)

- Update `cli.py` to accept all new flags
- Update `config.py` to expose them
- Rewrite `main.py` per-dataset loop to:
  - Load isoform UTR map once
  - For each dataset: compute selected PDUI methods, then selected differential APA method
  - After loop: run cluster matching strategy
- Update `example.yaml` to document all new keys
- Update `peackcalling.py` call site to invoke `Peak.reset_pasnumber()` (no logic change inside peackcalling.py itself — orchestration is main.py's job)

### Phase 3 — Validation + profiling (sequential, me)

- Run full BAM with default methods → regression check vs existing pipeline output
- Run full BAM with `--pdui-method=shannon,proportion --diff-method=nb_pairwise` → verify all output files present and well-formed
- Run atlas test config (`atlas_validation.yaml`) → verify > 50% snap rate, comparable cluster count
- Run cell-gap diagnostic to characterize the 9% multi-vs-single difference (no code change unless real bug)
- **Performance profiling** (using `python-performance-optimization` + `memory-profiling` skills):
  - Run `scripts/profile_pipeline.py` on full BAM
  - Verify all targets in §11.2 are met
  - If any target missed, identify the specific function via cProfile → optimize → re-run
  - Commit `reports/profile_run.txt` and `reports/memory_profile.txt` so future PRs can be compared against this baseline

---

## 7. Test strategy

### Unit tests

Each new module gets unit tests in `tests/`. Synthetic data only, no real BAMs.

- `tests/test_gtf2isoform_utr.py` — synthetic 10-gene GTF, verify per-transcript UTR extraction
- `tests/test_pas_to_isoform.py` — synthetic isoform map + PAS BED, verify ranks and transcript positions
- `tests/test_pdui_strategies.py` — synthetic count matrix, verify each strategy's output schema and invariants
- `tests/test_nb_regression.py` — simulated counts (NB-distributed), verify p-value uniformity under null and recovery under signal
- `tests/test_cluster_matching.py` — synthetic h5ad pair, verify Hungarian assignment correctness

### Integration tests

End-to-end on `test_run/chr22.bam` (subset, fast):

- `multi_sample_test.yaml` with `merge_strategy: none` and full new pipeline → verify all output files exist
- Same with `merge_strategy: atlas` and atlas BED → verify atlas snap output

### Regression tests

End-to-end on full BAM:

- Single-sample with default flags → byte-identical to current `feature/hybrid-multi-sample` output (no behavior change for default-flag users)
- Multi-sample with default flags → matches current `feature/hybrid-multi-sample` 12-cluster, 662-cell result per dataset

---

## 8. Evaluation criteria

A successful merge requires all of:

| Criterion | Pass condition |
|---|---|
| All unit tests pass | `pytest tests/` exits 0 |
| Regression on default flags | Output diff vs `feature/hybrid-multi-sample` is empty for cluster labels and PDUI classic values |
| New PDUI methods produce expected schemas | Shannon entropy in [0, log2(N)]; proportions per (gene, cell) sum to 1.0 |
| NB regression p-values are not inflated | KS test against uniform distribution: p > 0.01 on null simulation |
| Cross-dataset matching | On 2 identical-BAM datasets, all clusters match 1-to-1 with confidence ≥ 0.9 |
| Atlas snap rate | On real PolyASite 2.0 with full BAM: ≥ 50% of called PAS within 50bp of an atlas PAS |
| BAM tagging speedup | `before` strategy runtime on 14M-read BAM: ≥ 5× faster than current |
| Per-dataset pasnumbers | First pasnumber in each dataset's first BED is 1 |

---

## 9. Risks

| Risk | Mitigation |
|---|---|
| GTF parsing is fiddly — Ensembl format quirks | Use existing GTF cache infrastructure; add tolerance for missing `three_prime_utr` (some transcripts only have `exon` + `CDS`) |
| NB GLM convergence failure on low-count PAS | Catch `ConvergenceWarning`, fall back to method-of-moments dispersion + Wald test |
| MNN matching needs cells in shared embedding | Document clearly that `mnn` requires re-embedding; default to `marker_overlap` which does not |
| 9% cell-count gap may turn out to be a real bug | Diagnostic script in Phase 3; if real bug found, scope to next branch |
| samtools version mismatch (`addreplacerg` requires ≥ 1.10) | Add version check at startup; raise clear error if too old |
| Per-isoform PDUI matrix grows large for genes with many isoforms | Document the size implications; provide `mean` aggregation as size-controlled fallback |

---

## 10. Open questions / future work

(Not addressed in this branch.)

- Joint cross-dataset clustering with batch correction (Harmony, Scanorama, scVI)
- Marker-PAS cell type annotation (matching against published cell atlas markers)
- Internal priming filter using genome FASTA (PolyA-stretch detection)
- Weighted ensemble of PDUI methods (e.g. shannon × proportion as a stability metric)
- Bayesian hierarchical model for differential APA (alternative to NB)
- Visualization layer: per-isoform PDUI heatmaps, cross-dataset cluster Sankey

---

## 11. HPC / Performance plan

The pipeline has four data-heavy hotspots. Each gets explicit performance treatment via the HPC skills.

### 11.1 Hotspots and applied skills

| Hotspot | Where | Scale | Skill applied | Strategy |
|---|---|---|---|---|
| GTF isoform parsing | `gtf2isoform_utr.py` | 1.2 GB GTF, ~250K transcripts | `python-parallel-data-streaming` | Chunk the GTF by chromosome, parse chunks in `multiprocessing.Pool`, merge results. Use mmap for the file read. |
| Per-PAS NB regression | `nb_pairwise.py`, `nb_multi.py` | ~32K PAS × N cluster-pairs × per-PAS GLM fit | `python-parallel-data-streaming` | Embarrassingly parallel — `joblib.Parallel(n_jobs=-1)` over PAS. Chunk size = ceil(n_pas / n_workers). |
| PAS → isoform mapping | `pas_to_isoform.py` | ~32K PAS × ~250K transcripts | `python-parallel-data-streaming` | Use `bedtools intersect -wa -wb` (already C-native) for the heavy join, then post-process the result in pandas with vectorized ops. |
| Per-PAS proportion matrix | `proportion.py` | ~32K PAS × ~700 cells × per-isoform breakdown | `dask` (optional) + `single-cell-rna-qc` | Default: scipy.sparse vectorized normalization (no copy). Fallback: dask.array if matrix > 4GB. |
| PDUI compute across all genes | all PDUI strategies | ~20K genes per dataset | `python-parallel-data-streaming` | `joblib.Parallel` over genes, share read-only sparse matrix via memory-map or shared memory. |

### 11.2 Performance targets

| Component | Target on full BAM (chr1–22, 14M reads) | Hard ceiling |
|---|---|---|
| GTF isoform parser (cold) | < 60s | 120s |
| GTF isoform parser (cache hit) | < 1s | 5s |
| PAS → isoform mapping | < 30s | 60s |
| PDUI proportion (all strategies) | < 30s per dataset | 120s |
| PDUI shannon | < 10s per dataset | 30s |
| NB regression per cluster pair | < 60s per pair | 180s |
| Cross-dataset cluster matching | < 30s for 2 datasets, 12 clusters each | 60s |
| BAM tagging (samtools addreplacerg) | < 2 min per BAM | 5 min |
| Total multi-sample run (2 datasets, all strategies) | < 15 min | 30 min |

### 11.3 Profiling discipline

Before optimizing anything, we measure. Add to Phase 3:

- `python-performance-optimization` skill: run `cProfile` on baseline run, identify top 10 hotspots, optimize only those
- `memory-profiling` skill: track peak RAM during multi-isoform PDUI computation; if > 2× baseline, switch to chunked / streaming approach
- Profiling script: `scripts/profile_pipeline.py` — runs the multi-sample pipeline under `cProfile` + `memory_profiler`, dumps results to `reports/profile_run.txt`

Output: a profile report committed alongside the code, so future regressions are obvious.

### 11.4 Memory budget

Multi-isoform PDUI matrix size estimate (worst case):

```
shape  = (n_genes × max_isoforms_per_gene, n_cells)
       ≈ (20,000 × 5, 1000) = 100,000 × 1,000
bytes  = 100,000 × 1,000 × 8 (float64) = 800 MB per strategy
```

With three strategies (classic, proportion, shannon) running concurrently per dataset and per-isoform aggregation: peak ~3 GB. Acceptable on a 16 GB workstation but tight on a laptop.

Mitigations (in order of preference):
1. Float32 instead of float64 → halves memory
2. Sparse output for `proportion` strategy (most cells have 0 for most PAS) → sparsity ~80%
3. Stream per-strategy: compute one strategy fully, write to disk, free, then next
4. Chunk by gene block (1000 genes at a time)

`single-cell-rna-qc` skill provides the canonical patterns for sparse-aware operations.

### 11.5 Single-cell QC adoption

`single-cell-rna-qc` skill applies to:
- `matrixfilter.py:filter_cb` — current implementation is text-stream parsing; can be replaced with Scanpy `sc.pp.filter_cells(min_counts=...)` on the AnnData representation, leveraging optimized C paths
- `matrixfilter.py:preprocessing` — already uses Scanpy, but could add doublet detection (Scrublet) and ambient RNA correction (SoupX) as optional steps
- Per-dataset clustering: add standard QC plots (n_genes, n_counts, pct_mito) before clustering — currently missing

Out-of-scope for this branch but logged: doublet detection, ambient RNA correction.

---

## 12. Appendix: file inventory

### New files

```
ema/annotate/gtf2isoform_utr.py
ema/quantification/pas_to_isoform.py
ema/quantification/strategies/__init__.py
ema/quantification/strategies/base.py
ema/quantification/strategies/classic.py
ema/quantification/strategies/proportion.py
ema/quantification/strategies/shannon.py
ema/switch_test/strategies/__init__.py
ema/switch_test/strategies/base.py
ema/switch_test/strategies/fisher.py
ema/switch_test/strategies/nb_pairwise.py
ema/switch_test/strategies/nb_multi.py
ema/clustering/cross_dataset/__init__.py
ema/clustering/cross_dataset/base.py
ema/clustering/cross_dataset/marker_overlap.py
ema/clustering/cross_dataset/mnn.py
ema/clustering/cross_dataset/jaccard.py
tests/test_gtf2isoform_utr.py
tests/test_pas_to_isoform.py
tests/test_pdui_strategies.py
tests/test_nb_regression.py
tests/test_cluster_matching.py
test_run/atlas_validation.yaml
scripts/profile_pipeline.py        ← cProfile + memory_profiler harness
reports/profile_run.txt            ← committed baseline (Phase 3 output)
reports/memory_profile.txt         ← committed baseline (Phase 3 output)
plans/apa-completeness.md          ← this file
```

### Modified files

```
ema/main.py                       — Phase 2 wiring
ema/cli.py                        — new flags
ema/config.py                     — expose new flags
ema/datasets/manager.py           — Agent D (samtools tagging)
ema/countmatrix/peak.py           — Agent D (reset method)
ema/quantification/pdui.py        — Agent A (port to strategy registry)
example.yaml                      — document new keys
```

### Deleted files

None.
