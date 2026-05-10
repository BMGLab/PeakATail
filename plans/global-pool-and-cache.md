# Plan — Global Job Pool + RAM-adaptive Tiles + GTF Cache

**Branch**: `feature/global-pool` (off `feature/peak-calling-perf`)
**Author**: Trex
**Date**: 2026-05-10
**Status**: Pending approval

---

## 1. Summary

Three orthogonal extensions to the peak-calling performance work:

1. **Global job pool**: collapse the current nested loop (`for dataset in datasets: for tile in tiles:`) into a flat `JobSpec` list run through a single `multiprocessing.Pool`. This parallelises across datasets AND tiles in one go.

2. **RAM-adaptive tile sizing**: `ResourceManager` computes the tile size dynamically from free RAM, worker count, and BAM read density. Users with small machines get small tiles (more workers, less RAM/worker); users with big machines get large tiles (less merge overhead). No manual tuning needed.

3. **Reusable GTF cache + standalone parser**: GTF parsing already runs in a background thread inside `ema`, and we already have per-output-dir caching. Add (a) a standalone `ema_parse_gtf` CLI command, (b) a global cache fingerprinted by SHA256+mtime+size of the GTF, (c) a fast-path skip in main when cache hits.

A fourth, **cross-cutting** concern: all three extensions involve more parallelism, so the plan also formalises **global-state encapsulation** to make races structurally impossible — not just "safe by luck of spawn semantics".

---

## 2. Goals

1. **Parallelise across datasets AND tiles in one pool** — eliminate the sequential `for dataset` loop in `main.py`
2. **Eliminate samtools-view materialization** in workers — use `pysam.AlignmentFile.fetch(chrom, start, end)` directly (one disk read pass instead of three)
3. **Auto-size tiles by free RAM** — no OOM under load, no manual tuning required
4. **GTF cache survives across runs** — first parse populates `~/.cache/peakatail/gtf/`, subsequent runs skip
5. **Standalone `ema_parse_gtf` command** — pre-warm cache before any analysis
6. **Encapsulate global state** so workers can't accidentally share BarcodeIndex / Peak.pasnumber / default_sample_id between jobs
7. **Backwards compatibility** — existing `--tiles` and monolithic paths still work; new behavior opt-in via `--global-pool` flag (or default-on after validation)

Non-goals (out of scope for this branch):
- Distributed across multiple machines (Nextflow/Snakemake territory)
- Parallelising downstream stages (clustering, atlas snap)
- GPU offload
- Replacing the streaming bisect/SortedList core algorithm

---

## 3. Architecture

### 3.1 Single global job pool

```
For each dataset d in YAML:
    For each chromosome c in d.bam:
        For each tile t in tiles_for(c, ram_adaptive=True):
            For each direction in (pos, neg):
                jobs.append(JobSpec(d, c, t.start, t.end, t.fetch_start, t.fetch_end, direction, ...))

with multiprocessing.Pool(spawn, n_workers) as pool:
    results = pool.imap_unordered(worker_run_one, jobs, chunksize=1)
    for r in tqdm(results, total=len(jobs)):
        ...

# After all jobs done — group by (dataset, direction) and merge
for dataset_id in dataset_ids:
    merge_tile_outputs(dataset_id, jobs_for_dataset, output_paths)
```

`imap_unordered(chunksize=1)` so fast workers steal jobs from slow ones — best load balancing for variable-cost tiles.

### 3.2 Direct pysam fetch in workers

Today each tile worker calls `samtools view -bh chrom:start-end` to materialize a region BAM (extra I/O write + extra read pass). New worker calls `pysam.AlignmentFile(bam, threads=1).fetch(chrom, fetch_start, fetch_end)` directly. Saves 2/3 of the BAM I/O.

This requires refactoring `peak_calling()` to accept a `region=(chrom, start, end)` parameter (or a `read_iter` callable). Backward compatible: `region=None` means full BAM (current behavior).

### 3.3 RAM-adaptive tile sizing

```python
class ResourceManager:
    def get_tile_size(
        self,
        bam_path: str,
        n_workers: int,
        target_per_worker_mb: int = 300,
    ) -> int:
        """Compute safe tile size for current memory budget + BAM density."""
        free_ram_mb = self.free_ram_mb()
        budget_per_worker = (free_ram_mb * 0.7) / n_workers
        budget_per_worker = min(budget_per_worker, target_per_worker_mb)

        # Cheap density estimate: sample first 1Mb of BAM
        density = self._estimate_bam_density(bam_path)  # reads per Mb

        # Each read in cb_dict ~ 50 bytes; Peak overhead ~10MB; pysam buffer ~32MB
        base_overhead_mb = 60
        bytes_per_read = 50
        budget_for_reads = (budget_per_worker - base_overhead_mb) * 1024 * 1024
        max_reads_per_tile = budget_for_reads / bytes_per_read

        tile_size_bp = int(max_reads_per_tile / density)
        # Floor to keep merge overhead reasonable
        return max(5_000_000, min(tile_size_bp, 100_000_000))
```

The 1Mb density sample is fast (<100ms) and a reasonable proxy for the whole BAM.

### 3.4 Standalone `ema_parse_gtf` + cache

New CLI command in `ema/annotate/cli.py`:

```bash
ema_parse_gtf --gtf GTF [--cache-dir DIR]
```

Produces in cache dir:
- `gene_end.bed` — gene endpoint BED for find_close
- `utr_lengths.tsv` — gene_id → 3'UTR length
- `isoform_utrs.pkl` — per-transcript UTR structure (pickled dict for fast load)
- `cache_meta.json` — `{gtf_sha256_first4kb, gtf_mtime, gtf_size, parsed_at, parser_version}`

In `ema` main: before spawning the GTF parse thread, probe cache dir for matching meta. If hit, load directly from cache files (fast). If miss, run the parse thread as today.

Cache key: SHA256 of first 4KB + mtime + size. Cheap to compute (<10ms for 1.2GB GTF) and robust enough — accidental collision would require two distinct 1.2GB GTFs to match in the first 4KB AND have identical mtime AND size, which is essentially zero probability for real-world files.

### 3.5 Encapsulated global state (CRITICAL — see §4)

Workers access `BarcodeIndex` as a passed instance, not a module singleton. Same for `default_sample_id`. `Peak.pasnumber` becomes a per-instance counter. Module-level singletons remain for backward compatibility but the new code path is fully instance-based.

---

## 4. Global state management

This is the most important section because it's where parallelism most often goes wrong.

### 4.1 Inventory of process-global state today

| Variable | Location | Mutation | Reset mechanism |
|---|---|---|---|
| `Peak.pasnumber` | `peak.py` class attr | `+= 1` on each peak emit | `Peak.reset_pasnumber()` |
| `_index` (BarcodeIndex) | `indexing.py` module singleton | `[cb] = next_id` per CB | `reset_index()` |
| `_default_sample_id` | `read.py` module var | `set_default_sample_id(id)` | overwrite |
| `start_time` | `peackcalling.py` module | written at import | none (read-only after) |
| `_REGISTRY` (×3) | strategy `__init__.py` | populated at import | none (read-only after) |
| `filtered_cb_list` | `matrixfilter.py` module | overwritten in `filter_cb()` | reset in function |

### 4.2 Risks under different parallelism models

| Model | `Peak.pasnumber` | `_index` | `_default_sample_id` | Notes |
|---|---|---|---|---|
| **multiprocessing.Pool(spawn)** (today, tile workers) | Each worker has own copy ✓ | Each worker has own copy ✓ | Each worker has own copy ✓ | Spawn copies on fork; mutations isolated |
| **multiprocessing.Pool(fork)** | Inherits parent state — risk of stale counter | Inherits parent dict — risk of stale entries | Inherits parent value — risk of wrong sample | NEVER use fork for our workers |
| **threading.Thread** (today, GTF parse) | Shared mutable state — RACE | Shared mutable state — RACE | Shared mutable state — RACE | Currently the GTF thread doesn't touch any of these — safe by accident |
| **asyncio** | Single event loop, single thread — no race | No race | No race | Not applicable to our CPU-bound work |

The **structural risk**: today everything is "safe by spawn" — if anyone ever switches to `fork` or introduces threading on the wrong code path, races appear silently.

### 4.3 Encapsulation plan

**Phase 1 (in this branch)**: encapsulate `BarcodeIndex` as the highest-mutation-rate state. Keep others as documented-isolated module globals.

```python
# Today:
from ema.countmatrix.indexing import indexing  # uses module singleton
col = indexing(cb)

# After:
from ema.countmatrix.indexing import BarcodeIndex
index = BarcodeIndex()      # fresh instance per worker
col = index.get_index(cb)
```

Affected files (all touched):
- `peak.py` — `Peak.cb_counting(cb, index)` accepts index parameter
- `paswrite.py` — `matrix_write(cb_dict, pasnumber, output, index)` accepts index
- `peackcalling.py` — creates a `BarcodeIndex` at top, passes to `peak.cb_counting` and `matrix_write`
- `matrixfilter.py` — `filter_cb(input_matrix_paths, cb_list, ...)` already accepts cb_list explicitly; remove the `get_mapping()` fallback path
- `main.py` — orchestration unchanged (each peak_calling call still gets isolated state)

Backward compatibility shim: keep module-level `indexing()` / `reset_index()` / `get_mapping()` functions wrapping a global singleton, so any caller we miss still works (just not optimally).

**Phase 2 (this branch)**: encapsulate `Peak.pasnumber` as instance-level.

```python
# Today: class attribute shared by all Peak instances
Peak.pasnumber = 0
peak = Peak(...)
Peak.pasnumber += 1  # global counter

# After: instance attribute on each peak_calling run
class Peak:
    def __init__(self, ..., pasnumber_offset=0):
        self.pasnumber_counter = pasnumber_offset
        self.next_pasnumber = pasnumber_offset
    
    def emit_pas(self):
        self.next_pasnumber += 1
        return self.next_pasnumber

# Caller:
peak_state = PeakCallingState(initial_pasnumber=0)  # per-worker
peak_state.emit_pas()  # increments instance counter
```

Workers each have their own `PeakCallingState`. Main process renumbers globally during merge.

**Phase 3 (this branch)**: `default_sample_id` becomes a parameter on `read_check`, not a module variable.

```python
def read_check(read, direction, sample_id, ...):  # explicit param
    try:
        rg = read.get_tag('RG')
    except (KeyError, ValueError):
        rg = sample_id  # fallback uses arg, not module global
```

Hot path consideration: this is called per-read (14M times). Adding one positional argument has zero perf cost.

### 4.4 Locking?

**Not needed for spawn-based multiprocessing** (each worker is a separate OS process, no shared memory). The encapsulation above eliminates the need for locks because there's literally no shared state to race over.

If we ever need genuine shared state (e.g., shared progress counter across workers): use `multiprocessing.Value(typecode, lock=True)` — atomic, lock-free for primitive types.

If we ever need threading (e.g., async I/O for BAM reads): every encapsulated instance becomes per-thread via `threading.local()`.

### 4.5 Defensive assertions at worker boundaries

```python
def worker_run_one(job_spec):
    # Sanity: this worker should start with empty state
    assert Peak.pasnumber == 0, "Peak.pasnumber leaked between workers"
    # ... rest of work ...
```

Catches accidental state leaks during refactoring. Removed in production builds.

---

## 5. Component requirements

### 5.1 `JobSpec` dataclass (new file: `ema/jobs/job_spec.py`)

```python
@dataclass(frozen=True)
class JobSpec:
    job_id: str           # unique e.g. "sampleA_chr1_tile0_pos"
    dataset_id: str
    bam_path: str
    chrom: str
    tile_start: int       # core range (peaks emitted only if start ∈ [tile_start, tile_end))
    tile_end: int
    fetch_start: int      # extended range for context (10kb overlap)
    fetch_end: int
    direction: bool       # False = pos strand, True = neg strand
    strategy_name: str
    peak_kwargs: dict     # default_threshold, merge_len, etc.
    bam_threads: int = 1
    output_bed: Path
    output_mtx: Path
    output_cb: Path
```

Frozen + pickle-safe. All paths absolute.

### 5.2 `build_job_list(datasets, ram_adaptive_size=True)`

Walks YAML datasets, builds JobSpec list. Sorts by estimated cost (longest tiles first) for better load balancing on `imap_unordered`. Estimated cost = tile_length × dataset_density.

### 5.3 `worker_run_one(job_spec)` (new file: `ema/jobs/worker.py`)

```python
def worker_run_one(job: JobSpec) -> dict:
    # Open BAM with region restriction — no temp file
    bam = pysam.AlignmentFile(job.bam_path, 'rb', threads=job.bam_threads)
    read_iter = bam.fetch(job.chrom, job.fetch_start, job.fetch_end)

    # Per-worker state — fully encapsulated
    index = BarcodeIndex()
    state = PeakCallingState()

    # Run the streaming algorithm with the region iterator
    peaks_emitted = run_streaming_peak_calling(
        read_iter,
        direction=job.direction,
        sample_id=job.dataset_id,
        index=index,
        state=state,
        **job.peak_kwargs,
    )

    # Filter: emit only peaks whose start is in [tile_start, tile_end)
    write_peaks(
        [p for p in peaks_emitted if job.tile_start <= p.start < job.tile_end],
        job.output_bed, job.output_mtx,
    )
    write_cb_list(index.mapping, job.output_cb)

    return {'job_id': job.job_id, 'n_peaks': len(peaks_emitted), ...}
```

### 5.4 `global_pool_runner(jobs, n_workers)` (in `ema/jobs/runner.py`)

`multiprocessing.Pool` with `spawn` context, `imap_unordered(chunksize=1)`, tqdm progress bar.

### 5.5 Refactored `peak_calling(region=...)`

Accept optional `region` tuple. When given, use `bam.fetch(*region)` instead of full iteration. Backward compatible default.

### 5.6 `ResourceManager.get_tile_size(bam_path, n_workers, target_per_worker_mb=300)`

Adaptive sizing per §3.3.

### 5.7 `ema_parse_gtf` standalone CLI (new entry point in `pyproject.toml`)

```bash
ema_parse_gtf --gtf GTF [--cache-dir ~/.cache/peakatail/gtf] [--force]
```

Runs both `process_gtf_cached` (gene-end + UTR lengths) and `parse_isoform_utrs` (isoform UTR map). Writes all outputs + meta to cache dir.

### 5.8 GTF cache fingerprint check (in `main.py`)

```python
def _gtf_cache_hit(gtf_path, cache_dir):
    meta_path = Path(cache_dir) / "cache_meta.json"
    if not meta_path.exists(): return False
    meta = json.loads(meta_path.read_text())
    return (
        meta['gtf_size'] == os.path.getsize(gtf_path)
        and meta['gtf_mtime'] == os.path.getmtime(gtf_path)
        and meta['gtf_sha256_first4kb'] == _sha256_first_4kb(gtf_path)
    )
```

Hot path: <10ms even for 1.2GB GTF.

---

## 6. Phase plan

### Phase 1 — Encapsulate global state (sequential, blocking everything else)

This goes FIRST because all subsequent parallelism work depends on safe state isolation.

| Order | Component | Files | Owner |
|---|---|---|---|
| 1.1 | `BarcodeIndex` becomes parameter-passed instance | `peak.py`, `paswrite.py`, `peackcalling.py`, `matrixfilter.py` | Sonnet agent A |
| 1.2 | `Peak.pasnumber` → instance counter `PeakCallingState` | `peak.py`, `peackcalling.py` | Same agent |
| 1.3 | `default_sample_id` → explicit `read_check(...)` parameter | `read.py`, `peackcalling.py` | Same agent |
| 1.4 | Backward-compat shims for module-level functions | `indexing.py`, `read.py` | Same agent |
| 1.5 | Defensive assertions at worker entry points | `tile_runner.py`, future `worker.py` | Same agent |

After Phase 1: full BAM regression — `diff_peak_outputs.py` must still PASS (byte-identical to baseline).

### Phase 2 — Refactor `peak_calling(region=...)` + remove samtools view

| Order | Component | Files | Owner |
|---|---|---|---|
| 2.1 | Add `region` param to `peak_calling` | `peackcalling.py` | Sonnet agent B |
| 2.2 | Use `bam.fetch(*region)` instead of `for read in bamfile` when region given | `peackcalling.py` | Same |
| 2.3 | Update tile_runner to pass `region` instead of materializing temp BAM | `tile_runner.py` | Same |

After Phase 2: tile-mode regression test — must still PASS, expect ~30% speedup on tiles.

### Phase 3 — Global job pool

| Order | Component | Files | Owner |
|---|---|---|---|
| 3.1 | `JobSpec` dataclass | `ema/jobs/job_spec.py` (new) | Sonnet agent C |
| 3.2 | `build_job_list(datasets)` | `ema/jobs/builder.py` (new) | Same |
| 3.3 | `worker_run_one(job)` | `ema/jobs/worker.py` (new) | Same |
| 3.4 | `global_pool_runner(jobs)` | `ema/jobs/runner.py` (new) | Same |
| 3.5 | `merge_jobs_per_dataset` (extends current tile merge) | `ema/jobs/merger.py` (new) | Same |
| 3.6 | Integrate in `main.py` behind `--global-pool` flag | `main.py` | me |
| 3.7 | tqdm progress bar | `ema/jobs/runner.py` | Same agent |

After Phase 3: full BAM with multi-dataset config + `--global-pool` — must still PASS regression, expect 3-5× speedup on multi-dataset.

### Phase 4 — RAM-adaptive tile sizing

| Order | Component | Files | Owner |
|---|---|---|---|
| 4.1 | `ResourceManager.get_tile_size(...)` | `ema/utils/resource_manager.py` | Sonnet agent D |
| 4.2 | `_estimate_bam_density(bam_path)` helper | Same file | Same |
| 4.3 | Use in `build_job_list` | `ema/jobs/builder.py` | Same |
| 4.4 | CLI override: `--tile-size auto|N` | `cli.py` | Same |

After Phase 4: stress test on small (8GB RAM) and large (64GB) machines — must produce same output, just different worker/tile counts.

### Phase 5 — GTF cache + standalone parser

| Order | Component | Files | Owner |
|---|---|---|---|
| 5.1 | Cache fingerprint helper (SHA256 + mtime + size) | `ema/annotate/gtf_cache.py` | Sonnet agent E |
| 5.2 | `cache_meta.json` write/read | Same | Same |
| 5.3 | `ema_parse_gtf` CLI entry point | `ema/annotate/cli.py` (new) | Same |
| 5.4 | Add to `pyproject.toml` scripts | `pyproject.toml` | Same |
| 5.5 | Cache-hit fast-path in main.py | `main.py` | me |
| 5.6 | `--gtf-cache-dir` CLI arg | `cli.py` | Same agent |

After Phase 5: run twice — second run should skip GTF parse (60s saved).

### Phase 6 — Validation + bench

- Full BAM, baseline mode → diff vs `reports/baseline_apa_completeness/` PASS
- Full BAM, `--global-pool` → diff PASS, measure wall-clock
- Multi-dataset config (2-4 datasets) → diff PASS per dataset, measure wall-clock vs current sequential
- GTF cache: timed comparison cold vs warm

---

## 7. HPC skills usage

| Skill | Where applied | Pattern |
|---|---|---|
| `python-parallel-data-streaming` | Phase 3 global pool | `multiprocessing.Pool(spawn)` + `imap_unordered(chunksize=1)` for work-stealing load balance |
| `python-resource-management` | Phase 4 RAM-adaptive | Free RAM probe via psutil, density estimate via small BAM read, budget-based sizing |
| `python-performance-optimization` | Phase 6 bench | cProfile to verify peak calling is no longer the bottleneck |
| `python-background-jobs` | Phase 3 worker pool | tqdm progress, retry-on-failure for individual workers |
| `python-resource-management` | Phase 1 state encapsulation | Instance-per-worker pattern — eliminates need for locks |
| `python-anti-patterns` | Phase 1 review | Check for hidden module-level mutable state we missed |

### 7.1 Patterns explicitly applied

**Work-stealing pool**: `imap_unordered(chunksize=1)` — workers grab next available job immediately. Better than `chunksize=N` which can leave one worker with the last big chunk.

**Instance-per-worker state**: each worker creates its own `BarcodeIndex`, `PeakCallingState` at startup. No shared state, no locks, no races.

**Defensive assertions**: assert clean state at worker entry — catches state-leak bugs during refactoring instead of silently producing wrong output.

**Cheap-fast caching**: SHA256 of first 4KB (not full file) + mtime + size — robust enough, fast enough.

**Backward-compat shims**: module-level `indexing()` / `reset_index()` / `get_mapping()` keep working but wrap a global singleton — anything we miss in refactor still runs.

### 7.2 What we do NOT do

- **No `multiprocessing.Manager` shared dicts** — they're slow and we don't need shared state
- **No `threading.Lock`** — spawn-based isolation eliminates the need
- **No `dask` / `ray`** — overkill for single-machine scale
- **No `multiprocessing.fork`** — explicitly use `spawn` to avoid inherited state issues
- **No global mutable singletons** in new code — only kept as backward-compat shims around encapsulated instances

---

## 8. Test strategy

### Unit tests
- `tests/test_barcode_index_encapsulation.py` — verify two `BarcodeIndex` instances don't share state
- `tests/test_peak_state_encapsulation.py` — verify two `PeakCallingState` instances have independent counters
- `tests/test_job_spec_serialization.py` — verify JobSpec pickles cleanly (mp.Queue compatible)
- `tests/test_resource_manager_tile_size.py` — verify adaptive sizing under different RAM/CPU mocks
- `tests/test_gtf_cache_fingerprint.py` — verify fingerprint changes when GTF mtime/size/content changes

### Integration tests
- `tests/test_global_pool_smoke.py` — small synthetic BAM through global pool, verify per-dataset outputs match monolithic
- Full BAM regression: `scripts/diff_peak_outputs.py reports/baseline_apa_completeness/ emaout/` must PASS for every mode

### End-to-end on full BAM
1. Baseline mode (no flags) → must match baseline byte-for-byte
2. `--tiles` mode → must match
3. `--pipeline` mode → must match
4. **NEW**: `--global-pool` mode → must match
5. Multi-dataset YAML (2 datasets) with `--global-pool` → must match per-dataset baselines AND show speedup

### GTF cache test
1. Delete cache dir
2. Run `ema_parse_gtf --gtf big.gtf` → measure time T1 (~60s)
3. Run again → should skip with message, T2 < 5s
4. `touch big.gtf` (mtime change) → next run re-parses
5. Modify first 4KB → fingerprint catches it, re-parses

---

## 9. Evaluation criteria

A successful merge requires all of:

| Criterion | Pass condition |
|---|---|
| All unit tests pass | `pytest tests/` exits 0 |
| Full BAM regression | `diff_peak_outputs.py` ALL CHECKS PASS for baseline / tiles / pipeline / global-pool / multi-dataset |
| State encapsulation | `pytest tests/test_*_encapsulation.py` PASS |
| Multi-dataset speedup | 4 datasets on 16+ core box: `--global-pool` ≥ 2.5× faster than sequential |
| RAM safety | On 8GB-RAM machine: `--global-pool` does not exceed 6 GB peak RSS, completes successfully |
| GTF cache | Second run on same GTF: <5s GTF stage (was ~60s) |
| Backwards compatibility | All existing single-sample, multi-sample-merge, multi-sample-atlas configs still produce identical output without new flags |
| No fork-based MP | `git grep -E "Pool\(\)" ema/` returns 0 hits — only `get_context('spawn')` allowed |

---

## 10. Risks

| Risk | Mitigation |
|---|---|
| `BarcodeIndex` encapsulation breaks a caller we missed | Keep module-level shim functions; integration tests catch it |
| `Peak.pasnumber` instance-counter changes pasnumbering semantics | Pasnumbers are renumbered during merge anyway — order is what matters |
| `pysam.fetch` skips reads we currently iterate | Add fetch-mode regression test on small BAM with known peak set |
| RAM-adaptive sizing produces too-many tiny tiles | Floor of 5Mb tile; total tiles capped at 10× n_workers |
| BAM density sample (first 1Mb) misleading on uneven coverage | Acceptable — workers re-balance via work-stealing pool |
| GTF cache collision via SHA256 + mtime + size | Negligible probability for real-world files; `--force` flag to override |
| `ema_parse_gtf` running concurrently with `ema --gtf-cache-dir` (race on cache write) | File lock on cache dir; second writer waits or skips |
| Per-worker pysam handles leak file descriptors | Use `with pysam.AlignmentFile(...) as bam:` context manager |
| tqdm progress bar interferes with logging | Use `tqdm.write()` for log output |
| `multiprocessing.Pool` zombies on Ctrl+C | Use `pool.terminate()` in signal handler |

---

## 11. Open questions / future work

(Not addressed in this branch.)

- Distributed across multiple machines (Nextflow / Snakemake)
- Parallel downstream stages (per-dataset clustering in parallel)
- Cython/Rust rewrite of the streaming bisect inner loop
- Smart cache pre-warming based on usage patterns
- Cross-pipeline cache sharing (e.g., share GTF cache across Snakemake workflows)

---

## 12. Appendix: file inventory

### New files

```
ema/jobs/__init__.py
ema/jobs/job_spec.py             — JobSpec dataclass
ema/jobs/builder.py              — build_job_list(...)
ema/jobs/worker.py               — worker_run_one(...)
ema/jobs/runner.py               — global_pool_runner(...)
ema/jobs/merger.py               — merge_jobs_per_dataset(...)
ema/annotate/cli.py              — ema_parse_gtf entry point
ema/countmatrix/peak_state.py    — PeakCallingState class (replaces Peak.pasnumber class attr)
tests/test_barcode_index_encapsulation.py
tests/test_peak_state_encapsulation.py
tests/test_job_spec_serialization.py
tests/test_resource_manager_tile_size.py
tests/test_gtf_cache_fingerprint.py
tests/test_global_pool_smoke.py
plans/global-pool-and-cache.md   — this file
```

### Modified files

```
ema/countmatrix/peak.py          — Phase 1: instance-based pasnumber + cb_dict; index parameter on cb_counting
ema/countmatrix/indexing.py      — Phase 1: BarcodeIndex still exists; module-level functions become shims
ema/countmatrix/paswrite.py      — Phase 1: matrix_write accepts index parameter
ema/countmatrix/read.py          — Phase 1: read_check accepts sample_id parameter
ema/countmatrix/peackcalling.py  — Phase 1: instance state; Phase 2: region parameter
ema/countmatrix/tile_runner.py   — Phase 2: pass region instead of materializing
ema/matrixfilter.py              — Phase 1: remove get_mapping fallback path
ema/utils/resource_manager.py    — Phase 4: get_tile_size + density estimate
ema/annotate/gtf_cache.py        — Phase 5: fingerprint helper + cache_meta.json
ema/main.py                      — Phase 3 + 5: global-pool dispatch + GTF cache fast-path
ema/cli.py                       — Phases 3, 4, 5: --global-pool, --tile-size auto, --gtf-cache-dir
pyproject.toml                   — Phase 5: add ema_parse_gtf entry point
```

### Deleted files

None.
