# Plan — Peak Calling Performance Optimization

**Branch**: `feature/peak-calling-perf` (off `feature/apa-completeness` after that branch lands)
**Author**: Trex
**Date**: 2026-05-10
**Status**: Pending approval

---

## 1. Summary

Peak calling on a 14.7M-read BAM currently takes ~2 minutes single-threaded. Profiling-style measurements show two main cost centers:

- **Per-read pysam attribute access**: ~25–30% of total time
- **`bisect.insort` into a Python list**: ~25% of total time (O(n) per insert)
- **CB string hashing in dict lookups**: ~10–15% of total time

Disk I/O is **not** the bottleneck — pysam already buffers at the BGZF block level. We're leaving free wins on the table by not using `pysam.AlignmentFile(..., threads=N)` for parallel block decompression.

This plan delivers a 5–10× speedup on peak calling phase via five layered optimizations, applied in lowest-risk-first order. The streaming bisect endpoint algorithm remains intact — we improve the data structures around it, not the algorithm itself.

---

## 2. Goals

1. **Turn on parallel BGZF decompression** — `threads=4` on every `AlignmentFile` open (one-line change, ~2× free)
2. **Replace `bisect.insort` with `sortedcontainers.SortedList`** — O(log n) inserts vs O(n)
3. **2-bit encode CB barcodes to int64 once at read time** — int-keyed dicts replace string-keyed dicts everywhere
4. **Pipeline architecture (reader → finder → writer)** — overlap I/O with CPU
5. **Tile-based parallelism within a chromosome** — N-core scaling for the finder stage
6. **Explicit per-worker isolation + final merge** for `Peak.pasnumber` and `BarcodeIndex` so parallelism doesn't corrupt the count matrix

Non-goals (out of scope):
- Replacing pysam with custom Cython/htslib bindings
- Replacing the streaming bisect endpoint algorithm itself
- New peak-calling strategy (existing strategy registry stays unchanged)
- GPU acceleration

---

## 3. Architectural decisions

### 3.1 The pysam reality

pysam has **no public batched-read API**. Internally:
- BGZF blocks (~64KB compressed) are read from disk one at a time
- The `threads=N` parameter on `AlignmentFile` parallelizes block decompression — currently unused
- Python iteration over reads happens at one read per `next()` call but the I/O underneath is already block-batched

Implication: "batched reading" in our code does not mean a new pysam API call. It means **buffering parsed reads in our reader stage** before crossing process boundaries (IPC) or invoking downstream stages.

### 3.2 Strategy + registry stays untouched

The existing `ema/strategies/` registry (`original`, `lambda_poisson`, `lambda_gradient`, `sierra_iterative`) is unaffected. Each strategy operates on a `Peak` object and returns PAS positions. The optimizations here improve how reads flow into the `Peak` object and how the count matrix is keyed — not how the strategy decides where peaks are.

### 3.3 Per-worker isolation + main-process merge

The two pieces of process-global state (`Peak.pasnumber` class attribute and `BarcodeIndex` module singleton) cannot be shared across workers safely without locking, which destroys parallelism. Instead, **each worker has its own complete copy**, and the main process performs a deterministic merge after all workers complete:
- `pasnumber` is reassigned sequentially in genomic order during BED+MTX merging
- `BarcodeIndex` is rebuilt by union-ing per-worker CB sets and assigning canonical column indexes; per-worker MTX files have their column numbers remapped during stitching

This is the same pattern we already use in `concat_matrices` for multi-sample mode — applied at a finer (intra-chromosome) granularity.

### 3.4 No new dependencies that aren't actively maintained

Considered and **rejected**:
- `cykhash` (Cython int-int hash table): partial wheel coverage on Python 3.12+, last release 2023. Plain Python dict on int64 keys gets us ~80% of the win with zero deps.
- Custom Cython BAM parser: too much code to maintain.

Adopted:
- `sortedcontainers` — pure Python, supports 2.7+, rock solid, zero compatibility risk
- `joblib` — already installed transitively via scanpy
- `multiprocessing`, `mmap`, `queue` — stdlib

---

## 4. Bottleneck catalog and measurements

Measured on `data/SRR8325947_Aligned.sortedByCoord.out.bam` (14.7M reads, chr22 subset for benchmarks):

| Operation | Throughput today | Cost share |
|---|---|---|
| pysam iterate + 1 attribute | 2.58M reads/sec | baseline |
| pysam iterate + full read_check fields (5 attrs + 2 tags) | 1.30M reads/sec | 2× slowdown vs baseline = per-read attribute cost |
| samtools view subprocess + line iter | 1.47M reads/sec | comparable to pysam — no real win from leaving pysam |
| `bisect.insort` into 5000-element Python list | ~50K ops/sec | dominant in high-coverage regions |
| `dict[str, int]` lookup on 16-char CB | ~5M ops/sec | adds ~15% to per-peak emission |

For a full 14.7M-read run at 1.3M reads/sec → 11.3 seconds just for read iteration. The remaining ~110 seconds of total peak-calling time goes to bisect, CB indexing, peak state machine, and write.

---

## 5. Component requirements

### 5.1 Issue P-1 — Enable parallel BGZF decompression

**What we need**:
- Every `pysam.AlignmentFile(...)` open in the codebase passes `threads=N` where N is configurable (default 4, override via CLI `--bam-threads`)
- Audit current code (`peackcalling.py`, `manager.py`, anywhere else that opens BAM)
- Document tradeoff: more threads → more memory (each block buffer ~64KB × threads)

**Evaluation**:
- Run `time ema --config single_sample_full.yaml` before and after
- Expected: ~2× speedup on peak calling phase (from ~2min to ~1min)

**Test**:
- Smoke test: opens BAM with `threads=4`, reads first 1000 reads, no errors
- Memory check: peak RSS does not exceed baseline + 2 MB (since we're only adding 64KB per thread)

### 5.2 Issue P-2 — Replace `bisect.insort` with `SortedList`

**What we need**:
- Add `sortedcontainers>=2.4` to `requirments.txt`
- Replace `data_array = []` and `bisect.insort(data_array, end1)` with `from sortedcontainers import SortedList; data_array = SortedList()` and `data_array.add(end1)`
- Replace `bisect.bisect_left(data_array, start1)` with `data_array.bisect_left(start1)`
- Replace `data_array = data_array[slice_loc:]` with `del data_array[:slice_loc]` (in-place, O(slice_loc) instead of full copy)

**Evaluation**:
- Run before/after on chr22 BAM
- Expected: 3–5× speedup specifically on data_array operations
- Validate: same peak count, same coordinates as today (regression: `diff` BED outputs)

**Test**:
- Unit test: insert 10000 random integers, slice from middle, verify behavior matches old bisect+list
- Regression test: full pipeline on chr22, BED output byte-identical to baseline

### 5.3 Issue P-3 — 2-bit CB encoding to int64

**What we need**:
- New module: `ema/countmatrix/cb_encode.py`
- Function `encode_cb(cb_str: str) -> int` — single CB to int64 via 2-bit encoding (A=00, C=01, G=10, T=11)
- Function `encode_cb_batch(cb_list: list[str]) -> np.ndarray` — vectorized batch encoding via `np.frombuffer` + bitshift over byte buffers
- Function `decode_cb(cb_int: int, length: int = 16) -> str` — for output (writing CB strings back to filtered_cb.tsv)
- Handle non-ATCG characters (Ns) by returning a sentinel (e.g. `-1`) and filtering upstream
- Update `read.py:read_check` to encode CB to int64 before returning
- Update `BarcodeIndex` in `indexing.py` to use `dict[int, int]` instead of `dict[str, int]`
- Update `Peak.cb_dict` and `Peak.cb_positions` to use int64 keys (no API change, just type)
- Update `matrix_write` to call `decode_cb` only when writing the final CB list TSV (one decode per cell, not per read)

**Evaluation**:
- Encoding throughput: > 50M reads/sec for batch encode (target via numpy bit ops)
- Memory: CB int64 dict uses ~10× less memory than string dict for 200K cells
- Speedup on peak emission: ~3× (no string hashing in hot path)
- Verification: every CB string round-trips through encode→decode unchanged (assuming valid ATCG)

**Test**:
- Unit test: encode all combinations of 4-mer ATCG, verify decode returns original
- Unit test: batch encode 10000 random CBs, verify each matches single-encode
- Unit test: handle N character — returns sentinel, doesn't crash
- Regression test: full pipeline produces identical filterdcb.tsv as baseline (decoded order matches encoded order)

### 5.4 Issue P-4 — Batch CB resolution at peak emit

**What we need**:
- Defer global `BarcodeIndex` updates from per-read to per-peak-emit
- During peak accumulation: only touch the local `peak.cb_dict[cb_int64] += 1` (still per read but local-only, fast)
- At peak emit (in `matrix_write` or strategy hand-off): bulk-resolve all unique CBs in the peak via a single batch call to `BarcodeIndex.get_indices_batch(cb_list) -> np.ndarray`
- This new batch method assigns indices in one go, reducing dict-update overhead

**Evaluation**:
- Measure global `BarcodeIndex` lock contention if running in threading mode (should be near-zero since each peak emit is short)
- Speedup specifically on per-peak-emit phase: ~1.5× over per-read indexing

**Test**:
- Unit test: assign batch of 100 CBs, verify indices are sequential and consistent with single-CB calls
- Regression test: matrix.mtx output byte-identical to baseline

### 5.5 Issue P-5 — Pipeline architecture (Reader → Finder → Writer)

**What we need**:
- Three processes connected by `multiprocessing.Queue`:
  - **Reader**: opens BAM with `threads=4`, iterates reads, calls `read_check` (with batch CB encoding), buffers `BATCH_SIZE=10000` validated reads into a list of compact tuples, puts the batch on the reader→finder queue
  - **Finder**: consumes batches, iterates internally (no IPC), feeds reads to the existing streaming/bisect peak detection state machine, emits completed peaks (with cb_dict, positions, strategy outputs) onto the finder→writer queue
  - **Writer**: consumes peaks, calls strategy.find_pas(peak), writes BED + MTX
- Bounded queues (`maxsize=100` batches) for backpressure
- Sentinel value (`None`) on each queue to signal end-of-stream
- Each process logs its rate; main process logs overall throughput

**Evaluation**:
- Stage rate measurement: reader, finder, writer each report reads/sec (or batches/sec)
- Slowest stage = pipeline throughput
- Expected speedup: ~2–3× over current single-thread (depends on which stage was the bottleneck before)
- Memory: queue buffers bounded, peak RSS should not exceed 2 GB

**Test**:
- Smoke test: 100K-read BAM through the pipeline, verify same output as monolithic version
- Stress test: full BAM, verify no queue deadlocks, no memory growth over time
- Regression test: BED + MTX output byte-identical to monolithic baseline

### 5.6 Issue P-6 — Tile-based parallelism

**What we need**:
- Split each chromosome into tiles of ~25 Mb each (configurable via `--tile-size`)
- For each tile, fetch reads from `[tile_start - OVERLAP, tile_end + OVERLAP]` where `OVERLAP = 10000` bp (>> max peak width)
- Process N tiles in parallel via `multiprocessing.Pool`
- Each tile worker has its **own** `Peak.pasnumber` counter starting at 0 (use `Peak.reset_pasnumber()` already added by Agent D)
- Each tile worker has its **own** `BarcodeIndex` instance (use `reset_index()` already added)
- Each tile emits peaks **only if their start coordinate falls within `[tile_start, tile_end)`** — overlap reads provide context but peaks in overlap region are suppressed (the next tile owns them)
- After all tiles complete, main process merges:
  - Concatenate all tile BED files in genomic order
  - Reassign `pasnumber` sequentially (1, 2, 3, ...) in genomic order
  - Union all per-tile CB sets, assign canonical column indices
  - Rewrite each tile's MTX file with new pasnumber rows and remapped CB columns
  - Output a single merged BED + MTX

**Evaluation**:
- N-tile parallelism on K-core machine: expected ~min(N, K) × speedup
- Verify identical peak count and coordinates vs single-threaded baseline
- Verify no peaks at tile boundaries are dropped or duplicated

**Test**:
- Boundary test: place a peak deliberately spanning a tile boundary in synthetic data, verify it's emitted exactly once with correct boundaries
- Regression test: full pipeline on chr22, BED output identical to baseline (modulo deterministic pasnumber renumbering)
- Determinism test: run twice, verify byte-identical output

---

## 6. Phase plan

### Phase 0 — Baseline benchmarks (me, before any code change)

- Run `time ema --config single_sample_full.yaml` 3× on full BAM
- Profile with `python-performance-optimization` skill (cProfile)
- Measure peak RSS with `memory-profiling` skill
- Commit `reports/peak_calling_baseline.txt` so we can compare each optimization

### Phase 1 — Quick wins (sequential, me, ~2 hours)

| Order | Issue | Files | Skill | Risk |
|---|---|---|---|---|
| 1a | P-1 (`threads=4`) | `peackcalling.py`, `manager.py`, `cli.py` | `python-resource-management` | trivial |
| 1b | P-2 (`SortedList`) | `peackcalling.py`, `requirments.txt` | `python-performance-optimization` | trivial |

After Phase 1: re-run benchmarks. Expected ~3× speedup combined.

### Phase 2 — CB encoding (sequential, me, ~1 day)

| Order | Issue | Files | Skill |
|---|---|---|---|
| 2a | P-3 (2-bit CB encode) | new `cb_encode.py`, `read.py`, `indexing.py`, `paswrite.py`, `matrixfilter.py` | `python-performance-optimization`, `python-parallel-data-streaming` (numpy batch ops) |
| 2b | P-4 (batch resolve at emit) | `peak.py`, `paswrite.py`, `indexing.py` | `python-performance-optimization` |

After Phase 2: re-run benchmarks. Expected ~2× additional speedup.

### Phase 3 — Pipeline (Sonnet agent, ~1 day)

| Order | Issue | Files | Skill |
|---|---|---|---|
| 3 | P-5 (3-stage pipeline) | new `peak_pipeline.py`, modify `peackcalling.py` to dispatch | `python-parallel-data-streaming` (multiprocessing Queue patterns), `async-python-patterns` |

After Phase 3: re-run benchmarks. Expected ~2× additional speedup.

### Phase 4 — Tile parallelism (Sonnet agent, ~2 days)

| Order | Issue | Files | Skill |
|---|---|---|---|
| 4 | P-6 (tile parallelism + per-worker isolation + merge) | new `tile_runner.py`, modify `peackcalling.py` and `main.py` | `python-parallel-data-streaming`, `python-resource-management` (process pool lifecycle) |

After Phase 4: re-run benchmarks. Expected ~4× additional speedup on multi-core.

### Phase 5 — Validation (sequential, me)

- Run full pipeline on full BAM with all optimizations enabled, compare to Phase 0 baseline
- Verify byte-identical BED output (modulo pasnumber renumbering, which is deterministic)
- Verify byte-identical MTX output (modulo CB column ordering, which we sort canonically)
- Run all existing tests + new perf tests
- Profile final state, commit `reports/peak_calling_optimized.txt`

---

## 7. HPC skills usage and patterns

### 7.1 Skills mapped to phases

| Phase | Skill | What pattern |
|---|---|---|
| 0 | `python-performance-optimization` | cProfile → identify top 10 hotspots before any code change |
| 0 | `memory-profiling` | Track peak RSS baseline; flag any optimization that grows memory > 2× |
| 1a | `python-resource-management` | Ensure pysam handles close on all exit paths; threads param doesn't leak |
| 1b | `python-performance-optimization` | Profile-guided verification: re-run cProfile, confirm SortedList ops drop out of top 10 |
| 2a | `python-parallel-data-streaming` | Vectorized numpy batch encoding via `np.frombuffer` + bitshift; avoid per-string Python loops |
| 2b | `python-performance-optimization` | Verify dict-update count drops from per-read to per-peak |
| 3 | `python-parallel-data-streaming` | Producer-consumer with bounded `multiprocessing.Queue`; backpressure semantics; sentinel termination |
| 3 | `python-resource-management` | Process pool lifecycle: explicit `__enter__`/`__exit__` on Pool, proper SIGTERM handling, no zombie processes |
| 4 | `python-parallel-data-streaming` | Per-worker isolation pattern; main-process merge after `pool.map()`; deterministic ordering |
| 4 | `python-resource-management` | Tile worker cleanup: per-tile temp files removed after merge, memory freed via `gc.collect()` between tiles if needed |
| 5 | `python-performance-optimization` | Final cProfile to verify target met; if not, identify what changed unexpectedly |
| 5 | `memory-profiling` | Verify peak RSS within 2× of baseline despite parallel workers |

### 7.2 Patterns explicitly applied

**Producer-consumer with bounded queues** (Phase 3):
- Reader produces batches, Finder consumes; queue bounded to prevent OOM if Finder lags
- Symmetric for Finder→Writer
- Sentinel value (`None`) signals stream end; receiver re-puts sentinel before exiting (multi-consumer compatibility)

**Per-worker isolation + main-process merge** (Phase 4):
- Each worker is fork+exec'd with no shared mutable state
- Worker writes its output to `output_dir/tile_{tile_id}.{bed,mtx,cb}` — no inter-worker communication
- Main process collects, merges, renumbers, deletes per-tile files
- Same pattern we already use in `concat_matrices` for multi-sample mode

**Vectorized batch encoding** (Phase 2):
- 16-char ASCII CB → 16 bytes → `np.frombuffer(buf, dtype=np.uint8)` → translate via lookup table to 2-bit values → bitshift accumulate to int64
- One numpy call per batch of 10K CBs vs 10K Python loops

**Bounded memory via streaming** (Phase 3, 4):
- No stage accumulates more than a single batch in memory at any time
- Tile workers process one tile fully then discard state

### 7.3 What we explicitly do NOT do

- **No shared-memory dict** for `BarcodeIndex` — would require locks, kills parallelism
- **No `dask.distributed`** — overkill for our scale (single machine, < 100 GB data)
- **No `ray`** — extra runtime dependency, no win at our scale
- **No GPU offload** — peak detection is branch-heavy and not amenable to SIMD/GPU
- **No reimplementation of pysam** — htslib via pysam is fast enough; the bottleneck is Python-level

---

## 8. Test strategy

### Unit tests

- `tests/test_cb_encode.py` — round-trip every 4-mer, batch encode = single encode, N-character handling
- `tests/test_sorted_list_swap.py` — `SortedList` semantics match `bisect+list` for our usage
- `tests/test_pipeline_stages.py` — synthetic 1000-read input through reader→finder→writer, verify output
- `tests/test_tile_merge.py` — synthetic 2-tile peaks at boundary, verify dedup and pasnumber renumbering

### Regression tests

For each phase, the full-BAM output (BED + MTX + filterdcb.tsv) must be **functionally identical** to baseline:
- BED: same number of peaks, same coordinates, same strands; pasnumbers may renumber but order in file is preserved
- MTX: same nnz, same total counts; pasnumber rows may renumber; CB column order may differ but contents identical
- filterdcb.tsv: same set of cell barcodes (order may differ)

Comparison script: `scripts/diff_peak_outputs.py` — normalizes pasnumbers and CB column order, then byte-diffs.

### Integration test on full BAM

Run `time ema --config test_run/single_sample_full.yaml` after each phase. Record:
- Wall-clock time
- Peak RSS
- Peak count
- Cell count
- Cluster count

All four metrics must match baseline (within deterministic renumbering).

### Performance tests

- `pytest tests/test_perf.py --benchmark-only` (using `pytest-benchmark`)
- Each optimization must show measured speedup vs prior phase
- If a phase produces no measurable speedup → revert and investigate

---

## 9. Evaluation criteria

A successful merge requires all of:

| Criterion | Pass condition |
|---|---|
| All unit tests pass | `pytest tests/` exits 0 |
| Regression on full BAM | Peak count, MTX nnz, total counts match baseline within 0.1% |
| Phase 1 speedup | ≥ 2× faster than baseline |
| Phase 2 cumulative speedup | ≥ 4× faster than baseline |
| Phase 3 cumulative speedup | ≥ 6× faster than baseline |
| Phase 4 cumulative speedup | ≥ 8× faster than baseline |
| Memory | Peak RSS ≤ 2× baseline |
| Determinism | Two runs of full pipeline produce identical output (after canonical sort) |
| Backward compat | Existing single-sample legacy path produces same output (no opt-out flag needed) |

---

## 10. Risks

| Risk | Mitigation |
|---|---|
| `SortedList` slower than `bisect+list` for very small arrays | Profile both; if `data_array` rarely exceeds 100 elements, keep bisect for that case (hybrid) |
| 2-bit encoding loses information for non-ATCG CBs | Reject reads with N/IUPAC chars in CB (already filtered by `cb_len` check); document this |
| Pipeline IPC overhead exceeds win on small BAMs | Add `--no-pipeline` flag for small jobs; pipeline only enabled if BAM > 100 MB |
| Tile parallelism breaks for very small BAMs (no tile boundary needed) | Single-tile fallback when chromosome < 2 × tile_size |
| `multiprocessing` interaction with pysam objects | Pysam BAM handles do not pickle — workers must open their own handle from the path |
| Tile boundaries cut a peak in unusual ways | 10kb overlap >> max peak width (typically < 1kb); document the assumption |
| Peak.pasnumber reassignment breaks downstream tools that cached old numbers | Only renumber after merge; document that pasnumbers are not stable across runs (already true today) |
| `multiprocessing.Queue` pickling cost dominates on small reads | Batch size 10000 amortizes; configurable via `--batch-size` |

---

## 11. Open questions / future work

(Not addressed in this branch.)

- Custom Cython BAM parser bypassing pysam — only if profiling shows pysam attribute access still dominates after all phases
- GPU acceleration for coverage profile computation — research project, not engineering task
- Distributed peak calling across multiple machines (dask/ray) — only relevant for cohort-scale data (100s of BAMs)
- Memory-mapped output of MTX for streaming write — only if RAM becomes a constraint
- Profile-guided JIT compilation via Numba on the inner loop — possible additional 2× win but complicates debugging

---

## 12. Appendix: file inventory

### New files

```
ema/countmatrix/cb_encode.py          — Phase 2 (encode/decode/batch encode)
ema/countmatrix/peak_pipeline.py      — Phase 3 (Reader, Finder, Writer process classes)
ema/countmatrix/tile_runner.py        — Phase 4 (tile splitter, worker, merger)
tests/test_cb_encode.py               — Phase 2
tests/test_sorted_list_swap.py        — Phase 1
tests/test_pipeline_stages.py         — Phase 3
tests/test_tile_merge.py              — Phase 4
tests/test_perf.py                    — pytest-benchmark suite (all phases)
scripts/diff_peak_outputs.py          — regression diff tool
reports/peak_calling_baseline.txt     — Phase 0 cProfile output
reports/peak_calling_optimized.txt    — Phase 5 cProfile output
plans/peak-calling-perf.md            — this file
```

### Modified files

```
ema/countmatrix/peackcalling.py       — Phase 1 (threads, SortedList), Phase 3 (pipeline dispatch), Phase 4 (tile dispatch)
ema/countmatrix/peak.py               — Phase 2 (int64 cb_dict keys)
ema/countmatrix/read.py               — Phase 2 (encode CB before return)
ema/countmatrix/indexing.py           — Phase 2 (int64 keys, batch resolve method)
ema/countmatrix/paswrite.py           — Phase 2, Phase 4 (batch write, merge support)
ema/datasets/manager.py               — Phase 1 (threads on AlignmentFile opens)
ema/cli.py                            — new flags: --bam-threads, --tile-size, --batch-size, --no-pipeline
ema/config.py                         — expose new flags
ema/main.py                           — Phase 4 (dispatch tile runner)
requirments.txt                       — add sortedcontainers, pytest-benchmark
example.yaml                          — document new keys
```

### Deleted files

None.
