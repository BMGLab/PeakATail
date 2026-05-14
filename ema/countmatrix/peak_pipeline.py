"""3-stage pipeline for peak calling: Reader → Finder → Writer.

Overlaps BAM I/O with peak detection CPU work and BED/MTX write I/O by
running each stage in a dedicated subprocess connected by bounded
``multiprocessing.Queue`` instances.

Stage summary
-------------
Reader
    Opens the BAM file with ``threads=bam_threads`` for parallel BGZF
    decompression, iterates reads, calls :func:`read_check` to validate and
    extract primitive fields, and puts batches of ``batch_size`` tuples onto
    the reader→finder queue.  Each item on the queue is a ``list`` of
    ``(chrom, start, end, strand, cb)`` tuples — one put per *batch*, not per
    read, so IPC overhead is amortized.

Finder
    Consumes read batches, feeds reads into the same streaming
    ``SortedList``/bisect state machine used by the monolithic
    :func:`peak_calling`, and puts completed :class:`~ema.countmatrix.peak.Peak`
    snapshots onto the finder→writer queue.

Writer
    Consumes ``Peak`` objects, calls ``strategy.find_pas`` + ``get_cb_dict_for_pas``,
    and writes BED and MTX output using the same :func:`pas_write` /
    :func:`matrix_write` helpers as the monolithic path.

Design notes
------------
- ``multiprocessing.get_context('spawn')`` avoids fork-state issues with
  downstream BLAS / NumPy thread pools that may have been initialised in the
  parent process.
- Both queues are bounded (``maxsize=20``) for backpressure — the reader will
  block once the finder falls behind, keeping mLOPMENT - 594729 - 2526Bemory bounded.
- Sentinel ``None`` propagates end-of-stream.  Because each queue has exactly
  one consumer, re-putting the sentinel before exiting is not strictly required
  but is done for defensive correctness in case callers add consumers later.
- pysam ``AlignmentFile`` objects are **not** picklable.  The reader extracts
  primitive tuples before crossing the process boundary.
- Strategy objects *may* hold un-picklable state.  We pass the strategy name
  (a plain string) to subprocesses; each subprocess re-instantiates the
  strategy locally via :func:`~ema.strategies.get_strategy`.
- :class:`~ema.countmatrix.indexing.BarcodeIndex` is process-global state.
  The writer owns it (via the module-level ``_index`` singleton in
  ``indexing.py``).  The finder does not touch it — cb strings travel inside
  ``Peak.cb_dict`` and ``Peak.cb_positions`` as-is and are resolved to column
  indices only inside ``matrix_write`` (writer side).
- ``Peak.pasnumber`` is a class-level counter.  The writer owns it; the finder
  never increments it.
"""

from __future__ import annotations

import multiprocessing as mp
import sys
from collections import deque
from typing import Optional

from sortedcontainers import SortedList


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _flush_peak(
    peak,
    chro: str,
    signal: bool,
    strategy,
    bedfile,
    matrix,
    min_pas_spacing: int = 0,
    min_pas_prominence: float = 0.0,
) -> None:
    """Emit all PAS for *peak* to BED and MTX files.

    Mirrors the flush logic in the monolithic ``peak_calling``.  Imported
    inside subprocess functions to avoid top-level import at module load time
    (spawn context re-imports this module in each worker).

    Args:
        peak: :class:`~ema.countmatrix.peak.Peak` instance ready for emission.
        chro: Chromosome name for this peak.
        signal: Whether the peak was in-progress (signal=True) or completed
            by exceeding merge_len (signal=False but peak_list non-empty).
            Both paths flush; signal is kept for semantic clarity.
        strategy: Instantiated :class:`~ema.strategies.base.PeakFinderStrategy`.
        bedfile: Writable file-like for BED output.
        matrix: Writable file-like for MTX output.
    """
    from ema.countmatrix.peak import Peak
    from ema.countmatrix.paswrite import pas_write, matrix_write
    from ema.strategies.utils import merge_close_or_low_prominence

    pas_results = strategy.find_pas(peak)
    pas_results = merge_close_or_low_prominence(
        pas_results, peak, strategy, min_pas_spacing, min_pas_prominence,
    )
    for pas_1, pas_2 in pas_results:
        Peak.pasnumber += 1
        pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
        pas_write(chro, pas_1, pas_2, peak.peak_strand, pasnumber=Peak.pasnumber, output=bedfile)
        matrix_write(pas_cb_dict, Peak.pasnumber, matrix)


# ---------------------------------------------------------------------------
# Stage functions — each runs in its own subprocess
# ---------------------------------------------------------------------------


def reader_loop(
    bam_path: str,
    direction: bool,
    batch_size: int,
    bam_threads: int,
    default_sample_id: str,
    out_queue: mp.Queue,
) -> None:
    """Reader stage: iterate BAM, validate, batch, emit on *out_queue*.

    Each batch placed on the queue is a ``list`` of
    ``(chrom_str, start, end, strand_bool, cb_str)`` tuples.  Only valid reads
    (those for which :func:`~ema.countmatrix.read.read_check` returns a non-zero
    chrom) are included.  Invalid reads are silently skipped.

    A sentinel ``None`` is placed on the queue after the last batch to signal
    end-of-stream.

    Args:
        bam_path: Absolute path to the sorted, indexed BAM file.
        direction: Strand filter passed through to :func:`read_check`
            (``True`` = reverse/negative, ``False`` = forward/positive).
        batch_size: Number of validated reads per batch list.
        bam_threads: Threads for pysam BGZF block decompression.
        default_sample_id: Fallback sample ID for reads without an RG tag.
        out_queue: Bounded :class:`multiprocessing.Queue` to place batches on.
    """
    import pysam

    from ema.countmatrix.read import read_check, set_default_sample_id

    set_default_sample_id(default_sample_id)

    bamfile = pysam.AlignmentFile(bam_path, "rb", threads=bam_threads)

    batch: list[tuple[str, int, int, bool, str]] = []
    try:
        for read in bamfile:
            chro1, start1, end1, strand, cb = read_check(read=read, direction=direction)
            if chro1 == 0:
                continue
            batch.append((chro1, start1, end1, strand, cb))
            if len(batch) >= batch_size:
                out_queue.put(batch)
                batch = []

        # Flush the final partial batch
        if batch:
            out_queue.put(batch)
    finally:
        bamfile.close()
        # Signal end-of-stream
        out_queue.put(None)


def finder_loop(
    in_queue: mp.Queue,
    out_queue: mp.Queue,
    direction: bool,
    default_threshold: int,
    merge_len: int,
    strategy_name: str,
    dynamic_threshold: bool,
    floor_threshold: int,
    lambda_fold_change: float,
    lambda_window: int,
) -> None:
    """Finder stage: consume read batches, run streaming peak detection.

    Replicates the state-machine logic from the monolithic :func:`peak_calling`
    exactly, emitting completed :class:`~ema.countmatrix.peak.Peak` objects
    (as ``(chrom, peak)`` tuples) onto *out_queue* for the writer stage.

    The finder does **not** touch :class:`~ema.countmatrix.indexing.BarcodeIndex`
    or :attr:`Peak.pasnumber` — those are writer-owned.

    A sentinel ``None`` is placed on *out_queue* after processing is complete.

    Args:
        in_queue: Queue delivering batches from the reader stage.
        out_queue: Queue to place ``(chrom_str, Peak)`` tuples on.
        direction: Strand filter (passed through to :class:`Peak` constructor).
        default_threshold: Fixed peak height threshold (used when
            *dynamic_threshold* is ``False``).
        merge_len: Maximum gap in bp between sub-peaks before splitting.
        strategy_name: Name key for :func:`~ema.strategies.get_strategy`.
            Used here to confirm the strategy is importable; the finder does
            not call strategy methods — that happens in the writer.
        dynamic_threshold: Enable window-based dynamic threshold.
        floor_threshold: Minimum threshold in dynamic mode.
        lambda_fold_change: Multiplier on local lambda for dynamic threshold.
        lambda_window: Window size in bp for local lambda estimation.
    """
    from collections import deque

    from sortedcontainers import SortedList

    from ema.countmatrix.peak import Peak

    current_threshold = default_threshold
    background_deque: deque[int] = deque()

    data_array: SortedList = SortedList()
    signal = False
    chro = "1"
    l_end: int = 0
    i_end: int = 0
    i: int = 0
    peak = Peak(peak_strand=direction)

    try:
        while True:
            batch = in_queue.get()
            if batch is None:
                # Propagate sentinel and exit
                in_queue.put(None)  # defensive re-put for multi-consumer safety
                break

            for chro1, start1, end1, strand, cb in batch:

                # --- chromosome change: flush pending peak ---
                if chro1 != chro:
                    if signal or len(peak.peak_list) != 0:
                        out_queue.put((chro, peak))

                    signal = False
                    peak = Peak(
                        peak_start=0,
                        peak_strand=direction,
                        peak_list=[],
                        cb_dict={},
                        last_peak_end=0,
                        cb_positions={},
                    )
                    i_end, l_end, i = 0, 0, 0
                    data_array = SortedList()

                    if dynamic_threshold:
                        background_deque.clear()
                        current_threshold = floor_threshold

                # --- dynamic threshold update ---
                if dynamic_threshold:
                    background_deque.append(end1)
                    while background_deque and (end1 - background_deque[0]) > lambda_window:
                        background_deque.popleft()
                    if len(background_deque) > 10:
                        local_lambda = len(background_deque) / lambda_window
                        current_threshold = max(
                            floor_threshold, int(local_lambda * lambda_fold_change)
                        )

                # --- peak accumulation logic (mirrors monolithic exactly) ---
                if signal:
                    l_end = data_array[-current_threshold]
                    if start1 <= l_end:
                        peak.cb_counting(cb=cb)
                        peak.cb_position_counting(end1, cb)

                if signal and start1 > l_end:
                    signal = False
                    peak.last_peak_end = l_end
                    peak.peak_start = 0

                data_array.add(end1)

                if start1 > i_end:
                    slice_loc = data_array.bisect_left(start1)
                    if signal:
                        peak.peak_add(data_array=data_array, slice_loc=slice_loc)
                    del data_array[:slice_loc]
                    i_end = data_array[0]

                height = len(data_array)

                if not signal and height >= current_threshold:
                    signal = True
                    i += 1

                    if start1 - peak.last_peak_end > merge_len and i != 1:
                        # Completed peak — emit it
                        out_queue.put((chro1, peak))
                        peak = Peak(
                            peak_start=start1,
                            peak_strand=direction,
                            peak_list=[],
                            cb_dict={},
                            last_peak_end=0,
                            cb_positions={},
                        )
                    else:
                        peak.peak_start = start1

                chro = chro1

        # --- final flush after all batches consumed ---
        if signal or len(peak.peak_list) != 0:
            out_queue.put((chro, peak))

    finally:
        out_queue.put(None)


def writer_loop(
    in_queue: mp.Queue,
    bedfilepath: str,
    matrixpath: str,
    strategy_name: str,
    min_pas_spacing: int = 0,
    min_pas_prominence: float = 0.0,
) -> None:
    """Writer stage: consume peaks, run strategy, write BED and MTX.

    Owns the :class:`~ema.countmatrix.indexing.BarcodeIndex` singleton and the
    :attr:`~ema.countmatrix.peak.Peak.pasnumber` class counter for this run.
    All BED and MTX output is written here.

    Args:
        in_queue: Queue delivering ``(chrom_str, Peak)`` tuples from the finder.
        bedfilepath: Output path for the BED file.
        matrixpath: Output path for the count matrix (MTX) file.
        strategy_name: Strategy name key for :func:`~ema.strategies.get_strategy`.
    """
    from ema.countmatrix.peak import Peak
    from ema.countmatrix.paswrite import pas_write, matrix_write
    from ema.strategies import get_strategy
    from ema.strategies.utils import merge_close_or_low_prominence

    # Instantiate strategy locally (strategy objects may not be picklable)
    strategy = get_strategy(strategy_name)

    # Reset process-local global state so this worker starts clean
    Peak.reset_pasnumber()

    with open(bedfilepath, "w") as bedfile, open(matrixpath, "w") as matrix:
        while True:
            item = in_queue.get()
            if item is None:
                in_queue.put(None)  # defensive re-put
                break

            chro, peak = item

            # Skip empty peaks (shouldn't happen but guard defensively)
            if not peak.peak_list and not peak.cb_dict:
                continue

            pas_results = strategy.find_pas(peak)
            pas_results = merge_close_or_low_prominence(
                pas_results, peak, strategy, min_pas_spacing, min_pas_prominence,
            )
            for pas_1, pas_2 in pas_results:
                Peak.pasnumber += 1
                pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
                pas_write(
                    chro,
                    pas_1,
                    pas_2,
                    peak.peak_strand,
                    pasnumber=Peak.pasnumber,
                    output=bedfile,
                )
                matrix_write(pas_cb_dict, Peak.pasnumber, matrix)


# ---------------------------------------------------------------------------
# Pipeline orchestrator
# ---------------------------------------------------------------------------


def run_pipeline(
    direction: bool,
    bedfilepath: str,
    matrixpath: str,
    bamfile_dir: str,
    default_threshold: int = 5,
    merge_len: int = 100,
    strategy_name: str = "original",
    dynamic_threshold: bool = False,
    floor_threshold: int = 3,
    lambda_fold_change: float = 2.0,
    lambda_window: int = 5000,
    bam_threads: int = 4,
    batch_size: int = 10000,
    default_sample_id: str = "default",
    progress_client=None,
    min_pas_spacing: int = -1,
    min_pas_prominence: float = 5.0,
) -> None:
    """Run the 3-stage Reader → Finder → Writer pipeline.

    Spawns three subprocesses connected by bounded ``multiprocessing.Queue``
    instances.  Waits for all processes to complete, then re-raises any
    subprocess exceptions so the caller sees a clean failure.

    The writer process is joined last; it terminates only after the finder
    has flushed all peaks, ensuring BED and MTX files are complete before
    control returns to the caller.

    Args:
        direction: Strand direction (``True`` = reverse/negative,
            ``False`` = forward/positive).
        bedfilepath: Output path for the BED file.
        matrixpath: Output path for the count matrix.
        bamfile_dir: Path to the sorted, indexed BAM file.
        default_threshold: Fixed peak height threshold (default 5).
        merge_len: Maximum gap between sub-peaks before splitting (default 100).
        strategy_name: Strategy key for :func:`~ema.strategies.get_strategy`
            (default ``"original"``).
        dynamic_threshold: Enable window-based dynamic threshold (default
            ``False``).
        floor_threshold: Minimum threshold in dynamic mode (default 3).
        lambda_fold_change: Multiplier on local lambda (default 2.0).
        lambda_window: Window size in bp for lambda estimation (default 5000).
        bam_threads: Threads for pysam BGZF block decompression (default 4).
        batch_size: Number of validated reads per reader batch (default 10000).
        default_sample_id: Fallback sample ID for reads without an RG tag
            (default ``"default"``).
        progress_client: Optional :class:`~ema.progress.ProgressClient` for
            progress reporting.  When provided, ``advance(1)`` is called after
            the pipeline completes (one tick per direction call).  Pass ``None``
            (default) to disable — all callers without a ``ProgressManager``
            remain unaffected.  This replaces the tqdm-based per-batch progress
            that was planned for this stage; a single tick per strand pass is
            the safe, no-overhead alternative for the spawn-subprocess design.

    Raises:
        RuntimeError: If any subprocess exits with a non-zero exit code.
    """
    # Resolve auto-detect sentinel up front so the spawned writer subprocess
    # receives a concrete value (subprocesses can't read variable_config's
    # cache the spawned process has its own copy).
    if min_pas_spacing < 0:
        from ema.countmatrix.bam_utils import infer_median_read_length
        min_pas_spacing = infer_median_read_length(bamfile_dir)

    ctx = mp.get_context("spawn")

    # Two bounded queues for backpressure
    # reader → finder: batches of (chrom, start, end, strand, cb) tuples
    reader_to_finder: mp.Queue = ctx.Queue(maxsize=20)
    # finder → writer: (chrom_str, Peak) tuples
    finder_to_writer: mp.Queue = ctx.Queue(maxsize=20)

    reader_proc = ctx.Process(
        target=reader_loop,
        name="pipeline-reader",
        args=(
            bamfile_dir,
            direction,
            batch_size,
            bam_threads,
            default_sample_id,
            reader_to_finder,
        ),
        daemon=True,
    )

    finder_proc = ctx.Process(
        target=finder_loop,
        name="pipeline-finder",
        args=(
            reader_to_finder,
            finder_to_writer,
            direction,
            default_threshold,
            merge_len,
            strategy_name,
            dynamic_threshold,
            floor_threshold,
            lambda_fold_change,
            lambda_window,
        ),
        daemon=True,
    )

    writer_proc = ctx.Process(
        target=writer_loop,
        name="pipeline-writer",
        args=(
            finder_to_writer,
            bedfilepath,
            matrixpath,
            strategy_name,
            min_pas_spacing,
            min_pas_prominence,
        ),
        daemon=True,
    )

    # Start in pipeline order (reader last so queue is ready before it writes)
    writer_proc.start()
    finder_proc.start()
    reader_proc.start()

    # Wait for completion in data-flow order: reader → finder → writer
    reader_proc.join()
    finder_proc.join()
    writer_proc.join()

    # Surface subprocess failures to the caller
    failed: list[str] = []
    for proc in (reader_proc, finder_proc, writer_proc):
        if proc.exitcode != 0:
            failed.append(f"{proc.name} exited with code {proc.exitcode}")

    if failed:
        raise RuntimeError(
            "Pipeline subprocess failure(s):\n" + "\n".join(failed)
        )

    # Advance the progress bar by one tick (one strand pass complete).
    if progress_client is not None:
        progress_client.advance(1)
