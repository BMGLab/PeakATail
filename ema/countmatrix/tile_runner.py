"""Tile-based parallel peak calling (P-6).

Splits each chromosome into tiles of configurable size and processes them in
parallel via ``multiprocessing.Pool``.  Each worker is fully isolated — it has
its own ``Peak.pasnumber`` counter and ``BarcodeIndex`` singleton — and writes
per-tile BED/MTX/CB files to temporary directories.  The main process merges
all tiles into a single global output, renumbering pasnumbers and remapping CB
column indices for global uniqueness.

Key design decisions
--------------------
* **Phase 2 region fetch**: each tile worker calls the monolithic
  :func:`~ema.countmatrix.peackcalling.peak_calling` with
  ``region=(chrom, fetch_start, fetch_end)``.  This uses
  ``pysam.AlignmentFile.fetch()`` to stream only the relevant reads directly
  from the indexed source BAM, eliminating the extra disk write and read pass
  that materialising a per-tile BAM required.  The source BAM must be indexed.
* **Overlap buffer**: each tile fetches reads from
  ``[tile_start - OVERLAP, tile_end + OVERLAP]`` so peaks that straddle tile
  boundaries are detected with full context.  Only peaks whose *start*
  coordinate falls within ``[tile_start, tile_end)`` are retained; the
  adjacent tile owns boundary peaks.
* **spawn context**: ``multiprocessing.get_context('spawn')`` avoids fork
  issues with pysam file descriptors and avoids copying the parent process's
  module-level singletons (``Peak.pasnumber``, ``BarcodeIndex``) into workers.
* **Single-tile short-circuit**: when every chromosome fits inside one tile the
  pool is never created — we call ``peak_calling`` directly.  This avoids
  spawning overhead for small BAMs.
* **ResourceManager**: ``n_workers`` defaults to
  ``ResourceManager().get_n_jobs(per_worker_mb=500)``; callers can override.

Temp-file layout
----------------
Each worker writes into its own ``tempfile.mkdtemp()`` directory::

    /tmp/tile_<tile_id>_<random>/
        peaks.bed
        peaks.mtx
        barcodes.cb.tsv

No region BAM is materialised — reads are fetched directly from the indexed
source BAM via ``pysam.AlignmentFile.fetch(chrom, start, end)``.

The directory is removed in a ``try/finally`` block so cleanup happens even on
worker error.
"""

from __future__ import annotations

import gc
import json
import logging
import multiprocessing
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from ema.utils import ResourceManager

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Timing helpers
# ---------------------------------------------------------------------------


def _write_timings(timings: list[dict[str, Any]], out_path: Path) -> None:
    """Write per-tile timing records to a JSON file.

    Each record in *timings* must be JSON-serialisable.  The canonical schema
    produced by :func:`run_all_jobs` includes at minimum::

        {
            "dataset_id": str,
            "chrom": str,
            "tile_idx": int,
            "tile_start": int,
            "tile_end": int,
            "wall_seconds": float,
        }

    The file is written atomically: data is serialised to a string first; if
    serialisation fails the file is not touched.  Missing parent directories
    are created automatically.

    Args:
        timings: List of timing record dicts (one per tile job).
        out_path: Destination path for the JSON file.
    """
    try:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(timings, indent=2)
        out_path.write_text(payload)
        logger.debug("Tile timings written to %s (%d records)", out_path, len(timings))
    except Exception as exc:
        logger.warning("Could not write tile timings to %s: %s", out_path, exc)


# ---------------------------------------------------------------------------
# JobSpec — flat, pickle-safe descriptor for one tile × direction job
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class JobSpec:
    """Immutable descriptor for a single tile × direction peak-calling job.

    Frozen so it is hashable and safe to pickle across spawn boundaries.

    Attributes:
        job_id: Unique zero-based index across the entire multi-dataset job list.
        dataset_id: Logical dataset identifier (e.g. ``"sample1"``).
        bam_path: Path to the indexed source BAM file.
        chrom: Chromosome name (e.g. ``"chr22"``).
        tile_start: Start of the *core* tile window (inclusive, 0-based).
        tile_end: End of the core tile window (exclusive).
        fetch_start: Extended fetch start including overlap buffer.
        fetch_end: Extended fetch end including overlap buffer.
        direction: Strand — ``True`` = reverse/negative, ``False`` = forward/positive.
        default_threshold: Fixed peak height threshold.
        merge_len: Max gap before splitting peaks (bp).
        strategy_name: Registered peak-finding strategy name.
        dynamic_threshold: Enable window-based dynamic threshold.
        floor_threshold: Minimum threshold in dynamic mode.
        lambda_fold_change: Fold over local lambda for dynamic threshold.
        lambda_window: Background estimation window in bp.
        bam_threads: pysam BGZF decompression threads.
        default_sample_id: Fallback RG tag value for reads without an RG tag.
    """

    job_id: int
    dataset_id: str
    bam_path: str
    chrom: str
    tile_start: int
    tile_end: int
    fetch_start: int
    fetch_end: int
    direction: bool
    default_threshold: int = 5
    merge_len: int = 100
    strategy_name: str = "original"
    dynamic_threshold: bool = False
    floor_threshold: int = 3
    lambda_fold_change: float = 2.0
    lambda_window: int = 5000
    bam_threads: int = 4
    default_sample_id: str = "default"
    # Post-detection PAS merger (strategy-agnostic).  -1 for spacing triggers
    # auto-detect inside the spawned peak_calling; resolve before building
    # JobSpec for best efficiency (single header scan instead of N workers).
    min_pas_spacing: int = -1
    min_pas_prominence: float = 5.0
    # Stage 1: read-level poly(A) evidence (see ema/countmatrix/polya.py).
    polya_enabled: bool = True
    polya_min_clip: int = 6
    polya_min_purity: float = 0.8
    polya_window: int = 100
    polya_seed_window: int = 25
    polya_min_reads: int = 1
    polya_count_window: tuple = (-1, 25)


# ---------------------------------------------------------------------------
# Tile geometry helpers
# ---------------------------------------------------------------------------


def split_chromosome_into_tiles(
    chrom: str,
    chrom_length: int,
    tile_size: int = 25_000_000,
    overlap: int = 10_000,
) -> list[tuple[str, int, int, int, int]]:
    """Partition a chromosome into overlapping tiles.

    Args:
        chrom: Chromosome name (e.g. ``"chr22"`` or ``"22"``).
        chrom_length: Length of the chromosome in bp.
        tile_size: Size of each core tile window in bp (default 25 Mb).
        overlap: Extra context fetched on each side of the tile.  Peaks whose
            start coordinate falls inside the overlap zone are suppressed by
            the tile worker (the adjacent tile owns them).  Should be
            comfortably larger than the widest expected peak (default 10 kb).

    Returns:
        List of ``(chrom, tile_start, tile_end, fetch_start, fetch_end)``
        tuples, one per tile.  ``tile_start``/``tile_end`` define the *core*
        region (exclusive end); ``fetch_start``/``fetch_end`` are the extended
        region used for ``samtools view``.  Coordinates are 0-based and clamped
        to ``[0, chrom_length]``.
    """
    tiles: list[tuple[str, int, int, int, int]] = []
    pos = 0
    while pos < chrom_length:
        tile_start = pos
        tile_end = min(pos + tile_size, chrom_length)
        fetch_start = max(0, tile_start - overlap)
        fetch_end = min(chrom_length, tile_end + overlap)
        tiles.append((chrom, tile_start, tile_end, fetch_start, fetch_end))
        pos = tile_end
    return tiles


def get_chromosomes(bam_path: str) -> list[tuple[str, int]]:
    """Extract chromosome names and lengths from BAM header.

    Args:
        bam_path: Path to indexed BAM file.

    Returns:
        List of ``(chrom_name, length)`` tuples in header order.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return [
            (bam.get_reference_name(i), bam.get_reference_length(bam.get_reference_name(i)))
            for i in range(bam.nreferences)
        ]


# ---------------------------------------------------------------------------
# Worker function (runs in a spawned child process)
# ---------------------------------------------------------------------------


def tile_worker(args: "dict[str, Any] | JobSpec") -> dict[str, Any]:
    """Run peak calling on a single genomic tile.

    This function is executed inside a spawned worker process.  All process-
    global state is therefore independent from every other worker.

    Args:
        args: Either a :class:`JobSpec` (global-pool path, Phase 3) or a legacy
            ``dict`` (single-BAM :func:`run_tiled` path).

            When a ``JobSpec`` is passed the following fields are used:

            ``job_id`` (int)
                Unique zero-based index across the whole multi-dataset job list.
            ``dataset_id`` (str)
                Logical dataset identifier.
            ``bam_path`` (str)
                Path to the source indexed BAM.
            ``chrom`` (str)
                Chromosome name.
            ``tile_start`` / ``tile_end`` (int)
                Core tile window — peaks are only kept if
                ``peak_start in [tile_start, tile_end)``.
            ``fetch_start`` / ``fetch_end`` (int)
                Extended region fetched from the BAM (includes overlap buffer).
            ``direction`` (bool)
                Strand: ``True`` = reverse, ``False`` = forward.
            ``default_threshold`` (int)
            ``merge_len`` (int)
            ``strategy_name`` (str)
            ``dynamic_threshold`` (bool)
            ``floor_threshold`` (int)
            ``lambda_fold_change`` (float)
            ``lambda_window`` (int)
            ``bam_threads`` (int)
            ``default_sample_id`` (str)

    Returns:
        Dict with keys:

            ``tile_id`` (int) — equals ``job_id`` for JobSpec path
            ``dataset_id`` (str) — dataset this tile belongs to
            ``direction`` (bool) — strand
            ``chrom`` (str)
            ``tile_start`` / ``tile_end`` (int)
            ``bed_path`` (str) — per-tile BED (filtered to core region)
            ``mtx_path`` (str) — per-tile MTX (only rows for retained peaks)
            ``cb_path`` (str)  — ordered CB list (TSV, one CB per line)
            ``workdir`` (str)  — temp directory; main process deletes after merge

    Note:
        Phase 2: reads are fetched via ``pysam.AlignmentFile.fetch()`` using
        the ``region=(chrom, fetch_start, fetch_end)`` parameter of
        :func:`~ema.countmatrix.peackcalling.peak_calling`.  No temp BAM is
        written.  The source BAM must be indexed (``.bai`` must exist).
    """
    # --- per-worker imports (in spawned process, all modules are fresh) ---
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.indexing import reset_index, get_mapping
    from ema.countmatrix.peak import Peak
    from ema.countmatrix.read import set_default_sample_id
    from ema.strategies import get_strategy

    # Normalise: accept both JobSpec (Phase 3) and legacy dict (run_tiled).
    if isinstance(args, JobSpec):
        tile_id: int = args.job_id
        dataset_id: str = args.dataset_id
        chrom: str = args.chrom
        tile_start: int = args.tile_start
        tile_end: int = args.tile_end
        fetch_start: int = args.fetch_start
        fetch_end: int = args.fetch_end
        bam_path: str = args.bam_path
        direction: bool = args.direction
        _args_default_threshold = args.default_threshold
        _args_merge_len = args.merge_len
        _args_strategy_name = args.strategy_name
        _args_dynamic_threshold = args.dynamic_threshold
        _args_floor_threshold = args.floor_threshold
        _args_lambda_fold_change = args.lambda_fold_change
        _args_lambda_window = args.lambda_window
        _args_bam_threads = args.bam_threads
        _args_default_sample_id = args.default_sample_id
        _args_min_pas_spacing = args.min_pas_spacing
        _args_min_pas_prominence = args.min_pas_prominence
        _args_polya = dict(
            polya_enabled=args.polya_enabled,
            polya_min_clip=args.polya_min_clip,
            polya_min_purity=args.polya_min_purity,
            polya_window=args.polya_window,
            polya_seed_window=args.polya_seed_window,
            polya_min_reads=args.polya_min_reads,
            polya_count_window=args.polya_count_window,
        )
    else:
        tile_id = args["tile_id"]
        dataset_id = args.get("dataset_id", "default")
        chrom = args["chrom"]
        tile_start = args["tile_start"]
        tile_end = args["tile_end"]
        fetch_start = args["fetch_start"]
        fetch_end = args["fetch_end"]
        bam_path = args["bam_path"]
        direction = args["direction"]
        _args_default_threshold = args["default_threshold"]
        _args_merge_len = args["merge_len"]
        _args_strategy_name = args["strategy_name"]
        _args_dynamic_threshold = args["dynamic_threshold"]
        _args_floor_threshold = args["floor_threshold"]
        _args_lambda_fold_change = args["lambda_fold_change"]
        _args_lambda_window = args["lambda_window"]
        _args_bam_threads = args["bam_threads"]
        _args_default_sample_id = args["default_sample_id"]
        _args_min_pas_spacing = args.get("min_pas_spacing", -1)
        _args_min_pas_prominence = args.get("min_pas_prominence", 5.0)
        _args_polya = dict(
            polya_enabled=args.get("polya_enabled", True),
            polya_min_clip=args.get("polya_min_clip", 6),
            polya_min_purity=args.get("polya_min_purity", 0.8),
            polya_window=args.get("polya_window", 100),
            polya_seed_window=args.get("polya_seed_window", 25),
            polya_min_reads=args.get("polya_min_reads", 1),
            polya_count_window=args.get("polya_count_window", (-1, 25)),
        )

    # Per-worker isolation: reset all process-global mutable state
    reset_index()
    Peak.reset_pasnumber()
    set_default_sample_id(_args_default_sample_id)

    # Use a unique temporary directory so concurrent workers never collide
    workdir = Path(tempfile.mkdtemp(prefix=f"tile_{tile_id}_"))
    raw_bed = workdir / "raw_peaks.bed"
    raw_mtx = workdir / "raw_peaks.mtx"
    final_bed = workdir / "peaks.bed"
    final_mtx = workdir / "peaks.mtx"
    cb_path = workdir / "barcodes.cb.tsv"

    try:
        # ----------------------------------------------------------------
        # 1. Run peak_calling directly via region fetch (Phase 2).
        #    The 10 kb overlap buffer [fetch_start, fetch_end) is passed as
        #    the region so reads are fetched without materialising a temp BAM.
        # ----------------------------------------------------------------
        strategy = get_strategy(_args_strategy_name)
        peak_calling(
            direction,
            str(raw_bed),
            str(raw_mtx),
            bamfile_dir=bam_path,
            default_threshold=_args_default_threshold,
            merge_len=_args_merge_len,
            strategy=strategy,
            dynamic_threshold=_args_dynamic_threshold,
            floor_threshold=_args_floor_threshold,
            lambda_fold_change=_args_lambda_fold_change,
            lambda_window=_args_lambda_window,
            bam_threads=_args_bam_threads,
            region=(chrom, fetch_start, fetch_end),
            min_pas_spacing=_args_min_pas_spacing,
            min_pas_prominence=_args_min_pas_prominence,
            **_args_polya,
        )

        # ----------------------------------------------------------------
        # 2. Filter BED to core region [tile_start, tile_end)
        #    Collect the set of local pasnumbers that survive the filter
        #    so we can apply the same filter to the MTX rows.
        # ----------------------------------------------------------------
        kept_pasnumbers: set[int] = set()
        with open(raw_bed, "r") as src, open(final_bed, "w") as dst:
            for line in src:
                if not line.strip():
                    continue
                fields = line.split("\t")
                # BED cols: chrom, start, end, pasnumber, score, strand
                peak_start_coord = int(fields[1])
                pasnum = int(fields[3])
                if tile_start <= peak_start_coord < tile_end:
                    dst.write(line)
                    kept_pasnumbers.add(pasnum)

        # ----------------------------------------------------------------
        # 3. Filter MTX rows to kept pasnumbers
        # ----------------------------------------------------------------
        with open(raw_mtx, "r") as src, open(final_mtx, "w") as dst:
            for line in src:
                if not line.strip():
                    continue
                row_pas = int(line.split()[0])
                if row_pas in kept_pasnumbers:
                    dst.write(line)

        # ----------------------------------------------------------------
        # 4. Write ordered CB list (index → CB string)
        # ----------------------------------------------------------------
        mapping: dict[str, int] = get_mapping()
        # Sort by assigned column index to get canonical order
        ordered_cbs: list[str] = [
            cb for cb, _ in sorted(mapping.items(), key=lambda x: x[1])
        ]
        with open(cb_path, "w") as f:
            f.write("\n".join(ordered_cbs))
            if ordered_cbs:
                f.write("\n")

        # ----------------------------------------------------------------
        # 5. Remove raw temp files to free disk space early
        # ----------------------------------------------------------------
        raw_bed.unlink(missing_ok=True)
        raw_mtx.unlink(missing_ok=True)

        logger.info(
            "Tile %d done: dataset=%s chrom=%s core=[%d,%d) kept_peaks=%d",
            tile_id, dataset_id, chrom, tile_start, tile_end, len(kept_pasnumbers),
        )

        return {
            "tile_id": tile_id,
            "dataset_id": dataset_id,
            "direction": direction,
            "chrom": chrom,
            "tile_start": tile_start,
            "tile_end": tile_end,
            "bed_path": str(final_bed),
            "mtx_path": str(final_mtx),
            "cb_path": str(cb_path),
            "workdir": str(workdir),
        }

    except Exception:
        # Best-effort cleanup on error
        try:
            shutil.rmtree(workdir, ignore_errors=True)
        except Exception:
            pass
        raise


# ---------------------------------------------------------------------------
# Merge
# ---------------------------------------------------------------------------


def merge_tiles(
    tile_results: list[dict[str, Any]],
    final_bed: str,
    final_mtx: str,
    final_cb: str,
) -> None:
    """Merge per-tile BED/MTX/CB outputs into globally consistent files.

    The merge algorithm:

    1. Sort tile results by ``(chrom, tile_start)`` for deterministic genomic
       ordering.
    2. Concatenate BED files; assign new global pasnumbers 1..N sequentially.
       Build a ``tile_id → {local_pasnum → global_pasnum}`` remapping table.
    3. Union per-tile CB lists in encounter order; build a
       ``tile_id → {local_col → global_col}`` remapping table.
    4. Rewrite per-tile MTX files applying both remappings; append to final MTX.
    5. Write final CB list.
    6. Delete all per-tile temp directories.

    Args:
        tile_results: List of dicts returned by :func:`tile_worker`, one per
            tile.  May include tiles that produced zero peaks (empty files are
            handled gracefully).
        final_bed: Output path for the merged BED file.
        final_mtx: Output path for the merged MTX file.
        final_cb: Output path for the merged CB list.
    """
    # Sort tiles in genomic order
    sorted_tiles = sorted(tile_results, key=lambda r: (r["chrom"], r["tile_start"]))

    global_pasnum = 0
    # tile_id → { local_pas → global_pas }
    pas_remap: dict[int, dict[int, int]] = {}

    # Canonical CB list (encounter order across tiles in genomic order)
    canonical_cbs: list[str] = []
    cb_to_global_col: dict[str, int] = {}

    # tile_id → { local_col → global_col }
    col_remap: dict[int, dict[int, int]] = {}

    # ----------------------------------------------------------------
    # Pass 1: BED — assign global pasnumbers; build CB union
    # ----------------------------------------------------------------
    with open(final_bed, "w") as bed_out:
        for tile in sorted_tiles:
            tile_id = tile["tile_id"]
            pas_remap[tile_id] = {}

            bed_path = Path(tile["bed_path"])
            if not bed_path.exists() or bed_path.stat().st_size == 0:
                continue

            with open(bed_path, "r") as src:
                for line in src:
                    if not line.strip():
                        continue
                    fields = line.split("\t")
                    local_pas = int(fields[3])
                    global_pasnum += 1
                    pas_remap[tile_id][local_pas] = global_pasnum
                    fields[3] = str(global_pasnum)
                    bed_out.write("\t".join(fields))

            # Read per-tile CB list and extend canonical list
            col_remap[tile_id] = {}
            cb_path = Path(tile["cb_path"])
            if cb_path.exists():
                with open(cb_path, "r") as cf:
                    # local_col is 1-based (index in file = 1-based position)
                    for local_col_zero, cb in enumerate(cf):
                        cb = cb.strip()
                        if not cb:
                            continue
                        local_col = local_col_zero + 1  # 1-based
                        if cb not in cb_to_global_col:
                            canonical_cbs.append(cb)
                            cb_to_global_col[cb] = len(canonical_cbs)  # 1-based
                        col_remap[tile_id][local_col] = cb_to_global_col[cb]

    # ----------------------------------------------------------------
    # Pass 2: MTX — remap pasnumber rows and CB columns
    # ----------------------------------------------------------------
    with open(final_mtx, "w") as mtx_out:
        for tile in sorted_tiles:
            tile_id = tile["tile_id"]
            mtx_path = Path(tile["mtx_path"])
            if not mtx_path.exists() or mtx_path.stat().st_size == 0:
                continue

            p_map = pas_remap.get(tile_id, {})
            c_map = col_remap.get(tile_id, {})

            with open(mtx_path, "r") as src:
                for line in src:
                    if not line.strip():
                        continue
                    parts = line.split()
                    local_pas = int(parts[0])
                    local_col = int(parts[1])
                    count = parts[2]

                    new_pas = p_map.get(local_pas)
                    new_col = c_map.get(local_col)

                    if new_pas is None or new_col is None:
                        # Peak/CB was filtered out; skip
                        continue

                    mtx_out.write(f"{new_pas} {new_col} {count}\n")

    # ----------------------------------------------------------------
    # Pass 3: Write canonical CB list
    # ----------------------------------------------------------------
    with open(final_cb, "w") as cb_out:
        for cb in canonical_cbs:
            cb_out.write(cb + "\n")

    # ----------------------------------------------------------------
    # Pass 4: Cleanup temp directories
    # ----------------------------------------------------------------
    for tile in sorted_tiles:
        workdir = tile.get("workdir")
        if workdir:
            try:
                shutil.rmtree(workdir, ignore_errors=True)
            except Exception:
                logger.warning("Failed to remove temp dir %s", workdir)

    logger.info(
        "Merge complete: %d peaks, %d unique CBs",
        global_pasnum,
        len(canonical_cbs),
    )

    gc.collect()


# ---------------------------------------------------------------------------
# Global-pool helpers (Phase 3)
# ---------------------------------------------------------------------------


def build_job_specs(
    bam_list: list[tuple[str, str]],
    directions: list[bool] | None = None,
    tile_size: int = 25_000_000,
    tile_overlap: int = 10_000,
    default_threshold: int = 5,
    merge_len: int = 100,
    strategy_name: str = "original",
    dynamic_threshold: bool = False,
    floor_threshold: int = 3,
    lambda_fold_change: float = 2.0,
    lambda_window: int = 5000,
    bam_threads: int = 4,
    per_bam_tile_sizes: dict[str, int] | None = None,
    min_pas_spacing: int = -1,
    min_pas_prominence: float = 5.0,
    polya_enabled: bool = True,
    polya_min_clip: int = 6,
    polya_min_purity: float = 0.8,
    polya_window: int = 100,
    polya_seed_window: int = 25,
    polya_min_reads: int = 1,
    polya_count_window: tuple = (-1, 25),
) -> list[JobSpec]:
    """Build a flat list of :class:`JobSpec` across all datasets × chroms × tiles × directions.

    Args:
        bam_list: Ordered list of ``(dataset_id, bam_path)`` pairs.
        directions: Strand flags to process. Defaults to ``[False, True]``
            (forward then reverse).
        tile_size: Default core tile width in bp. Overridden per-BAM by
            *per_bam_tile_sizes* when provided.
        tile_overlap: Overlap buffer on each side of a tile in bp.
        default_threshold: Fixed peak height threshold.
        merge_len: Max gap before splitting peaks (bp).
        strategy_name: Registered peak-finding strategy name.
        dynamic_threshold: Enable window-based dynamic threshold.
        floor_threshold: Minimum threshold in dynamic mode.
        lambda_fold_change: Fold over local lambda for dynamic threshold.
        lambda_window: Background estimation window in bp.
        bam_threads: pysam BGZF decompression threads per worker.
        per_bam_tile_sizes: Optional ``{bam_path: tile_size_bp}`` dict for
            RAM-adaptive tile sizing (produced by Phase 4 logic in callers).

    Returns:
        Flat list of :class:`JobSpec` instances, one per
        (dataset × chrom × tile × direction).  The ``job_id`` field is a
        monotonically increasing index across the entire list.
    """
    if directions is None:
        directions = [False, True]

    specs: list[JobSpec] = []
    job_id = 0

    # Auto-detect resolution for ``min_pas_spacing < 0`` is performed lazily
    # inside each spawned worker's :func:`peak_calling` call (which caches
    # per-BAM in its own subprocess).  Doing it here would require opening
    # every BAM at job-build time, which is wrong for callers that may build
    # job specs before BAMs exist (e.g. unit tests with fake paths).

    for dataset_id, bam_path in bam_list:
        bam_tile_size = (
            per_bam_tile_sizes.get(str(bam_path), tile_size)
            if per_bam_tile_sizes
            else tile_size
        )
        per_bam_spacing = min_pas_spacing
        chromosomes = get_chromosomes(str(bam_path))
        for chrom, length in chromosomes:
            tiles = split_chromosome_into_tiles(chrom, length, bam_tile_size, tile_overlap)
            for (_, t_start, t_end, f_start, f_end) in tiles:
                for direction in directions:
                    specs.append(JobSpec(
                        job_id=job_id,
                        dataset_id=dataset_id,
                        bam_path=str(bam_path),
                        chrom=chrom,
                        tile_start=t_start,
                        tile_end=t_end,
                        fetch_start=f_start,
                        fetch_end=f_end,
                        direction=direction,
                        default_threshold=default_threshold,
                        merge_len=merge_len,
                        strategy_name=strategy_name,
                        dynamic_threshold=dynamic_threshold,
                        floor_threshold=floor_threshold,
                        lambda_fold_change=lambda_fold_change,
                        lambda_window=lambda_window,
                        bam_threads=bam_threads,
                        default_sample_id=dataset_id,
                        min_pas_spacing=per_bam_spacing,
                        min_pas_prominence=min_pas_prominence,
                        polya_enabled=polya_enabled,
                        polya_min_clip=polya_min_clip,
                        polya_min_purity=polya_min_purity,
                        polya_window=polya_window,
                        polya_seed_window=polya_seed_window,
                        polya_min_reads=polya_min_reads,
                        polya_count_window=tuple(polya_count_window),
                    ))
                    job_id += 1

    logger.info(
        "build_job_specs: %d datasets × chroms/tiles/dirs → %d total jobs",
        len(bam_list),
        len(specs),
    )
    return specs


def run_all_jobs(
    jobs: list[JobSpec],
    n_workers: int,
    progress_client=None,
    timings_path: Path | None = None,
) -> dict[tuple[str, bool], list[dict[str, Any]]]:
    """Dispatch all jobs through a single global ``multiprocessing.Pool``.

    Opens ONE spawn-context pool, dispatches all jobs via
    ``imap_unordered(chunksize=1)`` for work-stealing, and collects results
    grouped by ``(dataset_id, direction)``.

    Args:
        jobs: Flat list of :class:`JobSpec` instances produced by
            :func:`build_job_specs`.
        n_workers: Number of parallel worker processes.
        progress_client: Optional :class:`~ema.progress.ProgressClient`.
            When supplied, ``advance(1)`` is called after each tile result is
            collected.  Pass ``None`` (default) to disable progress ticking —
            all callers without a ``ProgressManager`` remain unaffected.
        timings_path: Optional path where per-tile timing JSON is written after
            all jobs complete.  When ``None`` (default) no file is written.
            When provided, :func:`_write_timings` is called with a list of
            records keyed by ``dataset_id``, ``chrom``, ``tile_idx``,
            ``tile_start``, ``tile_end``, and ``wall_seconds``.

    Returns:
        Dict mapping ``(dataset_id, direction)`` → list of tile-result dicts
        (same format returned by :func:`tile_worker`), in arrival order.
        Each value list can be passed directly to :func:`merge_tiles`.
    """
    import time

    if not jobs:
        return {}

    actual_workers = min(n_workers, len(jobs))
    logger.info(
        "run_all_jobs: dispatching %d jobs through Pool(%d) workers",
        len(jobs),
        actual_workers,
    )

    # Build a lookup from job_id → (chrom, tile_start, tile_end, tile_idx) for
    # the timing records.  tile_idx is the 0-based position within chrom for
    # each dataset×direction group.
    _job_meta: dict[int, tuple[str, int, int, int, str, bool]] = {}
    # Count per (dataset_id, chrom, direction) to derive tile_idx
    _tile_counters: dict[tuple[str, str, bool], int] = {}
    for j in jobs:
        _key = (j.dataset_id, j.chrom, j.direction)
        _idx = _tile_counters.get(_key, 0)
        _tile_counters[_key] = _idx + 1
        _job_meta[j.job_id] = (j.chrom, j.tile_start, j.tile_end, _idx, j.dataset_id, j.direction)

    ctx = multiprocessing.get_context("spawn")
    grouped: dict[tuple[str, bool], list[dict[str, Any]]] = {}
    timing_records: list[dict[str, Any]] = []

    _t_dispatch = time.monotonic()
    with ctx.Pool(processes=actual_workers) as pool:
        try:
            for result in pool.imap_unordered(tile_worker, jobs, chunksize=1):
                t_now = time.monotonic()
                key = (result["dataset_id"], result["direction"])
                grouped.setdefault(key, []).append(result)
                if progress_client is not None:
                    progress_client.advance(1)
                # Record timing if requested
                if timings_path is not None:
                    _jid = result.get("tile_id")
                    if _jid is not None and _jid in _job_meta:
                        _ch, _ts, _te, _tidx, _ds, _dir = _job_meta[_jid]
                        timing_records.append(
                            {
                                "dataset_id": _ds,
                                "chrom": _ch,
                                "tile_idx": _tidx,
                                "tile_start": _ts,
                                "tile_end": _te,
                                "direction": _dir,
                                "wall_seconds": round(t_now - _t_dispatch, 4),
                            }
                        )
        except Exception:
            pool.terminate()
            raise

    if timings_path is not None:
        _write_timings(timing_records, timings_path)

    total_results = sum(len(v) for v in grouped.values())
    logger.info(
        "run_all_jobs complete: %d results across %d (dataset, direction) groups",
        total_results,
        len(grouped),
    )
    return grouped


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------


def run_tiled(
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
    tile_size: int = 25_000_000,
    tile_overlap: int = 10_000,
    n_workers: int | None = None,
    default_sample_id: str = "default",
    min_pas_spacing: int = -1,
    min_pas_prominence: float = 5.0,
    polya_enabled: bool = True,
    polya_min_clip: int = 6,
    polya_min_purity: float = 0.8,
    polya_window: int = 100,
    polya_seed_window: int = 25,
    polya_min_reads: int = 1,
    polya_count_window: tuple = (-1, 25),
) -> None:
    """Run tile-parallel peak calling and merge results into final output files.

    Chromosomes are split into tiles of *tile_size* bp.  Each tile is processed
    by an independent worker process (spawned, not forked) with its own
    ``Peak.pasnumber`` counter and ``BarcodeIndex`` singleton.  After all
    workers complete the main process merges BED/MTX/CB in genomic order.

    Single-tile short-circuit: when there is only one tile across all
    chromosomes the function bypasses the pool entirely and calls
    :func:`~ema.countmatrix.peackcalling.peak_calling` directly.

    Args:
        direction: Strand — ``True`` = reverse/negative, ``False`` =
            forward/positive.
        bedfilepath: Output path for the final merged BED file.
        matrixpath: Output path for the final merged MTX file.
        bamfile_dir: Path to the indexed source BAM file.
        default_threshold: Fixed peak-height threshold (default 5).
        merge_len: Max gap before splitting peaks (default 100 bp).
        strategy_name: Name of the peak-finding strategy to use
            (passed to :func:`~ema.strategies.get_strategy`).
        dynamic_threshold: Enable window-based dynamic threshold.
        floor_threshold: Minimum threshold in dynamic mode (default 3).
        lambda_fold_change: Fold over local lambda for dynamic threshold.
        lambda_window: Background estimation window in bp (default 5000).
        bam_threads: pysam BGZF decompression threads per worker (default 4).
        tile_size: Core tile width in bp (default 25 Mb).
        tile_overlap: Overlap buffer on each side of a tile in bp (default
            10 kb).  Must be larger than the widest expected peak.
        n_workers: Number of parallel worker processes.  Defaults to
            ``ResourceManager().get_n_jobs(per_worker_mb=500)``.
        default_sample_id: Fallback RG tag value for reads without an RG
            tag.  Passed to :func:`~ema.countmatrix.read.set_default_sample_id`
            inside each worker.

    Note:
        The CB list file is written adjacent to *matrixpath* with the suffix
        ``_cb.tsv`` replaced by the conventional ``filterdcb.tsv`` name used
        elsewhere in the pipeline.  The caller is responsible for passing
        a *matrixpath* that is co-located with the intended CB output.
        The CB output path is derived as::

            <matrixpath_stem>_cb.tsv

        But since the rest of the pipeline writes ``filterdcb.tsv`` via
        :mod:`ema.output`, this function writes to a caller-supplied path
        derived at the call site.  For backward compat, :func:`run_tiled`
        infers the CB path as::

            Path(matrixpath).with_suffix('').with_suffix('') + '_cb.tsv'
            # i.e. /path/to/matrix_cb.tsv

        Callers that need a different CB path should move the file after
        this call.
    """
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.indexing import reset_index
    from ema.countmatrix.peak import Peak
    from ema.countmatrix.read import set_default_sample_id as _set_sid

    if n_workers is None:
        n_workers = ResourceManager().get_n_jobs(per_worker_mb=500)

    # Stage 1: once-per-BAM clip-rate QC in the dispatcher (tile workers use
    # region fetches and inherit this check — see peak_calling).
    _polya_kwargs = dict(
        polya_enabled=polya_enabled,
        polya_min_clip=polya_min_clip,
        polya_min_purity=polya_min_purity,
        polya_window=polya_window,
        polya_seed_window=polya_seed_window,
        polya_min_reads=polya_min_reads,
        polya_count_window=tuple(polya_count_window),
    )
    if polya_enabled:
        from ema.config import variable_config as _vc
        from ema.countmatrix.polya import check_clip_rate
        check_clip_rate(
            str(bamfile_dir),
            min_clip=polya_min_clip,
            min_purity=polya_min_purity,
            barcode_tag=_vc.barcode_tag or "CB",
        )

    # Derive CB output path (sits next to the MTX)
    mtx_path = Path(matrixpath)
    cb_output = str(mtx_path.parent / (mtx_path.stem + "_cb.tsv"))

    logger.info(
        "run_tiled: direction=%s bam=%s tile_size=%d overlap=%d n_workers=%d",
        direction, bamfile_dir, tile_size, tile_overlap, n_workers,
    )

    # Build full tile list across all chromosomes
    chromosomes = get_chromosomes(bamfile_dir)
    all_tiles: list[tuple[str, int, int, int, int]] = []
    for chrom, length in chromosomes:
        tiles = split_chromosome_into_tiles(chrom, length, tile_size, tile_overlap)
        all_tiles.extend(tiles)

    logger.info("Total tiles: %d across %d chromosomes", len(all_tiles), len(chromosomes))

    # ----------------------------------------------------------------
    # Single-tile short-circuit — avoids process pool overhead
    # ----------------------------------------------------------------
    if len(all_tiles) == 1:
        logger.info("Single tile — running monolithic peak_calling directly")
        reset_index()
        Peak.reset_pasnumber()
        _set_sid(default_sample_id)
        from ema.strategies import get_strategy
        strategy = get_strategy(strategy_name)
        peak_calling(
            direction,
            bedfilepath,
            matrixpath,
            bamfile_dir=bamfile_dir,
            default_threshold=default_threshold,
            merge_len=merge_len,
            strategy=strategy,
            dynamic_threshold=dynamic_threshold,
            floor_threshold=floor_threshold,
            lambda_fold_change=lambda_fold_change,
            lambda_window=lambda_window,
            bam_threads=bam_threads,
            min_pas_spacing=min_pas_spacing,
            min_pas_prominence=min_pas_prominence,
            **_polya_kwargs,
        )
        # Write CB list for consistency
        from ema.countmatrix.indexing import get_mapping
        mapping = get_mapping()
        ordered = [cb for cb, _ in sorted(mapping.items(), key=lambda x: x[1])]
        with open(cb_output, "w") as f:
            f.write("\n".join(ordered))
            if ordered:
                f.write("\n")
        return

    # ----------------------------------------------------------------
    # Build worker arg dicts
    # ----------------------------------------------------------------
    # Resolve auto-detect sentinel once for the whole pool before spawning
    # workers — avoids N workers each scanning the BAM header for read length.
    _resolved_spacing = min_pas_spacing
    if _resolved_spacing < 0:
        from ema.countmatrix.bam_utils import infer_median_read_length
        _resolved_spacing = infer_median_read_length(bamfile_dir)
        logger.info(
            "run_tiled: auto-detected median read length %d bp for %s",
            _resolved_spacing, bamfile_dir,
        )

    worker_args: list[dict[str, Any]] = []
    for tile_id, (chrom, tile_start, tile_end, fetch_start, fetch_end) in enumerate(all_tiles):
        worker_args.append({
            "tile_id": tile_id,
            "dataset_id": default_sample_id,
            "chrom": chrom,
            "tile_start": tile_start,
            "tile_end": tile_end,
            "fetch_start": fetch_start,
            "fetch_end": fetch_end,
            "bam_path": bamfile_dir,
            "direction": direction,
            "default_threshold": default_threshold,
            "merge_len": merge_len,
            "strategy_name": strategy_name,
            "dynamic_threshold": dynamic_threshold,
            "floor_threshold": floor_threshold,
            "lambda_fold_change": lambda_fold_change,
            "lambda_window": lambda_window,
            "bam_threads": bam_threads,
            "default_sample_id": default_sample_id,
            "min_pas_spacing": _resolved_spacing,
            "min_pas_prominence": min_pas_prominence,
            **_polya_kwargs,
        })

    # ----------------------------------------------------------------
    # Run workers — spawn context for fork safety
    # ----------------------------------------------------------------
    ctx = multiprocessing.get_context("spawn")
    actual_workers = min(n_workers, len(all_tiles))
    logger.info("Spawning pool of %d workers for %d tiles", actual_workers, len(all_tiles))

    tile_results: list[dict[str, Any]] = []
    with ctx.Pool(processes=actual_workers) as pool:
        try:
            tile_results = pool.map(tile_worker, worker_args)
        except Exception:
            pool.terminate()
            raise

    # ----------------------------------------------------------------
    # Merge
    # ----------------------------------------------------------------
    merge_tiles(tile_results, bedfilepath, matrixpath, cb_output)
    logger.info("run_tiled complete. BED=%s MTX=%s CB=%s", bedfilepath, matrixpath, cb_output)
