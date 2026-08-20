import logging
import pysam as ps
import bisect
import time
from collections import deque
from sortedcontainers import SortedList
from ema.countmatrix.peak import Peak
from ema.countmatrix.peak_state import PeakCallingState
from ema.countmatrix.read import read_check
from ema.countmatrix.paswrite import matrix_write, pas_write
from ema.countmatrix.polya import ClipSeeder, check_clip_rate, clip_site, read_umi
from ema.config import directory_config, variable_config
from ema.strategies.utils import merge_close_or_low_prominence
from typing import TYPE_CHECKING

log = logging.getLogger(__name__)

if TYPE_CHECKING:
    from ema.countmatrix.indexing import BarcodeIndex

start_time = time.time()


def peak_calling(
                    direction: bool, bedfilepath: str,
                    matrixpath: str, bamfile_dir=directory_config.bam_dir,
                    default_threshold=variable_config.default_threshold, merge_len=variable_config.merge_len,
                    strategy=None,
                    dynamic_threshold: bool = False,
                    floor_threshold: int = 3,
                    lambda_fold_change: float = 2.0,
                    lambda_window: int = 5000,
                    bam_threads: int = 4,
                    use_pipeline: bool = False,
                    batch_size: int = 10000,
                    # --- tile parallelism (P-6) ---
                    use_tiles: bool = False,
                    tile_size: int = 25_000_000,
                    tile_overlap: int = 10_000,
                    n_workers: int | None = None,
                    default_sample_id: str = "default",
                    # --- encapsulated state (Phase 1) ---
                    index: "BarcodeIndex | None" = None,
                    state: "PeakCallingState | None" = None,
                    sample_id: str | None = None,
                    # --- Phase 2: region fetch (avoid temp-BAM materialisation) ---
                    region: "tuple[str, int, int] | None" = None,
                    # --- per-chromosome progress reporting ---
                    progress_client=None,
                    # --- Post-detection PAS merger (strategy-agnostic) ---
                    # -1 (default) auto-detects median read length per BAM.
                    # 0 disables the distance tier.  5.0 is the default
                    # prominence floor (matching lambda_gradient's per-summit
                    # default).  0.0 disables the prominence tier.
                    min_pas_spacing: int = -1,
                    min_pas_prominence: float = 5.0,
                    # --- Stage 1: read-level poly(A) evidence -----------------
                    # polya_enabled default ON (annotate-only: the clip-read
                    # support of each PAS lands in BED column 5, previously a
                    # hardcoded 0; no coordinates or counts change).  The
                    # clip_seeded strategy additionally SEEDS PAS candidates
                    # from clip-site clusters (see ema/strategies/clip_seeded).
                    polya_enabled: bool = True,
                    polya_min_clip: int = 6,
                    polya_min_purity: float = 0.8,
                    polya_window: int = 100,
                    polya_seed_window: int = 25,
                    polya_min_reads: int = 1,
                ):

    """Stream a BAM file and call peaks using a sliding coverage window.

    Args:
        direction: Strand direction (True = reverse/negative, False = forward/positive).
        bedfilepath: Output path for BED file.
        matrixpath: Output path for count matrix.
        bamfile_dir: Path to indexed BAM file.
        default_threshold: Fixed peak height threshold when dynamic_threshold is off (default: 5).
        merge_len: Max gap between sub-peaks before splitting (default: 100).
        strategy: PeakFinderStrategy instance. Defaults to 'original'.
        dynamic_threshold: When True, use window-based local lambda for threshold.
        floor_threshold: Absolute minimum threshold in dynamic mode (default: 3).
        lambda_fold_change: Multiplier on local lambda for dynamic threshold (default: 2.0).
        lambda_window: Window size in bp for local lambda estimation (default: 5000).
        bam_threads: Number of threads for pysam BGZF block decompression (default: 4).
        use_pipeline: When True, dispatch to the 3-stage Reader->Finder->Writer
            pipeline instead of the monolithic single-threaded loop.  Default
            False preserves existing behaviour.
        batch_size: Read batch size for pipeline mode (default 10000).
        use_tiles: When True, dispatch to tile-based parallel peak calling.
            If both use_tiles and use_pipeline are True, tiles take precedence
            (tiles provide N-core scaling; pipeline only overlaps I/O with CPU
            on a single chromosome at a time -- tiles are strictly more parallel).
            Default False preserves existing behaviour.
        tile_size: Core tile width in bp (default 25 Mb). Only used when
            use_tiles=True.
        tile_overlap: Overlap buffer in bp on each side of a tile (default
            10 kb). Only used when use_tiles=True.
        n_workers: Number of parallel worker processes for tiled mode.
            Defaults to ResourceManager().get_n_jobs(per_worker_mb=500).
            Only used when use_tiles=True.
        default_sample_id: Fallback RG tag value for reads without an RG tag.
            Only used when use_tiles=True.
        index: Optional :class:`~ema.countmatrix.indexing.BarcodeIndex`
            instance.  When ``None`` (default), the module-level singleton is
            used -- preserving backward compatibility for all existing callers.
            New code (e.g. parallel workers) should pass an explicit instance.
            Only applies to the monolithic (non-tile, non-pipeline) path.
        state: Optional :class:`~ema.countmatrix.peak_state.PeakCallingState`
            instance.  When ``None`` (default), a fresh instance is created
            seeded from ``Peak.pasnumber`` -- preserving the legacy reset
            semantics.  New code should pass an explicit instance per worker.
            Only applies to the monolithic path.
        sample_id: Explicit fallback sample identifier for reads without an
            RG tag.  When ``None`` (default), ``_default_sample_id`` from
            :mod:`ema.countmatrix.read` is used (backward compatibility).
            Only applies to the monolithic path.
        region: Optional ``(chrom, start, end)`` 3-tuple specifying a genomic
            region to restrict read fetching.  When set, the BAM is opened
            with ``pysam.AlignmentFile(bam, "rb", threads=1)`` and reads are
            iterated via ``bamfile.fetch(chrom, start, end)`` rather than
            iterating the entire file.  The BAM **must** be indexed (a
            ``.bai`` file must exist alongside it); an ``IOError`` is raised
            with an actionable message if the index is missing.  Use this
            instead of materialising a temporary region BAM to save 2/3 of
            BAM I/O.  Only applies to the monolithic path (ignored when
            ``use_tiles=True`` or ``use_pipeline=True``).
        progress_client: Optional :class:`~ema.progress.ProgressClient`.  When
            supplied, the total is set to the number of chromosomes in the BAM
            and the client is advanced by 1 each time a chromosome boundary is
            crossed.  ``None`` (default) disables per-chromosome reporting and
            preserves backward compatibility for all existing callers.  Only
            applies to the monolithic path.
    """
    # --- Tile dispatch (takes precedence over pipeline) -------------------
    if use_tiles:
        from ema.countmatrix.tile_runner import run_tiled

        # Resolve strategy name for serialisation across process boundary
        if strategy is None:
            strategy_name = "original"
        elif isinstance(strategy, str):
            strategy_name = strategy
        else:
            from ema.strategies import _REGISTRY
            strategy_name = next(
                (k for k, v in _REGISTRY.items() if isinstance(strategy, v)),
                "original",
            )

        return run_tiled(
            direction=direction,
            bedfilepath=bedfilepath,
            matrixpath=matrixpath,
            bamfile_dir=bamfile_dir,
            default_threshold=default_threshold,
            merge_len=merge_len,
            strategy_name=strategy_name,
            dynamic_threshold=dynamic_threshold,
            floor_threshold=floor_threshold,
            lambda_fold_change=lambda_fold_change,
            lambda_window=lambda_window,
            bam_threads=bam_threads,
            tile_size=tile_size,
            tile_overlap=tile_overlap,
            n_workers=n_workers,
            default_sample_id=default_sample_id,
            min_pas_spacing=min_pas_spacing,
            min_pas_prominence=min_pas_prominence,
            polya_enabled=polya_enabled,
            polya_min_clip=polya_min_clip,
            polya_min_purity=polya_min_purity,
            polya_window=polya_window,
            polya_seed_window=polya_seed_window,
            polya_min_reads=polya_min_reads,
        )
    # --- End tile dispatch ------------------------------------------------

    # --- Pipeline dispatch -------------------------------------------------
    if use_pipeline:
        from ema.countmatrix.peak_pipeline import run_pipeline
        from ema.countmatrix.read import _default_sample_id

        # Resolve strategy name: strategy may be a string already (from CLI)
        # or an instantiated strategy object.  Extract the name so the
        # subprocess can re-instantiate locally (strategy objects may not
        # be picklable across spawn context).
        if strategy is None:
            strategy_name = "original"
        elif isinstance(strategy, str):
            strategy_name = strategy
        else:
            # Instantiated object — derive name from class registry
            from ema.strategies import _REGISTRY
            strategy_name = next(
                (k for k, v in _REGISTRY.items() if isinstance(strategy, v)),
                "original",
            )

        return run_pipeline(
            direction=direction,
            bedfilepath=bedfilepath,
            matrixpath=matrixpath,
            bamfile_dir=bamfile_dir,
            default_threshold=default_threshold,
            merge_len=merge_len,
            strategy_name=strategy_name,
            dynamic_threshold=dynamic_threshold,
            floor_threshold=floor_threshold,
            lambda_fold_change=lambda_fold_change,
            lambda_window=lambda_window,
            bam_threads=bam_threads,
            batch_size=batch_size,
            default_sample_id=_default_sample_id,
            min_pas_spacing=min_pas_spacing,
            min_pas_prominence=min_pas_prominence,
            polya_enabled=polya_enabled,
            polya_min_clip=polya_min_clip,
            polya_min_purity=polya_min_purity,
            polya_window=polya_window,
            polya_seed_window=polya_seed_window,
            polya_min_reads=polya_min_reads,
        )
    # --- End pipeline dispatch --------------------------------------------

    if strategy is None:
        from ema.strategies import get_strategy
        strategy = get_strategy("original")

    # --- Stage 1: poly(A) evidence setup -----------------------------------
    # A clip-seeding strategy needs the clip measurements whatever the
    # polya_enabled flag says — refuse the contradictory combination loudly
    # instead of silently emitting a coverage-only call set.
    _polya_seeded = bool(getattr(strategy, "seeds_from_clips", False))
    if _polya_seeded and not polya_enabled:
        raise ValueError(
            "peak strategy 'clip_seeded' requires poly(A) clip evidence; "
            "do not combine it with --polya-evidence off."
        )
    _polya_collect = polya_enabled or _polya_seeded
    _seeder = (
        ClipSeeder(
            direction,
            seed_window=polya_seed_window,
            min_reads=polya_min_reads,
            window=polya_window,
        )
        if _polya_seeded
        else None
    )
    # --- End Stage 1 setup --------------------------------------------------

    # --- Phase 1: encapsulated state setup --------------------------------
    # Resolve BarcodeIndex: use caller-supplied instance if given, otherwise
    # fall back to the module-level singleton (backward-compat default).
    if index is None:
        from ema.countmatrix.indexing import _index as index  # type: ignore[assignment]

    # Resolve PeakCallingState: honour any prior Peak.reset_pasnumber() call
    # by seeding the new state from Peak.pasnumber (the legacy class attr).
    if state is None:
        state = PeakCallingState(pasnumber=Peak.pasnumber)

    # Resolve sample_id fallback for reads without an RG tag.
    # None means "use module singleton inside read_check" (backward compat).
    # --- End Phase 1 setup -----------------------------------------------

    current_threshold = default_threshold

    # Window-based local lambda tracking (deque of recent read positions)
    # Only used when dynamic_threshold=True
    background_deque = deque()

    log.debug("peak_calling: bamfile_dir=%s", bamfile_dir)

    # --- Phase 2: BAI index check for region-fetch path -------------------
    if region is not None:
        import os as _os
        bai_path = bamfile_dir + ".bai"
        if not _os.path.exists(bai_path):
            raise IOError(
                f"BAM index not found for region fetch: {bai_path!r}. "
                "Run `samtools index {bamfile_dir}` to create it."
            )
    # --- End Phase 2 BAI check -------------------------------------------

    # R2 mitigation: estimate the poly(A) clip rate once per BAM per process
    # and warn loudly when the evidence channel looks destroyed.  Full-scan
    # invocations only — tile workers (region != None) inherit the check from
    # their dispatcher, and re-sampling per tile would be pure overhead.
    if _polya_collect and region is None:
        check_clip_rate(
            str(bamfile_dir),
            min_clip=polya_min_clip,
            min_purity=polya_min_purity,
            barcode_tag=variable_config.barcode_tag or "CB",
        )

    # Open BAM: use threads=1 in region mode (fetch already limits I/O);
    # use caller-specified bam_threads in full-scan mode.
    _open_threads = 1 if region is not None else bam_threads
    bamfile = ps.AlignmentFile(bamfile_dir, 'rb', threads=_open_threads)

    # Resolve the auto-detect sentinel for the PAS merger distance tier:
    # -1 means "infer median read length from this BAM, once, then cache".
    # Cache lives in variable_config.dataset_read_lengths so every BAM is
    # only inspected on its first peak_calling() invocation.
    if min_pas_spacing < 0:
        from ema.countmatrix.bam_utils import infer_median_read_length
        _cache = variable_config.dataset_read_lengths
        _key = str(bamfile_dir)
        if _key not in _cache:
            _cache[_key] = infer_median_read_length(_key)
            log.info(
                "Auto-detected median read length %d bp for %s "
                "(min_pas_spacing distance tier)",
                _cache[_key], bamfile_dir,
            )
        min_pas_spacing = _cache[_key]
    # Report total chromosomes to the progress bar (monolithic full-scan only).
    # `len(bamfile.references)` includes every reference in the BAM header —
    # for GRCh38 that's ~194 names with alts/decoys/HLA, most of which carry
    # zero reads.  We only fire one advance per chrom that actually has reads,
    # so denominate by `references with mapped reads` × 2 (one peak_calling
    # invocation per strand).  `get_index_statistics()` is O(n_refs) and reads
    # from the BAI, so it's effectively free.
    if progress_client is not None and region is None:
        try:
            nonempty = sum(1 for s in bamfile.get_index_statistics() if s.mapped > 0)
            progress_client.set_total(max(1, nonempty) * 2)
        except Exception:
            pass
    matrix = open(matrixpath, "w")
    bedfile = open(bedfilepath, "w")
    data_array = SortedList()
    signal = False
    chro = "1"
    l_end, i_end = 0, 0
    timercount = 0
    peak = Peak(peak_strand=direction)
    i = 0 # I forgot what is this but use in peakstarting block

    # --- Stage 1: emission helpers ----------------------------------------
    # _emit_peak replaces the five copy-pasted find_pas/merge/write blocks.
    # With poly(A) evidence off it performs exactly the same calls in the
    # same order as before (score column stays 0) — byte-identical output.
    # With evidence on (default) it additionally writes each PAS's clip-read
    # support into BED column 5.  Under a clip-seeding strategy the PAS are
    # not written immediately: coverage candidates are buffered in the
    # ClipSeeder and written per chromosome by _flush_seeded so clip-site
    # clusters (tier 1) and coverage-only peaks (tier 2) come out merged in
    # coordinate order.
    def _emit_peak(peak_obj, chro_out):
        pas_results = strategy.find_pas(peak_obj)
        pas_results = merge_close_or_low_prominence(
            pas_results, peak_obj, strategy, min_pas_spacing, min_pas_prominence
        )
        for pas_1, pas_2 in pas_results:
            pas_cb_dict = strategy.get_cb_dict_for_pas(peak_obj, pas_1, pas_2)
            if _seeder is not None:
                _seeder.add_coverage_pas(pas_1, pas_2, pas_cb_dict)
                continue
            pasnumber = state.bump_pasnumber()
            if _polya_collect:
                _support, _ = peak_obj.polya_support(
                    pas_1, pas_2, direction, polya_window
                )
            else:
                _support = 0
            pas_write(
                chro_out, pas_1, pas_2, direction,
                pasnumber=pasnumber, output=bedfile, score=_support,
            )
            matrix_write(pas_cb_dict, pasnumber, matrix, index=index)

    def _flush_seeded(chro_out):
        if _seeder is None:
            return
        for _start, _end, _support, _cb_dict in _seeder.flush():
            pasnumber = state.bump_pasnumber()
            pas_write(
                chro_out, _start, _end, direction,
                pasnumber=pasnumber, output=bedfile, score=_support,
            )
            matrix_write(_cb_dict, pasnumber, matrix, index=index)
    # --- End Stage 1 helpers -----------------------------------------------

    # Build the read iterator: region fetch or full-BAM scan.
    if region is not None:
        _r_chrom, _r_start, _r_end = region
        _read_iter = bamfile.fetch(_r_chrom, _r_start, _r_end)
    else:
        _read_iter = bamfile

    for read in _read_iter:
        if timercount%1000000 == 0:#controling time
            endtime = time.time()
            log.debug("elapsed: %.1fs", endtime - start_time)
        timercount += 1

        # checking read validity if it is not countinue to next ittirate
        #chro1 is play role as condition check
        # `strand` is unused: read_check only accepts reads whose strand
        # equals `direction`, and the emit helpers take the strand from
        # `direction` directly rather than from this loop variable (which,
        # on a pass that accepted no reads at all, would still hold the
        # sentinel 0 from the last rejected read).
        chro1, start1, end1, _strand, cb = read_check(
            read=read, direction=direction, sample_id=sample_id
        )

        if chro1 == 0:
            continue

        # Stage 1: take the clip measurement from the AlignedSegment at the
        # call site, after the read_check success — read_check's 5-tuple is
        # shared by the monolithic / pipeline / tile paths and must NOT be
        # widened (Stage-0 stop signal).  The UMI lookup only happens for
        # the ~1% of reads that actually carry a qualifying clip.
        _clip = _umi = None
        if _polya_collect:
            _clip = clip_site(read, polya_min_clip, polya_min_purity)
            if _clip is not None:
                _umi = read_umi(read)

        if chro1 != chro:  # Chromosome changed — flush peaks from OLD chromosome
            if signal:
                _emit_peak(peak, chro)
            elif len(peak.peak_list) != 0:
                _emit_peak(peak, chro)

            # Seeded mode: the OLD chromosome's clip clusters + buffered
            # coverage candidates are written here, in coordinate order.
            _flush_seeded(chro)

            signal = False
            peak = Peak(peak_start=0, peak_strand=direction, peak_list=[], cb_dict={}, last_peak_end=0, cb_positions={})  # make new instance of Peak class
            i_end, l_end, i = 0, 0, 0
            data_array = SortedList()

            # Advance per-chromosome progress bar (monolithic path).
            if progress_client is not None:
                try:
                    progress_client.advance(1)
                except Exception:
                    pass

            # Clear background deque on chromosome change
            if dynamic_threshold:
                background_deque.clear()
                current_threshold = floor_threshold

        # Seeded mode: EVERY qualifying clip read feeds the chromosome-level
        # accumulator, whether or not a coverage peak ever fires here — ~60%
        # of the clip evidence lies outside every coverage-peak window, which
        # is exactly why Phase 2 seeds candidates instead of filtering.
        # (Placed after the chromosome-change flush so a new chromosome's
        # clips never leak into the old chromosome's clusters.)
        if _seeder is not None and _clip is not None:
            _seeder.add_clip(_clip, cb, _umi)

        # Update window-based local lambda from trailing deque of read positions
        if dynamic_threshold:
            background_deque.append(end1)
            while background_deque and (end1 - background_deque[0]) > lambda_window:
                background_deque.popleft()
            if len(background_deque) > 10:  # need minimum reads for stable estimate
                local_lambda = len(background_deque) / lambda_window
                current_threshold = max(floor_threshold, int(local_lambda * lambda_fold_change))

        if signal:
            l_end = data_array[-current_threshold]  # it takes -N from end, where N is the active threshold
            if start1 <= l_end:  # only count if read is still within peak
                peak.cb_counting(cb=cb)
                peak.cb_position_counting(end1, cb)
                # Phase 1: per-peak clip accumulation, mirroring
                # cb_position_counting (annotate mode scores each emitted
                # PAS from its own peak's clip evidence).
                if _clip is not None:
                    peak.polya_counting(_clip, cb, _umi)


        # If newread start_point is more than l_end(where height is more than threshold)
        # it means peak has been finished
        if signal and  start1 > l_end: #in peak
            signal = False
            peak.last_peak_end = l_end
            peak.peak_start = 0

        data_array.add(end1)

        if start1 > i_end:

            slice_loc = data_array.bisect_left(start1)
            if signal:
                peak.peak_add(data_array=data_array, slice_loc=slice_loc)
            del data_array[:slice_loc]  # in-place slice; faster than full copy
            i_end = data_array[0]

        height = len(data_array)

        if signal == False and height >= current_threshold:
            signal = True
            i += 1  #I forgot what is this but will fix TODO

            if start1 - peak.last_peak_end > merge_len and i!=1:
                #Peak has been completed so find pas
                #write pas
                #make new instance of peak
                _emit_peak(peak, chro1)

                peak = Peak(peak_start=start1, peak_strand=direction, peak_list=[], cb_dict={}, last_peak_end=0, cb_positions={}) #make new instance of Peak class
            else:
                peak.peak_start = start1 #Peak has not been complete so will continue to ass items to peak_list

        chro = chro1# will check for chro check

    bamfile.close()

    # Final flush for last chromosome — without this, the last peak is dropped
    if signal:
        _emit_peak(peak, chro)
    elif len(peak.peak_list) != 0:
        _emit_peak(peak, chro)

    # Seeded mode: final chromosome's clusters + coverage candidates.
    _flush_seeded(chro)

    matrix.close()
    bedfile.close()

    # Advance progress bar for the final chromosome (chromosome-change event
    # does not fire for the last chrom since there is no subsequent read).
    if progress_client is not None and region is None:
        try:
            progress_client.advance(1)
        except Exception:
            pass

    # --- Phase 1: write state back to Peak class attr for backward compat --
    # Callers that read Peak.pasnumber after peak_calling() returns (e.g.
    # main.py checking the count) still see the updated value.
    Peak.pasnumber = state.pasnumber
    # --- End Phase 1 write-back -------------------------------------------

if __name__ == "__main__":
    peak_calling()


