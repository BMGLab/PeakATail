import pysam as ps
import bisect
import time
from collections import deque
from sortedcontainers import SortedList
from ema.countmatrix.peak import Peak
from ema.countmatrix.read import read_check
from ema.countmatrix.paswrite import matrix_write, pas_write
from ema.config import directory_config, variable_config

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
                ):

    '''
    Stream a BAM file and call peaks using a sliding coverage window.

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
            on a single chromosome at a time — tiles are strictly more parallel).
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
    '''
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
        )
    # --- End pipeline dispatch --------------------------------------------

    if strategy is None:
        from ema.strategies import get_strategy
        strategy = get_strategy("original")

    current_threshold = default_threshold

    # Window-based local lambda tracking (deque of recent read positions)
    # Only used when dynamic_threshold=True
    background_deque = deque()

    print(bamfile_dir)
    bamfile = ps.AlignmentFile(bamfile_dir, 'rb', threads=bam_threads)
    matrix = open(matrixpath, "w")
    bedfile = open(bedfilepath, "w")
    data_array = SortedList()
    signal = False
    chro = "1"
    l_end, i_end = 0, 0
    timercount = 0
    peak = Peak(peak_strand=direction)
    i = 0 # I forgot what is this but use in peakstarting block

    for read in bamfile:
        if timercount%1000000 == 0:#controling time
            endtime = time.time()
            print(f"{endtime-start_time}")
        timercount += 1
        
        # checking read validity if it is not countinue to next ittirate
        #chro1 is play role as condition check
        chro1, start1, end1, strand, cb = read_check(read=read, direction=direction)

        if chro1 == 0:
            continue
        
        if chro1 != chro:  # Chromosome changed — flush peaks from OLD chromosome
            if signal:

                pas_results = strategy.find_pas(peak)
                for pas_1, pas_2 in pas_results:
                    Peak.pasnumber += 1
                    pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
                    pas_write(chro, pas_1, pas_2, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(pas_cb_dict, Peak.pasnumber, matrix)

            elif len(peak.peak_list) != 0:

                pas_results = strategy.find_pas(peak)
                for pas_1, pas_2 in pas_results:
                    Peak.pasnumber += 1
                    pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
                    pas_write(chro, pas_1, pas_2, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(pas_cb_dict, Peak.pasnumber, matrix)

            signal = False
            peak = Peak(peak_start=0, peak_strand=direction, peak_list=[], cb_dict={}, last_peak_end=0, cb_positions={})  # make new instance of Peak class
            i_end, l_end, i = 0, 0, 0
            data_array = SortedList()

            # Clear background deque on chromosome change
            if dynamic_threshold:
                background_deque.clear()
                current_threshold = floor_threshold

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


        '''
        If newread start_point is more than l_end(where heghit is more than threshold)
        it means peak has been finisfhed
        '''
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
                pas_results = strategy.find_pas(peak)
                for pas_1, pas_2 in pas_results:
                    Peak.pasnumber += 1
                    pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
                    pas_write(chro1, pas_1, pas_2, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(pas_cb_dict, Peak.pasnumber, matrix)

                peak = Peak(peak_start=start1, peak_strand=direction, peak_list=[], cb_dict={}, last_peak_end=0, cb_positions={}) #make new instance of Peak class
            else:
                peak.peak_start = start1 #Peak has not been complete so will continue to ass items to peak_list

        chro = chro1# will check for chro check

    bamfile.close()

    # Final flush for last chromosome — without this, the last peak is dropped
    if signal:
        pas_results = strategy.find_pas(peak)
        for pas_1, pas_2 in pas_results:
            Peak.pasnumber += 1
            pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
            pas_write(chro, pas_1, pas_2, strand, pasnumber=Peak.pasnumber, output=bedfile)
            matrix_write(pas_cb_dict, Peak.pasnumber, matrix)
    elif len(peak.peak_list) != 0:
        pas_results = strategy.find_pas(peak)
        for pas_1, pas_2 in pas_results:
            Peak.pasnumber += 1
            pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
            pas_write(chro, pas_1, pas_2, strand, pasnumber=Peak.pasnumber, output=bedfile)
            matrix_write(pas_cb_dict, Peak.pasnumber, matrix)

    matrix.close()
    bedfile.close()

if __name__ == "__main__":
    peak_calling()


