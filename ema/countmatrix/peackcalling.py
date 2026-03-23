import pysam as ps
import bisect
import time
from collections import deque
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
                    lambda_window: int = 10000,
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
        lambda_window: Window size in bp for local lambda estimation (default: 10000).
    '''
    if strategy is None:
        from ema.strategies import get_strategy
        strategy = get_strategy("original")

    current_threshold = default_threshold

    # Window-based local lambda tracking (deque of recent read positions)
    # Only used when dynamic_threshold=True
    background_deque = deque()

    print(bamfile_dir)
    bamfile = ps.AlignmentFile(bamfile_dir, 'rb')
    matrix = open(matrixpath, "w")
    bedfile = open(bedfilepath, "w")
    data_array = []
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
        
        if chro1 != chro:  # TODO COMPLETE
            if signal:

                pas_results = strategy.find_pas(peak)
                for pas_1, pas_2 in pas_results:
                    Peak.pasnumber += 1
                    pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
                    pas_write(chro1, pas_1, pas_2, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(pas_cb_dict, Peak.pasnumber, matrix)

            elif len(peak.peak_list) != 0:  # TODO this block has code repetition

                pas_results = strategy.find_pas(peak)
                for pas_1, pas_2 in pas_results:
                    Peak.pasnumber += 1
                    pas_cb_dict = strategy.get_cb_dict_for_pas(peak, pas_1, pas_2)
                    pas_write(chro1, pas_1, pas_2, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(pas_cb_dict, Peak.pasnumber, matrix)

            signal = False
            peak = Peak(peak_start=0, peak_strand=direction, peak_list=[], cb_dict={}, last_peak_end=0, cb_positions={})  # make new instance of Peak class
            i_end, l_end, data_array, i = 0, 0, [], 0

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
            
        bisect.insort(data_array, end1)

        if start1 > i_end:
        
            slice_loc = bisect.bisect_left(data_array , start1)
            if signal:
                peak.peak_add(data_array=data_array, slice_loc=slice_loc)
            data_array = data_array[slice_loc: ]# update the list so the [0] index will always be i_end
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

if __name__ == "__main__":
    peak_calling() 


