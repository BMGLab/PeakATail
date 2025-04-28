import pysam as ps
import bisect
import time
from ema.countmatrix.peak import Peak
from ema.countmatrix.read import read_check
from ema.countmatrix.paswrite import matrix_write, pas_write
from ema.config import directory_config, variable_config
from collections import defaultdict

start_time = time.time()

def load_bed_intervals(bedfile_path):
    """
    Load BED file and return interval trees by chromosome.
    """
    intervals = defaultdict(list)
    with open(bedfile_path, "r") as bedfile:
        for line in bedfile:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            intervals[chrom].append((start, end))
    # Sort intervals for binary search
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals

def is_read_in_intervals(chrom, read_start, read_end, intervals_dict):
    """
    Check if a read overlaps any interval in intervals_dict.
    """
    if chrom not in intervals_dict:
        return False

    intervals = intervals_dict[chrom]
    # binary search for any overlap
    idx = bisect.bisect_left(intervals, (read_start, read_start))
    if idx < len(intervals) and intervals[idx][0] <= read_end and intervals[idx][1] >= read_start:
        return True
    if idx > 0 and intervals[idx-1][0] <= read_end and intervals[idx-1][1] >= read_start:
        return True
    return False

def peak_calling(direction: bool, 
                bedfilepath: str,
                matrixpath: str,
                bamfile_dir = directory_config.bam_dir,
                utr_bedfile: str = directory_config.utrbed, # "data/regions_3utr.bed",  # <<< NEW input
                default_threshold = variable_config.default_threshold,
                merge_len = variable_config.merge_len
                ):
    """
    Peak calling restricted to reads overlapping 3'-UTRs.
    """

    print(f"Loading BAM from: {bamfile_dir}")
    print(f"Restricting to UTRs from: {utr_bedfile}")

    # Load UTR intervals
    utr_intervals = load_bed_intervals(utr_bedfile)

    bamfile = ps.AlignmentFile(bamfile_dir, 'rb')
    matrix = open(matrixpath, "w")
    bedfile = open(bedfilepath, "w")
    data_array = [] 
    signal = False
    chro = "1"
    l_end, i_end = 0, 0
    timercount = 0
    peak = Peak(peak_strand=direction)
    i = 0 

    for read in bamfile:
        if timercount % 1000000 == 0:
            endtime = time.time()
            print(f"{timercount} reads processed. Time elapsed: {endtime - start_time:.2f} seconds.")
        timercount += 1

        chro1, start1, end1, strand, cb = read_check(read=read, direction=direction)

        if chro1 == 0:
            continue

        # UTR filtering:
        if not is_read_in_intervals(chro1, start1, end1, utr_intervals):
            continue

        if chro1 != chro:
            if signal:
                pas_1, pas_2 = peak.pasfind()
                if pas_1 != 0:
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)
            elif len(peak.peak_list) != 0:
                pas_1, pas_2 = peak.pasfind()
                if pas_1 != 0:
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)
            signal = False
            peak = Peak(peak_start=0, peak_strand=direction, peak_list=[], cb_dict={},last_peak_end=0)
            i_end, l_end, data_array, i = 0, 0, [], 0

        if signal:
            l_end = data_array[-default_threshold]
            peak.cb_counting(cb=cb)

        if signal and start1 > l_end:
            signal = False
            peak.last_peak_end = l_end
            peak.peak_start = 0

        bisect.insort(data_array, end1)

        if start1 > i_end:
            slice_loc = bisect.bisect_left(data_array, start1)
            if signal:
                peak.peak_add(data_array=data_array, slice_loc=slice_loc)
            data_array = data_array[slice_loc:]
            i_end = data_array[0]

        height = len(data_array)

        if signal == False and height >= default_threshold:
            signal = True
            i += 1
            if start1 - peak.last_peak_end > merge_len and i != 1:
                pas_1, pas_2 = peak.pasfind()
                if pas_1 != 0:
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)
                peak = Peak(peak_start=start1, peak_strand=direction, peak_list=[], cb_dict={}, last_peak_end=0)
            else:
                peak.peak_start = start1

        chro = chro1

    bamfile.close()
    matrix.close()
    bedfile.close()

if __name__ == "__main__":
    peak_calling()
