import pysam as ps
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
    default_threshold=variable_config.default_threshold,
    merge_len=variable_config.merge_len
):
    print(f"Processing BAM: {bamfile_dir}")
    bamfile = ps.AlignmentFile(bamfile_dir, 'rb')
    matrix = open(matrixpath, "w")
    bedfile = open(bedfilepath, "w")
    
    matrix_buffer, bed_buffer = [], []
    data_array = deque()
    signal = False
    chro = "1"
    l_end, i_end = 0, 0
    timercount = 0
    peak = Peak(peak_strand=direction)
    i = 0

    for read in bamfile:
        if timercount % 1_000_000 == 0:
            print(f"{timercount:,} reads processed. Time: {time.time() - start_time:.2f}s")
        timercount += 1
        
        chro1, start1, end1, strand, cb = read_check(read=read, direction=direction)
        if chro1 == 0:
            continue
        
        if chro1 != chro:
            if signal or peak.peak_list:
                pas_1, pas_2 = peak.pasfind()
                if pas_1 != 0:
                    Peak.pasnumber += 1
                    bed_buffer.append(pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=None))
                    matrix_buffer.append(matrix_write(peak.cb_dict, Peak.pasnumber, None))
            signal = False
            peak = Peak(peak_start=0, peak_strand=direction)
            data_array.clear()
            l_end, i_end, i = 0, 0, 0

        if signal:
            if len(data_array) >= default_threshold:
                l_end = list(data_array)[-default_threshold]
            peak.cb_counting(cb=cb)

        if signal and start1 > l_end:
            signal = False
            peak.last_peak_end = l_end
            peak.peak_start = 0

        while data_array and data_array[0] < start1:
            data_array.popleft()
        data_array.append(end1)

        if data_array and start1 > i_end:
            i_end = data_array[0]
            if signal:
                peak.peak_add(data_array=list(data_array), slice_loc=0)

        height = len(data_array)
        if not signal and height >= default_threshold:
            signal = True
            i += 1
            if start1 - peak.last_peak_end > merge_len and i != 1:
                pas_1, pas_2 = peak.pasfind()
                if pas_1 != 0:
                    Peak.pasnumber += 1
                    bed_buffer.append(pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=None))
                    matrix_buffer.append(matrix_write(peak.cb_dict, Peak.pasnumber, None))
                peak = Peak(peak_start=start1, peak_strand=direction)
            else:
                peak.peak_start = start1

        chro = chro1

        if timercount % 1_000_000 == 0:
            matrix.write(''.join(matrix_buffer))
            bedfile.write(''.join(bed_buffer))
            matrix_buffer.clear()
            bed_buffer.clear()

    matrix.write(''.join(matrix_buffer))
    bedfile.write(''.join(bed_buffer))
    matrix.close()
    bedfile.close()
    bamfile.close()

if __name__ == "__main__":
    peak_calling()
