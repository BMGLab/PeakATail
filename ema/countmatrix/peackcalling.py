import pysam as ps
import bisect
import time
from ema.countmatrix.peak import Peak
from ema.countmatrix.read import read_check
from ema.countmatrix.paswrite import matrix_write, pas_write
from ema.config import directory_config, variable_config

start_time = time.time()

def peak_calling(   
                    direction:bool, bedfilepath:str,
                    matrixpath:str, bamfile_dir = directory_config.bam_dir ,
                    default_threshold=variable_config.default_threshold, merge_len=variable_config.merge_len
                ):

    '''
    
    data_array
    '''
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
        
        if chro1 != chro:#TODO COMPLETE
            if signal:
                signal = False
                pas_1, pas_2 = peak.pasfind()

                if pas_1 != 0:
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)
            
            elif len(peak.peak_list) != 0:# TODO thsi block has code reaptition
                signal = False
                pas_1, pas_2 = peak.pasfind()

                if pas_1 != 0: # if pasfind method return don't False mean peak have valid pas
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)

            peak = Peak(peak_strand=direction) #make new instance of Peak class
            i_end, l_end = 0, 0

        if signal:
            l_end = data_array[-default_threshold] # it takes -5 from end in default it means where heghit is more than threshold
            peak.cb_counting(cb=cb)

        
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

        if signal == False and height >= default_threshold:
            signal = True
            i += 1  #I forgot what is this but will fix TODO

            if start1 - peak.last_peak_end > merge_len and i!=1:
                #Peak has been completed so find pas 
                #write pas
                #make new instance of peak
                pas_1, pas_2 = peak.pasfind()

                if pas_1 != 0: # if pasfind method return don't  mean peak have valid pas
                    Peak.pasnumber += 1
                    pas_write(chro1, pas_2, pas_1, strand, pasnumber=Peak.pasnumber, output=bedfile)
                    matrix_write(peak.cb_dict, Peak.pasnumber, matrix)

                peak = Peak(peak_start=start1, peak_strand=direction, peak_list=[], cb_dict={}) #make new instance of Peak class
            else:
                peak.peak_start = start1 #Peak has not been complete so will continue to ass items to peak_list

        chro = chro1# will check for chro check

    bamfile.close()

if __name__ == "__main__":
    peak_calling() 


