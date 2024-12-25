from ema.config import variable_config

def read_check(read, direction:bool, barcode=variable_config.barcode_tag, barcode_len=variable_config.cb_len ,seq_len=variable_config.seqlen, ignore_chro=variable_config.ignore_chro):
    '''
    fuction check read useubility 
    if it is useable return values
    if it is not return 0
    '''
    try:# do not calculate reads don't have CB
        cb = read.get_tag(barcode)
        if len(cb) != barcode_len:
            return 0, 0, 0, 0, 0
    except:
        return 0, 0, 0, 0, 0
    
    if read.get_tag('RG') != None:
        sample_id = read.get_tag('RG')
    else:
        sample_id = "sample"

    read_chro, read_start, read_end, read_strand = read.reference_name, read.reference_start, read.reference_end, read.is_reverse

    #skip reverse directions
    if direction != read_strand:
        return 0, 0, 0, 0, 0
     
    #that chro user Want to ignore them
    if read_chro in ignore_chro:
        return 0, 0, 0, 0, 0
    '''
    this block just get reads that len are standart len
    if it is > just skip read 
    if it is < change to standart
    '''
    #TODO 
    #this block will be fixed 
    #6
    if read_end - read_start > seq_len:
        return 0, 0, 0, 0, 0
    elif read_end - read_start < seq_len:
        read_end = read_start + seq_len
    
    cb = f"{sample_id}_{cb}"
    return read_chro, read_start, read_end, read_strand, cb
    
if __name__ == "__main__":
    read_check()