from config import variable_config

def read_check(read, barcode:str, direction:bool, barcode_len=variable_config.cb_len() ,seq_len=variable_config.seqlen(), ignore_chro=variable_config.ignore_chro()):
    '''
    fuction check read useubility 
    if it is useable return values
    if it is not return False
    '''
    try:# do not calculate reads don't have CB
        cb = read.get_tag(barcode)
        if len(cb) != barcode_len:
            return None
    except:
        return None
    
    read_chro, read_start, read_end, read_strand = read.reference_name, read.reference_start, read.reference_end, read.is_reverse

    #skip reverse directions
    if direction != read_strand:
        return None
     
    #that chro user Want to ignore them
    if read_chro in ignore_chro:
        return None
    '''
    this block just get reads that len are standart len
    if it is > just skip read 
    if it is < change to standart
    '''
    #TODO 
    #this block will be fixed 
    #6
    if read_end - read_start > seq_len:
        return None
    elif read_end - read_start < seq_len:
        return None 
    
    return read_chro, read_start, read_end, read_strand, cb
    

