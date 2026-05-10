from ema.config import variable_config
from ema.countmatrix.cb_encode import encode_cb

# Fallback sample_id used when a read has no RG tag. main.py sets this to the
# current dataset_id before each peak_calling call so multi-sample BAMs without
# pre-existing RG tags still produce unique CB prefixes per dataset.
_default_sample_id = "default"


def set_default_sample_id(sample_id: str) -> None:
    global _default_sample_id
    _default_sample_id = sample_id


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
    except KeyError:
        return 0, 0, 0, 0, 0

    # NOTE: invalid-char filtering is NOT done here — encode_cb is too slow
    # (~3μs × 14M reads = 40+ sec). BarcodeIndex.get_index handles invalid
    # CBs by encoding them to -1 (tuple key (sample_id, -1)); these end up in
    # the same column but are rare (1-letter Ns are <0.1% of CB sequencing data).

    read_chro, read_start, read_end, read_strand = read.reference_name, read.reference_start, read.reference_end, read.is_reverse
    try:
        sample_id = read.get_tag('RG')
    except (KeyError, ValueError):
        sample_id = _default_sample_id

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