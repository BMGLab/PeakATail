from ema.countmatrix.indexing import indexing

strand_char = {True: "-", False: "+"}
score = 0
def pas_write(chro:int, peak_start:int, l_end:int, strand:bool, pasnumber, output):
    bed_start = min(peak_start, l_end)
    bed_end = max(peak_start, l_end)
    peak_bed = f"{chro}\t{bed_start}\t{bed_end}\t{pasnumber}\t{score}\t{strand_char[strand]}\n"
    output.write(peak_bed)


def matrix_write(cb_dict:dict, pasnumber:int, output): 
    '''
    cb_dict: Peak object attribute that contain CellBarcodes of peak class
    pasnumber: Peak class level attribute
    output: Path to out put file 

    Write count matrix in MatrixMarket Format
    '''
    for cb in cb_dict:
        col, count = indexing(cb), cb_dict[cb]
        matx_str = f"{pasnumber} {col} {count}\n"
        output.write(matx_str)