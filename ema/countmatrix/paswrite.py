from ema.countmatrix.indexing import indexing
from ema.countmatrix.indexing import cb_total

direction_dict = {True:[1,"-"], False:[-1, "+"]}
score = 0
def pas_write(chro:int, peak_start:int, l_end:int, strand:bool, pasnumber, output):
    pas_list = [peak_start, l_end][::-direction_dict[strand][0]]
    peak_bed = f"{chro}\t{pas_list[0]}\t{pas_list[1]}\t{pasnumber}\t{score}\t{direction_dict[strand][1]}\n"
    output.write(peak_bed)


def matrix_write(cb_dict:dict, pasnumber:int, output): 
    '''
    cb_dict: Peak object attribute that contain CellBarcodes of peak class
    pasnumber: Peak class level attribute
    output: Path to out put file 

    Write count matrix in MatrixMarket Format
    '''
    global cb_total
    cb_num = len(cb_total)
    new_cb_set = set(cb_dict.keys())
    cleaned_cb_set = new_cb_set - cb_total
    index_map = {cb: cb_num+i for i, cb in enumerate(dict.fromkeys(cleaned_cb_set))}
    cb_total = cb_total | index_map
    for cb in cb_dict:


        col, count = cb_total[cb], cb_dict[cb]
        matx_str = f"{pasnumber} {col} {count}\n"
        output.write(matx_str)