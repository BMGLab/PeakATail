from ema.countmatrix.indexing import indexing

direction_dict = {True:[1,"-"], False:[-1, "+"]}
score = 0
def pas_write(chro:int, peak_start:int, peak_end:int, pas_start:int, pas_end:int, strand:bool, read_count:int, pasnumber, output):
    peak_list = [peak_start, peak_end][::direction_dict[strand][0]]
    pas_list = [pas_start, pas_end][::-direction_dict[strand][0]]
    
    peak_bed = f"{chro}\t{peak_list[0]}\t{peak_list[1]}\t{pasnumber}\t{read_count}\t{direction_dict[strand][1]}\t{pas_list[0]}\t{pas_list[1]}\n"
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