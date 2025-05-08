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
    next_index = len(cb_total)

    # Create a set of existing barcodes for fast membership checking
    existing_barcodes = set(cb_total.keys())

    # Add new barcodes that don't already exist in the global dictionary
    for barcode in cb_dict:
        if barcode not in existing_barcodes:
            cb_total[barcode] = next_index
            next_index += 1
            existing_barcodes.add(barcode)

    for cb in cb_dict:


        col, count = cb_total[cb], cb_dict[cb]
        matx_str = f"{pasnumber} {col} {count}\n"
        output.write(matx_str)