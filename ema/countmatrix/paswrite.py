from ema.countmatrix.indexing import indexing

direction_dict = {True:[1,"-"], False:[-1, "+"]}
score = 0

def writer(queue, output_file):

    while True:
        line = queue.get()
        if line == "DONE":
            break
        output_file.write(line)

def pas_write(chro:int, peak_start:int, l_end:int, strand:bool,queue, pasnumber, output):
    pas_list = [peak_start, l_end][::-direction_dict[strand][0]]
    peak_bed = f"{chro}\t{pas_list[0]}\t{pas_list[1]}\t{pasnumber}\t{score}\t{direction_dict[strand][1]}\n"
    queue.put(peak_bed)


def matrix_write(cb_dict:dict, pasnumber:int, output, queue):
    '''
    cb_dict: Peak object attribute that contain CellBarcodes of peak class
    pasnumber: Peak class level attribute
    output: Path to out put file 

    Write count matrix in MatrixMarket Format
    '''
    for cb in cb_dict:
        col, count = indexing(cb), cb_dict[cb]
        matx_str = f"{pasnumber} {col} {count}\n"
        queue.put(matx_str)