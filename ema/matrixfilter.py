from collections import defaultdict
from .countmatrix.indexing import cb_total

def filter_cb(negativematrixpath:str, positivematrixpath:str, output_path:str, filtered_cb_list_path:str, sorted_corrected_sparse_path:str, min_cb_count=2000):
    cb_list = list(cb_total.keys()) # take view of keys covert them to list
    cb_counts = defaultdict(int)

    with open(negativematrixpath, "r") as negativematrix, open (positivematrixpath, "r") as positivematrix:
        for line in negativematrix:# TODO make more efficient with mapping
            columns = line.split() 
            if len(columns) >= 2:
                cb_counts[columns[1]] += int(columns[2]) # columns[1] = cb_index, columns[2] = cbcount

        for line in positivematrix:
            columns = line.split() 
            if len(columns) >= 2:
                cb_counts[columns[1]] += int(columns[2]) # columns[1] = cb_index, columns[2] = cbcount

    # Filter out cb with counts less than min_cb_count 
    keep_cb = {int(cb) for cb, count in cb_counts.items() if count >= min_cb_count}

    lines_to_keep = []
    with open(negativematrixpath, "r") as negativematrix, open (positivematrixpath, "r") as positivematrix, open(output_path, "w"):
        for line in negativematrix:
            item = [int(item) for item in line.split()]
            if len(item) >= 2 and item[1] in keep_cb:
                lines_to_keep.append(item)

        for line in positivematrix:
            item = [int(item) for item in line.split()]
            if len(item) >= 2 and item[1] in keep_cb:
                lines_to_keep.append(item)
    
    sorted_spars_matrix = sorted(lines_to_keep, key=lambda x: x[1])

    last_corrected_index, corrected_index = 0, 0
    with open(filtered_cb_list_path, "w") as filtered_cb_list, open(sorted_corrected_sparse_path, "w") as sorted_corrected_sparse_list:
        for item in sorted_spars_matrix:
            col_index = item[2]
            if item[2] != last_corrected_index:
                corrected_index += 1
                last_corrected_index = col_index
                filtered_cb_list.write(f"{cb_list[col_index - 1]}\n")
                item[1] = corrected_index
                sorted_corrected_sparse_list.write(f"{item[0]} {item[1]} {item[2]}\n")
                
            else:
                item[1] = corrected_index
                sorted_corrected_sparse_list.write(f"{item[0]} {item[1]} {item[2]}\n")
    
    





