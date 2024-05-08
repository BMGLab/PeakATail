from collections import defaultdict
from ema.countmatrix.indexing import cb_total
from ema.config import filter_config, directory_config, variable_config
import pandas as pd
import scipy.io as sci
import anndata as ad
import scanpy as sc

filtered_cb_list = []


def filter_cb(negativematrixpath=directory_config.negmatrixpath,
                positivematrixpath=directory_config.posmatrixpath,
                sorted_corrected_sparse_path=directory_config.filterd_matrix, 
                min_read=filter_config.min_read, 
                filter_cb_file=directory_config.filtered_cb):
    
    """TODO
    """
    global filtered_cb_list
    #cb_list = list(cb_total.keys()) # take view of keys covert them to list
    cb_list = list(cb_total.keys())
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

    # Filter out cb with counts less than min_read 
    keep_cb = {int(cb) for cb, count in cb_counts.items() if count >= min_read}
    
    lines_to_keep = []
    pas_set = set()# use for count pas so write to countmatrix hedaer
    #get lines to keep
    with open(negativematrixpath, "r") as negativematrix, open (positivematrixpath, "r") as positivematrix:
        for line in negativematrix:
            item = [int(item) for item in line.split()]
            if len(item) >= 2 and item[1] in keep_cb:
                lines_to_keep.append(item)
                pas_set.add(item[0])

        for line in positivematrix:
            item = [int(item) for item in line.split()]
            if len(item) >= 2 and item[1] in keep_cb:
                lines_to_keep.append(item)
                pas_set.add(item[0])
    
    #get headers for countmatrix
    #sort lines by cb index 
    #looping over list and write them
    #filter cb that are has been removed
    matrix_pas_header, matrix_cb_header, matrix_nzero_header = max(pas_set), len(keep_cb), len(lines_to_keep)
    sorted_spars_matrix = sorted(lines_to_keep, key=lambda x: x[1])
    last_corrected_index, corrected_index = 0, 0
    with open(sorted_corrected_sparse_path, "w") as sorted_corrected_sparse_list:
        #write header for matrixmarket format
        sorted_corrected_sparse_list.write(variable_config.matrixmarketheader)
        sorted_corrected_sparse_list.write(f"{matrix_pas_header} {matrix_cb_header} {matrix_nzero_header}\n")
        
        for item in sorted_spars_matrix:
            col_index = item[1]
            if item[1] != last_corrected_index:
                corrected_index += 1
                last_corrected_index = col_index
                filtered_cb_list.append(cb_list[col_index - 1])
                item[1] = corrected_index
                sorted_corrected_sparse_list.write(f"{item[0]} {item[1]} {item[2]}\n")
                
            else:
                item[1] = corrected_index
                sorted_corrected_sparse_list.write(f"{item[0]} {item[1]} {item[2]}\n")
    
    with open(filter_cb_file, "w") as file:
        for item in filtered_cb_list:
            file.write(f"{item}\n")
    
def make_dataframe(matrixpath=directory_config.filterd_matrix, collist=filtered_cb_list):
    print("in")
    sparsematrix = sci.mmread(matrixpath)
    print(type(sparsematrix))
    print(sparsematrix)
    df = pd.DataFrame(sparsematrix.toarray(), columns=collist)
    print("frame ok")
    return df

def preprocessing(df:pd.DataFrame, min_cells=filter_config.min_cells, 
                  min_genes=filter_config.min_genes):
    
    #each Cell pas count
    pas_count = df[df != 0].count(axis=0)
    df = df.loc[:, pas_count >= min_genes] #filtering Cells
    print("cb ok")
    
    # each pas CsellBarcode count after filtering cells
    cell_count = df[df != 0].count(axis=1)
    df = df[cell_count >= min_cells] #filtering pas
    print("pasok")

    adata = ad.AnnData(X=df.transpose())
    print("adataok")
    return adata

