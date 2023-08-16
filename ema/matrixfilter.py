from collections import defaultdict
from .countmatrix.indexing import cb_total
from .config import filter_config, directory_config

import pandas as pd
import scipy.io as sci
import anndata as ad
import scanpy as sc

filtered_cb_list = []


def filter_cb(negativematrixpath=directory_config.negmatrixpath(), positivematrixpath=directory_config.posmatrixpath(), sorted_corrected_sparse_path=directory_config.filterd_matrix(), min_read=filter_config.min_read()):
    """TODO
    """
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

    # Filter out cb with counts less than min_read 
    keep_cb = {int(cb) for cb, count in cb_counts.items() if count >= min_read}

    lines_to_keep = []
    #get lines to keep
    with open(negativematrixpath, "r") as negativematrix, open (positivematrixpath, "r") as positivematrix:
        for line in negativematrix:
            item = [int(item) for item in line.split()]
            if len(item) >= 2 and item[1] in keep_cb:
                lines_to_keep.append(item)

        for line in positivematrix:
            item = [int(item) for item in line.split()]
            if len(item) >= 2 and item[1] in keep_cb:
                lines_to_keep.append(item)
 
    #sort lines by cb index 
    #looping over list and write them
    #filter cb that are has been removed
    sorted_spars_matrix = sorted(lines_to_keep, key=lambda x: x[1])
    last_corrected_index, corrected_index = 0, 0
    with open(sorted_corrected_sparse_path, "w") as sorted_corrected_sparse_list:
        for item in sorted_spars_matrix:
            col_index = item[2]
            if item[2] != last_corrected_index:
                corrected_index += 1
                last_corrected_index = col_index
                filtered_cb_list.append(cb_list[col_index - 1])
                item[1] = corrected_index
                sorted_corrected_sparse_list.write(f"{item[0]} {item[1]} {item[2]}\n")
                
            else:
                item[1] = corrected_index
                sorted_corrected_sparse_list.write(f"{item[0]} {item[1]} {item[2]}\n")
    
    
def make_dataframe(matrixpath=directory_config.matrixpath(), collist=filtered_cb_list):
    sparsematrix = sci.mmread(matrixpath)
    
    matrix_frame = pd.DataFrame(sparsematrix.toarray(), columns=collist)
    return matrix_frame.transpose()

def preprocessing(df:pd.DataFrame, min_cells=filter_config.min_cells(), min_genes=filter_config.min_genes()):
    adata = ad.AnnData(X=df)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata

