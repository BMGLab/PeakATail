import os 

output_dir = os.path.join(os.getcwd(), "emaout")

class directory_config():
    posmatrixpath = os.path.join(output_dir, "posmatrix.mtx")
    negmatrixpath = os.path.join(output_dir, "negmatrix.mtx")
    matrixpath = os.path.join(output_dir, "pascountmatrix.mtx")
    posbed = os.path.join(output_dir, "posbed.bed")
    negbed = os.path.join(output_dir, "negbed.bed")
    pasbed = os.path.join(output_dir, "pasbed.bed")
    filterd_matrix = os.path.join(output_dir, "filterdmatrix.mtx")
    filtered_cb = os.path.join(output_dir, "filterdcb.tsv")
    cluster_output = os.path.join(output_dir, "write", "pbmc3k.h5ad")
    cluster_labels = os.path.join(output_dir, "clusterlabels.csv")

class variable_config():
    seqlen = 98
    cb_len = 16
    default_threshold = 5
    merge_len = 100
    ignore_chro = ["MT", "mt"]

class filter_config():
    min_read = 2000
    min_cells = 3
    min_genes = 200