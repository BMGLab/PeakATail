from collections import defaultdict
from ema.countmatrix.indexing import get_mapping
from ema.config import filter_config, directory_config, variable_config
import numpy as np
import scipy.io as sci
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

filtered_cb_list = []


def filter_cb(negativematrixpath=directory_config.negmatrixpath,
                positivematrixpath=directory_config.posmatrixpath,
                sorted_corrected_sparse_path=directory_config.filterd_matrix,
                min_read=filter_config.min_read,
                filter_cb_file=directory_config.filtered_cb):

    """Filter cell barcodes by minimum read count and write a corrected sparse matrix.

    Reads positive and negative matrices, sums counts per barcode, keeps only
    barcodes with total counts >= min_read, re-indexes columns contiguously,
    and writes the filtered MatrixMarket file.
    """
    global filtered_cb_list
    cb_list = list(get_mapping().keys())
    cb_counts = defaultdict(int)

    with open(negativematrixpath, "r") as negativematrix, open(positivematrixpath, "r") as positivematrix:
        for line in negativematrix:
            columns = line.split()
            if len(columns) >= 3:
                cb_counts[columns[1]] += int(columns[2])

        for line in positivematrix:
            columns = line.split()
            if len(columns) >= 3:
                cb_counts[columns[1]] += int(columns[2])

    # Filter out cb with counts less than min_read
    keep_cb = {int(cb) for cb, count in cb_counts.items() if count >= min_read}

    lines_to_keep = []
    pas_set = set()
    with open(negativematrixpath, "r") as negativematrix, open(positivematrixpath, "r") as positivematrix:
        for line in negativematrix:
            item = [int(item) for item in line.split()]
            if len(item) >= 3 and item[1] in keep_cb:
                lines_to_keep.append(item)
                pas_set.add(item[0])

        for line in positivematrix:
            item = [int(item) for item in line.split()]
            if len(item) >= 3 and item[1] in keep_cb:
                lines_to_keep.append(item)
                pas_set.add(item[0])

    matrix_pas_header, matrix_cb_header, matrix_nzero_header = max(pas_set), len(keep_cb), len(lines_to_keep)
    sorted_spars_matrix = sorted(lines_to_keep, key=lambda x: x[1])
    last_corrected_index, corrected_index = 0, 0
    with open(sorted_corrected_sparse_path, "w") as sorted_corrected_sparse_list:
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


def make_dataframe(matrixpath=directory_config.filterd_matrix, collist=None):
    """Read filtered MatrixMarket file and return sparse matrix with PAS IDs preserved.

    Returns:
        tuple: (sparse_matrix, pas_ids, collist)
            - sparse_matrix: scipy.sparse.csc_matrix (rows=PAS, cols=cells)
            - pas_ids: numpy array of original PAS row indices (1-based from MatrixMarket)
            - collist: list of cell barcode strings
    """
    if collist is None:
        collist = filtered_cb_list

    sparse_coo = sci.mmread(matrixpath)
    sparse_csc = sp.csc_matrix(sparse_coo)

    # Extract unique PAS IDs (row indices) from the COO matrix
    # MatrixMarket is 1-based; scipy mmread converts to 0-based
    # We need the original 1-based PAS IDs for downstream annotation
    if sp.issparse(sparse_coo):
        coo = sparse_coo.tocoo()
        pas_ids = np.arange(1, sparse_csc.shape[0] + 1)
    else:
        pas_ids = np.arange(1, sparse_csc.shape[0] + 1)

    return sparse_csc, pas_ids, collist


def preprocessing(sparse_matrix, pas_ids, collist,
                  min_cells=filter_config.min_cells,
                  min_genes=filter_config.min_genes):
    """Filter sparse PAS-by-cell matrix and return AnnData.

    Args:
        sparse_matrix: scipy.sparse matrix (rows=PAS, cols=cells)
        pas_ids: array of PAS IDs corresponding to rows
        collist: list of cell barcode strings for columns
        min_cells: minimum number of cells a PAS must appear in
        min_genes: minimum number of PAS a cell must have (default lowered for APA)

    Returns:
        ad.AnnData: filtered AnnData object (obs=cells, var=PAS)
    """
    # Build AnnData: scanpy expects cells-by-features (cells x PAS)
    # Our matrix is PAS x cells, so transpose
    adata = ad.AnnData(
        X=sparse_matrix.T.tocsr(),
        obs={'barcode': collist},
        var={'pas_id': pas_ids}
    )
    adata.obs_names = collist
    adata.var_names = [str(pid) for pid in pas_ids]

    # Filter cells: keep cells with at least min_genes PAS detected
    # peaks_per_cell = number of PAS per cell (scanpy n_genes = features per obs)
    sc.pp.filter_cells(adata, min_genes=min_genes)

    # Filter PAS: keep PAS detected in at least min_cells cells
    # cells_per_peak = number of cells per PAS (scanpy n_cells = obs per feature)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    return adata
