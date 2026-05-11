from collections import defaultdict
from ema.countmatrix.indexing import get_mapping
from ema.config import filter_config, directory_config, variable_config
from ema.outputs import MATRIX_MARKET_HEADER
import numpy as np
import scipy.io as sci
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

filtered_cb_list = []


def filter_cb(input_matrix_paths: list = None,
              cb_list: list = None,
              negativematrixpath=None,
              positivematrixpath=None,
              sorted_corrected_sparse_path=None,
              min_read=None,
              filter_cb_file=None):
    """Filter cell barcodes by minimum read count and write a corrected sparse matrix.

    Reads one or more matrices, sums counts per barcode column index, keeps only
    barcodes with total counts >= min_read, re-indexes columns contiguously, and
    writes the filtered MatrixMarket file.

    Args:
        input_matrix_paths: Optional list of MatrixMarket file paths to read.
            When provided, only these files are read (multi-sample path).
            When None, falls back to reading negativematrixpath + positivematrixpath
            (legacy single-sample behaviour).
        cb_list: Optional list of CB strings, indexed by (column_index - 1).
            When provided, this is used as the authoritative CB lookup instead
            of ``get_mapping()``.  Required when ``input_matrix_paths`` is given
            and the global barcode index no longer reflects the correct
            per-dataset CBs.

            .. note::

                For parallel use (e.g. Phase 3 global pool) always supply an
                explicit ``cb_list`` rather than relying on the module-level
                ``get_mapping()`` singleton.  The singleton reflects only the
                last dataset processed in the current process; its value is
                undefined under fork-based multiprocessing.

        negativematrixpath: Legacy negative-strand matrix path (single-sample).
        positivematrixpath: Legacy positive-strand matrix path (single-sample).
        sorted_corrected_sparse_path: Output path for the filtered MatrixMarket file.
        min_read: Minimum total read count for a barcode to be kept.
        filter_cb_file: Output path for the filtered CB list (one CB per line).
    """
    global filtered_cb_list
    filtered_cb_list = []  # reset on each call to avoid cross-call accumulation

    # Resolve config-dependent defaults at call time, not import time —
    # otherwise the captured paths point at the pre-set_directory_config
    # ``emaout/`` defaults and the filter reads the wrong files.
    if negativematrixpath is None:
        negativematrixpath = directory_config.negmatrixpath
    if positivematrixpath is None:
        positivematrixpath = directory_config.posmatrixpath
    if sorted_corrected_sparse_path is None:
        sorted_corrected_sparse_path = directory_config.filterd_matrix
    if filter_cb_file is None:
        filter_cb_file = directory_config.filtered_cb
    if min_read is None:
        min_read = filter_config.min_read

    # Determine the CB lookup list (indexed by column_index - 1)
    if cb_list is not None:
        _cb_lookup = list(cb_list)
    else:
        _cb_lookup = list(get_mapping().keys())

    # Determine which matrix files to read
    if input_matrix_paths is not None:
        matrix_paths = [str(p) for p in input_matrix_paths]
    else:
        matrix_paths = [negativematrixpath, positivematrixpath]

    # Pass 1: sum counts per column index (1-based in file)
    cb_counts: dict[str, int] = defaultdict(int)
    for path in matrix_paths:
        with open(path, "r") as fh:
            first_data_line = True
            for line in fh:
                if not line.strip() or line.startswith("%"):
                    continue
                columns = line.split()
                if len(columns) < 3:
                    continue
                # When reading MatrixMarket files (from concat_matrices output),
                # the first non-% line is the dimension header — skip it.
                # Legacy per-BAM matrices written by peackcalling have no header.
                # Detect header: all three tokens are numeric digits AND it is the
                # first such line.  We use the `first_data_line` flag per file.
                if first_data_line:
                    first_data_line = False
                    # Heuristic: if column[2] is suspiciously large compared to a
                    # normal count value AND columns[0] and [1] are also large, it
                    # is likely the dimension header.  Simpler: just check whether
                    # reading it as a data entry would produce a valid col index.
                    # Actually, we must be careful — for legacy files there is NO
                    # header.  Use: if the number of non-zero entries (col[2]) is
                    # much larger than any realistic single-cell count, it is a
                    # header.  More robust: always try to use the line as data;
                    # if col_index (columns[1]) would be out of range for _cb_lookup
                    # treat as header and skip.
                    try:
                        col_idx = int(columns[1])
                    except ValueError:
                        continue  # not numeric at all — skip
                    if col_idx > len(_cb_lookup) and len(_cb_lookup) > 0:
                        # Column index exceeds CB count — this is a dimension header.
                        continue
                    # Otherwise treat as data
                    cb_counts[columns[1]] += int(columns[2])
                    continue
                cb_counts[columns[1]] += int(columns[2])

    # Filter: keep column indices whose total count meets threshold
    keep_cb = {int(cb) for cb, count in cb_counts.items() if count >= min_read}

    # Pass 2: collect data rows for kept CBs
    lines_to_keep = []
    pas_set: set[int] = set()
    for path in matrix_paths:
        with open(path, "r") as fh:
            first_data_line = True
            for line in fh:
                if not line.strip() or line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                if first_data_line:
                    first_data_line = False
                    try:
                        col_idx = int(parts[1])
                    except ValueError:
                        continue
                    if col_idx > len(_cb_lookup) and len(_cb_lookup) > 0:
                        continue  # dimension header
                    # Fall through to process as data
                try:
                    item = [int(x) for x in parts[:3]]
                except ValueError:
                    continue
                if item[1] in keep_cb:
                    lines_to_keep.append(item)
                    pas_set.add(item[0])

    matrix_pas_header = max(pas_set) if pas_set else 0
    matrix_cb_header = len(keep_cb)
    matrix_nzero_header = len(lines_to_keep)
    sorted_spars_matrix = sorted(lines_to_keep, key=lambda x: x[1])
    last_corrected_index, corrected_index = 0, 0

    with open(sorted_corrected_sparse_path, "w") as sorted_corrected_sparse_list:
        sorted_corrected_sparse_list.write(MATRIX_MARKET_HEADER)
        sorted_corrected_sparse_list.write(f"{matrix_pas_header} {matrix_cb_header} {matrix_nzero_header}\n")

        for item in sorted_spars_matrix:
            col_index = item[1]
            if item[1] != last_corrected_index:
                corrected_index += 1
                last_corrected_index = col_index
                # Use provided cb_list or the get_mapping() fallback
                filtered_cb_list.append(_cb_lookup[col_index - 1])
                item[1] = corrected_index
                sorted_corrected_sparse_list.write(f"{item[0]} {item[1]} {item[2]}\n")
            else:
                item[1] = corrected_index
                sorted_corrected_sparse_list.write(f"{item[0]} {item[1]} {item[2]}\n")

    with open(filter_cb_file, "w") as file:
        for item in filtered_cb_list:
            file.write(f"{item}\n")


def make_dataframe(matrixpath=None, collist=None):
    """Read filtered MatrixMarket file and return sparse matrix with PAS IDs preserved.

    Returns:
        tuple: (sparse_matrix, pas_ids, collist)
            - sparse_matrix: scipy.sparse.csc_matrix (rows=PAS, cols=cells)
            - pas_ids: numpy array of original PAS row indices (1-based from MatrixMarket)
            - collist: list of cell barcode strings
    """
    if matrixpath is None:
        matrixpath = directory_config.filterd_matrix
    if collist is None:
        collist = filtered_cb_list

    sparse_coo = sci.mmread(matrixpath)
    sparse_csc = sp.csc_matrix(sparse_coo)

    # Extract actual PAS IDs (row indices) from the COO matrix
    # MatrixMarket is 1-based; scipy mmread converts to 0-based
    # We need the original 1-based PAS IDs for downstream annotation
    # CRITICAL: use actual non-zero row indices, NOT sequential arange
    coo = sparse_coo if not sp.issparse(sparse_coo) else sparse_coo.tocoo()
    if sp.issparse(coo):
        pas_ids = np.sort(np.unique(coo.row + 1))  # +1: mmread converts 1-based to 0-based
    else:
        # Dense fallback: rows with any non-zero value
        pas_ids = np.sort(np.where(np.any(coo != 0, axis=1))[0] + 1)

    return sparse_csc, pas_ids, collist


def preprocessing(sparse_matrix, pas_ids, collist,
                  min_cells=None,
                  min_genes=None,
                  gene_ids=None):
    """Filter sparse PAS-by-cell matrix and return AnnData.

    Args:
        sparse_matrix: scipy.sparse matrix (rows=PAS, cols=cells)
        pas_ids: array of PAS IDs corresponding to rows
        collist: list of cell barcode strings for columns
        min_cells: minimum number of cells a PAS must appear in
        min_genes: minimum number of PAS a cell must have (default lowered for APA)
        gene_ids: Optional array (aligned to ``pas_ids``) of gene assignments
            from :func:`ema.annotate.annotate`.  When supplied, it is stored as
            ``adata.var['gene_id']`` so downstream tools (``ema switch
            length --isoform-agg per_gene``) can aggregate PAS within a gene
            without re-reading the BED/GTF.  Passing ``None`` keeps the legacy
            behaviour (no gene column).

    Returns:
        ad.AnnData: filtered AnnData object (obs=cells, var=PAS)
    """
    # Resolve config-dependent defaults at call time (see filter_cb above).
    if min_cells is None:
        min_cells = filter_config.min_cells
    if min_genes is None:
        min_genes = filter_config.min_genes

    # Build AnnData: scanpy expects cells-by-features (cells x PAS)
    # Our matrix is PAS x cells, so transpose
    var_dict = {'pas_id': pas_ids}
    if gene_ids is not None:
        if len(gene_ids) != len(pas_ids):
            raise ValueError(
                f"gene_ids length ({len(gene_ids)}) must match pas_ids "
                f"({len(pas_ids)}) — they index the same matrix rows."
            )
        var_dict['gene_id'] = gene_ids
    adata = ad.AnnData(
        X=sparse_matrix.T.tocsr(),
        obs={'barcode': collist},
        var=var_dict,
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
