from collections import defaultdict
from ema.countmatrix.indexing import get_mapping
from ema.config import filter_config, directory_config, variable_config
from ema.outputs import MATRIX_MARKET_HEADER
import numpy as np
import pandas as pd
import scipy.io as sci
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

filtered_cb_list = []


class _RaggedMatrixData(Exception):
    """Internal signal: a matrix file's tokenization is ambiguous for the
    vectorized fast-path reader (e.g. a data row has more than 3
    whitespace-separated tokens, or a required column fails int parsing).

    Callers catch this and fall back to :func:`_filter_cb_legacy`, which is
    a byte-for-byte copy of the original row-by-row implementation — so
    correctness is guaranteed on any input, and the vectorized path is a
    pure speed optimization for the well-formed (real-world) case.
    """


def _read_matrix_file_vectorized(path: str, n_cb: int) -> "pd.DataFrame":
    """Vectorized read of one MatrixMarket-ish matrix file.

    Mirrors, exactly, the per-line semantics of the original ``filter_cb``
    passes for a single file:

      * blank lines and ``%``-comment lines are skipped (``comment="%"``,
        ``skip_blank_lines=True``);
      * a data row with fewer than 3 whitespace tokens is dropped
        (``len(columns) < 3: continue``);
      * a data row with MORE than 3 whitespace tokens is ambiguous for a
        fixed-width vectorized parse (the original only ever reads the
        first 3 tokens and ignores the rest) — raises :class:`_RaggedMatrixData`
        so the caller falls back to the reference implementation;
      * the FIRST surviving row of the file is the MatrixMarket dimension
        header IFF ``int(row[1]) > n_cb`` (and ``n_cb > 0``); if
        ``int(row[1])`` isn't parseable at all the row is also dropped
        (matches ``first_data_line`` handling in both original passes —
        either way the header slot is "consumed" exactly once per file);
      * any row whose pas/cb/count token doesn't parse as an int raises
        :class:`_RaggedMatrixData` (pass 1 of the original crashes
        uncaught on a bad ``columns[2]``; deferring to the reference
        implementation reproduces that crash exactly instead of trying to
        reconcile pass-1-crashes vs. pass-2-silently-skips in vectorized
        code).

    Returns:
        DataFrame with int64 columns ``pas``, ``cb``, ``count`` and a
        ``cb_str`` column carrying the ORIGINAL (pre-int-cast) barcode-index
        token — Pass 1 of the original code accumulates counts keyed by
        that raw string, not by its int value, so this is preserved to stay
        byte-identical on pathological inputs (e.g. a leading-zero token).
    """
    try:
        df = pd.read_csv(
            path, sep=r"\s+", comment="%", header=None, dtype=str,
            engine="c", skip_blank_lines=True, na_filter=False,
        )
    except pd.errors.ParserError as exc:
        raise _RaggedMatrixData(f"{path}: ragged whitespace tokenization") from exc
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()

    empty = pd.DataFrame({"pas": pd.Series(dtype="int64"),
                           "cb_str": pd.Series(dtype=object),
                           "cb": pd.Series(dtype="int64"),
                           "count": pd.Series(dtype="int64")})
    if df.shape[0] == 0:
        return empty
    if df.shape[1] < 3:
        # every surviving row had <3 tokens -> all dropped
        return empty
    if df.shape[1] > 3:
        raise _RaggedMatrixData(f"{path}: rows with >3 whitespace tokens")

    valid = df[2] != ""  # rows with <3 tokens were padded with "" in col 2
    df = df[valid]
    if df.empty:
        return empty

    # -- header detection: mirrors `first_data_line` handling exactly --
    first_pos = df.index[0]
    first_col1 = df.at[first_pos, 1]
    drop_first = False
    try:
        col_idx = int(first_col1)
        if col_idx > n_cb and n_cb > 0:
            drop_first = True
    except ValueError:
        drop_first = True
    if drop_first:
        df = df.drop(index=first_pos)
    if df.empty:
        return empty

    pas_num = pd.to_numeric(df[0], errors="coerce")
    cb_num = pd.to_numeric(df[1], errors="coerce")
    count_num = pd.to_numeric(df[2], errors="coerce")
    if pas_num.isna().any() or cb_num.isna().any() or count_num.isna().any():
        raise _RaggedMatrixData(f"{path}: non-integer token in a data row")

    return pd.DataFrame({
        "pas": pas_num.to_numpy(dtype="int64"),
        "cb_str": df[1].to_numpy(),
        "cb": cb_num.to_numpy(dtype="int64"),
        "count": count_num.to_numpy(dtype="int64"),
    })


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

    # Vectorized fast path: each file is read ONCE via pandas' C parser and
    # reused for both the count-summation pass and the row-collection pass
    # (the original re-read every file once per pass). Falls back to the
    # byte-for-byte-identical row-by-row implementation whenever a file's
    # tokenization is ambiguous for a fixed-width vectorized parse — see
    # _read_matrix_file_vectorized / _RaggedMatrixData. That guarantees
    # correctness on any input while making the common (well-formed,
    # multi-million-row) case fast.
    try:
        per_file = [_read_matrix_file_vectorized(p, len(_cb_lookup)) for p in matrix_paths]
    except _RaggedMatrixData:
        _filter_cb_legacy(
            matrix_paths=matrix_paths,
            cb_lookup=_cb_lookup,
            min_read=min_read,
            sorted_corrected_sparse_path=sorted_corrected_sparse_path,
            filter_cb_file=filter_cb_file,
        )
        return

    if per_file:
        combined = pd.concat(per_file, ignore_index=True)
    else:
        combined = pd.DataFrame({
            "pas": pd.Series(dtype="int64"), "cb_str": pd.Series(dtype=object),
            "cb": pd.Series(dtype="int64"), "count": pd.Series(dtype="int64"),
        })

    # Pass 1 (vectorized): sum counts per RAW cb token across ALL files,
    # grouped on the STRING token — exactly like the original dict keyed by
    # `columns[1]` (not its int value), so a pathological leading-zero token
    # ("007" vs "7") reproduces the original's split-bucket behaviour.
    if combined.empty:
        keep_cb: set[int] = set()
    else:
        sums = combined.groupby("cb_str")["count"].sum()
        keep_cb = {int(cb) for cb, total in sums.items() if total >= min_read}

    # Pass 2 (vectorized): keep rows whose (int) cb is in keep_cb.
    kept = combined[combined["cb"].isin(keep_cb)]

    matrix_pas_header = int(kept["pas"].max()) if len(kept) else 0
    matrix_cb_header = len(keep_cb)
    matrix_nzero_header = len(kept)

    # Stable sort by ORIGINAL cb — ties keep the original file/row traversal
    # order, matching Python's stable `sorted(lines_to_keep, key=lambda x: x[1])`.
    kept_sorted = kept.sort_values(by="cb", kind="stable")
    cb_sorted = kept_sorted["cb"].to_numpy()

    # Contiguous 1..K re-index: a new group starts wherever cb differs from
    # the previous (already cb-sorted) row — identical to the original's
    # adjacent-duplicate walk over the cb-sorted rows.
    if len(cb_sorted):
        first_of_group = np.concatenate(([True], cb_sorted[1:] != cb_sorted[:-1]))
        new_cb = np.cumsum(first_of_group)
    else:
        first_of_group = np.empty(0, dtype=bool)
        new_cb = np.empty(0, dtype=np.int64)

    with open(sorted_corrected_sparse_path, "w") as sorted_corrected_sparse_list:
        sorted_corrected_sparse_list.write(MATRIX_MARKET_HEADER)
        sorted_corrected_sparse_list.write(f"{matrix_pas_header} {matrix_cb_header} {matrix_nzero_header}\n")
        if len(kept_sorted):
            out_df = pd.DataFrame({
                0: kept_sorted["pas"].to_numpy(),
                1: new_cb,
                2: kept_sorted["count"].to_numpy(),
            })
            out_df.to_csv(
                sorted_corrected_sparse_list, sep=" ", header=False,
                index=False, lineterminator="\n",
            )

    # filtered_cb_list: cb_lookup value for each newly-encountered ORIGINAL
    # cb in ascending order == sorted-unique original cbs present in kept rows.
    if len(cb_sorted):
        unique_sorted_cbs = cb_sorted[first_of_group]
        filtered_cb_list.extend(_cb_lookup[int(c) - 1] for c in unique_sorted_cbs)

    with open(filter_cb_file, "w") as file:
        if filtered_cb_list:
            file.write("\n".join(str(item) for item in filtered_cb_list) + "\n")


def _filter_cb_legacy(*, matrix_paths, cb_lookup, min_read,
                       sorted_corrected_sparse_path, filter_cb_file) -> None:
    """Byte-for-byte copy of the original row-by-row ``filter_cb`` body.

    Used as the correctness fallback whenever a matrix file's tokenization
    is too ambiguous for the vectorized fast path (see _RaggedMatrixData) —
    pathological/malformed input that should essentially never occur on
    real MatrixMarket output, but which the original implementation
    tolerated (or, for a bad ``columns[2]``, crashed on) on a line-by-line
    basis. Mutates the module-level ``filtered_cb_list`` exactly like the
    original function did.
    """
    global filtered_cb_list
    _cb_lookup = cb_lookup

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
                if first_data_line:
                    first_data_line = False
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
