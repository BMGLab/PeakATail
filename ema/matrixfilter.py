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


class _NonCanonicalBarcodeToken(_RaggedMatrixData):
    """Internal signal: a barcode-index token is not the canonical decimal
    form of its integer value (``"007"``, ``"+7"``, ``"7.0"``).

    The original implementation groups Pass-1 counts by the raw STRING
    token, so ``"007"`` and ``"7"`` are two buckets that both filter rows
    with ``cb == 7``.  The streaming reader groups by the integer -- which
    is exactly equivalent while every token is canonical -- so any token
    that is not defers the whole call to :func:`_filter_cb_legacy`, the
    byte-for-byte copy of the reference algorithm.
    """


#: Rows parsed per chunk.  4M rows x 3 int64 columns ~ 100 MB of working set.
_CHUNK_ROWS = 4_000_000

#: Refuse to build a per-barcode accumulator larger than this (a corrupt
#: file with an absurd column index would otherwise allocate GBs); such a
#: file goes to the reference implementation instead.
_MAX_CB_INDEX = 100_000_000


def _chunk_to_int(series):
    """Integer array for one parsed column, or :class:`_RaggedMatrixData`.

    ``errors="raise"`` reproduces the "non-integer token in a data row"
    rejection of the previous whole-file reader; a float-valued token (e.g.
    ``"3.5"``) makes pandas return a float column, which the reference
    implementation skips row-by-row rather than truncating -- so that also
    defers to :func:`_filter_cb_legacy`.
    """
    try:
        vals = pd.to_numeric(series, errors="raise")
    except (ValueError, TypeError) as exc:
        raise _RaggedMatrixData(f"non-integer token in a data row: {exc}") from exc
    if vals.dtype.kind not in "iu":
        raise _RaggedMatrixData("non-integer (float or out-of-range) token in a data row")
    return vals.to_numpy(dtype=np.int64)


def _downcast(arr: "np.ndarray") -> "np.ndarray":
    """int32 view of *arr* when every value fits — halves the row cost of
    the kept-row buffers (180M rows on the PBMC 10k run) and prints
    identically."""
    if arr.size and (arr.min() < np.iinfo(np.int32).min or arr.max() > np.iinfo(np.int32).max):
        return arr
    return arr.astype(np.int32, copy=False)


def _read_matrices_streaming(paths: list, n_cb: int):
    """Parse every matrix file ONCE, in chunks, into three integer arrays.

    Mirrors, exactly, the per-line semantics of the original ``filter_cb``
    passes (see :class:`_RaggedMatrixData` for the deferral cases):

      * blank lines and ``%``-comment lines are skipped;
      * a data row with fewer than 3 whitespace tokens is dropped;
      * a data row with MORE than 3 whitespace tokens is ambiguous for a
        fixed-width vectorized parse -> :class:`_RaggedMatrixData`;
      * the FIRST surviving row of each file is the MatrixMarket dimension
        header iff ``int(row[1]) > n_cb`` (with ``n_cb > 0``), or if that
        token is not an integer at all;
      * every barcode token must be the canonical decimal form of its value
        -> :class:`_NonCanonicalBarcodeToken`.

    Keeping only the integers (int32 where they fit) instead of a pandas
    frame with an object ``cb_str`` column per row is the memory fix: the
    previous reader held ~95 B per non-zero, i.e. an estimated 20-27 GB on
    the PBMC 10k matrices (200M non-zeros), which was the run's peak RSS
    once the dense TF-IDF was gone.

    Returns:
        ``(pas, cb, count)`` integer arrays, concatenated across *paths*.
    """
    pas_parts: list = []
    cb_parts: list = []
    cnt_parts: list = []
    seen_tokens: dict[str, int] = {}

    for path in paths:
        first_row_pending = True
        try:
            reader = pd.read_csv(
                path, sep=r"\s+", comment="%", header=None, dtype=str,
                engine="c", skip_blank_lines=True, na_filter=False,
                chunksize=_CHUNK_ROWS,
            )
        except pd.errors.EmptyDataError:
            continue
        try:
            for chunk in reader:
                if chunk.shape[1] > 3:
                    raise _RaggedMatrixData(f"{path}: rows with >3 whitespace tokens")
                if chunk.shape[1] < 3:
                    continue  # every row had <3 tokens -> all dropped
                chunk = chunk[chunk[2] != ""]  # rows with <3 tokens (padded)
                if chunk.empty:
                    continue
                if first_row_pending:
                    first_row_pending = False
                    first_pos = chunk.index[0]
                    try:
                        col_idx = int(chunk.at[first_pos, 1])
                        drop_first = col_idx > n_cb and n_cb > 0
                    except ValueError:
                        drop_first = True
                    if drop_first:
                        chunk = chunk.drop(index=first_pos)
                        if chunk.empty:
                            continue
                cb_tokens = chunk[1]
                for tok in pd.unique(cb_tokens.to_numpy()):
                    if tok in seen_tokens:
                        continue
                    try:
                        val = int(tok)
                    except ValueError as exc:
                        raise _RaggedMatrixData(
                            f"{path}: non-integer barcode token {tok!r}"
                        ) from exc
                    if str(val) != tok:
                        raise _NonCanonicalBarcodeToken(
                            f"{path}: non-canonical barcode token {tok!r}"
                        )
                    seen_tokens[tok] = val
                pas_parts.append(_downcast(_chunk_to_int(chunk[0])))
                cb_parts.append(_downcast(_chunk_to_int(cb_tokens)))
                cnt_parts.append(_downcast(_chunk_to_int(chunk[2])))
        except pd.errors.ParserError as exc:
            raise _RaggedMatrixData(f"{path}: ragged whitespace tokenization") from exc
        except pd.errors.EmptyDataError:
            continue

    if not pas_parts:
        empty = np.empty(0, dtype=np.int32)
        return empty, empty.copy(), empty.copy()
    if len(pas_parts) == 1:
        return pas_parts[0], cb_parts[0], cnt_parts[0]
    # Concatenate one column at a time and drop each part list as we go so
    # only one column is ever duplicated.
    out = []
    for parts in (pas_parts, cb_parts, cnt_parts):
        out.append(np.concatenate(parts))
        parts.clear()
    return out[0], out[1], out[2]


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

    # Vectorized fast path: every file is parsed ONCE, in chunks, into three
    # integer arrays (see _read_matrices_streaming) and both passes run on
    # those arrays.  Falls back to the byte-for-byte-identical row-by-row
    # implementation whenever a file's tokenization is ambiguous for a
    # fixed-width vectorized parse or a barcode token is not canonical — see
    # _RaggedMatrixData / _NonCanonicalBarcodeToken.  That guarantees
    # correctness on any input while making the common (well-formed,
    # multi-million-row) case fast AND small: the previous reader kept a
    # pandas frame with a per-row object ``cb_str`` column alive for the whole
    # call (~95 B per non-zero, an estimated 20-27 GB on the PBMC 10k
    # matrices); the integer arrays cost 12 B per non-zero.
    try:
        pas, cb, count = _read_matrices_streaming(matrix_paths, len(_cb_lookup))
    except _RaggedMatrixData:
        _filter_cb_legacy(
            matrix_paths=matrix_paths,
            cb_lookup=_cb_lookup,
            min_read=min_read,
            sorted_corrected_sparse_path=sorted_corrected_sparse_path,
            filter_cb_file=filter_cb_file,
        )
        return

    # Pass 1: total counts per barcode column index across ALL files.  The
    # original groups by the raw token; _read_matrices_streaming guarantees
    # every token is the canonical form of its integer, so grouping by the
    # integer is the same partition.
    n_rows = pas.size
    cb_max = int(cb.max()) if n_rows else 0
    if cb_max > _MAX_CB_INDEX or (n_rows and int(cb.min()) < 0):
        # Absurd column index: an accumulator that size is not worth
        # allocating — hand the file to the reference implementation.
        _filter_cb_legacy(
            matrix_paths=matrix_paths,
            cb_lookup=_cb_lookup,
            min_read=min_read,
            sorted_corrected_sparse_path=sorted_corrected_sparse_path,
            filter_cb_file=filter_cb_file,
        )
        return

    if n_rows:
        cb_i = cb.astype(np.intp, copy=False)
        # float64 weights are exact here: a per-barcode total is a sum of
        # read counts, far below 2**53.
        totals = np.bincount(cb_i, weights=count, minlength=cb_max + 1)
        occurrences = np.bincount(cb_i, minlength=cb_max + 1)
        keep_flags = (occurrences > 0) & (totals >= min_read)
        del cb_i, totals, occurrences
    else:
        keep_flags = np.zeros(1, dtype=bool)

    matrix_cb_header = int(keep_flags.sum())

    # Pass 2: keep the rows whose barcode survived.
    if n_rows:
        row_mask = keep_flags[cb]
        pas, cb, count = pas[row_mask], cb[row_mask], count[row_mask]
        del row_mask

    matrix_pas_header = int(pas.max()) if pas.size else 0
    matrix_nzero_header = int(pas.size)

    # Stable sort by ORIGINAL cb — ties keep the original file/row traversal
    # order, matching Python's stable `sorted(lines_to_keep, key=lambda x: x[1])`.
    if pas.size:
        order = np.argsort(cb, kind="stable")
        pas, cb, count = pas[order], cb[order], count[order]
        del order
        # Contiguous 1..K re-index: a new group starts wherever cb differs
        # from the previous (already cb-sorted) row — identical to the
        # original's adjacent-duplicate walk over the cb-sorted rows.
        first_of_group = np.empty(cb.size, dtype=bool)
        first_of_group[0] = True
        np.not_equal(cb[1:], cb[:-1], out=first_of_group[1:])
        new_cb = np.cumsum(first_of_group, dtype=np.int64)
    else:
        first_of_group = np.empty(0, dtype=bool)
        new_cb = np.empty(0, dtype=np.int64)

    with open(sorted_corrected_sparse_path, "w") as sorted_corrected_sparse_list:
        sorted_corrected_sparse_list.write(MATRIX_MARKET_HEADER)
        sorted_corrected_sparse_list.write(f"{matrix_pas_header} {matrix_cb_header} {matrix_nzero_header}\n")
        for lo in range(0, int(pas.size), _CHUNK_ROWS):
            hi = min(lo + _CHUNK_ROWS, int(pas.size))
            pd.DataFrame({0: pas[lo:hi], 1: new_cb[lo:hi], 2: count[lo:hi]}).to_csv(
                sorted_corrected_sparse_list, sep=" ", header=False,
                index=False, lineterminator="\n",
            )

    # filtered_cb_list: cb_lookup value for each newly-encountered ORIGINAL
    # cb in ascending order == sorted-unique original cbs present in kept rows.
    if cb.size:
        filtered_cb_list.extend(_cb_lookup[int(c) - 1] for c in cb[first_of_group])

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

    The returned matrix is ROW-ALIGNED with ``pas_ids``: row ``i`` holds the
    counts of PAS ``pas_ids[i]``.  Downstream (:func:`ema.annotate.annotate`)
    resolves a PAS to a matrix row by its POSITION in ``pas_ids``, so this
    invariant is load-bearing — see the guard in ``annotate()``.

    Returns:
        tuple: (sparse_matrix, pas_ids, collist)
            - sparse_matrix: scipy.sparse.csc_matrix (rows=PAS, cols=cells),
              ``shape[0] == len(pas_ids)``
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

    # RE-KEY THE MATRIX TO pas_ids.
    #
    # mmread returns a FULL-height matrix: rows 1..<the row count in the
    # MatrixMarket dimension header>, which is ``max(pas_id)`` over the kept
    # rows.  ``pas_ids`` above lists only the rows that actually carry counts.
    # Every PAS whose reads all belonged to barcodes dropped by the min_read
    # CB filter leaves an EMPTY row behind, so the two indexings drift apart by
    # the number of empty rows above each PAS (observed: 0 -> 1540 on a real
    # mouse testis run).  ``annotate()`` looks rows up by POSITION in
    # ``pas_ids``, so without this subset every PAS below the first gap would
    # be handed a different PAS's counts.
    if len(pas_ids):
        sparse_csc = sp.csc_matrix(sparse_csc.tocsr()[pas_ids - 1, :])
    else:
        sparse_csc = sp.csc_matrix((0, sparse_csc.shape[1]), dtype=sparse_csc.dtype)

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
            ``adata.var['gene_id']`` so downstream tools (``peakatail switch
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
    # Our matrix is PAS x cells, so transpose.
    # ``pas_ids`` becomes var_names positionally, so row i of sparse_matrix must
    # be PAS pas_ids[i]; a height mismatch would silently relabel every PAS
    # (AnnData's own error names shapes, not the cause), so check it here.
    if sparse_matrix.shape[0] != len(pas_ids):
        raise ValueError(
            f"sparse_matrix has {sparse_matrix.shape[0]} rows but pas_ids has "
            f"{len(pas_ids)} entries — pas_ids labels the rows positionally, "
            "so a mismatch mis-keys every PAS."
        )
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
