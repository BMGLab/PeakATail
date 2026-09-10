from ema.config import directory_config, variable_config
from ema.outputs import MATRIX_MARKET_HEADER
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite


def annotate(sparse_matrix, pas_ids, collist, genes,
             annotated_matrix=None,
             pasbed_dir=None,
             ):
    """Annotate PAS with gene assignments and filter the count matrix.

    Joins the sparse PAS-by-cell count matrix with gene annotations from
    find_close, preserving PAS IDs throughout. Writes the filtered sparse
    matrix and updated PAS BED file.

    Args:
        sparse_matrix: scipy.sparse matrix (rows=PAS, cols=cells) from make_dataframe.
        pas_ids: numpy array of PAS IDs (1-based, from MatrixMarket).
        collist: list of cell barcode strings.
        genes: DataFrame from find_close with PAS IDs as index, gene_id column.
        annotated_matrix: Path to write the filtered MatrixMarket file.
        pasbed_dir: Path to the PAS BED file (read and rewritten filtered).

    Returns:
        AnnotatedResult: object with sparse_matrix, pas_ids, collist, gene_ids.
    """
    # B2: resolve the output paths at CALL time, not at function-definition
    # time. The old signature bound ``directory_config.annotated_matrix`` /
    # ``.pasbed`` as default argument values, which captured the pre-
    # ``set_directory_config`` ``emaout/`` defaults at import — so a run wrote
    # stray shadow files (emaout/pasbed.bed, emaout/annotated_matrix.mtx) instead
    # of into the timestamped run dir. Sentinels + call-time resolution fix this.
    if annotated_matrix is None:
        annotated_matrix = directory_config.annotated_matrix
    if pasbed_dir is None:
        pasbed_dir = directory_config.pasbed

    # Build a DataFrame index from PAS IDs -- DO NOT overwrite with range(1,N)
    # The PAS IDs from make_dataframe() are the actual MatrixMarket row indices

    # Rows are resolved by POSITION in ``pas_ids`` (see pas_id_to_row below),
    # so the matrix MUST be row-aligned with pas_ids.  make_dataframe()
    # guarantees this; a full-height matrix paired with a compacted pas_ids
    # list silently mis-keys every PAS below the first empty row, so fail loud.
    if sparse_matrix.shape[0] != len(pas_ids):
        raise ValueError(
            f"sparse_matrix has {sparse_matrix.shape[0]} rows but pas_ids has "
            f"{len(pas_ids)} entries — rows are looked up by position in "
            "pas_ids, so a mismatch mis-keys the counts. Pass the matrix "
            "make_dataframe() returned (it is subset to pas_ids)."
        )

    # Find which PAS IDs have gene annotations
    assert genes.index.is_unique, "find_close produced duplicate PAS IDs — gene_id alignment will break"
    annotated_pas_ids = genes.index.intersection(pas_ids)

    if len(annotated_pas_ids) == 0:
        raise ValueError("No PAS IDs matched between count matrix and gene annotations. "
                         "Check that find_close and make_dataframe used consistent PAS IDs.")

    # Build mask: which rows in sparse_matrix correspond to annotated PAS.
    # keep_rows and keep_pas_ids are built from the SAME iteration order over
    # ``annotated_pas_ids``, so ``keep_rows[k]`` is by construction the row of
    # ``keep_pas_ids[k]`` -- never re-derive one from the other's position.
    pas_id_to_row = {pid: i for i, pid in enumerate(pas_ids)}
    keep_pas_ids = np.asarray(list(annotated_pas_ids))
    keep_rows = np.asarray([pas_id_to_row[pid] for pid in keep_pas_ids], dtype=np.intp)
    assert len(keep_rows) == len(keep_pas_ids), (
        "keep_rows/keep_pas_ids desynchronised — the PAS→row mapping is broken"
    )

    # Subset the sparse matrix to only annotated PAS
    filtered_sparse = sparse_matrix.tocsr()[keep_rows, :]
    assert filtered_sparse.shape[0] == len(keep_pas_ids), (
        f"annotated matrix has {filtered_sparse.shape[0]} rows for "
        f"{len(keep_pas_ids)} PAS IDs — row/PAS-ID alignment lost"
    )

    # Get gene IDs aligned with kept PAS
    gene_ids = genes.loc[keep_pas_ids, "gene_id"]

    # Save PAS-to-gene mapping
    pas_gene_df = pd.DataFrame({"gene_id": gene_ids.values}, index=keep_pas_ids)
    pas_gene_df.to_csv(directory_config.pas_geneid, header=False, sep="\t")

    # Save filtered sparse matrix
    mmwrite(target=annotated_matrix, a=sp.csc_matrix(filtered_sparse),
            comment=MATRIX_MARKET_HEADER)

    # Filter PAS BED file to keep only annotated PAS
    pas_df = pd.read_csv(pasbed_dir, sep="\t", header=None)
    filtered_pas = pas_df[pas_df[3].isin(keep_pas_ids)]
    filtered_pas.to_csv(pasbed_dir, sep="\t", header=False, index=False)

    # Return structured result for downstream preprocessing
    return AnnotatedResult(
        sparse_matrix=sp.csc_matrix(filtered_sparse),
        pas_ids=keep_pas_ids,
        collist=collist,
        gene_ids=gene_ids.values,
    )


class AnnotatedResult:
    """Container for annotated matrix results passed to preprocessing."""

    __slots__ = ("sparse_matrix", "pas_ids", "collist", "gene_ids")

    def __init__(self, sparse_matrix, pas_ids, collist, gene_ids):
        self.sparse_matrix = sparse_matrix
        self.pas_ids = pas_ids
        self.collist = collist
        self.gene_ids = gene_ids


if __name__ == "__main__":
    annotate()
