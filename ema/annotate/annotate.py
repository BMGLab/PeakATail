from ema.config import directory_config, variable_config
from ema.outputs import MATRIX_MARKET_HEADER
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scipy.io import mmwrite


def annotate(sparse_matrix, pas_ids, collist, genes,
             annotated_matrix=directory_config.annotated_matrix,
             pasbed_dir=directory_config.pasbed,
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
    # Build a DataFrame index from PAS IDs -- DO NOT overwrite with range(1,N)
    # The PAS IDs from make_dataframe() are the actual MatrixMarket row indices

    # Find which PAS IDs have gene annotations
    assert genes.index.is_unique, "find_close produced duplicate PAS IDs — gene_id alignment will break"
    annotated_pas_ids = genes.index.intersection(pas_ids)

    if len(annotated_pas_ids) == 0:
        raise ValueError("No PAS IDs matched between count matrix and gene annotations. "
                         "Check that find_close and make_dataframe used consistent PAS IDs.")

    # Build mask: which rows in sparse_matrix correspond to annotated PAS
    pas_id_to_row = {pid: i for i, pid in enumerate(pas_ids)}
    keep_rows = np.array([pas_id_to_row[pid] for pid in annotated_pas_ids
                          if pid in pas_id_to_row])
    keep_pas_ids = np.array([pid for pid in annotated_pas_ids
                             if pid in pas_id_to_row])

    # Subset the sparse matrix to only annotated PAS
    filtered_sparse = sparse_matrix.tocsr()[keep_rows, :]

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
