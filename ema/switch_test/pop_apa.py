"""Load filtered count matrix and build a PAS×Cell DataFrame with gene annotation.

Handles sparse PAS IDs — the matrix may have more rows than the gene mapping
because some PAS IDs were filtered out during annotation. Only PAS with gene
assignments are included in the output.
"""

from scipy.io import mmread
import numpy as np
import pandas as pd
from ema.config import directory_config


def pop_apa(spars_matrix_dir=directory_config.filterd_matrix,
            cellbarcodes_dir=directory_config.filtered_cb,
            genes_dir=directory_config.pas_geneid) -> pd.DataFrame:
    """Load sparse matrix and gene mapping, return annotated count matrix.

    The matrix rows are PAS IDs (1-based, potentially sparse — not every
    row has data). The gene mapping file has PAS_ID → gene_id pairs for
    only the annotated PAS. This function joins them correctly regardless
    of dimension mismatches.

    Args:
        spars_matrix_dir: Path to filtered MatrixMarket matrix.
        cellbarcodes_dir: Path to filtered cell barcodes TSV.
        genes_dir: Path to PAS → gene mapping TSV (pas_id, gene_id).

    Returns:
        DataFrame with MultiIndex (Ensemble_ID, PAS) and cell barcode columns,
        filtered to genes with >1 PAS (required for APA analysis).
    """
    # Read sparse matrix
    sparse_mat = mmread(spars_matrix_dir).tocsr()
    cellbarcodes = pd.read_csv(cellbarcodes_dir, sep="\t", header=None)[0].tolist()
    genes = pd.read_csv(genes_dir, sep="\t", header=None, names=["pas_id", "gene_id"])

    n_matrix_rows = sparse_mat.shape[0]
    n_matrix_cols = sparse_mat.shape[1]

    # The matrix has shape (max_pas_id, n_cells)
    # genes has pas_id values that index into the matrix (1-based)
    # Only keep PAS IDs that are within matrix range
    genes = genes[genes["pas_id"] <= n_matrix_rows].copy()
    genes = genes[genes["pas_id"] >= 1].copy()

    # Extract only the rows from the matrix that have gene annotations
    pas_ids = genes["pas_id"].values
    gene_ids = genes["gene_id"].values

    # Build count matrix for annotated PAS only (avoids dense conversion of full matrix)
    rows_to_extract = pas_ids - 1  # 0-based indexing
    sub_matrix = sparse_mat[rows_to_extract, :].toarray()

    # Build DataFrame
    multi_index = pd.MultiIndex.from_arrays(
        [gene_ids, pas_ids],
        names=["Ensemble_ID", "PAS"]
    )

    count_matrix = pd.DataFrame(
        sub_matrix,
        index=multi_index,
        columns=cellbarcodes[:n_matrix_cols]
    )

    # Keep only genes with >1 PAS (APA requires alternative sites)
    filtered_matrix = count_matrix.groupby("Ensemble_ID").filter(lambda x: len(x) > 1)

    return filtered_matrix
