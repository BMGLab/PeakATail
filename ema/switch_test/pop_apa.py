from scipy.io import mmread
import pandas as pd
from ema.config import directory_config
import pandas as pd
def pop_apa(spars_matrix_dir=directory_config.filterd_matrix,
        cellbarcodes_dir=directory_config.filtered_cb,
        genes_dir=directory_config.pas_geneid
        )-> pd.DataFrame:
    """
    This function takes in the directory paths for sparse matrix, cell barcodes, and genes,
    and returns a filtered count matrix based on the given conditions.

    Parameters:
    - spars_matrix_dir (str): Directory path for the sparse matrix file.
    - cellbarcodes_dir (str): Directory path for the cell barcodes file.
    - genes_dir (str): Directory path for the genes file.

    Returns:
    - filtered_matrix (pd.DataFrame): Alternativ PAS .
    """
    spars_matrix = mmread(spars_matrix_dir)
    cellbarcodes = pd.read_csv(cellbarcodes_dir, sep="\t", header=None)
    genes = pd.read_csv(genes_dir, sep="\t", header=None)

    count_matrix = pd.DataFrame(spars_matrix.toarray(), index=genes, columns=cellbarcodes[0])

    count_matrix.index = pd.MultiIndex.from_tuples(count_matrix.index.map(lambda x: (x[1], x[0])),
                            names=['Ensemble_ID', "PAS"])

    filtered_matrix = count_matrix.groupby("Ensemble_ID").filter(lambda x: len(x) > 1)
    return filtered_matrix
    
