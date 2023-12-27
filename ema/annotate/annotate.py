from ema.config import directory_config, variable_config
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

def annotate(countmatrix:pd.DataFrame(), 
            genes:pd.DataFrame(), 
            annotated_matrix= directory_config.filterd_matrix,
            pasbed_dir = directory_config.pasbed,
            ) -> pd.DataFrame():

    # set row names as index number
    countmatrix.index = range(1, len(countmatrix) + 1)
    print(countmatrix)
    #get genes has annotated and has not been filter matrix filter module
    annotated_frame = pd.concat([genes, countmatrix], axis=1, join="outer")
    ann_gene = annotated_frame.iloc[:, 0]
    ann_gene.to_csv(directory_config.pas_geneid, index=True, header=False, sep="\t")
    # remove gene_id column 
    annotated_frame = annotated_frame.drop(annotated_frame.columns[0], axis=1)
    #save new sparse matrix
    sparse_matrix = sp.csc_matrix(annotated_frame.values)
    mmwrite(target=annotated_matrix, a=sparse_matrix, comment=variable_config.matrixmarketheader)
    #save new pasbed file
    pas_df = pd.read_csv(pasbed_dir, sep="\t", header=None)
    filtered_pas = pas_df[pas_df[3].isin(annotated_frame.index)]
    filtered_pas.to_csv(pasbed_dir, sep="\t", header=False, index=False)
    return annotated_frame

if __name__ == "__main__":
    annotate()