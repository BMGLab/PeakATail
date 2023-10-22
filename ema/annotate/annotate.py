from ema.config import directory_config
import pandas as pd

def annotate(countmatrix:pd.DataFrame(), genes:pd.DataFrame()) -> pd.DataFrame():

    # set row names as index number
    countmatrix.index = range(1, len(countmatrix) + 1)
    print(countmatrix)
    #get genes has annotated and has not been filter matrix filter module
    annotated_frame = pd.concat([genes, countmatrix], axis=1, join="inner")
    ann_gene = annotated_frame.iloc[:, 0]
    ann_gene.to_csv(directory_config.pas_geneid, index=True, header=False, sep="\t")
    # remove gene_id column 
    annotated_frame = annotated_frame.drop(annotated_frame.columns[0], axis=1)

    return annotated_frame

if __name__ == "__main__":
    annotate()