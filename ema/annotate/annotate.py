import pandas as pd

def annotate(countmatrix:pd.DataFrame(), genes:pd.DataFrame()) -> pd.DataFrame():

    # set row names as index number
    countmatrix.set_index(pd.Index(1, len(countmatrix)+1))

    annotated_frame = pd.concat([genes, countmatrix], axis=1, join="left")
    annotated_frame = annotated_frame.set_index(annotated_frame.columns[0])
    annotated_frame.index.name = "gene"

    return annotated_frame

if __name__ == "__main__":
    annotate()