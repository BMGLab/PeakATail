from ema.annotate.annotate import annotate
from ema.matrixfilter import filter_cb, make_dataframe
from test_find_close import find_close
from ema.config import directory_config
import pandas as pd
def test_annotate():
    genes = find_close(posbed_dir=directory_config.posbed, 
                negbe_dir=directory_config.negbed,
                genomebed_dir=directory_config.endbed,
                mergebed=directory_config.pasbed,
                annotatedbed_dir=directory_config.annotatedbed)

    filter_cb() 
    
    count_matrix = make_dataframe()

    annotatedframe:pd.DataFrame() = annotate(countmatrix=count_matrix, genes= genes)
    annotatedframe.to_csv("~/D/BAMdata/proje/frame.csv", index=True, header=True, sep =",")

    print(annotatedframe)
 

test_annotate()