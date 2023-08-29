from .countmatrix.peackcalling import peak_calling
from .config import directory_config, variable_config
from .matrixfilter import filter_cb, make_dataframe, preprocessing
from .clustering import clustering
import multiprocessing

def main():
    peak_calling(True, bedfilepath=directory_config.negbed, matrixpath=directory_config.negmatrixpath)
    peak_calling(False, bedfilepath=directory_config.posbed, matrixpath=directory_config.posmatrixpath)
    filter_cb()
    matrix_df = make_dataframe()
    adata = preprocessing(df=matrix_df)
    clustering(adata=adata)
    
if __name__ == "__main__":
    main()