from ema.countmatrix.peackcalling import peak_calling
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
from ema.clustering import clustering
from ema.config import directory_config
import multiprocessing

def main():
    peak_calling(True, bedfile=directory_config.negbed(), matrix=directory_config.negmatrixpath())
    peak_calling(False, bedfile=directory_config.posbed(), matrix=directory_config.posmatrixpath())
    filter_cb()
    matrix_df = make_dataframe()
    adata = preprocessing(df=matrix_df)
    clustering(adata=adata)
    