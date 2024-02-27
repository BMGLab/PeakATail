from ema.countmatrix.peackcalling import peak_calling
from ema.config import directory_config, variable_config
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
from ema.clustering import clustering
from ema.annotate.annotate import annotate
from ema.annotate.find_close import find_close
from ema.annotate.gtftobed import gtf_bed
import multiprocessing as mp

def main():
    peak_calling(True, bedfilepath=directory_config.negbed, matrixpath=directory_config.negmatrixpath)
    peak_calling(False, bedfilepath=directory_config.posbed, matrixpath=directory_config.posmatrixpath)
    filter_cb()
    gtf_bed()
    genes = find_close()
    matrix_df = make_dataframe()
    annotated_frame = annotate(countmatrix=matrix_df, genes=genes)
    adata = preprocessing(df=annotated_frame)
    clustering(adata=adata)
    
if __name__ == "__main__":
    main()