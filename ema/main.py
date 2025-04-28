from ema.config import directory_config, variable_config
from ema.countmatrix.peakcalling import peak_calling
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
from ema.clustering.clustering import clustering
from ema.annotate.annotate import annotate
from ema.annotate.find_close import find_close
from ema.annotate.gtftobed import gtf_bed
from ema.annotate.gtf2utr_parser import gtf_3utr_bed
import multiprocessing as mp

def main():
    gtf_3utr_bed(gtffile=directory_config.gtf_dir,bedfile=directory_config.utrbed)
    peak_calling(False, bedfilepath=directory_config.posbed, matrixpath=directory_config.posmatrixpath)
    peak_calling(True, bedfilepath=directory_config.negbed, matrixpath=directory_config.negmatrixpath)
    filter_cb()
    gtf_bed()
    genes = find_close()
    matrix_df = make_dataframe()
    annotated_frame = annotate(countmatrix=matrix_df, genes=genes)
    #adata = preprocessing(df=annotated_frame)
    #clustering(adata=adata)
    
if __name__ == "__main__":
    main()