from ema.countmatrix.peackcalling import peak_calling
from ema.config import directory_config, variable_config, args
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
from ema.clustering.clustering import clustering
from ema.annotate.annotate import annotate
from ema.annotate.find_close import find_close
from ema.annotate.gtftobed import gtf_bed
from ema.strategies import get_strategy
import multiprocessing as mp

def main():
    strategy = get_strategy(args.strategy)
    peak_kwargs = dict(
        strategy=strategy,
        dynamic_threshold=args.dynamic_threshold,
        floor_threshold=args.floor_threshold,
        lambda_fold_change=args.lambda_fold_change,
        lambda_window=args.lambda_window,
    )
    peak_calling(True, bedfilepath=directory_config.negbed, matrixpath=directory_config.negmatrixpath, **peak_kwargs)
    peak_calling(False, bedfilepath=directory_config.posbed, matrixpath=directory_config.posmatrixpath, **peak_kwargs)
    filter_cb()
    gtf_bed()
    genes = find_close()
    sparse_matrix, pas_ids, collist = make_dataframe()
    annotated_frame = annotate(sparse_matrix=sparse_matrix, pas_ids=pas_ids,
                               collist=collist, genes=genes)
    adata = preprocessing(sparse_matrix=annotated_frame.sparse_matrix,
                          pas_ids=annotated_frame.pas_ids,
                          collist=collist)
    clustering(adata=adata)
    
if __name__ == "__main__":
    main()