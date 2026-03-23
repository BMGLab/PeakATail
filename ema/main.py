import threading

from ema.countmatrix.peackcalling import peak_calling
from ema.config import directory_config, variable_config, args
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
from ema.clustering.clustering import clustering
from ema.annotate.annotate import annotate
from ema.annotate.find_close import find_close
from ema.annotate.gtf_cache import process_gtf_cached
from ema.strategies import get_strategy


def main():
    strategy = get_strategy(args.strategy)
    peak_kwargs = dict(
        strategy=strategy,
        dynamic_threshold=args.dynamic_threshold,
        floor_threshold=args.floor_threshold,
        lambda_fold_change=args.lambda_fold_change,
        lambda_window=args.lambda_window,
    )

    # Start GTF pre-processing in a background thread (with caching)
    # This runs in parallel with peak_calling so annotation data is ready
    # before find_close needs it.
    gtf_result = {}

    def _gtf_worker():
        gtf_result["utr_lengths"] = process_gtf_cached(
            gtf_path=directory_config.gtf_dir,
            output_dir=directory_config.output_dir,
            endbed_path=directory_config.endbed,
            features_path=directory_config.raw_features,
        )

    gtf_thread = threading.Thread(target=_gtf_worker, name="gtf-preprocess")
    gtf_thread.start()

    # Peak calling (runs in parallel with GTF processing)
    peak_calling(True, bedfilepath=directory_config.negbed,
                 matrixpath=directory_config.negmatrixpath, **peak_kwargs)
    peak_calling(False, bedfilepath=directory_config.posbed,
                 matrixpath=directory_config.posmatrixpath, **peak_kwargs)
    filter_cb()

    # Wait for GTF processing to complete before find_close
    gtf_thread.join()
    utr_lengths = gtf_result.get("utr_lengths", {})

    # Find closest gene for each PAS with adaptive distance thresholding
    genes = find_close(
        utr_lengths=utr_lengths,
        max_distance=getattr(args, 'max_gene_distance', 5000),
        utr_multiplier=getattr(args, 'utr_multiplier', 2.0),
        include_extended=getattr(args, 'include_extended', False),
    )

    # Build sparse matrix (PAS IDs preserved)
    sparse_matrix, pas_ids, collist = make_dataframe()

    # Annotate: join gene assignments with count matrix
    result = annotate(sparse_matrix=sparse_matrix, pas_ids=pas_ids,
                      collist=collist, genes=genes)

    # Preprocess and cluster
    adata = preprocessing(sparse_matrix=result.sparse_matrix,
                          pas_ids=result.pas_ids,
                          collist=collist)
    clustering(adata=adata)


if __name__ == "__main__":
    main()
