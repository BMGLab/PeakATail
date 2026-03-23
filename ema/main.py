import threading

from ema.countmatrix.peackcalling import peak_calling
from ema.config import directory_config, variable_config, args
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
from ema.clustering.clustering import clustering
from ema.annotate.annotate import annotate
from ema.annotate.find_close import find_close
from ema.annotate.gtf_cache import process_gtf_cached
from ema.strategies import get_strategy
from ema.output import OutputManager


def main():
    # Set up output directory structure
    output_mgr = OutputManager(base_dir=directory_config.output_dir)
    output_mgr.setup()

    # Save run configuration
    output_mgr.save_run_config(vars(args))

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
    gtf_error = {}

    def _gtf_worker():
        try:
            gtf_result["utr_lengths"] = process_gtf_cached(
                gtf_path=directory_config.gtf_dir,
                output_dir=directory_config.output_dir,
                endbed_path=directory_config.endbed,
                features_path=directory_config.raw_features,
            )
        except Exception as e:
            gtf_error["exception"] = e

    gtf_thread = threading.Thread(target=_gtf_worker, name="gtf-preprocess")
    gtf_thread.start()

    # Peak calling (runs in parallel with GTF processing)
    peak_calling(True, bedfilepath=directory_config.negbed,
                 matrixpath=directory_config.negmatrixpath, **peak_kwargs)
    peak_calling(False, bedfilepath=directory_config.posbed,
                 matrixpath=directory_config.posmatrixpath, **peak_kwargs)

    # Save peak calling stats
    output_mgr.save_stats("peak_calling", {
        "strategy": args.strategy,
        "dynamic_threshold": args.dynamic_threshold,
        "floor_threshold": args.floor_threshold,
        "lambda_fold_change": args.lambda_fold_change,
        "lambda_window": args.lambda_window,
    })

    filter_cb()

    # Save CB filter stats
    output_mgr.save_stats("cb_filter", {
        "min_read": args.min_pas_per_cell,
    })

    # Wait for GTF processing to complete before find_close
    gtf_thread.join()
    if "exception" in gtf_error:
        raise RuntimeError(
            f"GTF processing failed: {gtf_error['exception']}"
        ) from gtf_error["exception"]
    utr_lengths = gtf_result.get("utr_lengths", {})

    # Save GTF annotation stats
    output_mgr.save_stats("gtf_annotation", {
        "gtf_path": directory_config.gtf_dir,
        "utr_count": len(utr_lengths),
    })

    # Find closest gene for each PAS with adaptive distance thresholding
    genes = find_close(
        utr_lengths=utr_lengths,
        max_distance=getattr(args, 'max_gene_distance', 5000),
        utr_multiplier=getattr(args, 'utr_multiplier', 2.0),
        include_extended=getattr(args, 'include_extended', False),
    )

    # Save PAS-gene assignment stats
    output_mgr.save_stats("pas_gene", {
        "max_gene_distance": getattr(args, 'max_gene_distance', 5000),
        "utr_multiplier": getattr(args, 'utr_multiplier', 2.0),
        "assigned_pas_count": len(genes),
    })

    # Build sparse matrix (PAS IDs preserved)
    sparse_matrix, pas_ids, collist = make_dataframe()

    # Annotate: join gene assignments with count matrix
    result = annotate(sparse_matrix=sparse_matrix, pas_ids=pas_ids,
                      collist=collist, genes=genes)

    # Save annotation stats
    output_mgr.save_stats("annotated", {
        "input_pas_count": len(pas_ids),
        "annotated_pas_count": len(result.pas_ids),
        "cell_count": len(result.collist),
    })

    # Preprocess and cluster
    adata = preprocessing(sparse_matrix=result.sparse_matrix,
                          pas_ids=result.pas_ids,
                          collist=collist)

    # Save preprocessing stats
    output_mgr.save_stats("preprocessing", {
        "cells_after_filter": adata.n_obs,
        "pas_after_filter": adata.n_vars,
    })

    clustering(adata=adata)

    # Save clustering stats
    output_mgr.save_stats("clustering", {
        "final_cells": adata.n_obs,
        "final_pas": adata.n_vars,
    })


if __name__ == "__main__":
    main()
