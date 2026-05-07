import threading
import shutil

from ema.countmatrix.peackcalling import peak_calling
from ema.config import directory_config, variable_config, args
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
from ema.clustering.clustering import clustering
from ema.annotate.annotate import annotate
from ema.annotate.find_close import find_close
from ema.annotate.gtf_cache import process_gtf_cached
from ema.strategies import get_strategy
from ema.output import OutputManager

try:
    from ema.datasets.manager import DatasetManager
    from ema.datasets.pas_merge import merge_pas_beds, concat_matrices
except ImportError:
    DatasetManager = None  # type: ignore[assignment,misc]
    merge_pas_beds = None  # type: ignore[assignment]
    concat_matrices = None  # type: ignore[assignment]


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

    # Prepare BAM list via DatasetManager (multi-sample support)
    if DatasetManager is not None and directory_config.datasets:
        dm = DatasetManager(output_dir=directory_config.output_dir)
        bam_list = dm.prepare(directory_config.datasets)
    else:
        # Fallback: single BAM from legacy --bamDir
        bam_list = [("default", directory_config.bam_dir)]

    # Collect per-dataset peak-calling outputs
    pos_beds: list[str] = []
    neg_beds: list[str] = []
    pos_matrices: list[str] = []
    neg_matrices: list[str] = []

    for dataset_id, bam_path in bam_list:
        bam_str = str(bam_path)
        pos_bed = bam_str.replace('.bam', f'_{dataset_id}_pos.bed')
        neg_bed = bam_str.replace('.bam', f'_{dataset_id}_neg.bed')
        pos_mtx = bam_str.replace('.bam', f'_{dataset_id}_pos.mtx')
        neg_mtx = bam_str.replace('.bam', f'_{dataset_id}_neg.mtx')

        peak_calling(False, bedfilepath=pos_bed, matrixpath=pos_mtx,
                     bamfile_dir=bam_str, **peak_kwargs)
        peak_calling(True, bedfilepath=neg_bed, matrixpath=neg_mtx,
                     bamfile_dir=bam_str, **peak_kwargs)

        pos_beds.append(pos_bed)
        neg_beds.append(neg_bed)
        pos_matrices.append(pos_mtx)
        neg_matrices.append(neg_mtx)

    # Route outputs to the expected downstream paths
    if len(bam_list) == 1:
        shutil.copy(pos_beds[0], directory_config.posbed)
        shutil.copy(neg_beds[0], directory_config.negbed)
        shutil.copy(pos_matrices[0], directory_config.posmatrixpath)
        shutil.copy(neg_matrices[0], directory_config.negmatrixpath)
    else:
        # Multi-sample path: merge PAS coordinates then concatenate matrices
        if merge_pas_beds is None:
            raise RuntimeError("ema.datasets.pas_merge not available")

        # Determine which dataset uses merge_after strategy
        after_datasets = {
            ds['id'] for ds in directory_config.datasets
            if ds.get('merge_strategy') == 'after'
        }

        # Collect all pos+neg BED files for merging
        all_beds = pos_beds + neg_beds
        merged_bed, mapping_path = merge_pas_beds(
            all_beds,
            output_dir=directory_config.output_dir,
            gap=getattr(args, 'pas_gap', 100),
        )
        shutil.copy(merged_bed, directory_config.pasbed)

        # Concatenate matrices — requires per-dataset CB lists and PAS ID files
        # These are written by peak_calling alongside each MTX
        all_mtx = pos_matrices + neg_matrices
        cb_lists = [p.replace('.mtx', '.cb.tsv') for p in all_mtx]
        merged_mtx, merged_cb = concat_matrices(
            mtx_paths=all_mtx,
            cb_list_paths=cb_lists,
            mapping_path=mapping_path,
            output_path=directory_config.matrixpath,
            output_cb_path=directory_config.filtered_cb,
        )
        print(f"[multi-sample] Merged {len(bam_list)} datasets → {merged_bed} ({merged_mtx})")

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
