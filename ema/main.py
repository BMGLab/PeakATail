import threading
import shutil
from pathlib import Path

import scipy.io as sci
import scipy.sparse as sp

from ema.countmatrix.peackcalling import peak_calling
from ema.countmatrix.indexing import get_mapping, reset_index
from ema.countmatrix.read import set_default_sample_id
from ema.config import directory_config, variable_config, args, filter_config
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

try:
    from ema.datasets.atlas_snap import snap_beds_to_atlas
except ImportError:
    snap_beds_to_atlas = None  # type: ignore[assignment]


def extract_per_dataset_mtx(
    input_mtx: Path,
    keep_col_indices: list[int],
    output_mtx: Path,
) -> None:
    """Extract a column subset from a MatrixMarket COO file.

    Reads the unified concatenated MTX (which has a proper MatrixMarket header
    written by concat_matrices), keeps only entries whose 0-based column index is
    in keep_col_indices, and re-numbers columns 1..N contiguously based on the
    order in keep_col_indices.

    Args:
        input_mtx: Path to the input MatrixMarket file (with header).
        keep_col_indices: Ordered list of 0-based column indices to retain.
            The output column for keep_col_indices[i] will be i+1 (1-based).
        output_mtx: Destination path for the filtered MatrixMarket file.
    """
    M = sci.mmread(str(input_mtx)).tocsc()  # CSC makes column slicing efficient
    sub = M[:, keep_col_indices]
    sci.mmwrite(str(output_mtx), sub.astype(int), field="integer")


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

    # Start GTF pre-processing in a background thread (with caching).
    # Runs in parallel with peak_calling so annotation data is ready before
    # find_close needs it.
    gtf_result: dict = {}
    gtf_error: dict = {}

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

    # -------------------------------------------------------------------------
    # Peak-calling loop — one iteration per (dataset_id, bam_path) pair.
    # Each iteration gets its own isolated CB column space via reset_index().
    # -------------------------------------------------------------------------
    output_dir = Path(directory_config.output_dir)
    peakcalling_dir = output_dir / "peakcalling"
    peakcalling_dir.mkdir(parents=True, exist_ok=True)

    # Track per-BAM outputs
    all_pos_beds: list[str] = []
    all_neg_beds: list[str] = []
    all_pos_mtxs: list[str] = []
    all_neg_mtxs: list[str] = []
    all_pos_cbs: list[str] = []   # cb.tsv paths (one per BAM)
    all_neg_cbs: list[str] = []   # same cb.tsv path as pos (shared index)
    all_dataset_ids_for_pos: list[str] = []
    all_dataset_ids_for_neg: list[str] = []

    # Track next BAM index per dataset so filenames are unique
    dataset_bam_indices: dict[str, int] = {}

    for dataset_id, bam_path in bam_list:
        idx = dataset_bam_indices.get(dataset_id, 0)
        dataset_bam_indices[dataset_id] = idx + 1

        reset_index()  # CRITICAL: isolate CB column space per (dataset, bam)
        set_default_sample_id(dataset_id)  # fallback when BAM has no RG tag

        pos_bed = peakcalling_dir / f"{dataset_id}_{idx}.pos.bed"
        neg_bed = peakcalling_dir / f"{dataset_id}_{idx}.neg.bed"
        pos_mtx = peakcalling_dir / f"{dataset_id}_{idx}.pos.mtx"
        neg_mtx = peakcalling_dir / f"{dataset_id}_{idx}.neg.mtx"
        cb_tsv = peakcalling_dir / f"{dataset_id}_{idx}.cb.tsv"

        peak_calling(
            False,
            bedfilepath=str(pos_bed),
            matrixpath=str(pos_mtx),
            bamfile_dir=str(bam_path),
            **peak_kwargs,
        )
        peak_calling(
            True,
            bedfilepath=str(neg_bed),
            matrixpath=str(neg_mtx),
            bamfile_dir=str(bam_path),
            **peak_kwargs,
        )

        # Dump CB list — shared by pos+neg (both used the same _index instance)
        mapping = get_mapping()
        ordered_cbs = [cb for cb, _ in sorted(mapping.items(), key=lambda x: x[1])]
        with open(cb_tsv, "w") as f:
            for cb in ordered_cbs:
                f.write(cb + "\n")

        all_pos_beds.append(str(pos_bed))
        all_neg_beds.append(str(neg_bed))
        all_pos_mtxs.append(str(pos_mtx))
        all_neg_mtxs.append(str(neg_mtx))
        # Both pos and neg MTX share the SAME cb.tsv (same barcode index)
        all_pos_cbs.append(str(cb_tsv))
        all_neg_cbs.append(str(cb_tsv))
        all_dataset_ids_for_pos.append(dataset_id)
        all_dataset_ids_for_neg.append(dataset_id)

    # Save peak calling stats
    output_mgr.save_stats("peak_calling", {
        "strategy": args.strategy,
        "dynamic_threshold": args.dynamic_threshold,
        "floor_threshold": args.floor_threshold,
        "lambda_fold_change": args.lambda_fold_change,
        "lambda_window": args.lambda_window,
    })

    # =========================================================================
    # Single-sample path (len(bam_list) == 1)
    # Copy outputs to legacy paths and continue with the existing pipeline.
    # =========================================================================
    if len(bam_list) == 1:
        shutil.copy(all_pos_beds[0], directory_config.posbed)
        shutil.copy(all_neg_beds[0], directory_config.negbed)
        shutil.copy(all_pos_mtxs[0], directory_config.posmatrixpath)
        shutil.copy(all_neg_mtxs[0], directory_config.negmatrixpath)

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

        # Find closest gene for each PAS
        genes = find_close(
            utr_lengths=utr_lengths,
            max_distance=getattr(args, "max_gene_distance", 5000),
            utr_multiplier=getattr(args, "utr_multiplier", 2.0),
            include_extended=getattr(args, "include_extended", False),
        )

        # Save PAS-gene assignment stats
        output_mgr.save_stats("pas_gene", {
            "max_gene_distance": getattr(args, "max_gene_distance", 5000),
            "utr_multiplier": getattr(args, "utr_multiplier", 2.0),
            "assigned_pas_count": len(genes),
        })

        # Build sparse matrix (PAS IDs preserved)
        sparse_matrix, pas_ids, collist = make_dataframe()

        # Annotate: join gene assignments with count matrix
        result = annotate(
            sparse_matrix=sparse_matrix,
            pas_ids=pas_ids,
            collist=collist,
            genes=genes,
        )

        # Save annotation stats
        output_mgr.save_stats("annotated", {
            "input_pas_count": len(pas_ids),
            "annotated_pas_count": len(result.pas_ids),
            "cell_count": len(result.collist),
        })

        # Preprocess and cluster
        adata = preprocessing(
            sparse_matrix=result.sparse_matrix,
            pas_ids=result.pas_ids,
            collist=collist,
        )

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

        return  # done with single-sample path

    # =========================================================================
    # Multi-sample path (len(bam_list) > 1)
    # =========================================================================
    if merge_pas_beds is None or concat_matrices is None:
        raise RuntimeError("ema.datasets.pas_merge is not available for multi-sample mode")

    unified_dir = output_dir / "unified"
    unified_dir.mkdir(parents=True, exist_ok=True)

    # Combine pos and neg lists for BED/MTX unification
    all_beds = all_pos_beds + all_neg_beds
    all_mtxs = all_pos_mtxs + all_neg_mtxs
    all_cbs = all_pos_cbs + all_neg_cbs
    all_ds_ids = all_dataset_ids_for_pos + all_dataset_ids_for_neg

    # Dispatch atlas vs. coordinate-merge based on config
    if directory_config.atlas:
        if snap_beds_to_atlas is None:
            raise RuntimeError("ema.datasets.atlas_snap is not available")
        unified_bed, mapping_path = snap_beds_to_atlas(
            all_beds,
            all_ds_ids,
            atlas_bed=directory_config.atlas,
            output_dir=unified_dir,
            distance=directory_config.atlas_distance,
        )
    else:
        unified_bed, mapping_path = merge_pas_beds(
            all_beds,
            all_ds_ids,
            output_dir=unified_dir,
            gap=getattr(args, "pas_gap", 100),
        )

    unified_mtx = unified_dir / "concatenated.mtx"
    unified_cb = unified_dir / "concatenated_cbs.tsv"
    concat_matrices(
        mtx_paths=all_mtxs,
        dataset_ids=all_ds_ids,
        cb_paths=all_cbs,
        mapping_path=mapping_path,
        output_mtx=unified_mtx,
        output_cb=unified_cb,
    )

    print(
        f"[multi-sample] Unified {len(bam_list)} BAMs "
        f"({'atlas' if directory_config.atlas else 'merge'}) "
        f"-> {unified_bed}"
    )

    # Read the concatenated CB list once — used for per-dataset column selection
    with open(unified_cb) as f:
        all_cb_strings = [line.strip() for line in f if line.strip()]

    # Wait for GTF processing before per-dataset annotation
    gtf_thread.join()
    if "exception" in gtf_error:
        raise RuntimeError(
            f"GTF processing failed: {gtf_error['exception']}"
        ) from gtf_error["exception"]
    utr_lengths = gtf_result.get("utr_lengths", {})

    # Run find_close ONCE on the unified PAS coordinate set.
    # find_close expects pos+neg BED inputs (it cats them) — split unified BED by strand.
    shutil.copy(str(unified_bed), directory_config.pasbed)
    with open(unified_bed) as _u, \
         open(directory_config.posbed, "w") as _p, \
         open(directory_config.negbed, "w") as _n:
        for line in _u:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6 and parts[5] == "+":
                _p.write(line)
            elif len(parts) >= 6 and parts[5] == "-":
                _n.write(line)
    genes = find_close(
        utr_lengths=utr_lengths,
        max_distance=getattr(args, "max_gene_distance", 5000),
        utr_multiplier=getattr(args, "utr_multiplier", 2.0),
        include_extended=getattr(args, "include_extended", False),
    )

    output_mgr.save_stats("gtf_annotation", {
        "gtf_path": directory_config.gtf_dir,
        "utr_count": len(utr_lengths),
    })
    output_mgr.save_stats("pas_gene", {
        "max_gene_distance": getattr(args, "max_gene_distance", 5000),
        "utr_multiplier": getattr(args, "utr_multiplier", 2.0),
        "assigned_pas_count": len(genes),
    })

    # Per-dataset clustering (each dataset_id clusters independently)
    unique_ds_ids: list[str] = list(dict.fromkeys(all_ds_ids))  # dedupe, preserve order
    per_dataset_dir = output_dir / "per_dataset"
    per_dataset_dir.mkdir(parents=True, exist_ok=True)

    for ds_id in unique_ds_ids:
        ds_dir = per_dataset_dir / ds_id
        ds_dir.mkdir(parents=True, exist_ok=True)

        # Find which columns in the concatenated MTX belong to this dataset.
        # CB strings are prefixed "{dataset_id}_{cb_seq}" by read.py.
        sub_indices = [
            i for i, cb in enumerate(all_cb_strings)
            if cb.startswith(f"{ds_id}_")
        ]
        sub_cbs = [all_cb_strings[i] for i in sub_indices]

        if not sub_indices:
            print(f"[multi-sample] WARNING: no cells found for dataset '{ds_id}' — skipping")
            continue

        # Extract per-dataset sub-matrix from the unified concatenated MTX
        pre_filter_mtx = ds_dir / "pre_filter.mtx"
        extract_per_dataset_mtx(
            input_mtx=unified_mtx,
            keep_col_indices=sub_indices,
            output_mtx=pre_filter_mtx,
        )

        filtered_mtx = ds_dir / "filtered_matrix.mtx"
        filtered_cb_path = ds_dir / "filtered_cb.tsv"
        filter_cb(
            input_matrix_paths=[str(pre_filter_mtx)],
            cb_list=sub_cbs,
            sorted_corrected_sparse_path=str(filtered_mtx),
            filter_cb_file=str(filtered_cb_path),
            min_read=filter_config.min_read,
        )

        sparse_matrix, pas_ids, _ = make_dataframe(matrixpath=str(filtered_mtx))

        # Read filtered CBs from the file written by filter_cb
        with open(filtered_cb_path) as f:
            collist = [line.strip() for line in f if line.strip()]

        result = annotate(
            sparse_matrix=sparse_matrix,
            pas_ids=pas_ids,
            collist=collist,
            genes=genes,
        )

        adata = preprocessing(
            sparse_matrix=result.sparse_matrix,
            pas_ids=result.pas_ids,
            collist=collist,
        )

        # Save per-dataset cluster output
        cluster_h5ad = ds_dir / "clusters.h5ad"
        clustering(adata=adata, output_h5ad=str(cluster_h5ad))

        # Save per-dataset stats as a JSON next to the h5ad (OutputManager only
        # has predefined slots for the single-sample pipeline).
        import json
        with open(ds_dir / "clustering_stats.json", "w") as f:
            json.dump({
                "dataset_id": ds_id,
                "final_cells": int(adata.n_obs),
                "final_pas": int(adata.n_vars),
            }, f, indent=2)

        print(f"[multi-sample] Dataset '{ds_id}': {adata.n_obs} cells, {adata.n_vars} PAS -> {cluster_h5ad}")

    return  # done with multi-sample path


if __name__ == "__main__":
    main()
