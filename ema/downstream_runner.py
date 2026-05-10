"""Per-dataset downstream pipeline runner for parallel execution.

This module provides top-level, pickle-safe worker functions that run the
complete per-dataset downstream pipeline (extract → filter → make_dataframe →
annotate → preprocessing → clustering) for one dataset in isolation.

Designed for use with ``multiprocessing.get_context("spawn").Pool`` so that:
- No shared mutable state leaks between workers.
- Each worker carries all inputs as plain Python values (no closures / lambdas).
- The worker is importable at the top level (required for pickle under spawn).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import scipy.io as sci


def extract_per_dataset_mtx(
    input_mtx: Path,
    keep_col_indices: list[int],
    output_mtx: Path,
) -> None:
    """Extract a column subset from a MatrixMarket COO file.

    Reads the unified concatenated MTX (which has a proper MatrixMarket header
    written by ``concat_matrices``), keeps only entries whose 0-based column
    index is in ``keep_col_indices``, and re-numbers columns 1..N
    contiguously based on the order in ``keep_col_indices``.

    Args:
        input_mtx: Path to the input MatrixMarket file (with header).
        keep_col_indices: Ordered list of 0-based column indices to retain.
            The output column for ``keep_col_indices[i]`` will be ``i+1``
            (1-based).
        output_mtx: Destination path for the filtered MatrixMarket file.
    """
    M = sci.mmread(str(input_mtx)).tocsc()  # CSC makes column slicing efficient
    sub = M[:, keep_col_indices]
    sci.mmwrite(str(output_mtx), sub.astype(int), field="integer")


def run_one_dataset_downstream(
    ds_id: str,
    sub_indices: list[int],
    sub_cbs: list[str],
    unified_mtx: str,
    per_dataset_dir: str,
    genes_pkl: bytes,
    min_read: int,
    filter_min_cells: int,
    filter_min_genes: int,
    *,
    log_queue=None,
    progress_client=None,
) -> dict[str, Any]:
    """Run the downstream pipeline for a single dataset.

    Executes ``extract_per_dataset_mtx`` → ``filter_cb`` → ``make_dataframe``
    → ``annotate`` → ``preprocessing`` → ``clustering`` for *ds_id* in isolation.
    All inputs are plain serialisable Python values so the function is safe to
    pickle across process boundaries under the ``spawn`` start method.

    Args:
        ds_id: Dataset identifier string (e.g. ``"sampleA"``).
        sub_indices: Ordered list of 0-based column indices in the unified MTX
            that belong to this dataset.
        sub_cbs: Cell-barcode strings corresponding to ``sub_indices``
            (same length, same order).
        unified_mtx: Absolute path to the unified concatenated MatrixMarket
            file produced by ``concat_matrices``.
        per_dataset_dir: Absolute path to the parent directory where a
            ``{ds_id}/`` subdirectory will be created for all outputs.
        genes_pkl: ``pickle.dumps()`` of the ``genes`` DataFrame returned by
            ``find_close``.  Serialised so it can be passed through the process
            boundary without relying on shared memory.
        min_read: Minimum total read count for a barcode to be retained by
            ``filter_cb``.
        filter_min_cells: Minimum number of cells a PAS must appear in
            (passed to ``preprocessing``).
        filter_min_genes: Minimum number of PAS a cell must have
            (passed to ``preprocessing``).

    Returns:
        A stats dictionary with keys ``dataset_id``, ``final_cells``,
        ``final_pas``.  Also writes ``{per_dataset_dir}/{ds_id}/clusters.h5ad``
        and ``{per_dataset_dir}/{ds_id}/clustering_stats.json`` on disk.

    Raises:
        RuntimeError: If no sub_indices were supplied (caller should skip
            before dispatching).
    """
    import pickle

    if log_queue is not None:
        from ema.logging_config import setup_worker_logging
        setup_worker_logging(log_queue)
    log = logging.getLogger(__name__)

    prefix = f"[ds={ds_id}]"

    # Guard BEFORE heavy imports so tests can verify the error without
    # triggering ema.config initialisation via ema.matrixfilter.
    if not sub_indices:
        raise RuntimeError(f"{prefix} no sub_indices supplied — caller must skip empty datasets")

    # Late imports so the heavy stack is only loaded in the worker process.
    # NOTE: extract_per_dataset_mtx is defined in this same module (no ema.main import).
    from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
    from ema.annotate.annotate import annotate
    from ema.clustering.clustering import clustering

    genes = pickle.loads(genes_pkl)

    ds_dir = Path(per_dataset_dir) / ds_id
    ds_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Extract per-dataset sub-matrix from unified MTX                  #
    # ------------------------------------------------------------------ #
    pre_filter_mtx = ds_dir / "pre_filter.mtx"
    log.info("%s extracting %d columns from unified MTX", prefix, len(sub_indices))
    extract_per_dataset_mtx(
        input_mtx=Path(unified_mtx),
        keep_col_indices=sub_indices,
        output_mtx=pre_filter_mtx,
    )

    # ------------------------------------------------------------------ #
    # 2. Filter barcodes by minimum read count                            #
    # ------------------------------------------------------------------ #
    filtered_mtx = ds_dir / "filtered_matrix.mtx"
    filtered_cb_path = ds_dir / "filtered_cb.tsv"
    log.info("%s filtering barcodes (min_read=%d)", prefix, min_read)
    filter_cb(
        input_matrix_paths=[str(pre_filter_mtx)],
        cb_list=sub_cbs,
        sorted_corrected_sparse_path=str(filtered_mtx),
        filter_cb_file=str(filtered_cb_path),
        min_read=min_read,
    )

    # ------------------------------------------------------------------ #
    # 3. Build sparse matrix                                              #
    # ------------------------------------------------------------------ #
    log.info("%s building sparse matrix", prefix)
    sparse_matrix, pas_ids, _ = make_dataframe(matrixpath=str(filtered_mtx))

    with open(filtered_cb_path) as fh:
        collist = [line.strip() for line in fh if line.strip()]

    # ------------------------------------------------------------------ #
    # 4. Annotate PAS with gene assignments                               #
    # ------------------------------------------------------------------ #
    log.info("%s annotating PAS", prefix)
    result = annotate(
        sparse_matrix=sparse_matrix,
        pas_ids=pas_ids,
        collist=collist,
        genes=genes,
    )

    # ------------------------------------------------------------------ #
    # 5. Preprocess (filter cells/PAS)                                    #
    # ------------------------------------------------------------------ #
    log.info("%s preprocessing", prefix)
    adata = preprocessing(
        sparse_matrix=result.sparse_matrix,
        pas_ids=result.pas_ids,
        collist=collist,
        min_cells=filter_min_cells,
        min_genes=filter_min_genes,
    )

    # ------------------------------------------------------------------ #
    # 6. Cluster                                                          #
    # ------------------------------------------------------------------ #
    cluster_h5ad = ds_dir / "clusters.h5ad"
    log.info("%s clustering (%d cells x %d PAS)", prefix, adata.n_obs, adata.n_vars)
    clustering(adata=adata, output_h5ad=str(cluster_h5ad))

    # ------------------------------------------------------------------ #
    # 7. Persist stats                                                    #
    # ------------------------------------------------------------------ #
    stats: dict[str, Any] = {
        "dataset_id": ds_id,
        "final_cells": int(adata.n_obs),
        "final_pas": int(adata.n_vars),
    }
    with open(ds_dir / "clustering_stats.json", "w") as fh:
        json.dump(stats, fh, indent=2)

    log.info("%s done — %d cells, %d PAS -> %s", prefix, adata.n_obs, adata.n_vars, cluster_h5ad)
    return stats


def downstream_worker_star(args: tuple) -> dict:
    """Top-level pickle-safe shim for ``pool.imap_unordered``.

    ``imap_unordered`` requires a callable that accepts a *single* argument.
    This shim unpacks the pre-built argument tuple and delegates to
    ``run_one_dataset_downstream``.  Must be defined at module level (not as a
    lambda or nested function) so it is picklable under the ``spawn`` context.

    Args:
        args: Tuple of positional arguments for ``run_one_dataset_downstream``,
            optionally followed by a ``log_queue`` value as the 10th element
            (or as a ``{"log_queue": ...}`` dict appended to the tuple).

    Returns:
        The stats dict returned by ``run_one_dataset_downstream``.
    """
    # Support optional log_queue appended as last element of the tuple.
    if len(args) == 10:
        *pos_args, log_queue = args
        return run_one_dataset_downstream(*pos_args, log_queue=log_queue)
    return run_one_dataset_downstream(*args)
