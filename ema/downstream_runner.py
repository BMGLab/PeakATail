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
    output_dir: str,
    genes_pkl: bytes,
    min_read: int,
    filter_min_cells: int,
    filter_min_genes: int,
    *,
    log_queue=None,
    progress_client=None,
    cluster_method: str = "leiden_tfidf",
    cluster_resolution: float = 1.0,
    cluster_n_pcs: int = 40,
    cluster_random_seed: int = 42,
    cluster_external_clusters: str | None = None,
    plot_engines: list[str] | None = None,
    cluster_n_neighbors: int | None = None,
    cluster_tfidf_scale_factor: float = 1e4,
    cluster_depth_corr_threshold: float = 0.75,
    cluster_n_svd_components: int = 50,
    cluster_n_top_hvg: int = 2000,
    atlas_of: dict | None = None,
    ip_of: dict | None = None,
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
        output_dir: Absolute path to the pipeline run root.  The worker
            initialises ``directory_config.output_dir`` with this value so all
            ``directory_config.*_for(ds_id)`` accessors resolve correctly in
            the spawned process.
        genes_pkl: ``pickle.dumps()`` of the ``genes`` DataFrame returned by
            ``find_close``.  Serialised so it can be passed through the process
            boundary without relying on shared memory.
        min_read: Minimum total read count for a barcode to be retained by
            ``filter_cb``.
        filter_min_cells: Minimum number of cells a PAS must appear in
            (passed to ``preprocessing``).
        filter_min_genes: Minimum number of PAS a cell must have
            (passed to ``preprocessing``).
        atlas_of: (D9) Optional ``{pas_id: (atlas_match, atlas_distance_bp)}``
            map, forwarded to ``write_pas_gene_artifacts`` /
            ``record_pas_drops`` to populate the atlas-snap status columns.
            ``{}`` (default) when atlas snapping is disabled.
        ip_of: (D9) Optional ``{pas_id: internal_priming_bool}`` map, same
            plumbing as ``atlas_of``. ``{}`` (default) when the
            internal-priming filter is disabled.

    Returns:
        A stats dictionary with keys ``dataset_id``, ``final_cells``,
        ``final_pas``.  Also writes ``07_clustering/{ds_id}/clusters.h5ad``
        and a ``clustering_stats.json`` alongside it on disk.

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

    # Initialise directory_config in this (possibly spawned) process so all
    # directory_config.*_for() accessors resolve to the correct run root.
    from ema.config import set_directory_config
    set_directory_config(output_dir=Path(output_dir))

    # Late imports so the heavy stack is only loaded in the worker process.
    # NOTE: extract_per_dataset_mtx is defined in this same module (no ema.main import).
    from ema.config import directory_config
    from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
    from ema.annotate.annotate import annotate
    from ema.clustering.clustering import clustering

    genes = pickle.loads(genes_pkl)

    # Scratch dir: use the clustering stage dir for this dataset as the working
    # area (pre_filter.mtx and filtered_matrix.mtx are transient).
    ds_scratch = directory_config.clustering_dir / ds_id
    ds_scratch.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # 1. Extract per-dataset sub-matrix from unified MTX                  #
    # ------------------------------------------------------------------ #
    pre_filter_mtx = ds_scratch / "pre_filter.mtx"
    log.info("%s extracting %d columns from unified MTX", prefix, len(sub_indices))
    extract_per_dataset_mtx(
        input_mtx=Path(unified_mtx),
        keep_col_indices=sub_indices,
        output_mtx=pre_filter_mtx,
    )
    if progress_client is not None:
        progress_client.advance(1)  # tick 1/6: extract

    # ------------------------------------------------------------------ #
    # 2. Filter barcodes by minimum read count                            #
    # ------------------------------------------------------------------ #
    filtered_mtx = ds_scratch / "filtered_matrix.mtx"
    filtered_cb_path = ds_scratch / "filtered_cb.tsv"
    log.info("%s filtering barcodes (min_read=%d)", prefix, min_read)
    filter_cb(
        input_matrix_paths=[str(pre_filter_mtx)],
        cb_list=sub_cbs,
        sorted_corrected_sparse_path=str(filtered_mtx),
        filter_cb_file=str(filtered_cb_path),
        min_read=min_read,
    )

    # Persist the kept barcode list for this dataset (canonical name +
    # min_read header for traceability).  See ema/outputs.py.
    # B3: ALWAYS write it — never gate on .exists(), which silently left
    # 02_cb_filter/<ds>/ with only the stats JSON.
    from ema.outputs import write_filtered_cb
    if filtered_cb_path.exists():
        _kept_cbs = [
            b.strip()
            for b in filtered_cb_path.read_text().splitlines()
            if b.strip()
        ]
    else:
        _kept_cbs = []
        log.warning(
            "%s cb_filter produced no filtered_cb.tsv at %s; writing an empty "
            "per-dataset filtered_cb.tsv (header only) for traceability.",
            prefix, filtered_cb_path,
        )
    write_filtered_cb(Path(output_dir), ds_id, _kept_cbs, min_read)

    if progress_client is not None:
        progress_client.advance(1)  # tick 2/6: filter

    # ------------------------------------------------------------------ #
    # 3. Build sparse matrix                                              #
    # ------------------------------------------------------------------ #
    log.info("%s building sparse matrix", prefix)
    sparse_matrix, pas_ids, _ = make_dataframe(matrixpath=str(filtered_mtx))

    with open(filtered_cb_path) as fh:
        collist = [line.strip() for line in fh if line.strip()]
    if progress_client is not None:
        progress_client.advance(1)  # tick 3/6: build_sparse

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

    # Persist canonical PAS->gene mapping + annotatedpas.bed + annotated
    # count matrix for this dataset.  See ema/outputs.py for the layout.
    from ema.outputs import write_pas_gene_artifacts, write_annotated_matrix
    write_pas_gene_artifacts(
        Path(output_dir), ds_id,
        result.pas_ids, result.gene_ids,
        atlas_of=atlas_of, ip_of=ip_of,
    )
    write_annotated_matrix(
        Path(output_dir), ds_id,
        result.sparse_matrix, result.pas_ids, result.collist,
    )

    if progress_client is not None:
        progress_client.advance(1)  # tick 4/6: annotate

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
        gene_ids=result.gene_ids,
    )

    # Persist the post-filter AnnData snapshot (pre-clustering).
    from ema.outputs import write_preprocessed_h5ad
    write_preprocessed_h5ad(Path(output_dir), ds_id, adata)

    if progress_client is not None:
        progress_client.advance(1)  # tick 5/6: preprocess

    # ------------------------------------------------------------------ #
    # 6. Cluster                                                          #
    # ------------------------------------------------------------------ #
    cluster_h5ad = directory_config.clusters_h5ad_for(ds_id)
    cluster_h5ad.parent.mkdir(parents=True, exist_ok=True)
    log.info(
        "%s clustering (%d cells x %d PAS) — method=%s resolution=%s",
        prefix, adata.n_obs, adata.n_vars,
        cluster_method, cluster_resolution,
    )
    clustering(
        adata=adata,
        method=cluster_method,
        resolution=cluster_resolution,
        n_pcs=cluster_n_pcs,
        random_seed=cluster_random_seed,
        external_clusters=cluster_external_clusters,
        output_h5ad=str(cluster_h5ad),
        n_neighbors=cluster_n_neighbors,
        tfidf_scale_factor=cluster_tfidf_scale_factor,
        depth_corr_threshold=cluster_depth_corr_threshold,
        n_svd_components=cluster_n_svd_components,
        n_top_hvg=cluster_n_top_hvg,
    )
    if progress_client is not None:
        progress_client.advance(1)  # tick 6/6: cluster

    # ------------------------------------------------------------------ #
    # 7. Persist stats                                                    #
    # All visualizations happen post-pipeline in render_run_outputs()     #
    # (see ema/viz/pipeline_hooks.py) so workers stay focused on data.    #
    # ------------------------------------------------------------------ #
    stats: dict[str, Any] = {
        "dataset_id": ds_id,
        "final_cells": int(adata.n_obs),
        "final_pas": int(adata.n_vars),
    }
    with open(cluster_h5ad.parent / "clustering_stats.json", "w") as fh:
        json.dump(stats, fh, indent=2)

    # B5: multi-sample runs previously emitted NO per-dataset stage funnel (the
    # single-sample branch of main.py did, via output_mgr.save_stats). Write the
    # per-dataset cell/PAS drop funnel so multi-sample cohort runs are auditable.
    stage_stats: dict[str, Any] = {
        "dataset_id": ds_id,
        "cb_filter": {
            "min_read": int(min_read),
            "cells_kept": len(_kept_cbs),
        },
        "input_matrix": {
            "input_pas": int(len(pas_ids)),
            "cells": int(len(collist)),
        },
        "annotated": {
            "annotated_pas": int(len(result.pas_ids)),
            "cells": int(len(result.collist)),
        },
        "preprocessing": {
            "final_cells": int(adata.n_obs),
            "final_pas": int(adata.n_vars),
            "min_cells": int(filter_min_cells),
            "min_genes": int(filter_min_genes),
        },
    }
    stage_stats_path = cluster_h5ad.parent / "stage_stats.json"
    with open(stage_stats_path, "w") as fh:
        json.dump(stage_stats, fh, indent=2)
    log.info("%s per-dataset stage funnel -> %s", prefix, stage_stats_path)

    # ------------------------------------------------------------------ #
    # 8. E3 provenance: record PAS/cell drops at the per-dataset drop     #
    #    sites (pas_gene, preprocess, cb_filter) and self-check the       #
    #    integrity invariant  surviving PAS == n_vars(clusters.h5ad).     #
    #    Each worker owns its OWN per-dataset ledger dir (no cross-process #
    #    sharing); atlas-snap drops live in the run-level ledger.         #
    # ------------------------------------------------------------------ #
    try:
        from ema.provenance import (
            ProvenanceLedger,
            check_survivor_invariant,
            record_cell_drops,
            record_pas_drops,
        )

        final_pas = [str(v) for v in adata.var_names]
        input_pas = [str(p) for p in pas_ids]
        annotated_pas = [str(p) for p in result.pas_ids]
        gene_of = {str(p): str(g) for p, g in zip(result.pas_ids, result.gene_ids)}

        led = ProvenanceLedger(cluster_h5ad.parent)
        surviving_pas = record_pas_drops(
            led,
            input_pas_ids=input_pas,
            annotated_pas_ids=annotated_pas,
            final_pas_ids=final_pas,
            dataset_id=ds_id,
            gene_of=gene_of,
            atlas_of=atlas_of,
            ip_of=ip_of,
        )
        record_cell_drops(
            led,
            input_cbs=[str(c) for c in sub_cbs],
            kept_cbs=[str(c) for c in _kept_cbs],
            final_cbs=[str(c) for c in adata.obs_names],
            dataset_id=ds_id,
        )
        led.flush()
        # Self-check: surviving PAS must equal the clustered matrix var count.
        # raise_on_fail=False → log a warning rather than kill a cohort worker;
        # a violation means the drop accounting missed a site (bug), not bad data.
        if not check_survivor_invariant(surviving_pas, int(adata.n_vars), raise_on_fail=False):
            log.warning(
                "%s provenance invariant OFF: surviving PAS=%d != n_vars=%d "
                "(pas_ledger accounting incomplete)",
                prefix, surviving_pas, int(adata.n_vars),
            )
        else:
            log.info("%s provenance invariant OK: %d surviving PAS == n_vars",
                     prefix, surviving_pas)
    except Exception:  # provenance is auxiliary — never fail the run over it
        log.warning("%s provenance ledger step failed (non-fatal)", prefix, exc_info=True)

    log.info("%s done — %d cells, %d PAS -> %s", prefix, adata.n_obs, adata.n_vars, cluster_h5ad)
    return stats


def downstream_worker_star(args: tuple) -> dict:
    """Top-level pickle-safe shim for ``pool.imap_unordered``.

    ``imap_unordered`` requires a callable that accepts a *single* argument.
    This shim unpacks the pre-built argument tuple and delegates to
    ``run_one_dataset_downstream``.  Must be defined at module level (not as a
    lambda or nested function) so it is picklable under the ``spawn`` context.

    Args:
        args: Tuple of positional arguments for ``run_one_dataset_downstream``.
            Supported arities:

            - 10 elements: ``(*pos_args[9], log_queue)`` — no progress client
            - 11 elements: ``(*pos_args[9], log_queue, progress_client)`` —
              progress client added by Phase 10 wiring
            - 12 elements: ``(*pos_args[9], log_queue, progress_client,
              cluster_kwargs_dict)`` — clustering hyperparameters added by
              YAML/CLI wiring fix.  ``cluster_kwargs_dict`` keys map to
              the ``cluster_*`` keyword args of
              :func:`run_one_dataset_downstream`.
            - 14 elements: ``(*pos_args[9], log_queue, progress_client,
              cluster_kwargs_dict, plot_engines, prov_kwargs_dict)`` (D9) —
              ``prov_kwargs_dict`` carries ``atlas_of``/``ip_of``, the
              per-unified-PAS atlas-match/internal-priming status maps.

    Returns:
        The stats dict returned by ``run_one_dataset_downstream``.
    """
    if len(args) == 14:
        # D9: atlas_of/ip_of (per-unified-PAS status maps) added as a single
        # trailing dict, mirroring the cluster_kwargs pattern, so
        # record_pas_drops() / write_pas_gene_artifacts() can populate the
        # ledger + annotatedpas.bed status columns per worker.
        *pos_args, log_queue, progress_client, cluster_kwargs, plot_engines, prov_kwargs = args
        return run_one_dataset_downstream(
            *pos_args,
            log_queue=log_queue,
            progress_client=progress_client,
            plot_engines=plot_engines,
            **(cluster_kwargs or {}),
            **(prov_kwargs or {}),
        )
    if len(args) == 13:
        *pos_args, log_queue, progress_client, cluster_kwargs, plot_engines = args
        return run_one_dataset_downstream(
            *pos_args,
            log_queue=log_queue,
            progress_client=progress_client,
            plot_engines=plot_engines,
            **(cluster_kwargs or {}),
        )
    if len(args) == 12:
        *pos_args, log_queue, progress_client, cluster_kwargs = args
        return run_one_dataset_downstream(
            *pos_args,
            log_queue=log_queue,
            progress_client=progress_client,
            **(cluster_kwargs or {}),
        )
    if len(args) == 11:
        *pos_args, log_queue, progress_client = args
        return run_one_dataset_downstream(
            *pos_args, log_queue=log_queue, progress_client=progress_client
        )
    if len(args) == 10:
        *pos_args, log_queue = args
        return run_one_dataset_downstream(*pos_args, log_queue=log_queue)
    return run_one_dataset_downstream(*args)
