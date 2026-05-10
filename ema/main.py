import logging
import threading
import shutil
from pathlib import Path

log = logging.getLogger(__name__)

from ema.countmatrix.peackcalling import peak_calling
from ema.countmatrix.indexing import get_mapping, reset_index
from ema.countmatrix.read import set_default_sample_id
from ema.config import directory_config, variable_config, args, filter_config
from ema.utils import get_resource_manager
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

# New strategy registries (Phase 1 of feature/apa-completeness)
# PDUI + diff strategies are imported by ema_switch (separate command), not here.
from ema.clustering.cross_dataset import get_match_strategy

try:
    from ema.datasets.atlas_snap import snap_beds_to_atlas
except ImportError:
    snap_beds_to_atlas = None  # type: ignore[assignment]


# extract_per_dataset_mtx is defined in downstream_runner to keep it
# pickle-safe and importable without triggering ema.config initialisation.
from ema.downstream_runner import extract_per_dataset_mtx  # noqa: E402


def run(
    cfg: dict | None = None,
    out_dir=None,
    progress=None,
    **kwargs,
) -> int:
    """Entry point called by ``ema run`` (Click CLI).

    Bridges the Click-driven ``cfg`` dict and keyword overrides into the
    module-level config dataclasses (``directory_config``, ``variable_config``,
    ``filter_config``, ``args``) that the existing pipeline body reads, then
    calls :func:`main` to execute the full pipeline.

    When ``cfg`` is ``None`` this function behaves identically to the legacy
    ``main()`` call — it consumes whatever the argparse-populated globals
    already contain.

    Args:
        cfg: Resolved YAML/CLI config dict from ``ema.cli.run``.  Must contain
            at minimum a ``"datasets"`` key.  Other recognised keys mirror the
            YAML schema (``gtf``, ``atlas``, ``atlas_distance``, ``seqlen``,
            ``cb_len``, ``barcode_tag``, ``min_read``, ``min_cells``,
            ``min_pas_per_cell``, ``pas_gap``).
        out_dir: Output directory (``pathlib.Path`` or ``str``).  When provided
            the pipeline writes all outputs here instead of the default
            ``emaout/`` directory.
        progress: A :class:`~ema.progress.ProgressManager` instance (already
            entered as a context manager).  When *None* progress wiring is a
            no-op — no bars are rendered.
        **kwargs: Additional per-run overrides forwarded from the Click command
            (e.g. ``peak_strategy``, ``dynamic_threshold``, ``tiles``, etc.).

    Returns:
        Exit code (0 on success).
    """
    # ------------------------------------------------------------------ #
    # 1. Bridge cfg + kwargs into module-level config dataclasses          #
    # ------------------------------------------------------------------ #
    if cfg is not None:
        # Output directory
        if out_dir is not None:
            directory_config.output_dir = str(out_dir)

        # Datasets list
        if "datasets" in cfg:
            directory_config.datasets = cfg["datasets"]

        # Optional path overrides from YAML
        if "gtf" in cfg:
            directory_config.gtf_dir = cfg["gtf"]
        if "atlas" in cfg:
            directory_config.atlas = cfg["atlas"]
        if "atlas_distance" in cfg:
            directory_config.atlas_distance = cfg["atlas_distance"]

        # variable_config scalars
        if "seqlen" in cfg:
            variable_config.seqlen = cfg["seqlen"]
        if "cb_len" in cfg:
            variable_config.cb_len = cfg["cb_len"]
        if "barcode_tag" in cfg:
            variable_config.barcode_tag = cfg["barcode_tag"]

        # filter_config scalars.
        #
        # The legacy `ema/cli/__init__.py::cli()` shim that runs at module
        # import time (via `ema/config.py`) only bridged `min_pas_per_cell`
        # from YAML — `min_read`, `min_cells`, `min_genes` silently kept
        # their class defaults (2000, 3, 50) regardless of YAML values.
        # The new `ema run` correctly applies all four. The atlas baseline
        # in reports/baseline_apa_completeness/ was regenerated to reflect
        # this — see CHANGELOG.md (0.2.0) for details.
        if "min_read" in cfg:
            filter_config.min_read = cfg["min_read"]
        if "min_cells" in cfg:
            filter_config.min_cells = cfg["min_cells"]
        if "min_pas_per_cell" in cfg:
            # YAML key `min_pas_per_cell` is the per-cell PAS threshold —
            # `preprocessing()` reads it as `filter_config.min_genes`.
            # Set both so the YAML key has its intended effect.
            filter_config.min_pas_per_cell = cfg["min_pas_per_cell"]
            filter_config.min_genes = cfg["min_pas_per_cell"]
        if "min_genes" in cfg:
            filter_config.min_genes = cfg["min_genes"]

        # Bridge kwargs into the argparse-style args namespace so the
        # pipeline body can read them via the existing `args.<attr>` pattern.
        _kwarg_to_args_map = {
            "peak_strategy": "strategy",
            "dynamic_threshold": "dynamic_threshold",
            "floor_threshold": "floor_threshold",
            "lambda_fold_change": "lambda_fold_change",
            "lambda_window": "lambda_window",
            "bam_threads": "bam_threads",
            "tiles": "tiles",
            "tile_size": "tile_size",
            "tile_overlap": "tile_overlap",
            "pipeline": "pipeline",
            "batch_size": "batch_size",
            "max_gene_distance": "max_gene_distance",
            "utr_multiplier": "utr_multiplier",
            "include_extended": "include_extended",
            "min_pas_per_cell": "min_pas_per_cell",
            "pas_gap": "pas_gap",
            "cluster_method": "clustering_method",
            "resolution": "resolution",
            "n_pcs": "n_pcs",
            "external_clusters": "external_clusters",
            "random_seed": "random_seed",
            "match_method": "cluster_match_method",
            "n_top_markers": "n_top_markers",
            # Strategy-tunable hyperparameters (forwarded into get_strategy())
            "lambda_method": "lambda_method",
            "max_pas": "max_pas",
            "smoothing_window": "smoothing_window",
            "min_prominence": "min_prominence",
        }
        for kw_key, args_attr in _kwarg_to_args_map.items():
            if kw_key in kwargs and kwargs[kw_key] is not None:
                setattr(args, args_attr, kwargs[kw_key])

    # ------------------------------------------------------------------ #
    # 2. Delegate to the pipeline body (progress wired inside)             #
    # ------------------------------------------------------------------ #
    _run_pipeline_body(progress=progress)
    return 0


def _run_pipeline_body(progress=None) -> None:
    """Execute the full pipeline using the current module-level config state.

    Args:
        progress: Optional :class:`~ema.progress.ProgressManager`.  Bars are
            registered here and :class:`~ema.progress.ProgressClient` handles
            are passed into workers.  All wiring is no-op when *progress* is
            ``None``.
    """
    # Helper: safely call add_stage even when progress is None
    def _add_stage(name: str, total=None) -> int | None:
        if progress is None:
            return None
        return progress.add_stage(name, total=total)

    def _add_subtask(parent, name: str, total=None) -> int | None:
        if progress is None or parent is None:
            return None
        return progress.add_subtask(parent, name, total=total)

    def _client(task_id) -> object | None:
        if progress is None or task_id is None:
            return None
        return progress.client(task_id)

    def _advance(task_id, n: int = 1) -> None:
        if progress is not None and task_id is not None:
            progress.client(task_id).advance(n)

    # ------------------------------------------------------------------
    # From this point on the code is identical to the original main() body.
    # ------------------------------------------------------------------

    # Set up output directory structure
    output_mgr = OutputManager(base_dir=directory_config.output_dir)
    output_mgr.setup()

    # Save run configuration
    output_mgr.save_run_config(vars(args))

    # Forward strategy-tunable hyperparameters that the user supplied via
    # CLI / YAML.  Only pass kwargs the strategy actually accepts (introspect
    # the constructor) so unrelated strategies aren't broken.
    import inspect
    _strategy_cls_kwargs: dict = {}
    for _attr in ("lambda_method", "max_pas", "smoothing_window", "min_prominence"):
        if hasattr(args, _attr) and getattr(args, _attr) is not None:
            _strategy_cls_kwargs[_attr] = getattr(args, _attr)
    try:
        from ema.strategies import _REGISTRY as _STRAT_REG
        _strat_cls = _STRAT_REG.get(args.strategy)
        if _strat_cls is not None:
            _accepted = set(inspect.signature(_strat_cls.__init__).parameters)
            _filtered_kwargs = {k: v for k, v in _strategy_cls_kwargs.items() if k in _accepted}
        else:
            _filtered_kwargs = {}
    except Exception:
        _filtered_kwargs = {}
    strategy = get_strategy(args.strategy, **_filtered_kwargs)
    peak_kwargs = dict(
        strategy=strategy,
        dynamic_threshold=args.dynamic_threshold,
        floor_threshold=args.floor_threshold,
        lambda_fold_change=args.lambda_fold_change,
        lambda_window=args.lambda_window,
        bam_threads=getattr(args, 'bam_threads', 4),
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
        dm = DatasetManager(
            output_dir=directory_config.output_dir,
            threads=getattr(args, 'bam_threads', 4),
        )
        bam_list = dm.prepare(directory_config.datasets)
    else:
        # Fallback: single BAM from legacy --bamDir
        bam_list = [("default", directory_config.bam_dir)]

    # Register top-level progress stages now that we know bam_list length.
    _peak_stage = _add_stage("Peak calling", total=len(bam_list))

    # -------------------------------------------------------------------------
    # Peak-calling — either global tile pool (--tiles) or sequential per-BAM.
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

    use_tiles: bool = getattr(args, "tiles", False)

    if use_tiles:
        # -----------------------------------------------------------------
        # TILE MODE (Phase 3): build a single flat JobSpec list across ALL
        # datasets × chroms × tiles × directions and dispatch through ONE
        # global multiprocessing.Pool with imap_unordered work-stealing.
        # -----------------------------------------------------------------
        from ema.countmatrix.tile_runner import (
            build_job_specs,
            run_all_jobs,
            merge_tiles,
        )
        from ema.utils.resource_manager import ResourceManager

        rm = get_resource_manager()
        n_workers = rm.get_n_jobs(per_worker_mb=300, stage="peak_tile")

        # Phase 4: RAM-adaptive tile sizing (unless --tile-size was explicitly
        # set by the user, i.e. differs from the CLI default of 25_000_000).
        _cli_tile_default = 25_000_000
        user_tile_size: int = getattr(args, "tile_size", _cli_tile_default)
        _tile_size_is_auto = (user_tile_size == _cli_tile_default)

        per_bam_tile_sizes: dict[str, int] = {}
        if _tile_size_is_auto:
            for _ds_id, _bam_path in bam_list:
                _bam_key = str(_bam_path)
                if _bam_key not in per_bam_tile_sizes:
                    per_bam_tile_sizes[_bam_key] = rm.get_tile_size(
                        bam_path=_bam_key,
                        n_workers=n_workers,
                        target_per_worker_mb=300,
                    )
            effective_tile_size = user_tile_size  # fallback if dict lookup misses
        else:
            effective_tile_size = user_tile_size

        tile_overlap: int = getattr(args, "tile_overlap", 10_000)

        # Resolve strategy name for cross-process serialisation
        _strategy_name: str
        if strategy is None:
            _strategy_name = "original"
        elif isinstance(strategy, str):
            _strategy_name = strategy
        else:
            from ema.strategies import _REGISTRY
            _strategy_name = next(
                (k for k, v in _REGISTRY.items() if isinstance(strategy, v)),
                "original",
            )

        jobs = build_job_specs(
            bam_list=[(ds_id, str(bp)) for ds_id, bp in bam_list],
            directions=[False, True],  # False=pos, True=neg
            tile_size=effective_tile_size,
            tile_overlap=tile_overlap,
            default_threshold=variable_config.default_threshold,
            merge_len=variable_config.merge_len,
            strategy_name=_strategy_name,
            dynamic_threshold=peak_kwargs.get("dynamic_threshold", False),
            floor_threshold=peak_kwargs.get("floor_threshold", 3),
            lambda_fold_change=peak_kwargs.get("lambda_fold_change", 2.0),
            lambda_window=peak_kwargs.get("lambda_window", 5000),
            bam_threads=peak_kwargs.get("bam_threads", 4),
            per_bam_tile_sizes=per_bam_tile_sizes if _tile_size_is_auto else None,
        )

        log.info(
            "%d dataset(s) -> %d jobs through Pool(%d)",
            len(bam_list), len(jobs), n_workers,
        )

        # Dispatch all jobs through one global pool
        # Pass a progress client so run_all_jobs can advance the peak bar per tile.
        _peak_client = _client(_peak_stage)
        grouped = run_all_jobs(jobs, n_workers=n_workers, progress_client=_peak_client)

        # Merge per-(dataset_id, direction) group and reconstruct file paths
        _ds_bam_indices: dict[str, int] = {}
        # Collect unique (dataset_id, bam_path) pairs in original order
        _seen_ds_bam: list[tuple[str, str]] = []
        for _ds_id, _bp in bam_list:
            _pair = (_ds_id, str(_bp))
            if _pair not in _seen_ds_bam:
                _seen_ds_bam.append(_pair)

        for _ds_id, _bam_path in _seen_ds_bam:
            _bam_path_str = str(_bam_path)
            _idx = _ds_bam_indices.get(_ds_id, 0)
            _ds_bam_indices[_ds_id] = _idx + 1

            pos_bed = peakcalling_dir / f"{_ds_id}_{_idx}.pos.bed"
            neg_bed = peakcalling_dir / f"{_ds_id}_{_idx}.neg.bed"
            pos_mtx = peakcalling_dir / f"{_ds_id}_{_idx}.pos.mtx"
            neg_mtx = peakcalling_dir / f"{_ds_id}_{_idx}.neg.mtx"
            cb_tsv = peakcalling_dir / f"{_ds_id}_{_idx}.cb.tsv"

            pos_cb = peakcalling_dir / f"{_ds_id}_{_idx}.pos_cb.tsv"
            neg_cb = peakcalling_dir / f"{_ds_id}_{_idx}.neg_cb.tsv"

            pos_tiles = grouped.get((_ds_id, False), [])
            neg_tiles = grouped.get((_ds_id, True), [])

            merge_tiles(pos_tiles, str(pos_bed), str(pos_mtx), str(pos_cb))
            merge_tiles(neg_tiles, str(neg_bed), str(neg_mtx), str(neg_cb))

            # Unify the two strand CB lists into one shared cb.tsv
            _all_cbs: list[str] = []
            _seen_cbs: set[str] = set()
            for _cb_file in (pos_cb, neg_cb):
                _p = Path(_cb_file)
                if _p.exists():
                    with open(_p) as _f:
                        for _line in _f:
                            _cb = _line.strip()
                            if _cb and _cb not in _seen_cbs:
                                _all_cbs.append(_cb)
                                _seen_cbs.add(_cb)
            with open(cb_tsv, "w") as _f:
                for _cb in _all_cbs:
                    _f.write(_cb + "\n")
            # Remove strand-split CB temp files
            for _p in (pos_cb, neg_cb):
                try:
                    Path(_p).unlink(missing_ok=True)
                except Exception:
                    pass

            all_pos_beds.append(str(pos_bed))
            all_neg_beds.append(str(neg_bed))
            all_pos_mtxs.append(str(pos_mtx))
            all_neg_mtxs.append(str(neg_mtx))
            all_pos_cbs.append(str(cb_tsv))
            all_neg_cbs.append(str(cb_tsv))
            all_dataset_ids_for_pos.append(_ds_id)
            all_dataset_ids_for_neg.append(_ds_id)

    else:
        # -----------------------------------------------------------------
        # SEQUENTIAL MODE: one dataset at a time (non-tile, existing path).
        # -----------------------------------------------------------------
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

            # Advance peak-calling bar once per dataset (sequential mode).
            _advance(_peak_stage)

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

        # Save CB filter stats — record what filter_cb actually used.
        output_mgr.save_stats("cb_filter", {
            "min_read": filter_config.min_read,
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

        # Forward all CLI/YAML clustering hyperparameters into clustering().
        # Without this, --cluster-method / --resolution / --n-pcs /
        # --random-seed / --external-clusters were silently ignored.
        clustering(
            adata=adata,
            method=getattr(args, "clustering_method", "leiden_tfidf"),
            resolution=getattr(args, "resolution", 1.0),
            n_pcs=getattr(args, "n_pcs", 40),
            random_seed=getattr(args, "random_seed", 42),
            external_clusters=getattr(args, "external_clusters", None),
        )

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
    _atlas_stage = _add_stage(
        "Atlas snap" if directory_config.atlas else "PAS merge",
        total=None,  # indeterminate spinner — we don't know peak count yet
    )
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
    # Mark atlas/merge stage complete
    _advance(_atlas_stage)

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

    log.info(
        "Unified %d BAMs (%s) -> %s",
        len(bam_list),
        "atlas" if directory_config.atlas else "merge",
        unified_bed,
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

    # ------------------------------------------------------------------ #
    # Phase 6: parallel per-dataset downstream pipeline                   #
    # ------------------------------------------------------------------ #
    # Pre-build (sub_indices, sub_cbs) for every dataset in the parent
    # process — O(n_cells) each, cheap.  Empty datasets are filtered out.
    import json
    import multiprocessing
    import pickle

    from ema.downstream_runner import (
        run_one_dataset_downstream,
        downstream_worker_star,
    )
    from ema.logging_config import get_log_queue

    # Serialise the genes DataFrame once; each worker deserialises its own copy.
    genes_pkl: bytes = pickle.dumps(genes)

    # Grab the parent's log queue so spawn workers can route records back.
    _log_queue = get_log_queue()

    # Register downstream stage (total = number of non-empty datasets).
    # We build worker_args first so we know n_datasets before adding the stage.
    _raw_worker_args: list[tuple] = []
    for ds_id in unique_ds_ids:
        sub_indices = [
            i for i, cb in enumerate(all_cb_strings)
            if cb.startswith(f"{ds_id}_")
        ]
        if not sub_indices:
            log.warning("no cells found for dataset '%s' — skipping", ds_id)
            continue
        sub_cbs = [all_cb_strings[i] for i in sub_indices]
        _raw_worker_args.append((
            ds_id,
            sub_indices,
            sub_cbs,
            str(unified_mtx),
            str(per_dataset_dir),
            genes_pkl,
            filter_config.min_read,
            filter_config.min_cells,
            filter_config.min_genes,
            _log_queue,
        ))

    # Register downstream progress stage now that we know how many datasets.
    _downstream_stage = _add_stage(
        "Per-dataset downstream", total=len(_raw_worker_args)
    )

    # Build final worker_args — for inline path we add per-dataset progress
    # clients (6 ticks per dataset: extract/filter/build/annotate/preprocess/cluster).
    # For pool path we append the client to the tuple so downstream_worker_star can unpack it.
    # Forward all clustering hyperparameters from CLI/YAML.  Without this,
    # --cluster-method / --resolution / --n-pcs / --random-seed / --external-clusters
    # were silently ignored in the multi-sample path.
    _cluster_kwargs: dict = {
        "cluster_method": getattr(args, "clustering_method", "leiden_tfidf"),
        "cluster_resolution": getattr(args, "resolution", 1.0),
        "cluster_n_pcs": getattr(args, "n_pcs", 40),
        "cluster_random_seed": getattr(args, "random_seed", 42),
        "cluster_external_clusters": getattr(args, "external_clusters", None),
    }
    worker_args: list[tuple] = []
    for _warg in _raw_worker_args:
        _ds_id_for_sub = _warg[0]
        _sub_stage = _add_subtask(_downstream_stage, _ds_id_for_sub, total=6)
        _sub_client = _client(_sub_stage)
        worker_args.append(_warg + (_sub_client, _cluster_kwargs))

    # Decide worker count: cap by RAM budget (each AnnData ~ 200-500 MB).
    n_datasets = len(worker_args)
    n_workers = min(
        n_datasets,
        get_resource_manager().get_n_jobs(per_worker_mb=500, stage="downstream"),
    )

    log.info(
        "Per-dataset downstream: %d datasets, n_workers=%d",
        n_datasets, n_workers,
    )

    if n_workers <= 1 or n_datasets == 1:
        # Inline path: no spawn overhead, backward-compatible.
        for arg_tuple in worker_args:
            # arg_tuple has: ds_id, sub_indices, sub_cbs, unified_mtx,
            #   per_dataset_dir, genes_pkl, min_read, min_cells, min_genes,
            #   log_queue, progress_client, cluster_kwargs  (12 elements)
            *pos_args, lq, pc, ck = arg_tuple
            run_one_dataset_downstream(
                *pos_args,
                log_queue=lq,
                progress_client=pc,
                **(ck or {}),
            )
            _advance(_downstream_stage)
    else:
        ctx = multiprocessing.get_context("spawn")
        with ctx.Pool(processes=n_workers) as pool:
            for _stats in pool.imap_unordered(
                downstream_worker_star, worker_args, chunksize=1
            ):
                _advance(_downstream_stage)  # one dataset complete

    # NOTE: PDUI and differential APA are NOT run here — they belong to the
    # separate `ema_switch` command. That command lets the user select which
    # cluster pairs to test and which marker-PAS subset to use, instead of
    # running all 20K PAS x all C(K,2) pairs unconditionally.

    # Cross-dataset cluster matching (only meaningful if >1 dataset)
    if len(unique_ds_ids) > 1:
        try:
            h5ad_paths = [per_dataset_dir / ds / "clusters.h5ad" for ds in unique_ds_ids]
            existing = [(p, ds) for p, ds in zip(h5ad_paths, unique_ds_ids) if p.exists()]
            if len(existing) > 1:
                match_strategy = get_match_strategy(args.cluster_match_method)
                match_df = match_strategy.match(
                    h5ad_paths=[p for p, _ in existing],
                    dataset_ids=[ds for _, ds in existing],
                    n_top_markers=args.n_top_markers,
                )
                cross_dir = output_dir / "cross_dataset"
                cross_dir.mkdir(exist_ok=True)
                match_df.to_csv(cross_dir / "canonical_cluster_map.tsv", sep="\t", index=False)
                n_canonical = match_df['canonical_cluster'].nunique() if 'canonical_cluster' in match_df.columns else 0
                log.info(
                    "Cross-dataset matching (%s): %d cluster entries -> %d canonical clusters",
                    args.cluster_match_method, len(match_df), n_canonical,
                )
        except Exception as e:
            log.warning("Cross-dataset matching failed: %s", e)

    return  # done with multi-sample path


def main() -> None:
    """Legacy entry point: runs the pipeline using the argparse-populated globals.

    Delegates to :func:`_run_pipeline_body` with no progress bars.
    Retained for backward compatibility (``if __name__ == '__main__'`` and any
    direct callers that have not yet migrated to :func:`run`).
    """
    _run_pipeline_body(progress=None)


if __name__ == "__main__":
    main()
