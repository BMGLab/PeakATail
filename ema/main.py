import json
import logging
import threading
import shutil
import time
from pathlib import Path

log = logging.getLogger(__name__)

from ema.countmatrix.peackcalling import peak_calling
from ema.countmatrix.indexing import get_mapping, reset_index
from ema.countmatrix.read import set_default_sample_id
from ema.config import directory_config, variable_config, args, filter_config, set_directory_config
from ema.utils import get_resource_manager
from ema.matrixfilter import filter_cb, make_dataframe, preprocessing
from ema.clustering.clustering import clustering
from ema.annotate.annotate import annotate
from ema.annotate.find_close import find_close
from ema.annotate.gtf_cache import process_gtf_cached
from ema.strategies import get_strategy
from ema.outputs import OutputManager

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


# ---------------------------------------------------------------------------
# Resource sampler (Tier 4 viz support)
# ---------------------------------------------------------------------------


class _ResourceSampler(threading.Thread):
    """Background thread that periodically samples process memory and CPU usage.

    Samples are appended as newline-delimited JSON records to *out_path*.
    Each record has the schema::

        {"elapsed_s": float, "rss_gb": float, "cpu_pct": float}

    The thread is a daemon thread so it never blocks process exit.  Call
    :meth:`stop` and then :meth:`join` to request a clean shutdown.

    The sampler degrades gracefully when ``psutil`` is not installed: it
    still runs (so callers need not check for its presence) but writes no
    records and logs a one-time warning.

    Args:
        out_path: Path where the JSONL file is written.
        interval_s: Sampling interval in seconds (default 5.0).
    """

    def __init__(self, out_path: Path, interval_s: float = 5.0) -> None:
        super().__init__(name="resource-sampler", daemon=True)
        self._out_path = Path(out_path)
        self._interval_s = interval_s
        self._stop_event = threading.Event()
        self._t0 = time.monotonic()

    def stop(self) -> None:
        """Signal the sampler to stop at the next interval boundary."""
        self._stop_event.set()

    def join(self, timeout: float | None = None) -> None:
        """Join the thread only if it was ever started."""
        if self.is_alive() or self._started.is_set():  # type: ignore[attr-defined]
            super().join(timeout=timeout)

    def run(self) -> None:  # noqa: D102 — threading.Thread override
        try:
            import psutil  # type: ignore[import-untyped]
        except ImportError:
            log.warning(
                "psutil not installed — resource_timeline sampling disabled. "
                "Install it with: pip install psutil"
            )
            return

        self._out_path.parent.mkdir(parents=True, exist_ok=True)
        proc = psutil.Process()
        # Prime the CPU measurement (first call always returns 0.0)
        try:
            proc.cpu_percent(interval=None)
        except Exception:
            pass

        with open(self._out_path, "a") as fh:
            while not self._stop_event.is_set():
                try:
                    rss_gb = proc.memory_info().rss / (1024 ** 3)
                    cpu_pct = proc.cpu_percent(interval=None)
                    elapsed = time.monotonic() - self._t0
                    record = {
                        "elapsed_s": round(elapsed, 3),
                        "rss_gb": round(rss_gb, 4),
                        "cpu_pct": round(cpu_pct, 2),
                    }
                    fh.write(json.dumps(record) + "\n")
                    fh.flush()
                except Exception as exc:
                    log.debug("Resource sampler error: %s", exc)

                self._stop_event.wait(timeout=self._interval_s)


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

    Implementation -- single bridge via :class:`ema.cli.config_schema.RunConfig`.
    The two hand-maintained maps (``_kwarg_to_args_map`` and
    ``_cfg_to_args_map``) that used to live here have been replaced by
    ``RunConfig.apply_to_legacy_globals()``.  See ema/cli/config_schema.py
    for the field-by-field bridge metadata.

    Args:
        cfg: Resolved YAML/CLI config dict from ``ema.cli.run``.  Must contain
            at minimum a ``"datasets"`` key.  Other recognised keys mirror the
            schema's yaml_key field (gtf, atlas, atlas_distance, seqlen,
            cb_len, barcode_tag, min_read, min_cells, min_pas_per_cell,
            pas_gap, ...).
        out_dir: Output directory (``pathlib.Path`` or ``str``).
        progress: A :class:`~ema.progress.ProgressManager` instance.
        **kwargs: Additional per-run overrides from the Click command
            (e.g. ``peak_strategy``, ``dynamic_threshold``, ``tiles``).

    Returns:
        Exit code (0 on success).
    """
    if cfg is not None:
        from ema.cli.config_schema import (
            RunConfig,
            field_specs,
            yaml_key_to_field_name,
            legacy_alias_to_field_name,
        )

        # ------------------------------------------------------------------
        # Step 1: build a RunConfig instance starting from YAML, then layer
        # kwargs on top.  ``run.py::_pipeline_kwargs()`` filters kwargs to
        # only those the user explicitly typed -- when the user didn't set
        # a flag, kwargs[k] is absent, so the YAML value (or schema default)
        # wins automatically.
        # ------------------------------------------------------------------
        rc = RunConfig.from_yaml_dict(cfg)

        # Map kwargs (snake_case Click variable names) onto RunConfig fields.
        field_names = {f.name for f in __import__("dataclasses").fields(RunConfig)}
        for kw_key, val in kwargs.items():
            if kw_key in field_names:
                setattr(rc, kw_key, val)

        # ------------------------------------------------------------------
        # Step 2: directory + datasets that don't fit the schema's scalar
        # bridge model (datasets is a list of dicts; output_dir we already
        # know from the resolved Path passed in).
        # ------------------------------------------------------------------
        # Build the kwargs for set_directory_config; only pass keys that were
        # actually supplied so we don't overwrite defaults with None.
        _dc_kwargs: dict = {}
        if out_dir is not None:
            from pathlib import Path as _Path
            _dc_kwargs["output_dir"] = _Path(out_dir)
        if "datasets" in cfg:
            _dc_kwargs["datasets"] = cfg["datasets"]
        if cfg.get("gtf"):
            _dc_kwargs["gtf_dir"] = cfg["gtf"]
        if cfg.get("bam_dir"):
            _dc_kwargs["bam_dir"] = cfg["bam_dir"]
        if cfg.get("atlas"):
            _dc_kwargs["atlas"] = cfg["atlas"]
        if cfg.get("atlas_distance") is not None:
            _dc_kwargs["atlas_distance"] = cfg["atlas_distance"]
        # Filename overrides from RunConfig.filenames (YAML key "filenames")
        if rc.filenames:
            _dc_kwargs["filenames"] = rc.filenames
        if _dc_kwargs:
            set_directory_config(**_dc_kwargs)

        # ------------------------------------------------------------------
        # Step 3: single bridge -- schema-driven mutation of legacy globals.
        # The bridge respects each field's ``legacy_dataclass_attr`` /
        # ``legacy_args_attr`` instructions; new fields land on the right
        # legacy attribute automatically.
        # ------------------------------------------------------------------
        rc.apply_to_legacy_globals()

    # ------------------------------------------------------------------ #
    # Delegate to the pipeline body (progress wired inside)               #
    # ------------------------------------------------------------------ #
    _plot_engines: list[str] | None = kwargs.pop("plot_engines", None)

    # Tier 4 viz: start the resource sampler iff plotting is enabled.
    # The sampler writes a JSONL of (elapsed, rss_gb, cpu_pct) records
    # next to the run's outputs; the resource_timeline viz reads it later.
    _sampler: _ResourceSampler | None = None
    if _plot_engines:
        try:
            _resource_jsonl = Path(directory_config.output_dir) / "resources.jsonl"
            _sampler = _ResourceSampler(_resource_jsonl, interval_s=5.0)
            _sampler.start()
        except Exception as e:
            log.warning("resource sampler failed to start: %s", e)
            _sampler = None

    _pipeline_result: dict | None = None
    try:
        _pipeline_result = _run_pipeline_body(progress=progress, plot_engines=_plot_engines)
    finally:
        if _sampler is not None:
            try:
                _sampler.stop()
                _sampler.join(timeout=10)
            except Exception as e:
                log.warning("resource sampler clean shutdown failed: %s", e)

    # All viz happens in a single post-pipeline orchestrator call (see
    # ema/viz/pipeline_hooks.py). The pipeline body has written every
    # artifact the orchestrator needs (h5ads, BEDs, atlas mapping, JSONL).
    if _pipeline_result is not None:
        from ema.viz.pipeline_hooks import render_run_outputs
        render_run_outputs(
            output_dir=Path(directory_config.output_dir),
            engines=_plot_engines,
            is_single_sample=_pipeline_result["is_single_sample"],
            unique_ds_ids=_pipeline_result["unique_ds_ids"],
            bed_paths=_pipeline_result["bed_paths"],
            atlas_enabled=_pipeline_result["atlas_enabled"],
            clustering_dir=_pipeline_result.get("clustering_dir"),
            single_sample_h5ad=_pipeline_result.get("single_sample_h5ad"),
        )

    # E2: write the run manifest — the contract artifact the hub indexes. Built
    # from the RESOLVED config (B0) + auto-discovered artifacts + id grammar.
    try:
        from ema.outputs import OutputManager, build_resolved_run_config
        _mgr = OutputManager(base_dir=str(directory_config.output_dir))
        _mgr.write_manifest(build_resolved_run_config())
    except Exception as _e:  # never fail a run over the manifest
        log.warning("run_manifest.json write failed: %s", _e)

    return 0


def _run_pipeline_body(progress=None, plot_engines: list[str] | None = None) -> dict:
    """Execute the full pipeline using the current module-level config state.

    Returns a metadata dict consumed by :func:`ema.viz.pipeline_hooks.render_run_outputs`
    so all visualisation can happen post-pipeline from on-disk artifacts.

    Args:
        progress: Optional :class:`~ema.progress.ProgressManager`.  Bars are
            registered here and :class:`~ema.progress.ProgressClient` handles
            are passed into workers.  All wiring is no-op when *progress* is
            ``None``.
        plot_engines: List of engine names to use for visualizations
            (e.g. ``["matplotlib", "plotly"]``).  ``None`` uses the default
            ``["matplotlib", "plotly"]``.  Empty list disables all plots.
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

    # Save run configuration (B0: serialize the RESOLVED config actually in
    # effect — directory_config/variable_config/filter_config — not the raw
    # argparse defaults, which record atlas=null on runs that snapped).
    from ema.outputs import build_resolved_run_config
    output_mgr.save_run_config(build_resolved_run_config())

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
    except (ImportError, TypeError, AttributeError) as exc:
        log.warning("strategy %r kwarg filter failed: %s — using empty kwargs", args.strategy, exc)
        _filtered_kwargs = {}
    strategy = get_strategy(args.strategy, **_filtered_kwargs)
    peak_kwargs = dict(
        strategy=strategy,
        dynamic_threshold=args.dynamic_threshold,
        floor_threshold=args.floor_threshold,
        lambda_fold_change=args.lambda_fold_change,
        lambda_window=args.lambda_window,
        bam_threads=getattr(args, 'bam_threads', 4),
        # Post-detection PAS merger (strategy-agnostic).
        # -1 spacing triggers auto-detect (median read length per BAM) inside
        # peak_calling / run_tiled / run_pipeline.  Prominence is the static
        # fallback used by non-lambda strategies; lambda strategies override
        # to compute_lambda(heights).
        min_pas_spacing=variable_config.min_pas_spacing,
        min_pas_prominence=variable_config.min_pas_prominence,
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
    # Total is set to 0 here; the per-chromosome progress_client will call
    # set_total() once the BAM is opened (sequential mode) or the tile pool
    # drives the bar directly via run_all_jobs (tile mode).
    _peak_stage = _add_stage("Peak calling", total=0)

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
            min_pas_spacing=variable_config.min_pas_spacing,
            min_pas_prominence=variable_config.min_pas_prominence,
        )

        log.info(
            "%d dataset(s) -> %d jobs through Pool(%d)",
            len(bam_list), len(jobs), n_workers,
        )

        # Dispatch all jobs through one global pool
        # Pass a progress client so run_all_jobs can advance the peak bar per tile.
        # Pass a timings_path so per-tile wall-seconds get persisted for the
        # tile_timing viz strategy (see Tier 4 viz block at end of run()).
        _peak_client = _client(_peak_stage)
        _timings_path = Path(directory_config.output_dir) / "tile_timings.json"
        grouped = run_all_jobs(
            jobs,
            n_workers=n_workers,
            progress_client=_peak_client,
            timings_path=_timings_path,
        )

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

            # Pass the progress client to the pos-strand call so the bar
            # shows per-chromosome progress.  set_total() fires once the BAM
            # is opened; advance() fires on each chromosome boundary.  Both
            # strand passes share the same client — peak_calling() sets total
            # to `nonempty_refs * 2` so the advances from both passes fill
            # the bar.
            _peak_client = _client(_peak_stage)
            peak_calling(
                False,
                bedfilepath=str(pos_bed),
                matrixpath=str(pos_mtx),
                bamfile_dir=str(bam_path),
                progress_client=_peak_client,
                **peak_kwargs,
            )
            peak_calling(
                True,
                bedfilepath=str(neg_bed),
                matrixpath=str(neg_mtx),
                bamfile_dir=str(bam_path),
                progress_client=_peak_client,
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

    # Force the peak-calling bar to 100 % — inner read filtering can drop
    # chroms below the chrom-change trigger threshold, so the advance count
    # is sometimes a few short of the apriori `nonempty * 2` estimate.
    if progress is not None and _peak_stage is not None:
        progress.finish(_peak_stage)

    # Save peak calling stats
    output_mgr.save_stats("peak_calling", {
        "strategy": args.strategy,
        "dynamic_threshold": args.dynamic_threshold,
        "floor_threshold": args.floor_threshold,
        "lambda_fold_change": args.lambda_fold_change,
        "lambda_window": args.lambda_window,
    })

    # Per-stage data snapshots — every step that mutates the data gets a
    # canonical file on disk.  See ema/outputs.py for the layout.
    from ema.outputs import write_per_dataset_beds, write_raw_peak_outputs
    output_dir = Path(directory_config.output_dir)
    _ds_pos: dict[str, list[str]] = {}
    _ds_neg: dict[str, list[str]] = {}
    _ds_pos_mtx: dict[str, list[str]] = {}
    _ds_neg_mtx: dict[str, list[str]] = {}
    _ds_cb: dict[str, str] = {}
    for _ds_id, _bed in zip(all_dataset_ids_for_pos, all_pos_beds):
        _ds_pos.setdefault(_ds_id, []).append(_bed)
    for _ds_id, _bed in zip(all_dataset_ids_for_neg, all_neg_beds):
        _ds_neg.setdefault(_ds_id, []).append(_bed)
    for _ds_id, _mtx in zip(all_dataset_ids_for_pos, all_pos_mtxs):
        _ds_pos_mtx.setdefault(_ds_id, []).append(_mtx)
    for _ds_id, _mtx in zip(all_dataset_ids_for_neg, all_neg_mtxs):
        _ds_neg_mtx.setdefault(_ds_id, []).append(_mtx)
    for _ds_id, _cb in zip(all_dataset_ids_for_pos, all_pos_cbs):
        _ds_cb.setdefault(_ds_id, _cb)
    # 1. RAW peak-calling snapshot (pre-filter): per_dataset/<ds>/raw/
    write_raw_peak_outputs(output_dir, _ds_pos, _ds_neg, _ds_pos_mtx, _ds_neg_mtx, _ds_cb)
    # 2. Combined post-merge BEDs:          per_dataset/<ds>/{posbed,negbed,pasbed}.bed
    write_per_dataset_beds(output_dir, _ds_pos, _ds_neg)

    # =========================================================================
    # Single-sample path (len(bam_list) == 1)
    # Copy outputs to legacy paths and continue with the existing pipeline.
    # =========================================================================
    if len(bam_list) == 1:
        shutil.copy(all_pos_beds[0], directory_config.posbed)
        shutil.copy(all_neg_beds[0], directory_config.negbed)
        shutil.copy(all_pos_mtxs[0], directory_config.posmatrixpath)
        shutil.copy(all_neg_mtxs[0], directory_config.negmatrixpath)

        _cb_filter_stage = _add_stage("CB filter", total=1)
        filter_cb()
        _advance(_cb_filter_stage)

        # Save CB filter stats + persist the kept barcode list for this dataset.
        output_mgr.save_stats("cb_filter", {
            "min_read": filter_config.min_read,
        })
        # B3: ALWAYS write the per-dataset filtered_cb.tsv. The old
        # ``if _filtered_cb_path.exists()`` guard silently skipped the write, so
        # real runs left 02_cb_filter/ with only the stats JSON — making the
        # min_read survivor set unrecoverable and the documented samtools -D
        # workflow impossible. Prefer the authoritative in-memory list that
        # filter_cb() just populated; fall back to the on-disk run-root file.
        from ema.outputs import write_filtered_cb
        import ema.matrixfilter as _mf
        _filtered_cb_path = Path(directory_config.filtered_cb)
        if _mf.filtered_cb_list:
            _kept_cbs = list(_mf.filtered_cb_list)
        elif _filtered_cb_path.exists():
            _kept_cbs = [
                b.strip()
                for b in _filtered_cb_path.read_text().splitlines()
                if b.strip()
            ]
        else:
            _kept_cbs = []
            log.warning(
                "cb_filter produced no filtered barcode list for %r; writing an "
                "empty filtered_cb.tsv (header only) for traceability.",
                bam_list[0][0],
            )
        write_filtered_cb(
            output_dir, bam_list[0][0], _kept_cbs, filter_config.min_read,
        )

        # Wait for GTF processing to complete before find_close
        _gtf_stage = _add_stage("GTF annotation", total=1)
        gtf_thread.join()
        if "exception" in gtf_error:
            raise RuntimeError(
                f"GTF processing failed: {gtf_error['exception']}"
            ) from gtf_error["exception"]
        utr_lengths = gtf_result.get("utr_lengths", {})
        _advance(_gtf_stage)

        # Save GTF annotation stats
        output_mgr.save_stats("gtf_annotation", {
            "gtf_path": directory_config.gtf_dir,
            "utr_count": len(utr_lengths),
        })

        # Find closest gene for each PAS
        _pas_gene_stage = _add_stage("PAS→gene assignment", total=1)
        genes = find_close(
            utr_lengths=utr_lengths,
            max_distance=getattr(args, "max_gene_distance", 5000),
            utr_multiplier=getattr(args, "utr_multiplier", 2.0),
            include_extended=getattr(args, "include_extended", False),
        )
        _advance(_pas_gene_stage)

        # Save PAS-gene assignment stats
        output_mgr.save_stats("pas_gene", {
            "max_gene_distance": getattr(args, "max_gene_distance", 5000),
            "utr_multiplier": getattr(args, "utr_multiplier", 2.0),
            "assigned_pas_count": len(genes),
        })

        # Build sparse matrix (PAS IDs preserved)
        _annot_stage = _add_stage("Annotated matrix", total=1)
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

        # Persist canonical PAS->gene mapping + annotatedpas.bed +
        # the annotated count matrix (post PAS-gene join, pre cell/PAS
        # filter).  See ema/outputs.py for the file layout.
        from ema.outputs import write_pas_gene_artifacts, write_annotated_matrix
        write_pas_gene_artifacts(
            output_dir, bam_list[0][0],
            result.pas_ids, result.gene_ids,
        )
        write_annotated_matrix(
            output_dir, bam_list[0][0],
            result.sparse_matrix, result.pas_ids, result.collist,
        )
        _advance(_annot_stage)

        # Preprocess and cluster.
        # Pass min_cells/min_genes EXPLICITLY: matrixfilter.preprocessing's
        # default kwargs are evaluated at function-def time, so they snapshot
        # filter_config.min_cells/min_genes from import.  Without this explicit
        # forward, YAML overrides applied to filter_config above would be
        # silently ignored on the single-sample path.
        _preproc_stage = _add_stage("Preprocessing", total=1)
        adata = preprocessing(
            sparse_matrix=result.sparse_matrix,
            pas_ids=result.pas_ids,
            collist=collist,
            min_cells=filter_config.min_cells,
            min_genes=filter_config.min_genes,
            gene_ids=result.gene_ids,
        )

        # Save preprocessing stats + persist the post-filter AnnData
        # (cluster labels haven't been added yet — clusters.h5ad will
        # supersede this once clustering() runs, but having the snapshot
        # lets users inspect the filter step in isolation).
        output_mgr.save_stats("preprocessing", {
            "cells_after_filter": adata.n_obs,
            "pas_after_filter": adata.n_vars,
        })
        from ema.outputs import write_preprocessed_h5ad
        write_preprocessed_h5ad(output_dir, bam_list[0][0], adata)
        _advance(_preproc_stage)

        # Forward all CLI/YAML clustering hyperparameters into clustering().
        # Persist the clustered AnnData to the canonical on-disk path via
        # directory_config (07_clustering/<ds>/clusters.h5ad).
        _cluster_stage = _add_stage("Clustering", total=1)
        _ss_h5ad = directory_config.clusters_h5ad_for(bam_list[0][0])
        _ss_h5ad.parent.mkdir(parents=True, exist_ok=True)
        clustering(
            adata=adata,
            method=getattr(args, "clustering_method", "leiden_tfidf"),
            resolution=getattr(args, "resolution", 1.0),
            n_pcs=getattr(args, "n_pcs", 40),
            random_seed=getattr(args, "random_seed", 42),
            external_clusters=getattr(args, "external_clusters", None),
            output_h5ad=str(_ss_h5ad),
            # New tunable hyperparameters from Wave 2A.  Multi-sample path
            # forwards these via _cluster_kwargs; single-sample path must
            # too or the flags are silently ignored for 1-dataset runs.
            n_neighbors=getattr(args, "n_neighbors", None),
            tfidf_scale_factor=getattr(args, "tfidf_scale_factor", 1e4),
            depth_corr_threshold=getattr(args, "depth_corr_threshold", 0.75),
            n_svd_components=getattr(args, "n_svd_components", 50),
            n_top_hvg=getattr(args, "n_top_hvg", 2000),
        )
        _advance(_cluster_stage)

        # Save clustering stats
        output_mgr.save_stats("clustering", {
            "final_cells": adata.n_obs,
            "final_pas": adata.n_vars,
        })

        return {  # done with single-sample path; viz happens in run() wrapper
            "is_single_sample": True,
            "unique_ds_ids": [bam_list[0][0]],
            "bed_paths": all_pos_beds + all_neg_beds,
            "atlas_enabled": bool(directory_config.atlas),
            "single_sample_h5ad": _ss_h5ad,
            "clustering_dir": directory_config.clustering_dir,
        }

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
    # B1: the SAME dataset_id appears once per strand in all_ds_ids; carry an
    # explicit strand per entry so the merge/atlas count-routing key becomes
    # (dataset_id, strand, pasnumber) and pos/neg PAS #N cannot collide.
    all_strands = (
        ["+"] * len(all_dataset_ids_for_pos) + ["-"] * len(all_dataset_ids_for_neg)
    )

    # Dispatch atlas vs. coordinate-merge based on config
    _atlas_stage = _add_stage(
        "Atlas snap" if directory_config.atlas else "PAS merge",
        total=None,  # indeterminate spinner — we don't know peak count yet
    )
    # E3: provenance ledger for the biggest silent drop site (atlas snap).
    _prov_ledger = None
    if directory_config.atlas:
        if snap_beds_to_atlas is None:
            raise RuntimeError("ema.datasets.atlas_snap is not available")
        try:
            from ema.provenance import ProvenanceLedger
            _prov_ledger = ProvenanceLedger(output_dir)
        except Exception as _e:  # provenance must never break a run
            log.warning("provenance ledger init failed: %s", _e)
            _prov_ledger = None
        unified_bed, mapping_path = snap_beds_to_atlas(
            all_beds,
            all_ds_ids,
            atlas_bed=directory_config.atlas,
            output_dir=unified_dir,
            distance=directory_config.atlas_distance,
            strands=all_strands,
            ledger=_prov_ledger,
        )
        if _prov_ledger is not None:
            try:
                _prov_ledger.flush()
                log.info(
                    "provenance: atlas-snap pas_ledger written (%d surviving PAS)",
                    _prov_ledger.count_surviving("pas"),
                )
            except Exception as _e:
                log.warning("provenance ledger flush failed: %s", _e)
    else:
        unified_bed, mapping_path = merge_pas_beds(
            all_beds,
            all_ds_ids,
            output_dir=unified_dir,
            gap=getattr(args, "pas_gap", 100),
            strands=all_strands,
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
        strands=all_strands,
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
            str(output_dir),
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
        # New tunable hyperparameters (default values match strategy constructors).
        "cluster_n_neighbors": getattr(args, "n_neighbors", None),
        "cluster_tfidf_scale_factor": getattr(args, "tfidf_scale_factor", 1e4),
        "cluster_depth_corr_threshold": getattr(args, "depth_corr_threshold", 0.75),
        "cluster_n_svd_components": getattr(args, "n_svd_components", 50),
        "cluster_n_top_hvg": getattr(args, "n_top_hvg", 2000),
    }
    # Resolve plot_engines: default to matplotlib-only (plotly is opt-in via
    # --plot-engine plotly|both).
    _resolved_plot_engines: list[str] = (
        plot_engines if plot_engines is not None else ["matplotlib"]
    )
    worker_args: list[tuple] = []
    for _warg in _raw_worker_args:
        _ds_id_for_sub = _warg[0]
        _sub_stage = _add_subtask(_downstream_stage, _ds_id_for_sub, total=6)
        _sub_client = _client(_sub_stage)
        worker_args.append(_warg + (_sub_client, _cluster_kwargs, _resolved_plot_engines))

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
            #   log_queue, progress_client, cluster_kwargs, plot_engines  (13 elements)
            *pos_args, lq, pc, ck, pe = arg_tuple
            run_one_dataset_downstream(
                *pos_args,
                log_queue=lq,
                progress_client=pc,
                plot_engines=pe,
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

    # Cross-dataset cluster matching (only meaningful if >1 dataset).
    #
    # Pre-fix this block had a bare ``except Exception`` that swallowed every
    # error including ImportErrors, KeyErrors, and pipeline-level
    # programming bugs.  We now narrow the catch to:
    #
    #   * FileNotFoundError    — h5ad missing on disk (data issue)
    #   * KeyError             — strategy/registry key missing (config issue)
    #   * ValueError           — match strategy rejected the inputs
    #   * pandas/IO errors during the to_csv() write
    #
    # Anything else (TypeError, ImportError, AttributeError, ...) is a true
    # programming bug and re-raises so it is visible in the run log instead
    # of being silently demoted to a one-line WARNING.
    if len(unique_ds_ids) > 1:
        h5ad_paths = [directory_config.clusters_h5ad_for(ds) for ds in unique_ds_ids]
        existing = [(p, ds) for p, ds in zip(h5ad_paths, unique_ds_ids) if p.exists()]
        if len(existing) > 1:
            try:
                match_strategy = get_match_strategy(args.cluster_match_method)
                match_df = match_strategy.match(
                    h5ad_paths=[p for p, _ in existing],
                    dataset_ids=[ds for _, ds in existing],
                    n_top_markers=args.n_top_markers,
                )
                cross_dir = output_dir / "cross_dataset"
                cross_dir.mkdir(exist_ok=True)
                match_df.to_csv(cross_dir / "canonical_cluster_map.tsv", sep="\t", index=False)
                n_canonical = (
                    match_df['canonical_cluster'].nunique()
                    if 'canonical_cluster' in match_df.columns else 0
                )
                log.info(
                    "Cross-dataset matching (%s): %d cluster entries -> %d canonical clusters",
                    args.cluster_match_method, len(match_df), n_canonical,
                )
                # B6: round-trip the canonical map back into each dataset's obs
                # so cross-sample comparisons join on a shared id space instead
                # of the per-sample (non-comparable) leiden labels.
                from ema.clustering.cross_dataset.roundtrip import write_canonical_clusters
                try:
                    write_canonical_clusters(match_df, existing)
                except (KeyError, OSError, ValueError) as e:
                    log.warning(
                        "canonical_cluster round-trip into obs skipped (%s): %s",
                        type(e).__name__, e,
                    )
            except (FileNotFoundError, KeyError, ValueError, OSError) as e:
                log.warning(
                    "Cross-dataset matching skipped (%s): %s",
                    type(e).__name__, e,
                )

    # Run report is generated by the post-pipeline viz orchestrator so it
    # indexes figures the orchestrator has already written.
    return {  # multi-sample metadata for the post-pipeline viz orchestrator
        "is_single_sample": False,
        "unique_ds_ids": unique_ds_ids,
        "bed_paths": all_pos_beds + all_neg_beds,
        "atlas_enabled": bool(directory_config.atlas),
        "clustering_dir": directory_config.clustering_dir,
        "single_sample_h5ad": None,
    }


def main() -> None:
    """Legacy entry point: runs the pipeline using the argparse-populated globals.

    Delegates to :func:`_run_pipeline_body` with no progress bars.
    Retained for backward compatibility (``if __name__ == '__main__'`` and any
    direct callers that have not yet migrated to :func:`run`).
    """
    _run_pipeline_body(progress=None)


if __name__ == "__main__":
    main()
