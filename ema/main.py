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

try:
    from ema.datasets.atlas_annotate import annotate_pas_against_atlas
except ImportError:
    annotate_pas_against_atlas = None  # type: ignore[assignment]


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
        # D9: fold atlas-match / internal-priming counts into entity_counts
        # so the hub can read the QC funnel straight off the manifest.
        _entity_counts: dict = {}
        if _pipeline_result is not None:
            _atlas_stats = _pipeline_result.get("atlas_stats") or {}
            for _k in ("n_atlas_matched", "n_atlas_unmatched", "atlas_match_rate"):
                if _k in _atlas_stats:
                    _entity_counts[_k] = _atlas_stats[_k]
            _ip_stats = _pipeline_result.get("ip_stats") or {}
            for _k in ("n_ip_flagged", "ip_flag_rate"):
                if _ip_stats.get(_k) is not None:
                    _entity_counts[_k] = _ip_stats[_k]
        _mgr.write_manifest(build_resolved_run_config(), entity_counts=_entity_counts)
    except Exception as _e:  # never fail a run over the manifest
        log.warning("run_manifest.json write failed: %s", _e)

    return 0


# ---------------------------------------------------------------------------
# D6: internal-priming (--ip-filter) / annotation (--annot-filter) filters.
#
# The filter LOGIC lives in ema.experimental.peak_filters / internal_priming
# (pre-existing, not modified here). These two helpers are the wiring seam
# that makes `ema run` actually invoke it:
#
#   * _validate_pas_filter_config() — fail loud, before any pipeline compute,
#     if a filter is enabled but its required input is missing.
#   * _apply_pas_filters(output_mgr) — apply the enabled filter(s) to
#     directory_config.posbed / .negbed *in place*, at the one point common
#     to both the single-sample and multi-sample paths: right after those
#     BEDs hold their final pre-annotation PAS set and right before
#     find_close() consumes them to build the gene-assignment DataFrame
#     that annotate() joins the count matrix against. Filtering here means
#     a dropped PAS is absent from `genes.index`, so annotate()'s
#     `genes.index.intersection(pas_ids)` naturally drops it from the
#     count matrix too — no other file needs touching.
#
# Both flags default to False, and this module is only imported when at
# least one is enabled, so the default-off pipeline never enters this code
# path at all — output is byte-identical to pre-D6 behaviour.
# ---------------------------------------------------------------------------
def _write_run_support(bed_paths, out_path) -> None:
    """Concatenate the caller BEDs' poly(A) support sidecars into one
    run-root ``pas_support.tsv`` (pas_id, clip_reads, clip_umis,
    clip_reads_f3844, clip_umis_f3844, window_reads, tier -- plus the
    ``--pas-features on`` call-time columns, whose presence is read off the
    first source file's own header rather than assumed).

    Best-effort and non-fatal: the sidecar is an annotation, never an input
    to the pipeline.  Only written when the pas_ids are still the caller's
    (single-BAM path); multi-BAM runs re-key their PAS in merge_pas_beds,
    so their sidecars stay next to each caller BED.
    """
    from ema.countmatrix.paswrite import SUPPORT_COLUMNS, support_path_for
    try:
        srcs = [Path(support_path_for(b)) for b in bed_paths]
        srcs = [s for s in srcs if s.exists()]
        if not srcs:
            return
        # The header must come from the CALLER's own sidecars, not from this
        # module's compile-time SUPPORT_COLUMNS: --pas-features on appends
        # columns to every row, and a hardcoded seven-column header over
        # nine-column rows silently mislabels every field after `tier`.
        # (Found on the real chr19+21 slice, not by a unit test -- the merge
        # had no coverage at all.  It now does: see
        # tests/test_pas_features.py::test_the_run_root_sidecar_keeps_the_callers_header.)
        header = "\t".join(SUPPORT_COLUMNS) + "\n"
        with open(srcs[0]) as fh:
            first = fh.readline()
        if first.startswith("pas_id"):
            header = first if first.endswith("\n") else first + "\n"
        with open(out_path, "w") as out:
            out.write(header)
            for s in srcs:
                with open(s) as fh:
                    first = fh.readline()
                    if not first.startswith("pas_id") and first.strip():
                        out.write(first)
                    for line in fh:
                        if line.strip():
                            out.write(line)
    except Exception as exc:  # pragma: no cover - annotation only
        log.warning("poly(A) support sidecar not written: %s", exc)


def _pas_features_enabled() -> bool:
    """``--pas-features on`` (peakAtail-prime).

    Read from ``variable_config`` and not from ``args``: the argparse
    namespace carries schema DEFAULTS for anything Click resolved (bug B0's
    class of defect), while ``variable_config`` is what the caller itself
    read when it wrote the sidecar's header.
    """
    from ema.config import variable_config as _vc

    return str(getattr(_vc, "pas_features", "off")).lower() == "on"


def _collect_pas_features_standalone(collector, genome_fasta) -> None:
    """Fill *collector* when the internal-priming filter is NOT running.

    With a genome FASTA this performs the SAME single pass the IP filter
    would have performed (``mode="annotate"``, ``output_path=None``: it drops
    nothing and writes no BED) -- so a run never makes two passes over the
    FASTA, whichever way ``--ip-filter`` is set.  Without one, only the
    BED-derived context columns are real and every sequence column is ``NA``.
    """
    import os

    from ema.countmatrix.pas_features import collect_from_bed

    beds = [str(directory_config.posbed), str(directory_config.negbed)]
    beds = [b for b in beds if os.path.exists(b) and os.path.getsize(b) > 0]
    if genome_fasta and os.path.exists(genome_fasta):
        from ema.experimental.internal_priming import filter_internal_priming

        for bed in beds:
            filter_internal_priming(
                bed, genome_fasta, None,
                window_left=getattr(args, "ip_window_left", 10),
                window_right=getattr(args, "ip_window_right", 30),
                a_stretch=getattr(args, "ip_a_stretch", 6),
                a_fraction=getattr(args, "ip_a_fraction", 0.7),
                mode="annotate", features=collector,
            )
    else:
        if beds:
            log.warning(
                "--pas-features on but no --genome-fasta: the sequence "
                "columns of pas_support.tsv (downstream A-content, hexamer, "
                "ip_tool_*) will be NA. Supply --genome-fasta to fill them."
            )
        for bed in beds:
            collect_from_bed(bed, collector)


def _pas_score_mode() -> str:
    """``--pas-score`` (peakAtail-prime), read from ``variable_config``.

    Not from ``args``: that namespace carries schema DEFAULTS for anything the
    Click layer resolved (bug B0's class of defect), while ``variable_config``
    is what this seam actually acts on.
    """
    from ema.config import variable_config as _vc

    return str(getattr(_vc, "pas_score", "none")).lower()


def _apply_pas_score(collector, stats: dict, output_mgr) -> None:
    """Compute the calibrated per-site score, and optionally SELECT on it.

    Runs at the same seam as the internal-priming veto, immediately after it,
    because that is where the per-site features exist and where a dropped PAS
    still disappears cleanly from ``genes.index`` (and therefore from the count
    matrix).  Three things happen here and nothing else:

    1. every candidate's probability is computed from the shipped constants
       (:mod:`ema.countmatrix.pas_score`; numpy only, no scikit-learn);
    2. the probabilities are handed to :func:`_write_pas_features`, which
       appends them to ``pas_support.tsv`` as one more column;
    3. **only** in ``--pas-score select``: a **tier-1** PAS whose probability is
       below the threshold is dropped from the pos/neg BEDs.

    THE HARD GATES STAY HARD.  The score is a re-ranker INSIDE tier-1 and
    inside the internal-priming veto -- it replaces the ">= 2 clip molecules"
    threshold and nothing else.  It can only ever REMOVE a tier-1 candidate:
    it cannot rescue an internally-primed one (the veto has already dropped it
    when ``--ip-filter --ip-filter-mode filter`` ran) and it cannot promote a
    coverage-only tier-2 one.  Every configuration in which a score was allowed
    to override the veto looked spectacular on the curated atlas and no better
    on long reads (`results/algo_headroom/VERIFY/` §1.4, §7.1).

    A PAS whose sequence window could not be read (``seq_ok != 1``: a contig
    missing from the FASTA, or a candidate too close to a contig edge) gets
    ``pas_score`` ``NA`` and is EXEMPT from the gate rather than silently
    dropped -- "we could not score it" must not be spelled the same way as
    "we scored it and it lost".
    """
    import os

    import numpy as np

    mode = _pas_score_mode()
    if mode == "none":
        return
    from ema.countmatrix.pas_features import SEAM_FEATURE_COLUMNS
    from ema.countmatrix.pas_score import load_model

    if collector is None:
        raise ValueError(
            "--pas-score requires --pas-features on: the score is computed "
            "from the per-site feature columns and cannot be evaluated "
            "without them.")

    support = Path(directory_config.output_dir) / "pas_support.tsv"
    if not support.exists():
        # Multi-BAM runs re-key their PAS ids in merge_pas_beds, so there is no
        # run-root sidecar to read the caller columns from.  Annotating is then
        # simply skipped -- but SELECTING must not be, because a gate that
        # silently does not apply produces a call set the run record claims was
        # filtered.  Fail instead.
        if mode == "select":
            raise ValueError(
                "--pas-score select needs the run-root pas_support.tsv for the "
                "caller's feature columns, and this run has none (multi-BAM "
                "runs re-key their PAS ids at merge time). Refusing to report a "
                "score-selected call set that was never selected.")
        log.warning(
            "--pas-score: no run-root pas_support.tsv (multi-BAM runs re-key "
            "their PAS ids in merge_pas_beds); no score computed.")
        return

    from ema.config import variable_config as _vc

    model = load_model(getattr(_vc, "pas_score_model", None) or "prime1")
    feats = collector.finish()

    import pandas as pd

    caller_cols = [c for c in model.features if c not in SEAM_FEATURE_COLUMNS]
    seam_cols = [c for c in model.features if c in SEAM_FEATURE_COLUMNS]
    sup = pd.read_csv(support, sep="\t", usecols=["pas_id"] + caller_cols,
                      dtype={"pas_id": str}, na_values=["NA"])
    n = len(sup)
    ids = sup.pas_id.values

    # The seam block is stored as tab-joined TEXT (one string per candidate, so
    # a genome-wide run does not hold a 22-key dict each).  Pull only the
    # columns this model names, straight into preallocated arrays -- building a
    # 650k x 22 object DataFrame first would cost most of a gigabyte for values
    # that are about to become floats.
    take = [(SEAM_FEATURE_COLUMNS.index(c), c) for c in seam_cols]
    cols = {c: np.full(n, np.nan) for _, c in take}
    ok = np.zeros(n, dtype=bool)
    seq_ok_at = SEAM_FEATURE_COLUMNS.index("seq_ok")
    n_missing = 0
    for i, pid in enumerate(ids):
        txt = feats.get(pid)
        if txt is None:
            n_missing += 1
            continue
        parts = txt.split("\t")
        ok[i] = parts[seq_ok_at] == "1"
        for j, c in take:
            v = parts[j]
            if v != "NA":
                cols[c][i] = float(v)
    for c in caller_cols:
        cols[c] = pd.to_numeric(sup[c], errors="coerce").values.astype(float)

    # Transforms are applied over the WHOLE candidate population and only then
    # subset to the scoreable rows.  It matters for any within-run statistic
    # (the `rankpct` transform): a percentile computed over a subset is a
    # different number from the one the offline fitter computed over the run.
    X = model.build_matrix(cols)
    prob = np.full(n, np.nan)
    if ok.any():
        prob[ok] = model.predict_proba_matrix(X[ok])
    scores = {str(pid): prob[i] for i, pid in enumerate(ids)}

    thr = getattr(_vc, "pas_score_min", -1.0)
    thr = model.threshold if thr is None or float(thr) < 0 else float(thr)
    st = {"mode": mode, "model": model.name, "threshold": float(thr),
          "n_scored": int(ok.sum()), "n_unscoreable": int((~ok).sum()),
          "n_missing_features": int(n_missing),
          "mean_prob": float(np.nanmean(prob)) if ok.any() else None}

    if mode == "select":
        dropped = kept = 0
        for bed_path in (str(directory_config.posbed), str(directory_config.negbed)):
            if not os.path.exists(bed_path) or os.path.getsize(bed_path) == 0:
                continue
            tmp = bed_path + ".pas_score_tmp"
            with open(bed_path) as src, open(tmp, "w") as dst:
                for line in src:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 6:
                        dst.write(line)
                        continue
                    try:
                        tier1 = float(parts[4]) > 0
                    except ValueError:
                        tier1 = False
                    p = scores.get(parts[3])
                    if tier1 and p is not None and p == p and p < thr:
                        dropped += 1
                        continue
                    dst.write(line)
                    kept += 1
            os.replace(tmp, bed_path)
        st.update({"n_dropped": dropped, "n_kept": kept})
        log.info("pas_score: select at p >= %.6f (model %s): kept %d, dropped %d "
                 "tier-1 candidates", thr, model.name, kept, dropped)
    else:
        log.info("pas_score: %d candidates scored with model %s "
                 "(threshold %.6f is RECORDED, not applied -- use "
                 "--pas-score select to apply it)", int(ok.sum()), model.name, thr)

    stats["pas_score"] = scores
    stats["pas_score_stats"] = st
    try:
        with open(output_mgr.path("pas_gene", "pas_score_stats.json"), "w") as fh:
            json.dump(st, fh, indent=2)
    except Exception as exc:  # pragma: no cover - bookkeeping only
        log.warning("pas_score: stats write failed: %s", exc)


def _write_pas_features(pas_filter_result) -> None:
    """Append the seam's feature columns to the run-root sidecar.

    APPEND ONLY: ``pas_support.tsv`` keeps every column it already had, in
    place, and gains :data:`~ema.countmatrix.pas_features.SEAM_FEATURE_COLUMNS`
    on the end.  Best-effort and non-fatal, exactly like the sidecar itself:
    the features are an annotation and nothing in the pipeline reads them.

    Multi-BAM runs re-key their PAS ids in ``merge_pas_beds``, so there is no
    run-root ``pas_support.tsv`` to extend; those runs get a standalone
    ``pas_features.tsv`` in the merged id space instead.
    """
    collector = (pas_filter_result or {}).get("pas_features")
    if collector is None:
        return
    try:
        from ema.countmatrix.pas_features import (
            append_columns, write_features_tsv,
        )

        from ema.countmatrix.pas_features import SEAM_FEATURE_COLUMNS

        feats = collector.finish()
        columns = SEAM_FEATURE_COLUMNS
        # --pas-score: the probability is one more column on the SAME append,
        # folded into the collector's own tab-joined text exactly the way
        # finish() folds the context block in -- so the sidecar is rewritten
        # once, not twice, whether or not the score was asked for.
        scores = (pas_filter_result or {}).get("pas_score")
        if scores:
            from ema.countmatrix.pas_score import ScoredFeatures

            columns = SEAM_FEATURE_COLUMNS + ("pas_score",)
            feats = ScoredFeatures(feats, scores)
        root = Path(directory_config.output_dir)
        support = root / "pas_support.tsv"
        if support.exists():
            n = append_columns(support, feats, columns)
            log.info(
                "pas_support.tsv: %d feature columns appended to %d rows (%s)",
                len(columns), n, collector.stats(),
            )
        else:
            out = root / "pas_features.tsv"
            n = write_features_tsv(out, feats, columns)
            log.info("pas_features.tsv: %d rows written (%s)", n, collector.stats())
    except Exception as exc:  # pragma: no cover - annotation only
        log.warning("pas_support.tsv feature columns not written: %s", exc)


def _resolve_annotation_bed() -> str:
    """Return the annotation BED source for --annot-filter.

    Uses the explicit ``--annotation-bed`` override when given, otherwise
    falls back to the GTF-derived gene BED that GTF preprocessing already
    writes to ``directory_config.endbed`` (see ema/annotate/gtftobed.py).
    """
    override = getattr(args, "annotation_bed", None)
    if override:
        return str(override)
    return str(directory_config.endbed)


def _validate_pas_score_config() -> None:
    """Fail BEFORE the run when --pas-score is misconfigured.

    A score computed from NA sequence columns would be a number that looks
    exactly like a real one in the sidecar, so the two things it cannot do
    without are checked here rather than warned about later.
    """
    import os

    mode = _pas_score_mode()
    if mode == "none":
        return
    if not _pas_features_enabled():
        raise ValueError(
            "--pas-score requires --pas-features on: the score is computed "
            "from the per-site feature columns.")
    genome_fasta = getattr(args, "genome_fasta", None)
    if not genome_fasta or not os.path.exists(genome_fasta):
        raise ValueError(
            "--pas-score requires --genome-fasta: the model uses the "
            "canonical-hexamer columns, which are NA without a genome "
            f"(got: {genome_fasta!r}).")
    from ema.countmatrix.pas_score import load_model

    model = load_model(getattr(variable_config, "pas_score_model", None) or "prime1")
    if mode == "select":
        thr = getattr(variable_config, "pas_score_min", -1.0)
        thr = model.threshold if thr is None or float(thr) < 0 else float(thr)
        if not (0.0 <= float(thr) <= 1.0):
            raise ValueError(
                f"--pas-score select needs a probability threshold in [0, 1]; "
                f"model {model.name!r} carries {model.threshold!r} and "
                f"--pas-score-min resolved to {thr!r}.")


def _validate_pas_filter_config() -> None:
    """Raise a clear error BEFORE running if --ip-filter / --annot-filter
    are enabled but misconfigured. Never silently skips a requested filter.
    """
    import os

    ip_filter = bool(getattr(args, "ip_filter", False))
    annot_filter = bool(getattr(args, "annot_filter", False))
    if not ip_filter and not annot_filter:
        return

    if ip_filter:
        genome_fasta = getattr(args, "genome_fasta", None)
        if not genome_fasta or not os.path.exists(genome_fasta):
            raise ValueError(
                "--ip-filter is enabled but --genome-fasta is missing or does "
                f"not exist (got: {genome_fasta!r}). Provide a genome FASTA "
                "indexed with pyfaidx (.fai) via --genome-fasta / YAML "
                "`genome_fasta:`."
            )

    if annot_filter:
        annotation_bed_override = getattr(args, "annotation_bed", None)
        if annotation_bed_override:
            if not os.path.exists(annotation_bed_override):
                raise ValueError(
                    "--annot-filter is enabled but --annotation-bed does not "
                    f"exist: {annotation_bed_override!r}."
                )
        elif not directory_config.gtf_dir or not os.path.exists(str(directory_config.gtf_dir)):
            # No explicit override and no GTF to derive one from at
            # find_close() time -- the GTF-derived endbed cannot be produced.
            raise ValueError(
                "--annot-filter is enabled but no --annotation-bed was "
                "given and no --gtf is set to derive one from. Provide "
                "--annotation-bed or --gtf."
            )


def _apply_polya_gate(output_mgr) -> dict | None:
    """Enforce ``--polya-mode filter|require`` on the pos/neg PAS BEDs, in place.

    Stage 1 seam: runs at the exact point --ip-filter runs (right before
    find_close() derives the gene assignment), so a dropped PAS disappears
    from ``genes.index`` and annotate() naturally drops it from the count
    matrix — no other file needs touching.

    ``filter`` drops every PAS whose BED score column (poly(A) clip
    support in distinct molecules, written by peak calling when
    --polya-evidence is on) is 0 —
    i.e. exactly the coverage-only tier.  ``require`` does the same and then
    FAILS the run when nothing clip-supported remains (the R2 chemistry
    failure mode: poly(A) trimmed upstream leaves support at 0 everywhere).

    No-op returning ``None`` when --polya-evidence is off or --polya-mode is
    'annotate' (the default), so default runs never enter this code path.
    """
    import os

    polya_on = str(getattr(args, "polya_evidence", "on")).lower() != "off"
    polya_mode = str(getattr(args, "polya_mode", "annotate"))
    if not polya_on or polya_mode not in ("filter", "require"):
        return None

    stats: dict = {"polya_mode": polya_mode}
    total_kept = 0
    total_seen = 0
    for strand_label, bed_path in (
        ("pos", directory_config.posbed),
        ("neg", directory_config.negbed),
    ):
        bed_path = str(bed_path)
        if not os.path.exists(bed_path) or os.path.getsize(bed_path) == 0:
            stats[strand_label] = {"total": 0, "kept": 0, "dropped": 0}
            continue
        tmp_out = bed_path + ".polya_gate_tmp"
        total = kept = 0
        with open(bed_path) as src, open(tmp_out, "w") as dst:
            for line in src:
                if not line.strip():
                    continue
                total += 1
                parts = line.rstrip("\n").split("\t")
                try:
                    support = int(float(parts[4]))
                except (IndexError, ValueError):
                    support = 0
                if support >= 1:
                    dst.write(line)
                    kept += 1
        os.replace(tmp_out, bed_path)
        stats[strand_label] = {
            "total": total, "kept": kept, "dropped": total - kept,
        }
        total_seen += total
        total_kept += kept

    if polya_mode == "require" and total_kept == 0:
        raise RuntimeError(
            "--polya-mode require: none of the "
            f"{total_seen} PAS carry any poly(A) clip support (BED score "
            "column all 0). Either the BAM's poly(A) evidence was destroyed "
            "upstream (trimming / chemistry — see the startup clip-rate "
            "warning) or peak calling ran with --polya-evidence off. "
            "Refusing to ship an unsupported call set."
        )

    log.info(
        "polya_gate: mode=%s pos(kept=%s/%s) neg(kept=%s/%s)",
        polya_mode,
        stats["pos"]["kept"], stats["pos"]["total"],
        stats["neg"]["kept"], stats["neg"]["total"],
    )
    try:
        _stats_path = output_mgr.path("pas_gene", "polya_gate_stats.json")
        with open(_stats_path, "w") as _f:
            json.dump(stats, _f, indent=2)
    except Exception as e:  # bookkeeping only — never fail the run
        log.warning("polya_gate: stats write failed: %s", e)
    return stats


def _apply_pas_filters(output_mgr) -> dict | None:
    """Apply the enabled PAS filter(s) to the pos/neg PAS BEDs, in place.

    No-op — ema.experimental.peak_filters is not even imported — when both
    --ip-filter and --annot-filter are off, so default-off runs take
    exactly the pre-D6 code path. Returns ``None`` in that case.

    D9: --ip-filter defaults to ``ip_filter_mode="annotate"`` (keep every
    PAS, flag it) rather than dropping it -- an internally-primed peak is
    candidate alternative-PAS signal, not noise, for a scientist hunting
    APA. ``--ip-filter-mode filter`` restores the pre-D9 drop behaviour.
    (``--annot-filter``, the gene-region membership filter, is a distinct
    concept and is unaffected -- it always drops non-overlapping peaks.)

    When run, the returned dict includes ``"ip_of"`` -- the merged
    ``{pas_id: internal_priming_bool}`` map across the pos+neg BEDs, for
    callers to thread into :func:`ema.outputs.write_pas_gene_artifacts` /
    :func:`ema.provenance.record_pas_drops`.
    """
    import os

    # Stage 1: the poly(A) support gate shares this seam (same point, same
    # in-place BED rewrite) and runs FIRST so --ip-filter statistics are
    # computed on the clip-supported set when both are enabled.
    polya_stats = _apply_polya_gate(output_mgr)

    ip_filter = bool(getattr(args, "ip_filter", False))
    annot_filter = bool(getattr(args, "annot_filter", False))

    # peakAtail-prime --pas-features: the per-site scoring covariates are
    # collected HERE, at the internal-priming seam, because that is where the
    # genome is already open and the full candidate set is still on disk.
    # When --ip-filter is on they ride inside its pass; when it is off the
    # collector's own pass is the only one.  Either way: one pass, never two.
    collector = None
    if _pas_features_enabled():
        from ema.countmatrix.pas_features import FeatureCollector

        collector = FeatureCollector()

    if not ip_filter and not annot_filter:
        _out: dict = {"polya_gate": polya_stats} if polya_stats else {}
        if collector is not None:
            _collect_pas_features_standalone(
                collector, getattr(args, "genome_fasta", None))
            _out["pas_features"] = collector
            _apply_pas_score(collector, _out, output_mgr)
        return _out or None

    ip_mode = str(getattr(args, "ip_filter_mode", "annotate"))

    from ema.experimental.peak_filters import apply_filters

    genome_fasta = getattr(args, "genome_fasta", None)
    annotation_bed = _resolve_annotation_bed() if annot_filter else None

    # Re-validate defensively: _validate_pas_filter_config() already ran at
    # the top of the pipeline, but this function may also be called directly
    # (e.g. from tests) without that pre-flight check having run first.
    if ip_filter and (not genome_fasta or not os.path.exists(genome_fasta)):
        raise ValueError(
            "--ip-filter is enabled but --genome-fasta is missing or does "
            f"not exist (got: {genome_fasta!r})."
        )
    if annot_filter and not os.path.exists(annotation_bed):
        raise ValueError(
            "--annot-filter is enabled but the resolved annotation BED does "
            f"not exist: {annotation_bed!r} (provide --annotation-bed or --gtf)."
        )

    combined_stats: dict = {
        "ip_filter": ip_filter,
        "ip_filter_mode": ip_mode if ip_filter else None,
        "annot_filter": annot_filter,
        "genome_fasta": genome_fasta if ip_filter else None,
        "annotation_bed": annotation_bed if annot_filter else None,
    }
    if polya_stats:
        combined_stats["polya_gate"] = polya_stats
    filtered_any = False
    ip_of: dict = {}  # D9: merged {pas_id: internal_priming_bool} across pos+neg
    for strand_label, bed_path in (
        ("pos", directory_config.posbed),
        ("neg", directory_config.negbed),
    ):
        bed_path = str(bed_path)
        if not os.path.exists(bed_path) or os.path.getsize(bed_path) == 0:
            # Nothing to filter for this strand (e.g. all-single-strand data).
            combined_stats[strand_label] = {"total": 0, "passed": 0, "filtered": 0}
            continue
        tmp_out = bed_path + ".pas_filter_tmp"
        stats = apply_filters(
            input_bed=bed_path,
            output_bed=tmp_out,
            genome_fasta=genome_fasta,
            annotation_bed=annotation_bed,
            enable_internal_priming=ip_filter,
            enable_annotation_filter=annot_filter,
            ip_window_left=getattr(args, "ip_window_left", 10),
            ip_window_right=getattr(args, "ip_window_right", 30),
            ip_a_stretch=getattr(args, "ip_a_stretch", 6),
            ip_a_fraction=getattr(args, "ip_a_fraction", 0.7),
            ip_mode=ip_mode,
            features=collector if ip_filter else None,
        )
        os.replace(tmp_out, bed_path)
        combined_stats[strand_label] = stats
        filtered_any = True
        ip_of.update(stats.get("internal_priming_flags") or {})

    if collector is not None:
        if not ip_filter:
            # --annot-filter alone does not open the genome, so the feature
            # pass still has to be made (or the BED-only fallback taken).
            _collect_pas_features_standalone(
                collector, getattr(args, "genome_fasta", None))
        combined_stats["pas_features"] = collector
        # The score runs AFTER the internal-priming veto, on the BEDs the veto
        # has already rewritten, so `--pas-score select` can only ever remove a
        # tier-1 candidate the veto let through -- never rescue one it dropped.
        _apply_pas_score(collector, combined_stats, output_mgr)

    if ip_filter:
        n_ip_flagged = sum(1 for v in ip_of.values() if v)
        combined_stats["n_ip_flagged"] = n_ip_flagged
        combined_stats["ip_flag_rate"] = (
            round(n_ip_flagged / len(ip_of), 4) if ip_of else 0.0
        )
        combined_stats["ip_of"] = ip_of

    # "peak_filters" is not one of OutputManager's pre-registered numbered
    # stage dirs (self.dirs), so save_stats("peak_filters", ...) would raise
    # KeyError. Use the existing pas_gene (04_pas_gene_assignment) stage dir
    # instead -- filtering happens immediately before that stage -- via a
    # distinctly-named file so it can't collide with pas_gene_stats.json.
    # Redact the per-PAS flag maps from the ON-DISK stats (they can be as
    # large as the PAS count and belong in the pasbed/ledger, not a JSON
    # blob) -- the full maps are still on the RETURNED dict for the caller.
    _disk_stats = {k: v for k, v in combined_stats.items()
                   if k not in ("ip_of", "pas_features")}
    for _strand_label in ("pos", "neg"):
        _sd = _disk_stats.get(_strand_label)
        if isinstance(_sd, dict):
            _sd = {k: v for k, v in _sd.items() if k != "internal_priming_flags"}
            if isinstance(_sd.get("internal_priming"), dict):
                _sd["internal_priming"] = {
                    k: v for k, v in _sd["internal_priming"].items() if k != "flags"
                }
            _disk_stats[_strand_label] = _sd
    _stats_path = output_mgr.path("pas_gene", "peak_filters_stats.json")
    with open(_stats_path, "w") as _f:
        json.dump(_disk_stats, _f, indent=2)
    # Register the rewritten pasbeds as filtered artifacts (best-effort --
    # never fail the run over manifest bookkeeping).
    if filtered_any:
        try:
            output_mgr.register_artifact(
                str(directory_config.posbed), stage="peak_filters",
                fmt="bed", schema_name="pas_bed_filtered",
            )
            output_mgr.register_artifact(
                str(directory_config.negbed), stage="peak_filters",
                fmt="bed", schema_name="pas_bed_filtered",
            )
        except Exception as e:  # pragma: no cover -- bookkeeping only
            log.warning("peak_filters: register_artifact failed: %s", e)
    log.info(
        "peak_filters: ip_filter=%s(mode=%s) annot_filter=%s pos(total=%s filtered=%s) "
        "neg(total=%s filtered=%s) n_ip_flagged=%s",
        ip_filter, ip_mode, annot_filter,
        combined_stats["pos"].get("total"), combined_stats["pos"].get("filtered"),
        combined_stats["neg"].get("total"), combined_stats["neg"].get("filtered"),
        combined_stats.get("n_ip_flagged"),
    )
    return combined_stats


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

    # D6: fail loud, before any peak-calling compute is spent, if --ip-filter
    # / --annot-filter are enabled but misconfigured (see _apply_pas_filters).
    _validate_pas_filter_config()
    _validate_pas_score_config()

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

    # Stage 1: read-level poly(A) evidence settings (see
    # ema/countmatrix/polya.py).  Threaded identically into the monolithic,
    # tile and pipeline paths so all three agree.
    from ema.countmatrix.polya import parse_count_window as _parse_count_window
    import sys as _sys
    if any(a.split("=", 1)[0] == "--polya-min-reads" for a in _sys.argv):
        log.warning(
            "--polya-min-reads is DEPRECATED and now spelled --polya-min-umis: "
            "the gate has always counted distinct (CB, UMI) molecules, and BED "
            "column 5 now reports the same unit (raw clip reads are in "
            "pas_support.tsv)."
        )
    _polya_kwargs = dict(
        polya_enabled=(str(getattr(args, "polya_evidence", "on")).lower() != "off"),
        polya_min_clip=int(getattr(args, "polya_min_clip", 6)),
        polya_min_purity=float(getattr(args, "polya_min_purity", 0.8)),
        polya_window=int(getattr(args, "polya_window", 100)),
        polya_seed_window=int(getattr(args, "polya_seed_window", 25)),
        # --polya-min-umis, with the deprecated --polya-min-reads spelling
        # mapping to the same (always molecule-based) gate.
        polya_min_umis=int(
            getattr(args, "polya_min_umis", None)
            if getattr(args, "polya_min_umis", None) is not None
            else getattr(args, "polya_min_reads", 1)
        ),
        polya_clip_filter=str(getattr(args, "polya_clip_filter", "none")),
        polya_count_window=_parse_count_window(
            getattr(args, "polya_count_window", "auto,25")
        ),
    )

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
        **_polya_kwargs,
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
            **_polya_kwargs,
        )

        # Stage 1: once-per-BAM clip-rate QC — tile workers use region
        # fetches and never full-scan, so the check runs here instead.
        if _polya_kwargs["polya_enabled"]:
            from ema.countmatrix.polya import check_clip_rate
            for _ds_id, _bam_path in {(d, str(b)) for d, b in bam_list}:
                check_clip_rate(
                    str(_bam_path),
                    min_clip=_polya_kwargs["polya_min_clip"],
                    min_purity=_polya_kwargs["polya_min_purity"],
                    barcode_tag=variable_config.barcode_tag or "CB",
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
        # PER-BAM MODE (default): one dataset/BAM at a time.  Peak calling
        # itself runs one spawned job per (contig, strand) with a
        # deterministic merge (ema/countmatrix/chrom_parallel.py) -- the
        # outputs are byte-identical to the legacy two-pass single-process
        # loop, which remains the fallback when the BAM has no index or
        # --peak-workers 1 is given.
        # -----------------------------------------------------------------
        from ema.countmatrix.chrom_parallel import (
            PER_WORKER_MB,
            can_run_parallel,
            load_index_from_cb_list,
            run_chrom_parallel,
        )

        _rm = get_resource_manager()
        _peak_workers_arg = getattr(args, "peak_workers", None)
        if _peak_workers_arg is not None and int(_peak_workers_arg) >= 1:
            _chrom_workers = int(_peak_workers_arg)
        else:
            _chrom_workers = _rm.get_n_jobs(
                per_worker_mb=PER_WORKER_MB, stage="peak_chrom"
            )
        # BGZF decompression threads are budgeted per worker so that
        # workers x threads never exceeds --threads (or the physical cores).
        _thread_budget = _rm.user_max_threads or _rm.physical_cpu_count()

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

            _peak_client = _client(_peak_stage)
            _parallel_ok = _chrom_workers > 1 and can_run_parallel(str(bam_path))
            if _chrom_workers > 1 and not _parallel_ok:
                log.warning(
                    "%s has no usable BAM index (.bai with mapped-read "
                    "statistics); peak calling falls back to the single-process "
                    "two-pass caller. Run `samtools index` to enable "
                    "per-chromosome parallelism.", bam_path,
                )

            if _parallel_ok:
                from ema.logging_config import get_log_queue
                _summary = run_chrom_parallel(
                    dataset_id=dataset_id,
                    bam_path=str(bam_path),
                    pos_bed=str(pos_bed), neg_bed=str(neg_bed),
                    pos_mtx=str(pos_mtx), neg_mtx=str(neg_mtx),
                    cb_tsv=str(cb_tsv),
                    n_workers=_chrom_workers,
                    strategy_name=args.strategy,
                    strategy_kwargs=_filtered_kwargs,
                    peak_kwargs=peak_kwargs,
                    workdir=str(peakcalling_dir / f"_jobs_{dataset_id}_{idx}"),
                    thread_budget=_thread_budget,
                    progress_client=_peak_client,
                    log_queue=get_log_queue(),
                )
                # The CB filter reads get_mapping() when no explicit list is
                # passed: make the singleton hold the merged column order.
                load_index_from_cb_list(_summary.pop("cb_list"))
                try:
                    with open(peakcalling_dir / f"{dataset_id}_{idx}.peak_jobs.json", "w") as _pj:
                        json.dump(_summary, _pj, indent=2)
                except Exception as _exc:  # pragma: no cover - bookkeeping only
                    log.warning("peak job timings not written: %s", _exc)
            else:
                # Pass the progress client to the pos-strand call so the bar
                # shows per-chromosome progress.  set_total() fires once the
                # BAM is opened; advance() fires on each chromosome boundary.
                # Both strand passes share the same client -- peak_calling()
                # sets total to `nonempty_refs * 2` so the advances from both
                # passes fill the bar.
                peak_calling(
                    False,
                    bedfilepath=str(pos_bed),
                    matrixpath=str(pos_mtx),
                    bamfile_dir=str(bam_path),
                    default_threshold=variable_config.default_threshold,
                    merge_len=variable_config.merge_len,
                    progress_client=_peak_client,
                    **peak_kwargs,
                )
                peak_calling(
                    True,
                    bedfilepath=str(neg_bed),
                    matrixpath=str(neg_mtx),
                    bamfile_dir=str(bam_path),
                    default_threshold=variable_config.default_threshold,
                    merge_len=variable_config.merge_len,
                    progress_client=_peak_client,
                    **peak_kwargs,
                )

                # Dump CB list -- shared by pos+neg (both used the same _index instance)
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
        **_polya_kwargs,
        "polya_mode": str(getattr(args, "polya_mode", "annotate")),
    })

    # ---- 3' cleavage-site offset correction (issue #72) ------------------
    # Called peak 3' ends stop ~90-105 nt short of the true cleavage site
    # (10x R2 coverage runs out before the poly(A) junction).  When the
    # opt-in --cleavage-offset flag is > 0 we shift each reported PAS 3' end
    # downstream, in place, on the per-(dataset,bam) strand BEDs produced by
    # BOTH the tiled and sequential paths.  Doing it here -- the single point
    # where all_pos_beds/all_neg_beds are finalised and before any snapshot,
    # legacy copy, find_close() or annotatedpas.bed derives from them -- keeps
    # the correction path-agnostic without threading a parameter through the
    # spawn-based peak-calling workers.  0 (default) is a no-op (legacy).
    _cleavage_offset = int(getattr(variable_config, "cleavage_offset", 0) or 0)
    _auto_offset = bool(getattr(variable_config, "auto_cleavage_offset", False))
    _offset_diag = None
    if _auto_offset:
        # Data-driven mode: infer the offset from the called peaks + FASTA.
        from ema.countmatrix.cleavage_offset import resolve_cleavage_offset

        _genome_fasta = getattr(args, "genome_fasta", None)
        _cleavage_offset, _offset_diag = resolve_cleavage_offset(
            _cleavage_offset,
            auto=True,
            bed_paths=list(all_pos_beds) + list(all_neg_beds),
            genome_fasta=_genome_fasta,
        )
    if _cleavage_offset > 0:
        from ema.countmatrix.cleavage_offset import rewrite_bed_3prime_offset

        # clip_seeded places tier-1 PAS at the observed poly(A) clip site (the
        # true cleavage coordinate, recorded as a >0 BED score); only its
        # coverage-only tier (score 0) carries the R2 read-length offset.
        _skip_supported = str(getattr(args, "strategy", "")) == "clip_seeded"
        _n_shifted = 0
        for _bed in list(all_pos_beds) + list(all_neg_beds):
            try:
                _n_shifted += rewrite_bed_3prime_offset(
                    _bed, _cleavage_offset, skip_supported=_skip_supported
                )
            except FileNotFoundError:
                # A strand may legitimately produce no BED for a dataset.
                continue
        log.info(
            "3' cleavage-offset correction: shifted %d PAS 3' ends downstream "
            "by %d bp (%s)%s", _n_shifted, _cleavage_offset,
            "auto-estimated" if _auto_offset else "--cleavage-offset",
            " [clip-supported tier exempt]" if _skip_supported else "",
        )
        _stats = {
            "cleavage_offset_bp": _cleavage_offset,
            "n_pas_shifted": _n_shifted,
            "auto": _auto_offset,
        }
        if _offset_diag is not None:
            _stats["diagnostics"] = _offset_diag
        # There is no dedicated "cleavage_offset" stage dir, so
        # save_stats("cleavage_offset", ...) KeyErrors on self.dirs (caught in
        # real-run validation, not unit tests). The offset is a peak-calling
        # correction -> persist alongside the peak-calling outputs. Use a local
        # alias: a later `import json` in this function makes bare `json`
        # function-local, so the module-level name is shadowed here.
        import json as _json_coff
        with open(
            output_mgr.path("peak_calling", "cleavage_offset_stats.json"), "w"
        ) as _coff:
            _json_coff.dump(_stats, _coff, indent=2, default=str)

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
        # poly(A) support sidecar: concatenate the two strand files into
        # <run>/pas_support.tsv, keyed by the same pas_id as pasbed.bed.
        _write_run_support(
            [all_pos_beds[0], all_neg_beds[0]],
            Path(directory_config.output_dir) / "pas_support.tsv",
        )
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

        # D6/D9: optional internal-priming / annotation filters, applied to the
        # pos/neg PAS BEDs in place, before find_close() derives the gene
        # assignment from them. No-op (not even imported) when both flags
        # are off — see _apply_pas_filters(). Note: atlas snapping never
        # runs on the single-sample path (only the multi-sample merge/atlas
        # stage does), so atlas_of stays empty here.
        _pas_filter_result = _apply_pas_filters(output_mgr)
        _ip_of: dict = (_pas_filter_result or {}).get("ip_of", {})
        # peakAtail-prime --pas-features: append the seam's columns to the
        # run-root sidecar written above.  Additive; nothing else moves.
        _write_pas_features(_pas_filter_result)

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
            "n_ip_flagged": (_pas_filter_result or {}).get("n_ip_flagged"),
            "ip_flag_rate": (_pas_filter_result or {}).get("ip_flag_rate"),
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
            atlas_of={}, ip_of=_ip_of,
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
            # D9: atlas never runs on the single-sample path.
            "atlas_stats": None,
            "ip_stats": {
                "n_ip_flagged": (_pas_filter_result or {}).get("n_ip_flagged"),
                "ip_flag_rate": (_pas_filter_result or {}).get("ip_flag_rate"),
            } if _pas_filter_result else None,
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

    # Dispatch atlas vs. coordinate-merge based on config.
    #
    # atlas_mode="annotate" (the default whenever --atlas is configured) is a
    # PURE OVERLAY: the unified PAS set is built EXACTLY like a no-atlas run
    # (merge_pas_beds proximity-merges ALL called PAS across datasets; atlas
    # plays no role in merging or id assignment), and atlas match/distance is
    # annotated onto that already-final set afterwards. This is what
    # guarantees a novel PAS seen in two datasets at one locus becomes ONE
    # unified feature regardless of whether an atlas is configured, and keeps
    # a no-atlas run byte-identical (that path never touches atlas code at
    # all).
    #
    # atlas_mode="filter" (opt-in) keeps the pre-existing snap-and-drop
    # behaviour unchanged: atlas snapping decides PAS identity and can drop
    # PAS with no atlas hit, via snap_beds_to_atlas.
    _atlas_stage = _add_stage(
        "Atlas snap" if directory_config.atlas else "PAS merge",
        total=None,  # indeterminate spinner — we don't know peak count yet
    )
    # E3: provenance ledger for the biggest silent drop site (atlas snap,
    # filter mode only — annotate mode never drops a PAS at this stage, so
    # there is nothing to ledger here; atlas_match/atlas_distance_bp still
    # reach the per-dataset ledgers via _atlas_of/record_pas_drops below).
    _prov_ledger = None
    # D9: {new_pas_id: (atlas_match, atlas_distance_bp)}, "" defaults when
    # atlas didn't run — see ema.datasets.atlas_snap's atlas_status.tsv
    # sidecar (filter mode) / ema.datasets.atlas_annotate's overlay sidecar
    # (annotate mode, the default).
    _atlas_of: dict = {}
    _atlas_stats: dict | None = None
    _atlas_mode = str(getattr(args, "atlas_mode", "annotate")) if directory_config.atlas else None

    if directory_config.atlas and _atlas_mode == "filter":
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
            mode=_atlas_mode,
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

        # D9: read the atlas_status.tsv / atlas_stats.json sidecars
        # snap_beds_to_atlas wrote, so downstream steps (per-dataset ledger,
        # annotatedpas.bed, manifest) can report match/no-match too.
        try:
            _status_path = unified_dir / "atlas_status.tsv"
            if _status_path.exists():
                with open(_status_path) as _f:
                    next(_f, None)  # header
                    for _line in _f:
                        _parts = _line.rstrip("\n").split("\t")
                        if len(_parts) != 3:
                            continue
                        _new_id, _match, _dist = _parts
                        _atlas_of[_new_id] = (
                            _match if _match != "" else "",
                            _dist if _dist != "" else "",
                        )
            _atlas_stats_path = unified_dir / "atlas_stats.json"
            if _atlas_stats_path.exists():
                _atlas_stats = json.loads(_atlas_stats_path.read_text())
        except Exception as _e:  # sidecars are auxiliary — never fail the run
            log.warning("atlas status sidecars unreadable (non-fatal): %s", _e)

        if _atlas_stats:
            output_mgr.save_stats("atlas_snap", _atlas_stats)
            log.info(
                "atlas snap (mode=%s): %d matched, %d unmatched (match_rate=%.4f)",
                _atlas_mode,
                _atlas_stats.get("n_atlas_matched", 0),
                _atlas_stats.get("n_atlas_unmatched", 0),
                _atlas_stats.get("atlas_match_rate", 0.0),
            )
    else:
        # Default-off (no atlas configured) AND atlas_mode="annotate" (the
        # default when an atlas IS configured) both build the unified PAS
        # set the same way — merge_pas_beds, with no atlas involvement.
        unified_bed, mapping_path = merge_pas_beds(
            all_beds,
            all_ds_ids,
            output_dir=unified_dir,
            gap=getattr(args, "pas_gap", 100),
            strands=all_strands,
        )
        if directory_config.atlas:
            # atlas_mode="annotate": overlay atlas match/distance onto the
            # ALREADY-FINAL unified PAS set. Pure annotation — never adds,
            # drops, renumbers, or moves a PAS.
            if annotate_pas_against_atlas is None:
                raise RuntimeError("ema.datasets.atlas_annotate is not available")
            try:
                _status_path, _stats_path = annotate_pas_against_atlas(
                    unified_pasbed_path=unified_bed,
                    atlas_bed_path=directory_config.atlas,
                    distance=directory_config.atlas_distance,
                    output_dir=unified_dir,
                )
                with open(_status_path) as _f:
                    next(_f, None)  # header
                    for _line in _f:
                        _parts = _line.rstrip("\n").split("\t")
                        if len(_parts) != 3:
                            continue
                        _new_id, _match, _dist = _parts
                        _atlas_of[_new_id] = (
                            _match if _match != "" else "",
                            _dist if _dist != "" else "",
                        )
                _atlas_stats = json.loads(_stats_path.read_text())
            except Exception as _e:  # overlay is auxiliary — never fail the run
                log.warning("atlas annotate-overlay failed (non-fatal): %s", _e)

            if _atlas_stats:
                output_mgr.save_stats("atlas_snap", _atlas_stats)
                log.info(
                    "atlas annotate-overlay (mode=%s): %d matched, %d unmatched (match_rate=%.4f)",
                    _atlas_mode,
                    _atlas_stats.get("n_atlas_matched", 0),
                    _atlas_stats.get("n_atlas_unmatched", 0),
                    _atlas_stats.get("atlas_match_rate", 0.0),
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

    # D6/D9: optional internal-priming / annotation filters, applied to the
    # unified pos/neg PAS BEDs in place, before find_close() derives the
    # gene assignment from them. No-op (not even imported) when both flags
    # are off — see _apply_pas_filters(). The unified BED name field is
    # already the rekeyed atlas/merge pas_id at this point (see the split
    # above), so the returned "ip_of" map is keyed correctly for
    # write_pas_gene_artifacts / record_pas_drops in each per-dataset worker.
    _pas_filter_result = _apply_pas_filters(output_mgr)
    _ip_of: dict = (_pas_filter_result or {}).get("ip_of", {})
    # peakAtail-prime --pas-features.  The multi-BAM path has no run-root
    # pas_support.tsv to extend (its PAS were re-keyed by merge_pas_beds), so
    # this writes a standalone pas_features.tsv in the merged id space.
    _write_pas_features(_pas_filter_result)

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
        "n_ip_flagged": (_pas_filter_result or {}).get("n_ip_flagged"),
        "ip_flag_rate": (_pas_filter_result or {}).get("ip_flag_rate"),
    })

    # Per-dataset clustering (each dataset_id clusters independently)
    unique_ds_ids: list[str] = list(dict.fromkeys(all_ds_ids))  # dedupe, preserve order

    # ------------------------------------------------------------------ #
    # Phase 6: parallel per-dataset downstream pipeline                   #
    # ------------------------------------------------------------------ #
    # Pre-build (sub_indices, sub_cbs) for every dataset in the parent
    # process — O(n_cells) each, cheap.  Empty datasets are filtered out.
    # NB: `json` is imported at module scope (top of file). A local `import
    # json` here would make the name function-local for the WHOLE function,
    # breaking earlier bare-`json` uses (cleavage-offset write, atlas-snap
    # reads) with UnboundLocalError — so do NOT re-import it locally.
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
        # Exact match on the sample half of "<sample_id>_<barcode>" (split on
        # the LAST underscore -- see ema.countmatrix.indexing.split_cb).
        # `startswith(f"{ds_id}_")` was ambiguous once sample ids may contain
        # "_": dataset "a" would also claim every cell of dataset "a_b".
        sub_indices = [
            i for i, cb in enumerate(all_cb_strings)
            if cb.rsplit("_", 1)[0] == ds_id
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
    # D9: same atlas_of/ip_of maps for every dataset (one unified pas_id
    # space) — built once, forwarded to each worker via _prov_kwargs so
    # record_pas_drops() / write_pas_gene_artifacts() in
    # run_one_dataset_downstream() can populate the ledger + annotatedpas.bed
    # status columns.
    _prov_kwargs: dict = {"atlas_of": _atlas_of, "ip_of": _ip_of}
    worker_args: list[tuple] = []
    for _warg in _raw_worker_args:
        _ds_id_for_sub = _warg[0]
        _sub_stage = _add_subtask(_downstream_stage, _ds_id_for_sub, total=6)
        _sub_client = _client(_sub_stage)
        worker_args.append(
            _warg + (_sub_client, _cluster_kwargs, _resolved_plot_engines, _prov_kwargs)
        )

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

    _all_stats: list = []  # E3: collect per-dataset stats (carry n_vars for reconcile)
    if n_workers <= 1 or n_datasets == 1:
        # Inline path: no spawn overhead, backward-compatible.
        for arg_tuple in worker_args:
            # arg_tuple has: ds_id, sub_indices, sub_cbs, unified_mtx,
            #   per_dataset_dir, genes_pkl, min_read, min_cells, min_genes,
            #   log_queue, progress_client, cluster_kwargs, plot_engines,
            #   prov_kwargs (atlas_of/ip_of, D9)  (14 elements)
            *pos_args, lq, pc, ck, pe, pk = arg_tuple
            _st = run_one_dataset_downstream(
                *pos_args,
                log_queue=lq,
                progress_client=pc,
                plot_engines=pe,
                **(ck or {}),
                **(pk or {}),
            )
            if isinstance(_st, dict):
                _all_stats.append(_st)
            _advance(_downstream_stage)
    else:
        ctx = multiprocessing.get_context("spawn")
        with ctx.Pool(processes=n_workers) as pool:
            for _stats in pool.imap_unordered(
                downstream_worker_star, worker_args, chunksize=1
            ):
                if isinstance(_stats, dict):
                    _all_stats.append(_stats)
                _advance(_downstream_stage)  # one dataset complete

    # E3 sidecar-then-reconcile: after the pool joins, the parent reconciles the
    # per-dataset ledger sidecars each worker wrote (no shared mutable state) and
    # self-checks the per-dataset invariant surviving PAS == n_vars(clusters.h5ad).
    # The run-level atlas-snap ledger is a DIFFERENT question and stays separate.
    try:
        from ema.provenance import reconcile_dataset_ledgers
        _recon_inputs = []
        for _st in _all_stats:
            _ds = _st.get("dataset_id")
            if _ds is None:
                continue
            _sidecar = directory_config.clusters_h5ad_for(_ds).parent / "provenance"
            _recon_inputs.append((_ds, _sidecar, int(_st.get("final_pas", 0))))
        if _recon_inputs:
            _cohort_dir = output_dir / "provenance" / "by_dataset"
            _recon = reconcile_dataset_ledgers(_recon_inputs, cohort_dir=_cohort_dir)
            import json as _json
            (output_dir / "provenance").mkdir(parents=True, exist_ok=True)
            (output_dir / "provenance" / "reconcile_summary.json").write_text(
                _json.dumps(_recon, indent=2)
            )
            _bad = [d["ds_id"] for d in _recon["datasets"] if not d["invariant_ok"]]
            if _bad:
                log.warning(
                    "provenance reconcile: invariant OFF for %d/%d dataset(s): %s "
                    "(pas_ledger accounting incomplete)",
                    len(_bad), _recon["n_datasets"], _bad,
                )
            else:
                log.info(
                    "provenance reconcile: invariant OK for all %d dataset(s) "
                    "(surviving PAS == n_vars)", _recon["n_datasets"],
                )
    except Exception as _e:  # provenance is auxiliary — never fail the run
        log.warning("provenance reconcile step failed (non-fatal): %s", _e)

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
        # D9: atlas/ip match-vs-no-match stats for the run manifest.
        "atlas_stats": _atlas_stats,
        "ip_stats": {
            "n_ip_flagged": (_pas_filter_result or {}).get("n_ip_flagged"),
            "ip_flag_rate": (_pas_filter_result or {}).get("ip_flag_rate"),
        } if _pas_filter_result else None,
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
