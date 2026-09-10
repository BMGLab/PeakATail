"""Shared body for ``peakatail reannotate`` — branch a completed ``peakatail run`` into a
new trim / clustering variant WITHOUT re-running peak calling.

Peak calling (streaming the BAMs) is the expensive stage.  The "trim" —
``find_close(max_gene_distance, utr_multiplier, include_extended)`` — and
everything downstream of it (annotate -> preprocess -> cluster) is cheap and
depends only on artifacts a base run already wrote to disk:

    <base_run>/posbed.bed              (unified +strand PAS)
    <base_run>/negbed.bed              (unified -strand PAS)
    <base_run>/unified/concatenated.mtx        (PAS x cell counts)
    <base_run>/unified/concatenated_cbs.tsv    (namespaced barcodes)

:func:`reannotate_run` reuses the *exact* tested internals the pipeline uses
(``find_close`` + the per-dataset worker ``run_one_dataset_downstream``), so a
branch is behaviourally identical to having run ``peakatail run`` with those trim /
clustering parameters — it just skips peak calling.  It mirrors the
multi-sample downstream section of ``ema/main.py`` (the code after "Run
find_close ONCE on the unified PAS coordinate set"), including the E3
provenance reconcile step and the E2 ``run_manifest.json`` write, so the
branch's ``--out`` directory is a COMPLETE, CHAINABLE run dir — structurally
indistinguishable from a base run's downstream output.  ``ema.data.Run.
from_dir()`` can load it directly and ``peakatail switch {diff,length,trend}`` can
consume its ``07_clustering/<ds>/clusters.h5ad`` files.

This module is the single implementation both callers delegate to:
  * ``ema/cli/reannotate.py``      — the ``peakatail reannotate`` Click subcommand.
  * ``scripts/reannotate_from_run.py`` — the original standalone script,
    kept working as a thin shim over this function.
"""
from __future__ import annotations

import fcntl
import functools
import json
import logging
import os
import pickle
import socket
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any

log = logging.getLogger(__name__)

__all__ = ["ReannotateError", "reannotate_run"]


class ReannotateError(RuntimeError):
    """Raised when the base run is missing a required artifact, or ``--out``
    collides with ``--base-run`` or with another live branch."""


def _require(path: Path, what: str) -> Path:
    if not path.exists():
        raise ReannotateError(f"base run missing {what}: {path}")
    return path


#: Whole-run exclusive lock file, at the branch's ``--out`` root.
OUT_DIR_LOCK_NAME = ".ema_reannotate.lock"


@contextmanager
def _claim_out_dir(out: Path):
    """Take an exclusive claim on ``out`` for the life of this branch.

    Two ``peakatail reannotate`` invocations pointed at the SAME ``--out`` are
    never a legitimate configuration -- they write the same
    ``04_pas_gene_assignment/<ds>/pas_gene.tsv``,
    ``05_annotated_matrix/<ds>/*`` and ``07_clustering/<ds>/clusters.h5ad``
    paths, so whichever artifacts survive are an arbitrary interleaving of
    two different parameter sets. That is not hypothetical: 5 branches of a
    sweep grid that shared a ``branch_name`` (hence a ``--out``) ran
    concurrently under Nextflow and produced 3 different ``pas_gene.tsv``
    row counts for IDENTICAL declared params, grouped by write time.
    Atomic writes (``ema.outputs.atomic_write``) stop a single file from
    being torn; only this guard stops the two runs from clobbering each
    other's artifacts wholesale. So: refuse, loudly and immediately, rather
    than silently interleave.

    The claim is an ``flock`` on ``<out>/.ema_reannotate.lock``, taken
    before any work starts. flock is held by the open file description, so
    the kernel drops it when this process exits for ANY reason -- a killed
    or crashed branch never leaves a stale lock that blocks the re-run
    (which a plain ``O_EXCL`` marker file would). It is advisory and
    process-scoped, so it does not protect against a run that ignores it,
    only against a second ``peakatail reannotate``.
    """
    out.mkdir(parents=True, exist_ok=True)
    lock_path = out / OUT_DIR_LOCK_NAME
    fd = os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o644)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            holder = ""
            try:
                holder = os.read(fd, 4096).decode("utf-8", "replace").strip()
            except OSError:
                pass
            raise ReannotateError(
                f"--out is already in use by another running `peakatail reannotate`: "
                f"{out}{f' (held by {holder})' if holder else ''} — two branches "
                "writing one output dir interleave their artifacts and silently "
                "corrupt both. Give every branch its OWN --out; a duplicate "
                "branch_name in a sweep grid is the usual cause. "
                f"(Lock file: {lock_path})"
            ) from None
        os.ftruncate(fd, 0)
        os.write(
            fd,
            f"pid={os.getpid()} host={socket.gethostname()} "
            f"started={time.strftime('%Y-%m-%dT%H:%M:%S')}\n".encode(),
        )
        yield lock_path
    finally:
        # Closing the fd releases the flock. The (dot-prefixed, never
        # manifested) lock file itself stays: unlinking it would race a
        # branch that already has the same inode open.
        os.close(fd)


def _guard_out_dir(fn):
    """Decorator: run ``fn`` holding an exclusive claim on its ``--out``.

    Wraps :func:`reannotate_run` (keyword-only, so ``out``/``base_run`` are
    always in ``kwargs``) instead of indenting its whole body under a
    ``with``.
    """
    @functools.wraps(fn)
    def _wrapper(**kwargs):
        out = Path(kwargs["out"]).resolve()
        if Path(kwargs["base_run"]).resolve() == out:
            # Let the body raise the more specific --out == --base-run error
            # rather than dropping a lock file into the base run dir.
            return fn(**kwargs)
        with _claim_out_dir(out):
            return fn(**kwargs)
    return _wrapper


def _load_or_compute_atlas_match(
    base: Path, atlas: str | Path | None, atlas_distance: int, out: Path,
) -> dict[str, tuple]:
    """``{pas_id: (atlas_match_bool, atlas_distance_bp)}`` for every PAS in the
    base run's unified set.

    Prefers ``<base>/unified/atlas_status.tsv`` -- already computed by the
    base run when it used ``--atlas-mode annotate`` (the default), so this
    is normally a plain read, no recompute. Falls back to a fresh
    :func:`ema.datasets.atlas_annotate.annotate_pas_against_atlas` pass
    (read-only, never drops/renumbers a PAS -- only used here to get the
    match status) when the base run has no cached sidecar, e.g. it never
    had an atlas configured. Requires *atlas* in that case.
    """
    cached = base / "unified" / "atlas_status.tsv"
    if cached.exists():
        log.info("atlas label: reusing cached %s", cached)
        status_path = cached
    else:
        if not atlas:
            raise ReannotateError(
                "--exclude-atlas-nonmatch is enabled but the base run has no "
                f"cached {cached} (it wasn't atlas-annotated) and no --atlas "
                "was given to recompute match status from."
            )
        from ema.datasets.atlas_annotate import annotate_pas_against_atlas

        unified_tmp = out / "_atlas_label_check" / "unified_pas.bed"
        unified_tmp.parent.mkdir(parents=True, exist_ok=True)
        with open(unified_tmp, "w") as f:
            for src in (base / "posbed.bed", base / "negbed.bed"):
                if src.exists():
                    f.write(src.read_text())
        status_path, _stats_path = annotate_pas_against_atlas(
            unified_pasbed_path=unified_tmp,
            atlas_bed_path=atlas,
            distance=atlas_distance,
            output_dir=unified_tmp.parent,
        )
        log.info("atlas label: recomputed match status -> %s", status_path)

    import pandas as pd

    df = pd.read_csv(status_path, sep="\t", dtype=str)
    return {
        row["unified_pas_id"]: (
            str(row["atlas_match"]).strip().lower() == "true",
            row.get("atlas_distance_bp", ""),
        )
        for _, row in df.iterrows()
    }


@_guard_out_dir
def reannotate_run(
    *,
    base_run: str | Path,
    out: str | Path,
    gtf: str | Path,
    max_gene_distance: int = 5000,
    utr_multiplier: float = 2.0,
    include_extended: bool = False,
    cluster_method: str = "leiden_tfidf",
    resolution: float = 1.0,
    n_neighbors: int | None = None,
    n_pcs: int = 40,
    n_svd_components: int = 50,
    n_top_hvg: int = 2000,
    random_seed: int = 42,
    tfidf_scale_factor: float = 1e4,
    depth_corr_threshold: float = 0.75,
    external_clusters: str | None = None,
    min_read: int = 1500,
    min_cells: int = 3,
    min_pas_per_cell: int = 50,
    threads: int | None = None,
    atlas: str | Path | None = None,
    atlas_distance: int = 50,
    genome_fasta: str | Path | None = None,
    ip_window_left: int = 10,
    ip_window_right: int = 30,
    ip_a_stretch: int = 6,
    ip_a_fraction: float = 0.7,
    annotation_bed: str | Path | None = None,
    exclude_atlas_nonmatch: bool = False,
    exclude_internal_priming: bool = False,
    exclude_not_in_3utr: bool = False,
) -> dict[str, Any]:
    """Branch ``base_run`` into ``out`` with new trim / label / mask / cluster params.

    Validates the base run's artifacts exist BEFORE importing anything heavy
    (mirrors ``scripts/reannotate_from_run.py``), then:

      0. (NEW, all optional/default-off) LABEL every PAS -- atlas_match,
         internal_priming, in_3utr -- from *atlas*/*genome_fasta*/
         *annotation_bed* when given. Labels always cover every PAS in the
         base run's unified set and are written into every dataset's
         ``annotatedpas.bed`` (never drops a row there). Independently, the
         ``exclude_*`` flags turn a label into a CLUSTERING-ONLY exclusion:
         the flagged PAS still get their full labeled ``annotatedpas.bed``
         row, but are masked out of the matrix ``preprocessing()``/
         ``clustering()`` (and therefore every downstream switch step) sees.
         See "PAS labels + clustering mask" below.
      1. Re-derives ``gene_end.bed`` / UTR lengths from ``gtf`` (two-tier
         cached — near-instant on repeat).
      2. Re-runs ``find_close`` (the trim) on the base run's unified PAS set
         -- ALWAYS the full, unfiltered set; labels/masking never touch this.
      3. Re-runs ``run_one_dataset_downstream`` per dataset (annotate ->
         [clustering-mask] -> preprocess -> cluster), producing a fresh
         ``07_clustering/<ds>/clusters.h5ad`` per dataset.
      4. Reconciles the per-dataset E3 provenance ledgers into
         ``provenance/by_dataset/{pas,cell}_ledger.tsv`` +
         ``reconcile_summary.json`` (same code path ``ema/main.py`` calls
         after its per-dataset workers finish).
      5. Writes ``run_manifest.json`` (E2) via
         ``ema.outputs.OutputManager.write_manifest``, exactly as ``peakatail run``
         does, so the branch is a complete, chainable run dir.

    Args:
        base_run: Completed ``peakatail run`` output dir to branch from.
        out: Fresh output dir for this branch (must differ from ``base_run``).
        gtf: Same GTF as the base run (or a different one, to re-annotate
            against a new annotation without re-calling peaks).
        max_gene_distance: TIER_3 distal cap (bp). Matters when
            ``include_extended``.
        utr_multiplier: TIER_2 boundary = gene UTR length x this.
        include_extended: Keep TIER_3 (distal/novel) PAS out to
            ``max_gene_distance``.
        cluster_method, resolution, n_neighbors, n_pcs, n_svd_components,
            n_top_hvg, random_seed, tfidf_scale_factor, depth_corr_threshold,
            external_clusters: Clustering hyperparameters forwarded verbatim
            to ``run_one_dataset_downstream`` / ``clustering()`` — every knob
            those accept is a parameter here too, none hardcoded, so an OFAT
            sweep can vary any of them purely via CLI flags.
            ``tfidf_scale_factor``/``depth_corr_threshold`` only affect
            ``leiden_tfidf``; ``n_top_hvg`` only affects ``leiden_libsize``;
            ``external_clusters`` is the label TSV path for
            ``--cluster-method external``.
        min_read, min_cells, min_pas_per_cell: Cell/PAS filters. Defaults
            MUST match ``peakatail run``'s schema defaults (``ema/cli/
            config_schema.py``) so a branch with unchanged params reproduces
            the base run's clustering.
        threads: Absolute worker ceiling wired into the ``ResourceManager``
            singleton (same mechanism ``ema/cli/run.py`` uses), so the
            downstream per-dataset ``Pool`` respects it. ``None`` leaves the
            existing/auto-detected ceiling untouched.
        atlas: Reference atlas BED. When given (or when the base run has a
            cached ``unified/atlas_status.tsv`` from its own ``--atlas-mode
            annotate``), every PAS is labeled ``atlas_match``/
            ``atlas_distance_bp`` in ``annotatedpas.bed`` -- pure
            annotation, nothing is ever dropped by this alone. Required
            (unless the cache exists) when *exclude_atlas_nonmatch* is set.
        atlas_distance: Max summit-to-atlas distance (bp) for a match, used
            only by the fallback recompute path (cache miss).
        genome_fasta, ip_window_left, ip_window_right, ip_a_stretch,
            ip_a_fraction: When *genome_fasta* is given, every PAS is
            labeled ``internal_priming`` (A-rich downstream stretch) in
            ``annotatedpas.bed`` via :func:`ema.experimental.
            internal_priming.filter_internal_priming` (mode="annotate" --
            label only, nothing dropped). Required when
            *exclude_internal_priming* is set. The other four are forwarded
            verbatim; same defaults as ``peakatail run``'s schema.
        annotation_bed: A region BED (e.g. a 3'UTR-only BED, for "PAS inside
            an annotated transcript 3'UTR"). When given, every PAS is
            labeled ``in_3utr`` in ``annotatedpas.bed`` via
            :func:`ema.experimental.peak_filters.label_pas_in_bed` (label
            only, nothing dropped). Required when *exclude_not_in_3utr* is
            set.
        exclude_atlas_nonmatch, exclude_internal_priming,
            exclude_not_in_3utr: All default False -- a branch with none of
            these set reproduces the base run's downstream exactly,
            byte-for-byte, same guarantee as the trim/cluster knobs above.
            When set, the corresponding label EXCLUDES that PAS from the
            CLUSTERING matrix only (see ``ema.downstream_runner.
            run_one_dataset_downstream``'s ``exclude_pas_ids``) -- every
            dataset's ``annotatedpas.bed`` still carries the FULL labeled
            PAS set regardless of which (if any) exclude flags are set.
            Multiple exclude flags OR together (a PAS matching any one is
            excluded). A PAS with no computed label for an active exclude
            flag (label source didn't cover it) is treated as excluded
            (conservative: can't confirm it belongs, so it doesn't survive
            a "keep only confirmed-X" mask) -- mirrors *atlas*'s existing
            "absent = unmatched" convention.

    Returns:
        The branch manifest dict written to ``<out>/branch_manifest.json``
        (also embeds ``manifest_path`` — the E2 ``run_manifest.json`` path —
        and ``reconcile_summary``).

    Raises:
        ReannotateError: ``out == base_run``; ``out`` is already claimed by
            another live ``peakatail reannotate`` (see :func:`_claim_out_dir`); or a
            required base-run artifact is missing.
    """
    base = Path(base_run).resolve()
    out = Path(out).resolve()
    if base == out:
        raise ReannotateError(
            "--out must differ from --base-run (would overwrite peaks)."
        )
    out.mkdir(parents=True, exist_ok=True)

    # Locate + validate the base artifacts before importing anything heavy.
    posbed = _require(base / "posbed.bed", "posbed.bed")
    negbed = _require(base / "negbed.bed", "negbed.bed")
    unified_mtx = _require(base / "unified" / "concatenated.mtx", "unified/concatenated.mtx")
    unified_cbs = _require(base / "unified" / "concatenated_cbs.tsv", "unified/concatenated_cbs.tsv")

    gtf = Path(gtf)

    # Heavy imports (need the peakatail venv on PYTHONPATH).
    from ema.config import set_directory_config, directory_config, filter_config
    from ema.annotate.gtf_cache import process_gtf_cached
    from ema.annotate.find_close import find_close
    from ema.downstream_runner import run_one_dataset_downstream, downstream_worker_star
    from ema.outputs import OutputManager, build_resolved_run_config
    from ema.provenance import reconcile_dataset_ledgers
    from ema.utils import get_resource_manager

    # Wire --threads into the ResourceManager the same way `peakatail run` does
    # (ema/cli/run.py), so the downstream per-dataset Pool respects it.
    if threads is not None:
        from ema.utils import reset_resource_manager
        from ema.utils.resource_manager import ResourceManager
        import ema.utils as _utils_mod
        reset_resource_manager()
        _utils_mod._RM_INSTANCE = ResourceManager(user_max_threads=threads)

    # Point the pipeline globals at the OUT dir (downstream writes go there).
    # output_dir MUST be a Path — DirectoryConfig stores it verbatim and every
    # path property does `self.output_dir / <name>`; a str would TypeError.
    set_directory_config(output_dir=out, gtf_dir=str(gtf))
    filter_config.min_read = min_read
    filter_config.min_cells = min_cells
    filter_config.min_genes = min_pas_per_cell  # min_pas_per_cell bridges to min_genes

    # 1. GTF -> gene_end.bed + utr_lengths (two-tier cached; fast on repeat).
    log.info("processing GTF (cached) -> gene_end.bed + utr_lengths")
    utr_lengths = process_gtf_cached(
        gtf_path=str(gtf),
        output_dir=str(out),
        endbed_path=str(directory_config.endbed),
        features_path=str(directory_config.raw_features),
    )
    log.info("utr_lengths: %d genes", len(utr_lengths))

    # 0. PAS LABELS (all optional, default off) — computed ONCE from the base
    #    run's unified posbed.bed/negbed.bed, keyed by the same pas_id
    #    namespace those beds (and find_close()'s output) use. NEVER mutates
    #    posbed/negbed or drops a row -- find_close() below always runs on
    #    the base run's ORIGINAL, full, unfiltered beds, so annotatedpas.bed
    #    for every dataset carries the FULL labeled PAS set regardless of
    #    which (if any) exclude_* flags are set.
    label_stats: dict = {}
    atlas_of: dict = {}
    ip_of: dict = {}
    in_3utr_of: dict = {}

    if atlas or (base / "unified" / "atlas_status.tsv").exists():
        atlas_of = _load_or_compute_atlas_match(base, atlas, atlas_distance, out)
        n_matched = sum(1 for v in atlas_of.values() if v[0])
        label_stats["atlas"] = {
            "n_labeled": len(atlas_of), "n_matched": n_matched,
            "n_unmatched": len(atlas_of) - n_matched,
        }
        log.info(
            "PAS label: atlas_match -- %d/%d matched", n_matched, len(atlas_of),
        )

    if genome_fasta:
        from ema.experimental.internal_priming import filter_internal_priming

        unified_tmp = out / "_ip_label_check" / "unified_pas.bed"
        unified_tmp.parent.mkdir(parents=True, exist_ok=True)
        with open(unified_tmp, "w") as f:
            for src in (posbed, negbed):
                if src.exists():
                    f.write(src.read_text())
        ip_discard = unified_tmp.parent / "unified_pas.ip_annotated.bed"
        ip_stats = filter_internal_priming(
            str(unified_tmp), str(genome_fasta), str(ip_discard),
            window_left=ip_window_left, window_right=ip_window_right,
            a_stretch=ip_a_stretch, a_fraction=ip_a_fraction,
            mode="annotate",
        )
        ip_of = ip_stats.get("flags") or {}
        n_flagged = sum(1 for v in ip_of.values() if v)
        label_stats["internal_priming"] = {
            "n_labeled": len(ip_of), "n_flagged": n_flagged,
        }
        log.info(
            "PAS label: internal_priming -- %d/%d flagged", n_flagged, len(ip_of),
        )

    if annotation_bed:
        from ema.experimental.peak_filters import label_pas_in_bed

        unified_tmp = out / "_3utr_label_check" / "unified_pas.bed"
        unified_tmp.parent.mkdir(parents=True, exist_ok=True)
        with open(unified_tmp, "w") as f:
            for src in (posbed, negbed):
                if src.exists():
                    f.write(src.read_text())
        label_result = label_pas_in_bed(str(unified_tmp), str(annotation_bed))
        in_3utr_of = label_result.get("flags") or {}
        label_stats["in_3utr"] = {
            "n_labeled": label_result["total"], "n_in_3utr": label_result["n_in_region"],
            "n_not_in_3utr": label_result["n_not_in_region"],
        }
        log.info(
            "PAS label: in_3utr -- %d/%d inside a 3'UTR",
            label_result["n_in_region"], label_result["total"],
        )

    # Independent of labeling: which labels ALSO exclude from clustering.
    if exclude_atlas_nonmatch and not atlas_of:
        raise ReannotateError(
            "--exclude-atlas-nonmatch is enabled but no atlas label was computed "
            "(no --atlas and no cached unified/atlas_status.tsv on the base run)."
        )
    if exclude_internal_priming and not ip_of:
        raise ReannotateError(
            "--exclude-internal-priming is enabled but no --genome-fasta was given "
            "to compute the internal_priming label."
        )
    if exclude_not_in_3utr and not in_3utr_of:
        raise ReannotateError(
            "--exclude-not-in-3utr is enabled but no --annotation-bed was given "
            "to compute the in_3utr label."
        )

    exclude_pas_ids: set = set()
    if exclude_atlas_nonmatch:
        exclude_pas_ids |= {pid for pid, (matched, _dist) in atlas_of.items() if not matched}
    if exclude_internal_priming:
        exclude_pas_ids |= {pid for pid, flagged in ip_of.items() if flagged}
    if exclude_not_in_3utr:
        exclude_pas_ids |= {pid for pid, in_utr in in_3utr_of.items() if not in_utr}
    if exclude_pas_ids:
        log.info(
            "clustering-exclusion mask: %d PAS excluded (atlas_nonmatch=%s "
            "internal_priming=%s not_in_3utr=%s)",
            len(exclude_pas_ids), exclude_atlas_nonmatch,
            exclude_internal_priming, exclude_not_in_3utr,
        )

    # 2. THE TRIM — re-run find_close on the base run's FULL, unfiltered
    #    unified PAS set (labels/exclusion never touch posbed/negbed).
    log.info(
        "find_close: max_distance=%d utr_multiplier=%.2f include_extended=%s",
        max_gene_distance, utr_multiplier, include_extended,
    )
    genes = find_close(
        posbed_dir=str(posbed),
        negbed_dir=str(negbed),
        genomebed_dir=str(directory_config.endbed),
        annotatedbed_dir=str(directory_config.annotatedbed),  # -> OUT
        mergebed=str(directory_config.pasbed),                # -> OUT
        utr_lengths=utr_lengths,
        max_distance=max_gene_distance,
        utr_multiplier=utr_multiplier,
        include_extended=include_extended,
    )
    log.info("find_close assigned %d PAS->gene rows", len(genes))
    genes_pkl = pickle.dumps(genes)

    # 3. Derive datasets from the namespaced barcodes.  cb = "<ds_id>_<barcode>";
    #    the barcode half is underscore-free but ds_id need not be (peakatail merge
    #    stamps RG = dataset_id), so split on the LAST '_' -- see
    #    ema.countmatrix.indexing.split_cb.  Splitting on the first '_' here
    #    truncated ids like "pbmc_10k_v3" to "pbmc", after which the
    #    `startswith(f"{ds_id}_")` selector below picked up the wrong columns.
    all_cbs = [ln.strip() for ln in unified_cbs.read_text().splitlines() if ln.strip()]
    unique_ds_ids: list[str] = list(dict.fromkeys(cb.rsplit("_", 1)[0] for cb in all_cbs))
    log.info("branching %d datasets: %s", len(unique_ds_ids), ", ".join(unique_ds_ids))

    # Record the dataset list on directory_config so the run manifest (E2)
    # carries DatasetRef entries for the branch, mirroring `peakatail run`.
    set_directory_config(datasets=[{"id": ds_id} for ds_id in unique_ds_ids])

    # write_pas_gene_artifacts() (called inside run_one_dataset_downstream)
    # extends <ds>/pasbed.bed with the gene_id column to build
    # annotatedpas.bed -- but that per-dataset pasbed.bed is normally written
    # by peak calling (ema.outputs.write_per_dataset_beds, called from
    # ema/main.py BEFORE the downstream section), the exact stage reannotate
    # skips.  Reuse the base run's own copy of that artifact (same peak-call
    # coordinates/pas_id namespace this branch's unified PAS set was built
    # from) so annotatedpas.bed still gets written for each dataset.
    for ds_id in unique_ds_ids:
        src_pasbed = base / "01_peak_calling" / ds_id / "pasbed.bed"
        if src_pasbed.exists():
            dst_pasbed = directory_config.pasbed_for(ds_id)
            dst_pasbed.parent.mkdir(parents=True, exist_ok=True)
            dst_pasbed.write_bytes(src_pasbed.read_bytes())
        else:
            log.warning(
                "base run has no 01_peak_calling/%s/pasbed.bed -- "
                "annotatedpas.bed will be skipped for this dataset", ds_id,
            )

    # Every run_one_dataset_downstream/clustering() knob is forwarded
    # verbatim from a reannotate_run() parameter -- nothing hardcoded here,
    # so an OFAT sweep can vary any of them purely via CLI flags.
    cluster_kwargs = dict(
        cluster_method=cluster_method,
        cluster_resolution=resolution,
        cluster_n_pcs=n_pcs,
        cluster_random_seed=random_seed,
        cluster_external_clusters=external_clusters,
        cluster_n_neighbors=n_neighbors,
        cluster_tfidf_scale_factor=tfidf_scale_factor,
        cluster_depth_corr_threshold=depth_corr_threshold,
        cluster_n_svd_components=n_svd_components,
        cluster_n_top_hvg=n_top_hvg,
    )
    # PAS labels + clustering-exclusion mask (step 0 above) -- forwarded
    # verbatim to every dataset's run_one_dataset_downstream() call. Full,
    # unified maps: write_pas_gene_artifacts() looks each PAS up by id, so
    # every dataset's annotatedpas.bed gets the same labels regardless of
    # which cells/PAS that dataset happens to carry.
    prov_kwargs = dict(
        atlas_of=atlas_of, ip_of=ip_of, in_3utr_of=in_3utr_of,
        exclude_pas_ids=exclude_pas_ids,
    )

    # 4. Run the tested per-dataset downstream worker for each dataset.
    worker_specs: list[tuple[str, list[int], list[str]]] = []
    for ds_id in unique_ds_ids:
        # Exact match on the sample half (see ema.countmatrix.indexing.split_cb);
        # `startswith(f"{ds_id}_")` let dataset "a" claim dataset "a_b"'s cells.
        sub_indices = [i for i, cb in enumerate(all_cbs) if cb.rsplit("_", 1)[0] == ds_id]
        if not sub_indices:
            log.warning("no cells for '%s' — skipping", ds_id)
            continue
        sub_cbs = [all_cbs[i] for i in sub_indices]
        worker_specs.append((ds_id, sub_indices, sub_cbs))

    n_datasets = len(worker_specs)
    n_workers = (
        min(n_datasets, get_resource_manager().get_n_jobs(per_worker_mb=500, stage="downstream"))
        if n_datasets else 0
    )

    results: list[dict] = []
    if n_workers <= 1 or n_datasets <= 1:
        # Inline path: no spawn overhead, matches scripts/reannotate_from_run.py.
        for ds_id, sub_indices, sub_cbs in worker_specs:
            log.info("  [%s] %d cells -> downstream", ds_id, len(sub_indices))
            stats = run_one_dataset_downstream(
                ds_id, sub_indices, sub_cbs,
                str(unified_mtx), str(out), genes_pkl,
                min_read, min_cells, min_pas_per_cell,
                **cluster_kwargs, **prov_kwargs,
            )
            results.append(stats)
    else:
        import multiprocessing

        log.info("Per-dataset downstream: %d datasets, n_workers=%d", n_datasets, n_workers)
        ctx = multiprocessing.get_context("spawn")
        worker_args = [
            (
                ds_id, sub_indices, sub_cbs,
                str(unified_mtx), str(out), genes_pkl,
                min_read, min_cells, min_pas_per_cell,
                # log_queue, progress_client, cluster_kwargs, plot_engines, prov_kwargs
                # (14-arg downstream_worker_star shape -- see its docstring)
                None, None, cluster_kwargs, None, prov_kwargs,
            )
            for ds_id, sub_indices, sub_cbs in worker_specs
        ]
        with ctx.Pool(processes=n_workers) as pool:
            for stats in pool.imap_unordered(downstream_worker_star, worker_args, chunksize=1):
                if isinstance(stats, dict):
                    results.append(stats)

    # 5. E3 sidecar-then-reconcile — same code path ema/main.py runs after its
    #    per-dataset pool joins: each worker already wrote its own ledger
    #    sidecar (no shared mutable state); reconcile concatenates them into
    #    the cohort ledger and self-checks surviving PAS == n_vars(h5ad).
    reconcile_summary: dict | None = None
    try:
        recon_inputs = []
        for st in results:
            ds = st.get("dataset_id")
            if ds is None:
                continue
            sidecar = directory_config.clusters_h5ad_for(ds).parent / "provenance"
            recon_inputs.append((ds, sidecar, int(st.get("final_pas", 0))))
        if recon_inputs:
            cohort_dir = out / "provenance" / "by_dataset"
            reconcile_summary = reconcile_dataset_ledgers(recon_inputs, cohort_dir=cohort_dir)
            (out / "provenance").mkdir(parents=True, exist_ok=True)
            (out / "provenance" / "reconcile_summary.json").write_text(
                json.dumps(reconcile_summary, indent=2)
            )
            bad = [d["ds_id"] for d in reconcile_summary["datasets"] if not d["invariant_ok"]]
            if bad:
                log.warning(
                    "provenance reconcile: invariant OFF for %d/%d dataset(s): %s "
                    "(pas_ledger accounting incomplete)",
                    len(bad), reconcile_summary["n_datasets"], bad,
                )
            else:
                log.info(
                    "provenance reconcile: invariant OK for all %d dataset(s) "
                    "(surviving PAS == n_vars)", reconcile_summary["n_datasets"],
                )
    except Exception as e:  # provenance is auxiliary — never fail the branch over it
        log.warning("provenance reconcile step failed (non-fatal): %s", e)

    # 5b. Cross-dataset cluster matching (only meaningful if >1 dataset) —
    #     mirrors ema/main.py's post-downstream block so a multi-dataset
    #     branch also gets `canonical_cluster` written into each dataset's
    #     obs, letting `peakatail switch diff --cluster-key canonical_cluster`
    #     compare clusters across datasets the same way it would for a base
    #     `peakatail run`. Best-effort / non-fatal, same narrow exception set as
    #     the code it mirrors.
    if len(unique_ds_ids) > 1:
        h5ad_paths = [directory_config.clusters_h5ad_for(ds) for ds in unique_ds_ids]
        existing = [(p, ds) for p, ds in zip(h5ad_paths, unique_ds_ids) if p.exists()]
        if len(existing) > 1:
            try:
                from ema.clustering.cross_dataset import get_match_strategy
                match_strategy = get_match_strategy("marker_overlap")
                match_df = match_strategy.match(
                    h5ad_paths=[p for p, _ in existing],
                    dataset_ids=[ds for _, ds in existing],
                )
                cross_dir = out / "cross_dataset"
                cross_dir.mkdir(exist_ok=True)
                match_df.to_csv(cross_dir / "canonical_cluster_map.tsv", sep="\t", index=False)
                n_canonical = (
                    match_df["canonical_cluster"].nunique()
                    if "canonical_cluster" in match_df.columns else 0
                )
                log.info(
                    "Cross-dataset matching: %d cluster entries -> %d canonical clusters",
                    len(match_df), n_canonical,
                )
                from ema.clustering.cross_dataset.roundtrip import write_canonical_clusters
                try:
                    write_canonical_clusters(match_df, existing)
                except (KeyError, OSError, ValueError) as e:
                    log.warning(
                        "canonical_cluster round-trip into obs skipped (%s): %s",
                        type(e).__name__, e,
                    )
            except (FileNotFoundError, KeyError, ValueError, OSError) as e:
                log.warning("Cross-dataset matching skipped (%s): %s", type(e).__name__, e)

    # 6. E2 run manifest — exactly as `peakatail run` does (ema/main.py), so the
    #    branch is a COMPLETE, CHAINABLE run dir: `ema.data.Run.from_dir(out)`
    #    loads it and `peakatail switch {diff,length,trend}` can consume its
    #    07_clustering/<ds>/clusters.h5ad files directly.
    entity_counts = {
        "n_datasets": len(results),
        "final_cells_total": sum(int(r.get("final_cells", 0)) for r in results),
        "final_pas_total": sum(int(r.get("final_pas", 0)) for r in results),
    }
    mgr = OutputManager(base_dir=str(directory_config.output_dir))
    manifest_path = mgr.write_manifest(build_resolved_run_config(), entity_counts=entity_counts)

    # 7. Branch manifest — human-readable trim/filter/cluster provenance for
    #    this specific branch (kept for scripts/reannotate_from_run.py
    #    backward compatibility; the E2 run_manifest.json above is what the
    #    hub / ema.data.Run actually read).
    branch_manifest = {
        "base_run": str(base),
        "gtf": str(gtf),
        "trim": {
            "max_gene_distance": max_gene_distance,
            "utr_multiplier": utr_multiplier,
            "include_extended": include_extended,
        },
        "clustering": {
            "method": cluster_method,
            "resolution": resolution,
            "n_neighbors": n_neighbors,
        },
        "filters": {
            "min_read": min_read,
            "min_cells": min_cells,
            "min_pas_per_cell": min_pas_per_cell,
        },
        "pas_labels_and_mask": {
            "exclude_atlas_nonmatch": exclude_atlas_nonmatch,
            "exclude_internal_priming": exclude_internal_priming,
            "exclude_not_in_3utr": exclude_not_in_3utr,
            # Aggregate across datasets -- every dataset's annotatedpas.bed
            # carries the FULL labeled set (n_kept_in_results); the
            # clustering matrix only sees n_used_for_clustering.
            "n_kept_in_results": sum(int(r.get("n_kept_in_results", 0)) for r in results),
            "n_used_for_clustering": sum(int(r.get("n_used_for_clustering", 0)) for r in results),
            "n_excluded_for_clustering": sum(
                int(r.get("n_excluded_for_clustering", 0)) for r in results
            ),
            "labels": label_stats,
        },
        "datasets": results,
        "manifest_path": manifest_path,
        "reconcile_summary": reconcile_summary,
    }
    (out / "branch_manifest.json").write_text(
        json.dumps(branch_manifest, indent=2, default=str)
    )
    log.info("DONE — %d datasets clustered under %s", len(results), out)
    return branch_manifest
