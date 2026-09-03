"""Pipeline output persistence.

Single module for everything the pipeline writes to disk other than viz
figures (which live in ``ema/viz/``).  Two concerns:

  1. :class:`OutputManager` — owns the numbered ``0X_<stage>/`` directory
     layout under the run root, plus the ``*_stats.json`` files each
     stage emits.  Was previously in ``ema/output.py`` (singular); merged
     here so callers have one place to look.

  2. ``write_*`` functions — persist actual data payloads (BEDs, count
     matrices, h5ads, PAS-gene mappings) under ``<run>/per_dataset/<ds>/``
     and ``<run>/per_dataset/<ds>/raw/`` so every pipeline step leaves
     inspectable artifacts, not just a stats JSON.

The pipeline body in ``ema/main.py`` and the per-dataset worker in
``ema/downstream_runner.py`` both delegate file IO to this module — see
also ``ema/viz/pipeline_hooks.py`` for the analogous viz separation.
"""
from __future__ import annotations

import json
import logging
import os
import threading
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Iterable

log = logging.getLogger(__name__)


def _unique_write_token() -> str:
    """A token unique enough to disambiguate concurrent atomic-write temp
    files, whether the concurrency is across OS processes (the real-world
    case -- two separate ``ema reannotate`` invocations) or threads within
    one process (e.g. test harnesses). PID alone only disambiguates the
    former; adding the thread ident closes the latter too.
    """
    return f"{os.getpid()}_{threading.get_ident()}"


@contextmanager
def atomic_write(dst: Path | str, mode: str = "w", **open_kwargs):
    """Open a temp file next to ``dst``, then ``os.replace()`` it into place.

    EVERY text artifact this module writes goes through here -- never write
    a final artifact path directly. The temp file lives in the SAME
    directory as ``dst`` (same filesystem, so ``os.replace()`` is a single
    atomic rename syscall), which closes a real corruption window: two
    pipeline invocations whose ``--out`` accidentally resolves to the same
    run dir (e.g. a duplicate ``ema reannotate`` branch name) run
    CONCURRENTLY, and two processes writing the same path with a plain
    ``open(..., "w")`` can interleave, leaving a reader with a torn file --
    the actual mechanism behind a real bug (5 A2/A3 trim/cluster branches
    with identical params producing 3 different ``pas_gene.tsv`` row counts,
    grouped by write time). Two atomic writers racing on the same path
    still can't corrupt the file -- the loser's complete temp file just gets
    rename()'d over; whichever wins is a complete, valid file, never a
    byte-level mix of both. A reader always sees either the complete prior
    content or the complete new content. If the body raises, the temp file
    is removed and ``dst`` is left untouched.

    (The underlying cause -- a duplicate ``--out`` -- is guarded separately
    and up front, see ``ema.reannotate.reannotate_run``; this makes the
    write itself safe regardless.)
    """
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f".{dst.name}.tmp{_unique_write_token()}")
    try:
        with open(tmp, mode, **open_kwargs) as fh:
            yield fh
        os.replace(tmp, dst)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


# MatrixMarket header written by matrixfilter and annotate.  Lives here
# (not in variable_config) because it is a fixed constant, not a user
# tunable.
MATRIX_MARKET_HEADER = "%%MatrixMarket matrix coordinate integer general\n"


class OutputManager:
    """Manages the numbered ``0X_<stage>/`` layout + per-stage stats JSONs.

    Stage names match the pipeline phases in ``ema/main.py``:

        peak_calling      -> 01_peak_calling/
        cb_filter         -> 02_cb_filter/
        gtf_annotation    -> 03_gtf_annotation/
        pas_gene          -> 04_pas_gene_assignment/
        annotated         -> 05_annotated_matrix/
        preprocessing     -> 06_preprocessing/
        clustering        -> 07_clustering/
        differential      -> 08_differential/
        gtf_cache         -> gtf_cache/   (shared)

    Use :meth:`path` to address a file within a stage directory and
    :meth:`save_stats` to drop the canonical ``<stage>_stats.json``.
    """

    def __init__(self, base_dir: str = "emaout") -> None:
        self.base_dir = base_dir
        # E2: registered artifacts for run_manifest.json.
        self._artifacts: list[dict] = []
        self.dirs = {
            "peak_calling": os.path.join(base_dir, "01_peak_calling"),
            "cb_filter": os.path.join(base_dir, "02_cb_filter"),
            "gtf_annotation": os.path.join(base_dir, "03_gtf_annotation"),
            "pas_gene": os.path.join(base_dir, "04_pas_gene_assignment"),
            "annotated": os.path.join(base_dir, "05_annotated_matrix"),
            "preprocessing": os.path.join(base_dir, "06_preprocessing"),
            "clustering": os.path.join(base_dir, "07_clustering"),
            "differential": os.path.join(base_dir, "08_differential"),
            "gtf_cache": os.path.join(base_dir, "gtf_cache"),
            # D9: atlas-snap stats (n_atlas_matched/n_atlas_unmatched/...)
            # land alongside the other atlas-snap outputs (atlas_mapping.tsv,
            # atlas_status.tsv, ...) under unified/, the multi-sample-path
            # merge/atlas output dir -- not a new numbered stage dir.
            "atlas_snap": os.path.join(base_dir, "unified"),
        }

    def setup(self) -> None:
        """Create all stage directories (idempotent)."""
        for d in self.dirs.values():
            os.makedirs(d, exist_ok=True)

    def path(self, stage: str, filename: str) -> str:
        """Resolve ``<run>/0X_<stage>/<filename>`` for a known stage."""
        return os.path.join(self.dirs[stage], filename)

    def save_stats(self, stage: str, stats: dict) -> None:
        """Write the canonical ``<stage>_stats.json`` for a pipeline stage."""
        path = self.path(stage, f"{stage}_stats.json")
        with atomic_write(path) as f:
            json.dump(stats, f, indent=2)

    def save_run_config(self, args_dict: dict) -> None:
        """Persist the resolved run configuration at the run root.

        ``args_dict`` should be the *resolved* configuration — see
        :func:`build_resolved_run_config`, which reconciles the argparse
        namespace with the resolved ``directory_config`` / ``variable_config`` /
        ``filter_config`` singletons the pipeline body actually reads. Passing a
        bare ``vars(args)`` here reintroduces bug B0 (records argparse defaults;
        e.g. ``atlas`` shows ``null`` on runs that snapped).
        """
        config = {"timestamp": datetime.now().isoformat(), **args_dict}
        path = os.path.join(self.base_dir, "run_config.json")
        with atomic_write(path) as f:
            json.dump(config, f, indent=2, default=str)

    # ------------------------------------------------------------------ #
    # E2: run_manifest.json — the contract artifact the hub indexes.     #
    # ------------------------------------------------------------------ #
    @staticmethod
    def _content_hash(abs_path) -> str | None:
        """sha256 of an artifact file as ``"sha256:<hex>"`` (None if unreadable).

        Lets any manifest consumer detect a stale index when a REFERENCED
        artifact (ledger/parquet/h5ad) changes but the manifest bytes don't —
        hashing only the manifest sha256 misses that. Streamed so large h5ads
        don't blow memory. The contract Artifact model doesn't yet field this,
        so it rides as an extra key (ignored on validate); consuming it needs a
        contract bump (flagged for hub-team in HANDOFF.md).
        """
        import hashlib

        try:
            p = Path(abs_path)
            if not p.is_file():
                return None
            h = hashlib.sha256()
            with p.open("rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            return f"sha256:{h.hexdigest()}"
        except OSError:
            return None

    def register_artifact(
        self,
        path: str,
        *,
        stage: str,
        fmt: str,
        schema_name: str,
        schema_version: str = "0.1.0",
        entity_counts: dict | None = None,
    ) -> None:
        """Register one produced artifact for the run manifest (E2).

        Fields conform to peakatail_contract.models.Artifact:
        ``path, stage, format, schema_name, schema_version, entity_counts``
        (plus an extra ``content_hash`` the contract ignores on validate).
        ``path`` is stored relative to the run root when possible so the
        manifest is relocatable.
        """
        try:
            rel = os.path.relpath(path, self.base_dir)
        except ValueError:
            rel = str(path)
        self._artifacts.append(
            {
                "path": rel,
                "stage": stage,
                "format": fmt,
                "schema_name": schema_name,
                "schema_version": schema_version,
                "entity_counts": entity_counts or {},
                "content_hash": self._content_hash(path),
            }
        )

    def _auto_discover_artifacts(self) -> list[dict]:
        """Best-effort scan of the run root for standard artifacts (E2).

        Complements explicit ``register_artifact`` calls so the manifest is
        useful even where registration isn't threaded through. Emits
        contract-conformant Artifact dicts (format is a Format enum value).
        """
        # (rel_path, stage, format, schema_name)
        patterns = [
            ("unified/pas_uid.tsv", "unified", "tsv", "pas_uid"),
            ("unified/atlas_mapping.tsv", "atlas_snap", "tsv", "atlas_mapping"),
            ("unified/multi_sample_pas_mapping.tsv", "pas_merge", "tsv", "pas_mapping"),
            ("provenance/pas_ledger.tsv", "provenance", "tsv", "PasLedgerRow"),
            ("provenance/cell_ledger.tsv", "provenance", "tsv", "CellLedgerRow"),
            # Multi-sample runs write the reconciled per-dataset ledgers under
            # by_dataset/ (the run-level provenance/pas_ledger.tsv only exists
            # when a run-level, e.g. atlas-snap, drop ledger was written).
            # Listed AFTER the run-level paths so, when both exist, the
            # run-level one is the first manifest entry for its schema_name;
            # on a no-atlas multi-sample run only these by_dataset entries are
            # emitted (the run-level files don't exist, so aren't registered).
            ("provenance/by_dataset/pas_ledger.tsv", "provenance", "tsv", "PasLedgerRow"),
            ("provenance/by_dataset/cell_ledger.tsv", "provenance", "tsv", "CellLedgerRow"),
            ("provenance/reconcile_summary.json", "provenance", "json", "reconcile_summary"),
            ("run_config.json", "run", "json", "run_config"),
        ]
        found: list[dict] = []
        base = Path(self.base_dir)
        seen = {a["path"] for a in self._artifacts}

        def _emit(rel, stage, fmt, schema_name):
            if rel in seen:
                return
            found.append({
                "path": rel, "stage": stage, "format": fmt,
                "schema_name": schema_name, "schema_version": "0.1.0",
                "entity_counts": {},
                "content_hash": self._content_hash(base / rel),
            })

        for rel, stage, fmt, schema_name in patterns:
            if (base / rel).exists():
                _emit(rel, stage, fmt, schema_name)
        # Per-dataset clustering h5ads + stage stats.
        clustering = base / "07_clustering"
        if clustering.exists():
            for h5ad in sorted(clustering.glob("*/clusters.h5ad")):
                _emit(os.path.relpath(h5ad, self.base_dir),
                      "clustering", "h5ad", "clusters.h5ad")
            for ss in sorted(clustering.glob("*/stage_stats.json")):
                _emit(os.path.relpath(ss, self.base_dir),
                      "qc", "json", "stage_stats")
        return found

    @staticmethod
    def _datasets_from_config(resolved_config: dict) -> list[dict]:
        """Build contract DatasetRef dicts from the resolved config (E2).

        Defensive against the several dataset dict shapes the config may use
        ({id/dataset_id/name}, {bams/bam_paths}).
        """
        raw = (resolved_config.get("directories", {}) or {}).get("datasets", []) or []
        out: list[dict] = []
        for ds in raw:
            if isinstance(ds, dict):
                ds_id = ds.get("id") or ds.get("dataset_id") or ds.get("name")
                bams = ds.get("bams") or ds.get("bam_paths") or []
                label = ds.get("label")
                if ds_id is None:
                    continue
                out.append({
                    "dataset_id": str(ds_id),
                    "bam_paths": [str(b) for b in bams] if isinstance(bams, (list, tuple)) else [str(bams)],
                    "label": label,
                })
            elif isinstance(ds, str):
                out.append({"dataset_id": ds, "bam_paths": [], "label": None})
        return out

    def write_manifest(
        self,
        resolved_config: dict,
        *,
        run_id: str | None = None,
        stratum_to_label: dict | None = None,
        entity_counts: dict | None = None,
        contract_version: str = "0.1.0",
    ) -> str:
        """Write ``run_manifest.json`` (E2) at the run root and return its path.

        The written object conforms to ``peakatail_contract.models.RunManifest``
        (frozen v0.1.0), i.e. it parses+validates cleanly there:
          * ``run_id`` / ``root`` — run identity (root = absolute run dir).
          * ``contract_version`` (semver) — the hub validates against this.
          * ``datasets`` — DatasetRef list derived from the resolved config.
          * ``artifacts`` — registered + auto-discovered Artifact dicts
            (path/stage/format/schema_name/schema_version/entity_counts).
          * ``resolved_config`` — the RESOLVED run config (B0), never argparse
            defaults; the ONLY config surface the hub reads.
          * ``stratum_to_label`` — full celltype labels for truncated stratum
            dir names (fixes D10's 48-char truncation join break).
          * ``entity_counts`` — run-level totals for the QC funnel.

        ``id_grammar`` is written as an extra, informational field (ignored by
        the contract model, which does not forbid extras).
        """
        base = Path(self.base_dir)
        artifacts = list(self._artifacts) + self._auto_discover_artifacts()
        datasets = self._datasets_from_config(resolved_config)
        counts = dict(entity_counts or {})
        counts.setdefault("n_datasets", len(datasets))
        manifest = {
            "run_id": run_id or base.name or str(base),
            "root": str(base.resolve()),
            "contract_version": contract_version,
            "timestamp": datetime.now().isoformat(),
            "datasets": datasets,
            "artifacts": artifacts,
            "resolved_config": resolved_config,
            "stratum_to_label": stratum_to_label or {},
            "entity_counts": counts,
            # informational only (not a contract field; extras are ignored)
            "id_grammar": {
                "pas_uid": "chrom:pos:strand (pos = end-1 on +, start on -)",
                "cell_uid": "{dataset_id}:{barcode}",
                "cluster_uid": "{dataset_id}:{leiden}",
                "canonical_cluster": "shared int across datasets (obs column)",
                "finding_uid": "{arm}:{celltype}:{gene_id}",
            },
        }
        path = os.path.join(self.base_dir, "run_manifest.json")
        with atomic_write(path) as f:
            json.dump(manifest, f, indent=2, default=str)
        return path


def build_resolved_run_config() -> dict:
    """Assemble the fully-resolved run configuration actually in effect.

    Bug B0: ``run_config.json`` used to serialize ``vars(args)`` — the argparse
    namespace — which holds *defaults* for everything supplied via YAML or
    ``set_directory_config`` (atlas path, gtf, datasets, atlas_distance, …). So
    a run that snapped to an atlas recorded ``atlas: null``, making the config
    non-reproducible and misleading the data controller.

    This reads the resolved module-level config singletons that the pipeline
    body reads from, so the persisted config matches what actually ran. Values
    are grouped by their source and JSON-safe (Paths → str via the caller's
    ``default=str``). Robust to missing attributes.
    """
    from ema.config import (
        directory_config as _dc,
        variable_config as _vc,
        filter_config as _fc,
        args as _args,
    )

    def _get(obj, name, default=None):
        try:
            return getattr(obj, name, default)
        except Exception:  # lazy proxies may raise on unset attrs
            return default

    resolved: dict = {}

    # --- directories / inputs (set via set_directory_config, NOT on args) ---
    resolved["directories"] = {
        "output_dir": str(_get(_dc, "output_dir", "")),
        "bam_dir": _get(_dc, "bam_dir"),
        "gtf_dir": _get(_dc, "gtf_dir"),
        "atlas": _get(_dc, "atlas"),
        "atlas_distance": _get(_dc, "atlas_distance"),
        "datasets": _get(_dc, "datasets", []),
        "filenames": _get(_dc, "filenames", {}),
    }

    # --- resolved scalar knobs the pipeline body actually reads ---
    resolved["variables"] = {
        k: _get(_vc, k)
        for k in (
            "seqlen", "cb_len", "barcode_tag", "default_threshold",
            "merge_len", "min_pas_spacing", "min_pas_prominence",
        )
    }
    resolved["filters"] = {
        k: _get(_fc, k)
        for k in ("min_read", "min_cells", "min_genes", "min_pas_per_cell")
    }

    # --- remaining argparse fields (strategy, thresholds, tiles, …) ---
    # Kept for completeness, but under a namespaced key so the resolved
    # directory/variable/filter values above are unambiguous. Filter out
    # private/callable entries.
    try:
        ns = vars(_args._get()) if hasattr(_args, "_get") else vars(_args)
    except Exception:
        ns = {}
    resolved["args"] = {
        k: v for k, v in ns.items()
        if not k.startswith("_") and not callable(v)
    }

    return resolved


def _concat_beds(srcs: Iterable[Path | str], dst: Path) -> None:
    """Concatenate BED files into ``dst``, skipping blank lines.

    Caller decides ordering (pos-then-neg, alphabetical, etc.) — this
    function preserves it.
    """
    with atomic_write(dst) as out:
        for src in srcs:
            with open(src) as f:
                # Bulk read + filter + a single writelines() call per file
                # instead of one out.write() per line — same bytes (each
                # surviving line, including its own original line ending,
                # is passed through unmodified), fewer Python-level I/O
                # calls on files with many peak rows.
                out.writelines(line for line in f if line.strip())


def write_raw_peak_outputs(
    output_dir: Path,
    pos_beds_by_ds: dict[str, list[str]],
    neg_beds_by_ds: dict[str, list[str]],
    pos_mtxs_by_ds: dict[str, list[str]] | None = None,
    neg_mtxs_by_ds: dict[str, list[str]] | None = None,
    cb_tsvs_by_ds: dict[str, str] | None = None,
) -> None:
    """Persist the **raw** peak-calling outputs per dataset, before any filter.

    Files written under ``<run>/01_peak_calling/<ds>/raw/``:

      * ``pos.bed`` / ``neg.bed`` — raw per-strand BEDs (concatenated across
        BAM replicates for that dataset)
      * ``pas.bed`` — union (pos+neg) raw PAS BED
      * ``pos.mtx`` / ``neg.mtx`` — raw count matrices (one column per
        unfiltered cell barcode)
      * ``cb.tsv`` — raw cell-barcode index aligned to the MTX columns
    """
    import shutil
    from ema.config import directory_config

    for ds_id in set(pos_beds_by_ds) | set(neg_beds_by_ds):
        raw_dir = directory_config.raw_dir_for(ds_id)
        raw_dir.mkdir(parents=True, exist_ok=True)
        if pos_beds_by_ds.get(ds_id):
            _concat_beds(pos_beds_by_ds[ds_id], raw_dir / "pos.bed")
        if neg_beds_by_ds.get(ds_id):
            _concat_beds(neg_beds_by_ds[ds_id], raw_dir / "neg.bed")
        _concat_beds(
            list(pos_beds_by_ds.get(ds_id, [])) + list(neg_beds_by_ds.get(ds_id, [])),
            raw_dir / "pas.bed",
        )
        if pos_mtxs_by_ds and pos_mtxs_by_ds.get(ds_id):
            # Replicate-level MTX files; we keep the first as the canonical
            # "raw" matrix for the dataset.  Multi-replicate merging happens
            # in the count matrix concat step later.
            shutil.copy(pos_mtxs_by_ds[ds_id][0], raw_dir / "pos.mtx")
        if neg_mtxs_by_ds and neg_mtxs_by_ds.get(ds_id):
            shutil.copy(neg_mtxs_by_ds[ds_id][0], raw_dir / "neg.mtx")
        if cb_tsvs_by_ds and cb_tsvs_by_ds.get(ds_id):
            shutil.copy(cb_tsvs_by_ds[ds_id], raw_dir / "cb.tsv")
        log.info("Persisted raw peak-calling outputs for %r at %s", ds_id, raw_dir)


def write_filtered_cb(
    output_dir: Path,
    dataset_id: str,
    filtered_barcodes: list[str] | "Iterable[str]",
    min_read: int,
) -> Path:
    """Persist the list of cell barcodes that passed the ``min_read`` filter.

    Writes ``<run>/02_cb_filter/<ds>/filtered_cb.tsv`` with a single column
    of barcodes plus a header ``barcode\tmin_read=<n>``.  Useful for
    downstream tools that need to subset other data to the same cells.
    """
    from ema.config import directory_config

    dst = directory_config.filtered_cb_for(dataset_id)
    with atomic_write(dst) as f:
        f.write(f"barcode\tmin_read={min_read}\n")
        for cb in filtered_barcodes:
            f.write(f"{cb}\n")
    log.info(
        "Persisted filtered cell barcodes for %r at %s",
        dataset_id, dst,
    )
    return dst


def write_annotated_matrix(
    output_dir: Path,
    dataset_id: str,
    sparse_matrix,
    pas_ids,
    cells: list[str],
) -> Path:
    """Persist the annotated count matrix (post PAS-gene join) for one dataset.

    Writes:
      * ``<run>/05_annotated_matrix/<ds>/annotated_matrix.mtx`` — MatrixMarket
        sparse matrix, rows = PAS (with gene assignment), cols = cells.
      * ``<run>/05_annotated_matrix/<ds>/annotated_pas_ids.tsv`` — row index.
      * ``<run>/05_annotated_matrix/<ds>/annotated_cells.tsv`` — column index.
    """
    import scipy.io as _sci
    from ema.config import directory_config

    # The .mtx and annotated_pas_ids.tsv written below are a POSITIONAL join:
    # line k of the id file names row k of the matrix.  Persisting a matrix
    # whose height differs from the id list silently mis-keys every downstream
    # consumer (this is exactly regression 04e0b3a), so fail loud here too.
    if sparse_matrix.shape[0] != len(pas_ids):
        raise ValueError(
            f"annotated matrix for {dataset_id!r} has {sparse_matrix.shape[0]} "
            f"rows but {len(pas_ids)} PAS IDs — annotated_pas_ids.tsv is the "
            "row index of annotated_matrix.mtx and must be the same length."
        )
    if sparse_matrix.shape[1] != len(cells):
        raise ValueError(
            f"annotated matrix for {dataset_id!r} has {sparse_matrix.shape[1]} "
            f"columns but {len(cells)} barcodes — annotated_cells.tsv is the "
            "column index of annotated_matrix.mtx and must be the same length."
        )

    mtx_path = directory_config.annotated_matrix_for(dataset_id)
    mtx_path.parent.mkdir(parents=True, exist_ok=True)
    _sci.mmwrite(str(mtx_path), sparse_matrix.astype(int), field="integer")
    with atomic_write(directory_config.annotated_pas_ids_for(dataset_id)) as f:
        f.write("pas_id\n" + "\n".join(str(p) for p in pas_ids) + "\n")
    with atomic_write(directory_config.annotated_cells_for(dataset_id)) as f:
        f.write("barcode\n" + "\n".join(cells) + "\n")
    log.info(
        "Persisted annotated count matrix for %r at %s (%d PAS × %d cells)",
        dataset_id, mtx_path, len(pas_ids), len(cells),
    )
    return mtx_path


def write_preprocessed_h5ad(
    output_dir: Path,
    dataset_id: str,
    adata,
) -> Path:
    """Persist the preprocessed (filtered) AnnData before clustering.

    Saves to ``<run>/06_preprocessing/<ds>/preprocessed.h5ad`` so users can
    inspect the filtered count matrix shape independently of cluster
    labels (which live in ``07_clustering/<ds>/clusters.h5ad``).
    """
    from ema.config import directory_config

    dst = directory_config.preprocessed_h5ad_for(dataset_id)
    dst.parent.mkdir(parents=True, exist_ok=True)
    adata.write(dst)
    log.info(
        "Persisted preprocessed AnnData for %r at %s (%d cells × %d PAS)",
        dataset_id, dst, adata.n_obs, adata.n_vars,
    )
    return dst


def write_per_dataset_beds(
    output_dir: Path,
    pos_beds_by_ds: dict[str, list[str]],
    neg_beds_by_ds: dict[str, list[str]],
) -> dict[str, Path]:
    """Write canonical posbed / negbed / pasbed per dataset.

    Files land under ``<run>/01_peak_calling/<ds>/``.

    Args:
        output_dir: Pipeline run-root (kept for API compat; not used directly).
        pos_beds_by_ds: ``{dataset_id: [path, ...]}`` positive-strand BEDs.
        neg_beds_by_ds: same for negative strand.

    Returns:
        ``{dataset_id: pasbed_path}`` for callers that need it (e.g. the
        annotation step needs the pasbed to build annotatedpas.bed).
    """
    from ema.config import directory_config

    pasbeds: dict[str, Path] = {}
    for ds_id in set(pos_beds_by_ds) | set(neg_beds_by_ds):
        posbed = directory_config.posbed_for(ds_id)
        posbed.parent.mkdir(parents=True, exist_ok=True)
        if pos_beds_by_ds.get(ds_id):
            _concat_beds(pos_beds_by_ds[ds_id], posbed)
        negbed = directory_config.negbed_for(ds_id)
        if neg_beds_by_ds.get(ds_id):
            _concat_beds(neg_beds_by_ds[ds_id], negbed)
        pasbed = directory_config.pasbed_for(ds_id)
        _concat_beds(
            list(pos_beds_by_ds.get(ds_id, [])) + list(neg_beds_by_ds.get(ds_id, [])),
            pasbed,
        )
        pasbeds[ds_id] = pasbed
        log.info("Persisted per-dataset BEDs for %r at %s", ds_id, posbed.parent)
    return pasbeds


def write_pas_gene_artifacts(
    output_dir: Path,
    dataset_id: str,
    pas_ids,
    gene_ids,
    *,
    atlas_of: dict | None = None,
    ip_of: dict | None = None,
    in_3utr_of: dict | None = None,
) -> tuple[Path, Path]:
    """Write ``pas_gene.tsv`` and ``annotatedpas.bed`` for one dataset.

    ``pas_gene.tsv`` is an explicit two-column mapping
    (``pas_id\tgene_id``) — useful for downstream tools that don't want
    to read the full AnnData.

    ``annotatedpas.bed`` extends ``pasbed.bed`` with a trailing gene_id
    column, and (D9) three further trailing status columns --
    ``atlas_match``, ``atlas_distance_bp``, ``internal_priming`` -- plus a
    fourth, ``in_3utr`` -- so the "keep everything, annotate with
    match/no-match" atlas-snap / internal-priming / 3'UTR-membership labels
    are queryable straight off the BED, not just the PAS ledger. This BED
    covers EVERY PAS this dataset's peak-calling produced (one row per
    ``pasbed.bed`` line), independent of whether a PAS made it into the
    clustering matrix -- a PAS excluded from clustering (see
    ``ema.downstream_runner.run_one_dataset_downstream``'s
    ``exclude_pas_ids``) still gets its full row here, just like every
    other PAS. The first 6 columns stay plain BED6 so existing BED
    consumers (bedtools, etc.) still parse the file; extra columns are
    appended, never inserted. It depends on ``pasbed.bed`` having been
    written first (see :func:`write_per_dataset_beds`); if the pasbed
    isn't on disk yet we skip the BED part and warn.

    Args:
        output_dir: Pipeline run-root (kept for API compat; not used directly).
        dataset_id: Dataset name.
        pas_ids: Array-like of PAS IDs (parallel to ``gene_ids``).
        gene_ids: Array-like of gene IDs aligned to ``pas_ids``.
        atlas_of: Optional ``{pas_id: (atlas_match, atlas_distance_bp)}`` --
            see ``ema.datasets.atlas_snap.snap_beds_to_atlas``'s
            ``atlas_status.tsv`` sidecar. PAS absent from the map (atlas
            didn't run, or this PAS predates atlas) get ``""`` for both
            columns.
        ip_of: Optional ``{pas_id: internal_priming_bool}`` -- see
            ``ema.experimental.peak_filters.apply_filters``'s
            ``"internal_priming_flags"`` stats key. PAS absent from the
            map get ``""``.
        in_3utr_of: Optional ``{pas_id: in_3utr_bool}`` -- whether the PAS
            overlaps an annotated transcript 3'UTR (see
            ``ema.experimental.peak_filters.label_pas_in_bed``). PAS absent
            from the map (label didn't run) get ``""``.

    Returns:
        ``(pas_gene_tsv_path, annotatedpas_bed_path)``.  The BED path
        may not exist if pasbed wasn't on disk.
    """
    import numpy as np  # local import — heavy module
    import pandas as pd  # local import — heavy module
    from ema.config import directory_config

    atlas_of = atlas_of or {}
    ip_of = ip_of or {}
    in_3utr_of = in_3utr_of or {}

    # pas_ids/gene_ids are a POSITIONAL pair (gene_ids[k] annotates pas_ids[k]).
    # Check the length explicitly and coerce a pandas Series to its values --
    # a Series would otherwise be re-aligned by ITS OWN index inside the
    # DataFrame constructor below, silently pairing PAS with the wrong gene.
    pas_ids = np.asarray(pas_ids)
    gene_ids = np.asarray(gene_ids)
    if len(pas_ids) != len(gene_ids):
        raise ValueError(
            f"pas_gene mapping for {dataset_id!r} has {len(pas_ids)} PAS IDs but "
            f"{len(gene_ids)} gene IDs — they index the same annotated rows."
        )

    # Both writes below go through atomic_write() (temp file in the SAME
    # directory + os.replace) -- never write the final path directly. See
    # that helper's docstring for the corruption window this closes.
    pas_gene_tsv = directory_config.pas_gene_for(dataset_id)
    with atomic_write(pas_gene_tsv) as dst:
        pd.DataFrame({"pas_id": pas_ids, "gene_id": gene_ids}).to_csv(
            dst, sep="\t", index=False,
        )

    annot_bed = directory_config.annotatedpas_for(dataset_id)
    annot_bed.parent.mkdir(parents=True, exist_ok=True)
    pasbed = directory_config.pasbed_for(dataset_id)
    if pasbed.exists():
        lookup = dict(zip([str(p) for p in pas_ids], gene_ids))
        with open(pasbed) as src, atomic_write(annot_bed) as dst:
            for line in src:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 4:
                    pas_id = parts[3]
                    gid = lookup.get(pas_id, "")
                    atlas_match, atlas_distance_bp = atlas_of.get(pas_id, ("", ""))
                    ip_flag = ip_of.get(pas_id, "")
                    in_3utr_flag = in_3utr_of.get(pas_id, "")
                    extra = [gid, str(atlas_match), str(atlas_distance_bp), str(ip_flag), str(in_3utr_flag)]
                    dst.write("\t".join(parts + extra) + "\n")
    else:
        log.warning(
            "annotatedpas.bed for %r skipped — pasbed.bed not on disk at %s",
            dataset_id, pasbed,
        )
    log.info(
        "Persisted PAS->gene mapping (%d rows) at %s",
        len(pas_ids), pas_gene_tsv.parent,
    )
    return pas_gene_tsv, annot_bed
