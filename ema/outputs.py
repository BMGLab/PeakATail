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
from datetime import datetime
from pathlib import Path
from typing import Iterable

log = logging.getLogger(__name__)


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
        with open(path, "w") as f:
            json.dump(stats, f, indent=2)

    def save_run_config(self, args_dict: dict) -> None:
        """Persist the resolved run configuration at the run root."""
        config = {"timestamp": datetime.now().isoformat(), **args_dict}
        path = os.path.join(self.base_dir, "run_config.json")
        with open(path, "w") as f:
            json.dump(config, f, indent=2)


def _concat_beds(srcs: Iterable[Path | str], dst: Path) -> None:
    """Concatenate BED files into ``dst``, skipping blank lines.

    Caller decides ordering (pos-then-neg, alphabetical, etc.) — this
    function preserves it.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w") as out:
        for src in srcs:
            with open(src) as f:
                for line in f:
                    if line.strip():
                        out.write(line)


def write_raw_peak_outputs(
    output_dir: Path,
    pos_beds_by_ds: dict[str, list[str]],
    neg_beds_by_ds: dict[str, list[str]],
    pos_mtxs_by_ds: dict[str, list[str]] | None = None,
    neg_mtxs_by_ds: dict[str, list[str]] | None = None,
    cb_tsvs_by_ds: dict[str, str] | None = None,
) -> None:
    """Persist the **raw** peak-calling outputs per dataset, before any filter.

    Files written under ``<run>/per_dataset/<ds>/raw/``:

      * ``pos.bed`` / ``neg.bed`` — raw per-strand BEDs (concatenated across
        BAM replicates for that dataset)
      * ``pas.bed`` — union (pos+neg) raw PAS BED
      * ``pos.mtx`` / ``neg.mtx`` — raw count matrices (one column per
        unfiltered cell barcode)
      * ``cb.tsv`` — raw cell-barcode index aligned to the MTX columns

    Downstream consumers reading e.g. ``per_dataset/<ds>/pas.bed`` always
    get the unfiltered call set, while ``per_dataset/<ds>/pasbed.bed``
    represents the filtered/post-merge view (see
    :func:`write_per_dataset_beds`).
    """
    import shutil

    for ds_id in set(pos_beds_by_ds) | set(neg_beds_by_ds):
        raw_dir = output_dir / "per_dataset" / ds_id / "raw"
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

    Writes ``<run>/per_dataset/<ds>/filtered_cb.tsv`` with a single column
    of barcodes plus a header ``barcode\tmin_read=<n>``.  Useful for
    downstream tools that need to subset other data to the same cells.
    """
    ds_root = output_dir / "per_dataset" / dataset_id
    ds_root.mkdir(parents=True, exist_ok=True)
    dst = ds_root / "filtered_cb.tsv"
    with open(dst, "w") as f:
        f.write(f"barcode\tmin_read={min_read}\n")
        for cb in filtered_barcodes:
            f.write(f"{cb}\n")
    log.info(
        "Persisted filtered cell barcodes (%d kept) for %r at %s",
        sum(1 for _ in filtered_barcodes if True), dataset_id, dst,
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
      * ``<run>/per_dataset/<ds>/annotated_matrix.mtx`` — MatrixMarket
        sparse matrix, rows = PAS (with gene assignment), cols = cells.
      * ``<run>/per_dataset/<ds>/annotated_pas_ids.tsv`` — row index.
      * ``<run>/per_dataset/<ds>/annotated_cells.tsv`` — column index.
    """
    import scipy.io as _sci

    ds_root = output_dir / "per_dataset" / dataset_id
    ds_root.mkdir(parents=True, exist_ok=True)
    mtx_path = ds_root / "annotated_matrix.mtx"
    _sci.mmwrite(str(mtx_path), sparse_matrix.astype(int), field="integer")
    (ds_root / "annotated_pas_ids.tsv").write_text(
        "pas_id\n" + "\n".join(str(p) for p in pas_ids) + "\n"
    )
    (ds_root / "annotated_cells.tsv").write_text(
        "barcode\n" + "\n".join(cells) + "\n"
    )
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

    Saves to ``<run>/per_dataset/<ds>/preprocessed.h5ad`` so users can
    inspect the filtered count matrix shape independently of cluster
    labels (which live in ``clusters.h5ad``).
    """
    ds_root = output_dir / "per_dataset" / dataset_id
    ds_root.mkdir(parents=True, exist_ok=True)
    dst = ds_root / "preprocessed.h5ad"
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

    Args:
        output_dir: Pipeline run-root.
        pos_beds_by_ds: ``{dataset_id: [path, ...]}`` positive-strand BEDs.
        neg_beds_by_ds: same for negative strand.

    Returns:
        ``{dataset_id: pasbed_path}`` for callers that need it (e.g. the
        annotation step needs the pasbed to build annotatedpas.bed).
    """
    pasbeds: dict[str, Path] = {}
    for ds_id in set(pos_beds_by_ds) | set(neg_beds_by_ds):
        ds_root = output_dir / "per_dataset" / ds_id
        ds_root.mkdir(parents=True, exist_ok=True)
        if pos_beds_by_ds.get(ds_id):
            _concat_beds(pos_beds_by_ds[ds_id], ds_root / "posbed.bed")
        if neg_beds_by_ds.get(ds_id):
            _concat_beds(neg_beds_by_ds[ds_id], ds_root / "negbed.bed")
        pasbed = ds_root / "pasbed.bed"
        _concat_beds(
            list(pos_beds_by_ds.get(ds_id, [])) + list(neg_beds_by_ds.get(ds_id, [])),
            pasbed,
        )
        pasbeds[ds_id] = pasbed
        log.info("Persisted per-dataset BEDs for %r at %s", ds_id, ds_root)
    return pasbeds


def write_pas_gene_artifacts(
    output_dir: Path,
    dataset_id: str,
    pas_ids,
    gene_ids,
) -> tuple[Path, Path]:
    """Write ``pas_gene.tsv`` and ``annotatedpas.bed`` for one dataset.

    ``pas_gene.tsv`` is an explicit two-column mapping
    (``pas_id\tgene_id``) — useful for downstream tools that don't want
    to read the full AnnData.

    ``annotatedpas.bed`` extends ``pasbed.bed`` with a trailing gene_id
    column.  It depends on ``pasbed.bed`` having been written first
    (see :func:`write_per_dataset_beds`); if the pasbed isn't on disk
    yet we skip the BED part and warn.

    Args:
        output_dir: Pipeline run-root.
        dataset_id: Dataset name (matches ``per_dataset/<ds>/``).
        pas_ids: Array-like of PAS IDs (parallel to ``gene_ids``).
        gene_ids: Array-like of gene IDs aligned to ``pas_ids``.

    Returns:
        ``(pas_gene_tsv_path, annotatedpas_bed_path)``.  The BED path
        may not exist if pasbed wasn't on disk.
    """
    import pandas as pd  # local import — heavy module

    ds_root = output_dir / "per_dataset" / dataset_id
    ds_root.mkdir(parents=True, exist_ok=True)

    pas_gene_tsv = ds_root / "pas_gene.tsv"
    pd.DataFrame({"pas_id": pas_ids, "gene_id": gene_ids}).to_csv(
        pas_gene_tsv, sep="\t", index=False,
    )

    annot_bed = ds_root / "annotatedpas.bed"
    pasbed = ds_root / "pasbed.bed"
    if pasbed.exists():
        lookup = dict(zip([str(p) for p in pas_ids], gene_ids))
        with open(pasbed) as src, open(annot_bed, "w") as dst:
            for line in src:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 4:
                    gid = lookup.get(parts[3], "")
                    dst.write("\t".join(parts) + "\t" + gid + "\n")
    else:
        log.warning(
            "annotatedpas.bed for %r skipped — pasbed.bed not on disk at %s",
            dataset_id, pasbed,
        )
    log.info(
        "Persisted PAS->gene mapping (%d rows) at %s",
        len(pas_ids), ds_root,
    )
    return pas_gene_tsv, annot_bed
