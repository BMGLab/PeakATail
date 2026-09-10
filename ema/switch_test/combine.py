"""A2: combine stage-labelled clustered h5ads for grouped switch analysis.

``peakatail switch {diff,length}`` compare cells grouped by an obs key (``--cluster-key``).
To ask a *cross-dataset* question — "does APA differ across disease stages /
timepoints, within each cell type?" — the per-dataset ``clusters.h5ad`` files
first have to be stitched into one AnnData whose ``obs[group_key]`` carries the
stage label (and, optionally, split by an existing cell-type obs column so each
cell type is contrasted independently).

That stitching was the only genuinely-new glue in
``scripts/b3_stage_celltype_switch.py`` (everything else just calls the real
``peakatail switch`` statistics). This module promotes it to a package feature. The
pure :func:`combine_labeled` carries the logic and is unit-testable on tiny
synthetic AnnData; :func:`combine_to_dir` is the on-disk wrapper used by
``peakatail switch combine``.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

__all__ = ["slug", "combine_labeled", "combine_to_dir"]


def slug(s: str) -> str:
    """Filesystem-safe short slug for a group/cell-type label."""
    return re.sub(r"[^0-9A-Za-z]+", "_", str(s)).strip("_")[:48] or "group"


def combine_labeled(
    adatas: list,
    labels: list[str],
    *,
    group_key: str = "stage",
    split_key: str | None = None,
    min_cells: int = 1,
) -> dict[str, Any]:
    """Stitch labelled AnnData into per-group combined AnnData.

    Args:
        adatas: list of AnnData (one per dataset/stage), same var space or
            outer-joinable.
        labels: parallel list of group labels (e.g. stage names), one per
            adata. Stamped onto ``obs[group_key]``.
        group_key: obs column that will carry ``labels`` (default ``"stage"``).
        split_key: optional existing obs column (e.g. ``"celltype"``). When
            given, cells are partitioned by its value and each partition is
            combined independently; when None a single ``"__all__"`` group is
            produced.
        min_cells: drop any ``group_key`` value with fewer than this many cells
            (within a split); a split with fewer than 2 surviving group values
            is omitted (can't contrast).

    Returns:
        ``{split_value: combined_adata}`` — each combined AnnData has
        ``obs[group_key]`` populated and only groups with >= ``min_cells``
        cells, and only splits spanning >= 2 group values. ``var["gene_id"]``
        is preserved (merge="first").

    Raises:
        ValueError: if ``len(adatas) != len(labels)`` or ``adatas`` is empty.
    """
    import anndata as ad

    if not adatas:
        raise ValueError("combine_labeled: no adatas supplied")
    if len(adatas) != len(labels):
        raise ValueError(
            f"combine_labeled: {len(adatas)} adatas but {len(labels)} labels"
        )

    stamped = []
    for adata, label in zip(adatas, labels):
        a = adata.copy()
        a.obs[group_key] = str(label)
        stamped.append(a)

    # Determine split values.
    if split_key is None:
        splits = {"__all__": stamped}
    else:
        splits = {}
        for a in stamped:
            if split_key not in a.obs.columns:
                continue
            for val in a.obs[split_key].astype(str).unique():
                if val in ("nan", "None", ""):
                    continue
                sub = a[a.obs[split_key].astype(str) == val]
                if sub.n_obs > 0:
                    splits.setdefault(val, []).append(sub.copy())

    out: dict[str, Any] = {}
    for split_val, subs in splits.items():
        if not subs:
            continue
        comb = ad.concat(subs, join="outer", merge="first", label=None)
        # Drop group values below the min-cells floor.
        counts = comb.obs[group_key].value_counts()
        keep_groups = set(counts[counts >= min_cells].index)
        if len(keep_groups) < 2:
            continue  # need >= 2 groups to contrast
        mask = comb.obs[group_key].isin(keep_groups).values
        comb = comb[mask].copy()
        if comb.obs[group_key].nunique() < 2:
            continue
        out[split_val] = comb
    return out


def combine_to_dir(
    inputs: list[tuple[str, str]],
    out_dir: str | Path,
    *,
    group_key: str = "stage",
    split_key: str | None = None,
    min_cells: int = 1,
) -> dict[str, dict]:
    """Read labelled h5ads, combine, and write one h5ad per group.

    Args:
        inputs: list of ``(label, h5ad_path)`` pairs.
        out_dir: directory to write ``<slug>.h5ad`` files into.
        group_key/split_key/min_cells: see :func:`combine_labeled`.

    Returns:
        manifest dict ``{split_value: {slug, path, n_cells, groups}}``.

    Raises:
        FileNotFoundError: if an input h5ad is missing.
    """
    import anndata as ad

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    adatas, labels = [], []
    for label, path in inputs:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"input h5ad not found: {p}")
        adatas.append(ad.read_h5ad(p))
        labels.append(label)

    combined = combine_labeled(
        adatas, labels, group_key=group_key, split_key=split_key, min_cells=min_cells,
    )
    manifest: dict[str, dict] = {}
    for split_val, comb in combined.items():
        s = slug(split_val)
        path = out_dir / f"{s}.h5ad"
        comb.write_h5ad(path)
        groups = {str(k): int(v) for k, v in comb.obs[group_key].value_counts().items()}
        manifest[split_val] = {
            "slug": s, "path": str(path),
            "n_cells": int(comb.n_obs), "groups": groups,
        }
    return manifest
