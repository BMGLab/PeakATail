"""Marker-PAS selection for downstream differential testing.

Standard scRNA-seq workflow: after clustering, identify the top differentially
expressed features per cluster (markers) via Wilcoxon rank-sum, then restrict
heavy downstream tests (NB regression, PDUI per-isoform expansion) to the
**union** of those markers.

Why: testing all 20K PAS against all C(K, 2) cluster pairs is wasteful — most
PAS are not informative for cluster discrimination. Markers are by construction
the ones with cluster-specific signal, so they are exactly the ones worth
formal differential testing.
"""

from __future__ import annotations
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc


def select_marker_pas(
    adata: ad.AnnData,
    cluster_key: str = "leiden",
    top_n_per_cluster: int = 200,
    method: str = "wilcoxon",
    min_pct_in_cluster: float = 0.10,
) -> list[int]:
    """Return the union of top-N marker PAS across all clusters.

    Args:
        adata: AnnData (cells x PAS) with cluster labels in ``adata.obs[cluster_key]``.
        cluster_key: Column in ``adata.obs`` holding cluster labels (default: "leiden").
        top_n_per_cluster: Number of top markers to keep from each cluster.
        method: ``rank_genes_groups`` method (default: "wilcoxon").
        min_pct_in_cluster: Drop markers detected in fewer than this fraction
            of cells in their own cluster (default: 10%).

    Returns:
        Sorted list of integer PAS IDs (the union across clusters).
    """
    if cluster_key not in adata.obs:
        raise ValueError(f"adata.obs has no column '{cluster_key}'")

    # scanpy's rank_genes_groups reads ``adata.obs[groupby].cat.categories``
    # internally, so the groupby column MUST be categorical. canonical_cluster
    # (and other real cluster keys) arrive as int/object from clusters.h5ad,
    # which crashes with "Can only use .cat accessor with a 'category' dtype"
    # on pandas >=2. Coerce here (idempotent — a no-op if already categorical).
    import pandas as pd
    if not isinstance(adata.obs[cluster_key].dtype, pd.CategoricalDtype):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method=method,
        n_genes=top_n_per_cluster,
        use_raw=False,
    )

    selected: set[str] = set()
    names = adata.uns["rank_genes_groups"]["names"]
    pcts = adata.uns["rank_genes_groups"].get("pts")  # detection rate per cluster

    for cluster_name in names.dtype.names:
        cluster_top = np.asarray(names[cluster_name])
        if pcts is not None:
            cluster_pct = np.asarray(pcts[cluster_name])
            keep = cluster_pct >= min_pct_in_cluster
            cluster_top = cluster_top[keep[: len(cluster_top)]]
        for pas in cluster_top[:top_n_per_cluster]:
            selected.add(str(pas))

    # Convert to int when possible (downstream uses int PAS IDs)
    out = []
    for s in selected:
        try:
            out.append(int(s))
        except ValueError:
            out.append(s)
    out.sort(key=lambda x: (isinstance(x, str), x))
    return out


def restrict_count_matrix(
    count_matrix,
    marker_pas: list,
    axis: str = "rows",
):
    """Restrict a (PAS x cells) or (cells x PAS) count matrix to marker_pas only.

    Args:
        count_matrix: pandas DataFrame.
        marker_pas: list of PAS IDs to keep.
        axis: "rows" if PAS are along rows, "cols" if along columns.

    Returns:
        Filtered DataFrame containing only the requested PAS, in the given order.
    """
    keep = [p for p in marker_pas if p in (count_matrix.index if axis == "rows" else count_matrix.columns)]
    if axis == "rows":
        return count_matrix.loc[keep]
    return count_matrix.loc[:, keep]


def save_markers(marker_pas: list, output_path: str | Path) -> None:
    """Write the marker PAS list to a TSV file (one per line)."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("pas_id\n")
        for p in marker_pas:
            f.write(f"{p}\n")
