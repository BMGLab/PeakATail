"""A1 core: marker-signature cell-type assignment + PAS-cluster concordance.

Given a per-cell signature-score matrix (cells x signatures, e.g. from
``scanpy.tl.score_genes`` over marker gene sets) and a cell->cluster labelling,
assign a cell type to each cluster by TWO cross-checked methods and report their
agreement:

    M1 (cluster-argmax): z-score each signature across cells, average per
        cluster, argmax the per-cluster mean → the cluster's cell type.
    M2 (cell-majority): z-score, argmax per cell → each cell's tentative type,
        then take the majority type within each cluster.

Two methods that disagree flag an ambiguous cluster. Concordance against the
PAS-based clustering is the standard ARI/AMI.

These functions are pure numpy/pandas/sklearn and unit-tested on synthetic
score matrices. The upstream GEX I/O (STARsolo pooling) + scanpy clustering /
``score_genes`` that PRODUCE the score matrix are intentionally not in the
package yet — they need real-data validation first (see the standalone
prototype scripts/gex_celltyping.py).
"""
from __future__ import annotations

from typing import Any

__all__ = ["zscore_across_cells", "assign_celltypes", "concordance"]


def zscore_across_cells(score_df):
    """Z-score each signature column across cells (mean 0, unit std).

    Assignment should reflect RELATIVE enrichment of a signature across cells,
    not its absolute magnitude (some signatures are intrinsically higher). A
    constant column (zero variance) is left at 0 rather than dividing by zero.
    """
    import numpy as np

    mean = score_df.mean(axis=0)
    std = score_df.std(axis=0, ddof=0).replace(0.0, np.nan)
    z = (score_df - mean) / std
    return z.fillna(0.0)


def assign_celltypes(score_df, cluster_labels) -> dict[str, Any]:
    """Assign a cell type per cluster via M1 (cluster-argmax) and M2 (majority).

    Args:
        score_df: DataFrame (cells x signatures); columns are cell-type names.
        cluster_labels: per-cell cluster label, aligned to ``score_df.index``
            (Series or array of length n_cells).

    Returns:
        dict with:
            ``m1_by_cluster``   {cluster: celltype} from cluster-mean argmax
            ``m2_by_cluster``   {cluster: celltype} from per-cell majority
            ``cell_argmax``     {cell: celltype} per-cell argmax (M2 basis)
            ``agreement``       fraction of clusters where M1 == M2 (0..1;
                                NaN if there are no clusters)
            ``n_clusters``      number of distinct clusters

    Raises:
        ValueError: if ``score_df`` has no columns, or lengths mismatch.
    """
    import numpy as np
    import pandas as pd

    if score_df.shape[1] == 0:
        raise ValueError("assign_celltypes: score_df has no signature columns")
    labels = pd.Series(list(cluster_labels), index=score_df.index)
    if len(labels) != len(score_df):
        raise ValueError(
            f"assign_celltypes: {len(labels)} cluster labels but {len(score_df)} cells"
        )

    z = zscore_across_cells(score_df)
    sigs = list(z.columns)

    # M1: per-cluster mean z, argmax.
    cluster_mean = z.groupby(labels).mean()
    m1_by_cluster = {
        clu: sigs[int(np.argmax(row.values))]
        for clu, row in cluster_mean.iterrows()
    }

    # M2: per-cell argmax → majority within cluster.
    cell_arg = z.values.argmax(axis=1)
    cell_argmax = {cell: sigs[int(a)] for cell, a in zip(z.index, cell_arg)}
    cell_type_series = pd.Series([sigs[int(a)] for a in cell_arg], index=z.index)
    m2_by_cluster = {}
    for clu, idx in cell_type_series.groupby(labels).groups.items():
        # majority vote; ties broken by first-seen signature order for determinism.
        counts = cell_type_series.loc[idx].value_counts()
        top = counts.max()
        winners = [s for s in sigs if counts.get(s, 0) == top]
        m2_by_cluster[clu] = winners[0]

    clusters = list(cluster_mean.index)
    n = len(clusters)
    agree = (sum(1 for c in clusters if m1_by_cluster[c] == m2_by_cluster[c]) / n
             if n else float("nan"))

    return {
        "m1_by_cluster": m1_by_cluster,
        "m2_by_cluster": m2_by_cluster,
        "cell_argmax": cell_argmax,
        "agreement": agree,
        "n_clusters": n,
    }


def concordance(labels_a, labels_b) -> dict[str, float]:
    """ARI + AMI between two same-length cell labellings (e.g. GEX vs PAS).

    Args:
        labels_a, labels_b: parallel per-cell label sequences (any hashable
            labels; only co-membership matters).

    Returns:
        ``{"ARI": float, "AMI": float, "n_cells": int}``.

    Raises:
        ValueError: if the two label sequences differ in length or are empty.
    """
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

    a = list(labels_a)
    b = list(labels_b)
    if len(a) != len(b):
        raise ValueError(f"concordance: length mismatch {len(a)} vs {len(b)}")
    if not a:
        raise ValueError("concordance: empty label sequences")
    return {
        "ARI": round(float(adjusted_rand_score(a, b)), 4),
        "AMI": round(float(adjusted_mutual_info_score(a, b)), 4),
        "n_cells": len(a),
    }
