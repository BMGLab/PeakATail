"""Clustering evaluation and comparison utilities.

Provides functions for:
- Comparing two clusterings (ARI, AMI)
- Computing cluster quality metrics (silhouette score)
- Finding marker peaks per cluster (Wilcoxon rank-sum)
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from typing import Dict, Optional, Tuple


def compare_clusterings(labels_a, labels_b) -> Dict[str, float]:
    """Compare two clustering results using standard metrics.

    Args:
        labels_a: First set of cluster labels (array-like).
        labels_b: Second set of cluster labels (array-like).
            Must have the same length as labels_a.

    Returns:
        Dictionary with:
            - ari: Adjusted Rand Index (-1 to 1, 1 = perfect agreement)
            - ami: Adjusted Mutual Information (0 to 1, 1 = perfect agreement)
    """
    labels_a = np.asarray(labels_a)
    labels_b = np.asarray(labels_b)

    if len(labels_a) != len(labels_b):
        raise ValueError(
            f"Label arrays must have same length. "
            f"Got {len(labels_a)} and {len(labels_b)}."
        )

    ari = adjusted_rand_score(labels_a, labels_b)
    ami = adjusted_mutual_info_score(labels_a, labels_b)

    return {"ari": ari, "ami": ami}


def compute_cluster_quality(adata: ad.AnnData,
                            labels_key: str = 'leiden',
                            use_rep: Optional[str] = None,
                            metric: str = 'euclidean',
                            sample_size: Optional[int] = None,
                            random_seed: int = 42) -> Dict[str, float]:
    """Compute silhouette score for a clustering.

    Args:
        adata: AnnData object with clustering labels.
        labels_key: Key in adata.obs for cluster labels (default: 'leiden').
        use_rep: Representation to use. If None, auto-detects
            'X_lsi' or 'X_pca'.
        metric: Distance metric for silhouette (default: 'euclidean').
            Use 'cosine' for TF-IDF/LSI embeddings.
        sample_size: If set, subsample this many cells for faster
            computation.
        random_seed: Random seed for subsampling.

    Returns:
        Dictionary with:
            - silhouette_mean: Mean silhouette score (-1 to 1)
            - n_clusters: Number of clusters
            - n_cells: Number of cells used
    """
    labels = adata.obs[labels_key].values

    # Auto-detect representation
    if use_rep is None:
        if 'X_lsi' in adata.obsm:
            use_rep = 'X_lsi'
            metric = 'cosine'
        elif 'X_pca' in adata.obsm:
            use_rep = 'X_pca'
        else:
            raise ValueError(
                "No dimensionality reduction found. "
                "Run reduce_dims() first or specify use_rep."
            )

    X = adata.obsm[use_rep]

    # Subsample for large datasets
    if sample_size is not None and sample_size < X.shape[0]:
        rng = np.random.RandomState(random_seed)
        indices = rng.choice(X.shape[0], size=sample_size, replace=False)
        X = X[indices]
        labels = labels[indices]

    # Need at least 2 clusters for silhouette
    n_clusters = len(np.unique(labels))
    if n_clusters < 2:
        return {
            "silhouette_mean": float('nan'),
            "n_clusters": n_clusters,
            "n_cells": len(labels),
        }

    sil = silhouette_score(X, labels, metric=metric)

    return {
        "silhouette_mean": sil,
        "n_clusters": n_clusters,
        "n_cells": len(labels),
    }


def find_marker_peaks(adata: ad.AnnData,
                      groupby: str = 'leiden',
                      method: str = 'wilcoxon',
                      n_genes: int = 50) -> pd.DataFrame:
    """Find marker PAS per cluster using rank-based test.

    Uses scanpy's rank_genes_groups, which works on any features
    (genes or PAS peaks).

    Args:
        adata: AnnData object with clustering labels. The .X matrix
            should contain normalized values (not raw counts) for
            proper statistical testing.
        groupby: Key in adata.obs for cluster labels (default: 'leiden').
        method: Statistical test method. One of 'wilcoxon', 't-test',
            't-test_overestim_var' (default: 'wilcoxon').
        n_genes: Number of top markers to report per cluster
            (default: 50).

    Returns:
        DataFrame with columns: group, names, scores, pvals, pvals_adj,
        logfoldchanges. Rows = top markers per cluster.
    """
    # Make a copy to avoid modifying original
    adata_test = adata.copy()

    sc.tl.rank_genes_groups(
        adata_test,
        groupby=groupby,
        method=method,
        n_genes=n_genes,
    )

    # Extract results into a DataFrame
    result = sc.get.rank_genes_groups_df(adata_test, group=None)
    return result
