"""Clustering module for PAS peak count data.

Supports multiple clustering strategies via the strategy registry:
    - leiden_tfidf (default): TF-IDF + LSI + Leiden (scATAC-seq style)
    - leiden_libsize: Library-size normalization + PCA + Leiden
    - external: Import pre-computed cluster labels from CSV

Usage:
    from ema.clustering.clustering import clustering

    # Using default TF-IDF strategy
    adata = clustering(adata)

    # Using library-size strategy
    adata = clustering(adata, method='leiden_libsize')

    # Using external labels
    adata = clustering(adata, method='external',
                       external_clusters='/path/to/labels.csv')
"""

import os
import anndata as ad
import pandas as pd
from ema.clustering.strategies import get_strategy


def clustering(adata: ad.AnnData,
               method='leiden_tfidf',
               resolution=1.0,
               n_pcs=40,
               random_seed=42,
               external_clusters=None,
               outputpath=None,
               output_h5ad=None):
    """Cluster cells based on PAS peak count data.

    Args:
        adata: AnnData object with raw PAS counts (cells x peaks).
        method: Clustering strategy name. One of 'leiden_tfidf',
            'leiden_libsize', or 'external'.
        resolution: Leiden clustering resolution (default: 1.0).
        n_pcs: Number of dimensions for neighbor computation (default: 40).
        random_seed: Random seed for reproducibility (default: 42).
        external_clusters: Path to CSV with pre-computed cluster labels.
            Required when method='external'.
        outputpath: Path to save cluster labels CSV. If None, uses
            default from directory_config.
        output_h5ad: Path to save AnnData as h5ad. If None, uses
            default from directory_config.

    Returns:
        AnnData object with cluster labels in .obs['leiden'],
        UMAP embedding in .obsm['X_umap'], and dimensionality
        reduction in .obsm['X_pca'] or .obsm['X_lsi'].
    """
    # Build strategy kwargs
    strategy_kwargs = {
        'resolution': resolution,
        'random_seed': random_seed,
    }

    if method == 'leiden_tfidf':
        strategy_kwargs['n_dims'] = n_pcs
    elif method == 'leiden_libsize':
        strategy_kwargs['n_pcs'] = n_pcs
    elif method == 'external':
        strategy_kwargs['labels_path'] = external_clusters

    strategy = get_strategy(method, **strategy_kwargs)

    print(f"Clustering with strategy: {method}")
    print(f"Parameters: {strategy.get_params()}")

    # Run the three-step pipeline
    adata = strategy.normalize(adata)
    adata = strategy.reduce_dims(adata)
    adata = strategy.cluster(adata)

    # Extract and save cluster labels
    cluster_labels = adata.obs['leiden'].copy()
    cluster_labels = cluster_labels.sort_index()

    if outputpath is not None:
        os.makedirs(os.path.dirname(outputpath), exist_ok=True)
        cluster_labels.to_csv(outputpath, header=True, index=True)
        print(f"Cluster labels saved to: {outputpath}")

    # Save full AnnData object (includes UMAP, dim reduction, etc.)
    if output_h5ad is not None:
        os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
        adata.write(output_h5ad)
        print(f"AnnData saved to: {output_h5ad}")

    # Report clustering summary
    n_clusters = cluster_labels.nunique()
    print(f"Found {n_clusters} clusters across {adata.n_obs} cells")
    print(f"Cluster sizes:\n{cluster_labels.value_counts().sort_index()}")

    return adata
