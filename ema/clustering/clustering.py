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


# ---------------------------------------------------------------------------
# Private helpers shared by all _do_* functions
# ---------------------------------------------------------------------------

def _save_and_report(adata: ad.AnnData, outputpath=None, output_h5ad=None):
    """Save cluster labels + AnnData and print summary; returns adata."""
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


# ---------------------------------------------------------------------------
# Module-level strategy implementations (called by the registry shims)
# ---------------------------------------------------------------------------

def _do_leiden_tfidf(adata: ad.AnnData,
                     resolution=1.0,
                     random_seed=42,
                     n_pcs=40,
                     external_clusters=None,
                     outputpath=None,
                     output_h5ad=None,
                     n_neighbors=30,
                     tfidf_scale_factor=1e4,
                     depth_corr_threshold=0.75,
                     n_svd_components=50,
                     **_ignored):  # n_top_hvg etc. — libsize-only knobs
    """Run the TF-IDF + LSI + Leiden pipeline."""
    strategy = get_strategy("leiden_tfidf",
                            resolution=resolution,
                            random_seed=random_seed,
                            n_dims=n_pcs,
                            n_neighbors=n_neighbors,
                            scale_factor=tfidf_scale_factor,
                            depth_corr_threshold=depth_corr_threshold,
                            n_components=n_svd_components)

    print(f"Clustering with strategy: leiden_tfidf")
    print(f"Parameters: {strategy.get_params()}")

    adata = strategy.normalize(adata)
    adata = strategy.reduce_dims(adata)
    adata = strategy.cluster(adata)

    return _save_and_report(adata, outputpath=outputpath, output_h5ad=output_h5ad)


def _do_leiden_libsize(adata: ad.AnnData,
                       resolution=1.0,
                       random_seed=42,
                       n_pcs=40,
                       external_clusters=None,
                       outputpath=None,
                       output_h5ad=None,
                       n_neighbors=10,
                       n_svd_components=50,
                       n_top_hvg=2000,
                       **_ignored):  # tfidf_scale_factor, depth_corr_threshold — tfidf-only
    """Run the library-size normalization + PCA + Leiden pipeline."""
    strategy = get_strategy("leiden_libsize",
                            resolution=resolution,
                            random_seed=random_seed,
                            n_pcs=n_pcs,
                            n_neighbors=n_neighbors,
                            n_comps=n_svd_components,
                            n_top_genes=n_top_hvg)

    print(f"Clustering with strategy: leiden_libsize")
    print(f"Parameters: {strategy.get_params()}")

    adata = strategy.normalize(adata)
    adata = strategy.reduce_dims(adata)
    adata = strategy.cluster(adata)

    return _save_and_report(adata, outputpath=outputpath, output_h5ad=output_h5ad)


def _do_external(adata: ad.AnnData,
                 resolution=1.0,
                 random_seed=42,
                 n_pcs=40,
                 external_clusters=None,
                 outputpath=None,
                 output_h5ad=None,
                 **_ignored):
    """Load pre-computed external cluster labels."""
    strategy = get_strategy("external",
                            resolution=resolution,
                            random_seed=random_seed,
                            labels_path=external_clusters)

    print(f"Clustering with strategy: external")
    print(f"Parameters: {strategy.get_params()}")

    adata = strategy.normalize(adata)
    adata = strategy.reduce_dims(adata)
    adata = strategy.cluster(adata)

    return _save_and_report(adata, outputpath=outputpath, output_h5ad=output_h5ad)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def clustering(adata: ad.AnnData,
               method='leiden_tfidf',
               resolution=1.0,
               n_pcs=40,
               random_seed=42,
               external_clusters=None,
               outputpath=None,
               output_h5ad=None,
               n_neighbors: int | None = None,
               tfidf_scale_factor: float = 1e4,
               depth_corr_threshold: float = 0.75,
               n_svd_components: int = 50,
               n_top_hvg: int = 2000):
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
        n_neighbors: kNN graph size. Defaults to 30 for leiden_tfidf and
            10 for leiden_libsize when None.
        tfidf_scale_factor: TF-IDF scale factor for leiden_tfidf (default: 10000).
        depth_corr_threshold: Correlation threshold for depth-component removal
            in leiden_tfidf (default: 0.75).
        n_svd_components: SVD/PCA components computed before filtering
            (default: 50).
        n_top_hvg: Highly variable genes selected by leiden_libsize (default: 2000).

    Returns:
        AnnData object with cluster labels in .obs['leiden'],
        UMAP embedding in .obsm['X_umap'], and dimensionality
        reduction in .obsm['X_pca'] or .obsm['X_lsi'].
    """
    # Resolve method-specific n_neighbors defaults when caller passes None.
    if n_neighbors is None:
        n_neighbors = 10 if method == "leiden_libsize" else 30

    from ema.clustering.registry import get_clustering_strategy
    return get_clustering_strategy(method)(
        adata,
        resolution=resolution,
        n_pcs=n_pcs,
        random_seed=random_seed,
        external_clusters=external_clusters,
        outputpath=outputpath,
        output_h5ad=output_h5ad,
        n_neighbors=n_neighbors,
        tfidf_scale_factor=tfidf_scale_factor,
        depth_corr_threshold=depth_corr_threshold,
        n_svd_components=n_svd_components,
        n_top_hvg=n_top_hvg,
    )
