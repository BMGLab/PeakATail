"""Library-size normalization + PCA + Leiden clustering.

This is the traditional gene-expression-style pipeline adapted for PAS peak data.
Pipeline: normalize_total -> log1p -> HVG -> scale -> PCA -> neighbors -> Leiden.
"""

import scanpy as sc
import anndata as ad
from typing import Dict, Any

from ema.clustering.strategies.base import ClusteringStrategy
from ema.clustering.strategies import register


@register("leiden_libsize")
class LeidenLibsizeStrategy(ClusteringStrategy):
    """Library-size normalization with PCA and Leiden clustering.

    Parameters:
        resolution: Leiden resolution parameter (default: 1.0).
        n_pcs: Number of principal components for neighbors (default: 40).
        n_comps: Number of PCA components to compute (default: 50).
        n_top_genes: Number of highly variable genes to select (default: 2000).
        n_neighbors: Number of neighbors for kNN graph (default: 10).
        random_seed: Random seed for reproducibility (default: 42).

    Tunable hyperparameters:
        resolution (default 1.0): Leiden resolution. CLI: ``--resolution`` /
            YAML: ``resolution``.
        n_comps (default 50): PCA components computed before neighbors.
            CLI: ``--n-svd-components`` / YAML: ``n_svd_components``.
        n_pcs (default 40): PCA components used for the kNN graph (must be
            ≤ n_comps). CLI: ``--n-pcs`` / YAML: ``n_pcs``.
        n_neighbors (default 10): kNN graph size. CLI: ``--n-neighbors`` /
            YAML: ``n_neighbors``.
        n_top_genes (default 2000): Highly variable PAS selected before PCA.
            CLI: ``--n-top-hvg`` / YAML: ``n_top_hvg``.
        random_seed (default 42): RNG seed. CLI: ``--random-seed`` /
            YAML: ``random_seed``.
    """

    def __init__(self, resolution=1.0, n_pcs=40, n_comps=50,
                 n_top_genes=2000, n_neighbors=10, random_seed=42):
        self.resolution = resolution
        self.n_pcs = n_pcs
        self.n_comps = n_comps
        self.n_top_genes = n_top_genes
        self.n_neighbors = n_neighbors
        self.random_seed = random_seed

    def normalize(self, adata: ad.AnnData) -> ad.AnnData:
        """Library-size normalization, log1p, HVG selection."""
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Adjusted HVG parameters for peak data (sparser, lower dynamic range)
        n_top = min(self.n_top_genes, adata.n_vars)
        sc.pp.highly_variable_genes(
            adata,
            min_mean=0.005,
            max_mean=5,
            min_disp=0.3,
            n_top_genes=n_top,
        )
        adata = adata[:, adata.var.highly_variable].copy()

        sc.pp.scale(adata, max_value=10)
        return adata

    def reduce_dims(self, adata: ad.AnnData) -> ad.AnnData:
        """PCA dimensionality reduction."""
        n_comps = min(self.n_comps, adata.n_vars - 1, adata.n_obs - 1)
        sc.tl.pca(adata, n_comps=n_comps, random_state=self.random_seed)
        return adata

    def cluster(self, adata: ad.AnnData) -> ad.AnnData:
        """kNN graph, UMAP, and Leiden clustering."""
        n_pcs = min(self.n_pcs, adata.obsm['X_pca'].shape[1])
        sc.pp.neighbors(
            adata,
            n_neighbors=self.n_neighbors,
            n_pcs=n_pcs,
            random_state=self.random_seed,
        )
        sc.tl.umap(adata, random_state=self.random_seed)
        sc.tl.leiden(
            adata,
            resolution=self.resolution,
            random_state=self.random_seed,
            flavor='igraph',
            n_iterations=2,
            directed=False,
        )
        return adata

    def get_params(self) -> Dict[str, Any]:
        return {
            "strategy": "leiden_libsize",
            "normalization": "library_size",
            "dim_reduction": "PCA",
            "clustering": "Leiden",
            "resolution": self.resolution,
            "n_pcs": self.n_pcs,
            "n_comps": self.n_comps,
            "n_top_genes": self.n_top_genes,
            "n_neighbors": self.n_neighbors,
            "random_seed": self.random_seed,
        }
