"""Import pre-computed cluster labels from an external CSV file.

This strategy skips all normalization, dimensionality reduction, and
clustering. It simply loads cluster labels from a user-provided CSV file
and assigns them to the AnnData object.

This enables using gene-expression-based clusters (from Seurat, scanpy, etc.)
for the APA switch test, matching the workflow of other APA tools.
"""

import pandas as pd
import scanpy as sc
import anndata as ad
from typing import Dict, Any

from ema.clustering.strategies.base import ClusteringStrategy
from ema.clustering.strategies import register


@register("external")
class ExternalStrategy(ClusteringStrategy):
    """Load pre-computed cluster labels from CSV.

    The CSV file must have two columns (no header):
        column 0: cell barcode
        column 1: cluster label

    Parameters:
        labels_path: Path to CSV file with cluster labels.
        random_seed: Random seed for UMAP computation (default: 42).
    """

    def __init__(self, labels_path=None, random_seed=42):
        if labels_path is None:
            raise ValueError(
                "ExternalStrategy requires labels_path. "
                "Provide --external-clusters on the CLI."
            )
        self.labels_path = labels_path
        self.random_seed = random_seed

    def normalize(self, adata: ad.AnnData) -> ad.AnnData:
        """No normalization for external labels. Apply log-norm for UMAP only."""
        # Light normalization just so UMAP is meaningful
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return adata

    def reduce_dims(self, adata: ad.AnnData) -> ad.AnnData:
        """PCA for UMAP visualization only (labels are pre-computed)."""
        n_comps = min(50, adata.n_vars - 1, adata.n_obs - 1)
        sc.tl.pca(adata, n_comps=n_comps, random_state=self.random_seed)
        return adata

    def cluster(self, adata: ad.AnnData) -> ad.AnnData:
        """Load external labels and compute UMAP for visualization."""
        # Load external labels
        labels_df = pd.read_csv(self.labels_path, header=None)
        labels_df.columns = ['barcode', 'cluster']
        labels_df['cluster'] = labels_df['cluster'].astype(str)
        labels_df = labels_df.set_index('barcode')

        # Match barcodes between adata and labels
        common_barcodes = adata.obs_names.intersection(labels_df.index)
        if len(common_barcodes) == 0:
            raise ValueError(
                f"No matching barcodes between AnnData ({adata.n_obs} cells) "
                f"and external labels ({len(labels_df)} entries). "
                "Check barcode format."
            )

        n_missing = adata.n_obs - len(common_barcodes)
        if n_missing > 0:
            print(
                f"Warning: {n_missing} cells in AnnData have no external label. "
                "They will be assigned cluster 'unassigned'."
            )

        # Assign labels
        adata.obs['leiden'] = 'unassigned'
        adata.obs.loc[common_barcodes, 'leiden'] = (
            labels_df.loc[common_barcodes, 'cluster'].values
        )
        adata.obs['leiden'] = pd.Categorical(adata.obs['leiden'])

        # Compute UMAP for visualization
        n_pcs = min(40, adata.obsm['X_pca'].shape[1])
        sc.pp.neighbors(
            adata,
            n_neighbors=10,
            n_pcs=n_pcs,
            random_state=self.random_seed,
        )
        sc.tl.umap(adata, random_state=self.random_seed)

        return adata

    def get_params(self) -> Dict[str, Any]:
        return {
            "strategy": "external",
            "labels_path": self.labels_path,
            "random_seed": self.random_seed,
        }
