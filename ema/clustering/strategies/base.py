"""Abstract base class for clustering strategies."""

from abc import ABC, abstractmethod
import anndata as ad
from typing import Dict, Any


class ClusteringStrategy(ABC):
    """Base class for clustering strategies.

    Each strategy handles the full pipeline from raw counts to cluster labels:
        normalize -> reduce_dims -> cluster

    Subclasses must implement all three steps plus get_params for logging.
    """

    @abstractmethod
    def normalize(self, adata: ad.AnnData) -> ad.AnnData:
        """Normalize the count matrix.

        Args:
            adata: AnnData object with raw counts.

        Returns:
            AnnData with normalized values in .X
        """
        pass

    @abstractmethod
    def reduce_dims(self, adata: ad.AnnData) -> ad.AnnData:
        """Perform dimensionality reduction.

        Args:
            adata: AnnData object with normalized data.

        Returns:
            AnnData with reduced representation stored
            (e.g., .obsm['X_pca'] or .obsm['X_lsi']).
        """
        pass

    @abstractmethod
    def cluster(self, adata: ad.AnnData) -> ad.AnnData:
        """Compute neighbors, UMAP, and cluster labels.

        Args:
            adata: AnnData object with dimensionality reduction computed.

        Returns:
            AnnData with .obs['leiden'] cluster labels and .obsm['X_umap'].
        """
        pass

    @abstractmethod
    def get_params(self) -> Dict[str, Any]:
        """Return current strategy parameters for logging/reproducibility."""
        pass
