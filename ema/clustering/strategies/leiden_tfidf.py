"""TF-IDF + LSI (TruncatedSVD) + Leiden clustering.

scATAC-seq style pipeline adapted for PAS peak data.
Uses Signac Method 1: log1p(TF * IDF * scale_factor).

PAS counts have higher dynamic range (0-100+) than scATAC (0-3),
so we use sublinear_tf=True to compress high counts.

Pipeline: TF-IDF -> TruncatedSVD (LSI) -> remove depth-correlated
components -> cosine neighbors -> Leiden.
"""

import numpy as np
import scanpy as sc
import anndata as ad
from scipy.sparse import issparse, csr_matrix
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize
from typing import Dict, Any

from ema.clustering.strategies.base import ClusteringStrategy
from ema.clustering.strategies import register


def _tfidf_signac_method1(X, scale_factor=1e4):
    """Signac Method 1 TF-IDF: log1p(TF * IDF * scale_factor).

    TF  = count / total_counts_per_cell
    IDF = n_cells / cells_with_peak

    Computed on the STORED entries only.  The previous implementation
    densified the cells x PAS matrix (``X.toarray()``) and held four
    float64 copies of it at once (``astype``, ``tf``, ``tf * idf``,
    ``log1p``), i.e. ``4 * n_cells * n_PAS * 8`` bytes -- 291 GB for the
    23,303 x 390,493 PBMC 10k matrix, which was >98 % of the whole run's
    peak RSS on every dataset measured.  The sparse form needs ~3 x nnz x 8
    bytes (~3 GB on that matrix) and is bit-identical to the dense result:

    * every zero of ``X`` maps to ``log1p(0 * idf * scale) == 0`` and was
      dropped by the trailing ``csr_matrix(tfidf)`` anyway, so the sparsity
      pattern is the same after :meth:`eliminate_zeros`;
    * the per-entry arithmetic is the same scalar sequence in the same
      order (``/ cell_total``, ``* idf``, ``* scale_factor``, ``log1p``);
    * the row sums and the per-PAS cell counts are sums of integers
      (counts), exact in float64 whatever the summation order.

    Verified identical (``X``, ``X_lsi`` and Leiden labels) against the
    dense implementation on two real runs (PBMC chr19+21 slice, mouse
    testis) and in ``tests/test_tfidf_sparse_identity.py``.

    Args:
        X: Count matrix (cells x peaks), sparse or dense.  Integer-valued
            counts are expected (bit-identity with the dense formula relies
            on the row sums being exact; non-integer input still gives the
            correct TF-IDF).
        scale_factor: Scaling factor (default: 10,000).

    Returns:
        Sparse CSR matrix (float64) with TF-IDF values.
    """
    X = csr_matrix(X, dtype=np.float64, copy=True)
    X.sum_duplicates()     # canonical: sorted, unique column indices per row
    X.eliminate_zeros()    # explicit zeros would be dropped by the dense path
    n_cells, n_peaks = X.shape

    # Term frequency: normalize each cell by total counts
    cell_totals = np.asarray(X.sum(axis=1), dtype=np.float64).ravel()
    cell_totals[cell_totals == 0] = 1  # avoid division by zero

    # Inverse document frequency: n_cells / cells_with_peak
    cells_with_peak = np.bincount(
        X.indices[X.data > 0], minlength=n_peaks
    ).astype(np.int64)
    cells_with_peak[cells_with_peak == 0] = 1  # avoid division by zero
    idf = n_cells / cells_with_peak

    # Signac Method 1: log1p(TF * IDF * scale_factor), same op order as the
    # dense expression ``np.log1p((X / cell_totals) * idf * scale_factor)``.
    rows = np.repeat(np.arange(n_cells), np.diff(X.indptr))
    data = X.data / cell_totals[rows]
    data *= idf[X.indices]
    data *= scale_factor
    np.log1p(data, out=data)

    tfidf = csr_matrix((data, X.indices, X.indptr), shape=(n_cells, n_peaks))
    tfidf.eliminate_zeros()
    return tfidf


def _remove_depth_correlated_components(X_lsi, total_counts, threshold=0.75):
    """Remove SVD components highly correlated with sequencing depth.

    ArchR-style filtering: remove components where |correlation| with
    total counts per cell exceeds the threshold.

    Args:
        X_lsi: LSI embedding (cells x components).
        total_counts: Total counts per cell.
        threshold: Correlation threshold for removal (default: 0.75).

    Returns:
        Filtered LSI embedding with correlated components removed,
        and list of kept component indices.
    """
    n_components = X_lsi.shape[1]
    keep_mask = np.ones(n_components, dtype=bool)

    for i in range(n_components):
        corr = np.abs(np.corrcoef(X_lsi[:, i], total_counts)[0, 1])
        if corr > threshold:
            keep_mask[i] = False

    return X_lsi[:, keep_mask], np.where(keep_mask)[0]


@register("leiden_tfidf")
class LeidenTfidfStrategy(ClusteringStrategy):
    """TF-IDF + LSI with cosine-distance Leiden clustering.

    This is the recommended strategy for PAS peak data, following the
    scATAC-seq paradigm (Signac/ArchR).

    Parameters:
        resolution: Leiden resolution parameter (default: 1.0).
        n_components: Number of SVD components to compute (default: 50).
        n_dims: Number of LSI dimensions to use for neighbors (default: 40).
        n_neighbors: Number of neighbors for kNN graph (default: 30).
        scale_factor: TF-IDF scale factor (default: 10,000).
        depth_corr_threshold: Threshold for removing depth-correlated
            components (default: 0.75).
        random_seed: Random seed for reproducibility (default: 42).

    Tunable hyperparameters:
        resolution (default 1.0): Leiden resolution. Higher values produce more,
            smaller clusters. CLI: ``--resolution`` / YAML: ``resolution``.
        n_components (default 50): SVD components computed before depth-
            correlation filtering. CLI: ``--n-svd-components`` /
            YAML: ``n_svd_components``.
        n_dims (default 40): LSI dimensions used for the kNN graph (applied
            after depth-correlation filtering). CLI: ``--n-pcs`` /
            YAML: ``n_pcs``.
        n_neighbors (default 30): kNN graph size. Larger values smooth cluster
            boundaries. CLI: ``--n-neighbors`` / YAML: ``n_neighbors``.
        scale_factor (default 10000): TF-IDF scale factor. Adjust if your
            counts have very different dynamic range. CLI: ``--tfidf-scale-factor``
            / YAML: ``tfidf_scale_factor``.
        depth_corr_threshold (default 0.75): Pearson |r| threshold for
            removing LSI components correlated with sequencing depth (ArchR-style).
            Set to 1.0 to disable. CLI: ``--depth-corr-threshold`` /
            YAML: ``depth_corr_threshold``.
        random_seed (default 42): RNG seed for SVD and Leiden.
            CLI: ``--random-seed`` / YAML: ``random_seed``.
    """

    def __init__(self, resolution=1.0, n_components=50, n_dims=40,
                 n_neighbors=30, scale_factor=1e4,
                 depth_corr_threshold=0.75, random_seed=42):
        self.resolution = resolution
        self.n_components = n_components
        self.n_dims = n_dims
        self.n_neighbors = n_neighbors
        self.scale_factor = scale_factor
        self.depth_corr_threshold = depth_corr_threshold
        self.random_seed = random_seed

    def normalize(self, adata: ad.AnnData) -> ad.AnnData:
        """Apply Signac Method 1 TF-IDF normalization."""
        # Store total counts before normalization (for depth correlation check)
        if issparse(adata.X):
            adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
        else:
            adata.obs['total_counts'] = adata.X.sum(axis=1)

        # Apply TF-IDF
        adata.X = _tfidf_signac_method1(adata.X, scale_factor=self.scale_factor)
        return adata

    def reduce_dims(self, adata: ad.AnnData) -> ad.AnnData:
        """LSI (TruncatedSVD) dimensionality reduction.

        Computes SVD, removes depth-correlated components (ArchR-style),
        and stores the result in .obsm['X_lsi'].
        """
        n_components = min(
            self.n_components,
            adata.n_vars - 1,
            adata.n_obs - 1,
        )
        if n_components < 2:
            raise ValueError(
                f"Dataset too small for LSI: n_obs={adata.n_obs}, n_vars={adata.n_vars}. "
                f"Need at least 3 cells and 3 features."
            )

        svd = TruncatedSVD(
            n_components=n_components,
            algorithm='randomized',
            random_state=self.random_seed,
        )
        X_lsi = svd.fit_transform(adata.X)

        # L2 normalize the embeddings (standard for cosine-based downstream)
        X_lsi = normalize(X_lsi, norm='l2')

        # Remove depth-correlated components (ArchR-style)
        total_counts = adata.obs['total_counts'].values
        X_lsi_filtered, kept_indices = _remove_depth_correlated_components(
            X_lsi, total_counts, threshold=self.depth_corr_threshold,
        )

        adata.obsm['X_lsi'] = X_lsi_filtered
        adata.uns['lsi_kept_components'] = kept_indices
        adata.uns['lsi_variance_ratio'] = svd.explained_variance_ratio_
        return adata

    def cluster(self, adata: ad.AnnData) -> ad.AnnData:
        """Cosine-distance kNN graph, UMAP, and Leiden clustering."""
        n_dims = min(self.n_dims, adata.obsm['X_lsi'].shape[1])

        sc.pp.neighbors(
            adata,
            use_rep='X_lsi',
            n_neighbors=self.n_neighbors,
            n_pcs=n_dims,
            metric='cosine',
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
            "strategy": "leiden_tfidf",
            "normalization": "TF-IDF (Signac Method 1)",
            "dim_reduction": "LSI (TruncatedSVD)",
            "clustering": "Leiden",
            "resolution": self.resolution,
            "n_components": self.n_components,
            "n_dims": self.n_dims,
            "n_neighbors": self.n_neighbors,
            "scale_factor": self.scale_factor,
            "depth_corr_threshold": self.depth_corr_threshold,
            "random_seed": self.random_seed,
        }
