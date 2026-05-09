"""Mutual Nearest Neighbor (MNN) cross-dataset cluster matching strategy.

Algorithm
---------
1. Load all h5ad files and concatenate their cell × feature matrices along
   the obs axis into a single combined AnnData.  A ``batch`` column records
   which dataset each cell comes from.

2. Re-embed all cells jointly using TruncatedSVD (LSI) on the concatenated
   matrix.  This places cells from all datasets in a shared latent space
   without explicit batch correction — useful as a fast approximation when
   batch effects are moderate.

3. For each pair of datasets (A, B), find Mutual Nearest Neighbors in the
   shared embedding.  A cell from A is an MNN partner of a cell from B if:
       - cell_a is among the k nearest neighbours of cell_b in dataset A's
         embedding, AND
       - cell_b is among the k nearest neighbours of cell_a in dataset B's
         embedding.

4. For each cluster_i in dataset A: count how many of its cells have MNN
   links to cells in each cluster_j of dataset B.  The dominant target
   cluster is the one with the most MNN votes.

5. Build a similarity score per (A_i, B_j) pair:
       sim(A_i, B_j) = (MNN votes from A_i to B_j) / |cells in A_i|
   This is the fraction of cluster A_i cells that have their dominant MNN
   partner in cluster B_j.

6. Hungarian assignment + transitive closure, identical to marker_overlap.

Limitations / Tradeoffs
-----------------------
- Requires re-embedding ALL datasets together, so memory scales with
  total cell count.  For >100K cells, consider subsampling first.
- The shared LSI is computed without batch correction (no Harmony/Scanorama).
  Strong batch effects may dominate the embedding, causing false negatives.
- k-NN search cost is O(n_cells² / n_datasets) per pair.  Slower than
  marker_overlap for large datasets (>10K cells per dataset).
- Requires all datasets to share the same feature space (same PAS BED).
- Use ``marker_overlap`` as the default; switch to ``mnn`` when marker gene
  sets are too dissimilar due to differing detection rates.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.optimize import linear_sum_assignment
from scipy.sparse import issparse, vstack
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize

from ema.clustering.cross_dataset.base import ClusterMatchStrategy
from ema.clustering.cross_dataset import register_match_strategy
from ema.clustering.cross_dataset.marker_overlap import (
    _components_to_dataframe,
    _hungarian_matches,
    _identity_mapping,
    _transitive_closure,
)


def _concat_adatas(
    h5ad_paths: list[Path],
    dataset_ids: list[str],
) -> ad.AnnData:
    """Load and concatenate h5ad files along obs axis.

    Intersects features (var_names) across datasets so the combined matrix
    has a consistent feature space.  Stores dataset origin in ``obs['batch']``
    and the original leiden label in ``obs['leiden_original']``.

    Args:
        h5ad_paths: Paths to h5ad files.
        dataset_ids: Corresponding dataset identifiers.

    Returns:
        Combined AnnData with ``obs['batch']`` and ``obs['leiden_original']``.

    Raises:
        ValueError: If any dataset is missing the ``leiden`` obs column.
        ValueError: If datasets share no common features.
    """
    adatas: list[ad.AnnData] = []
    for path, ds_id in zip(h5ad_paths, dataset_ids):
        adata = ad.read_h5ad(path)
        if "leiden" not in adata.obs.columns:
            raise ValueError(
                f"Dataset '{ds_id}' ({path}) has no 'leiden' column in obs."
            )
        adata.obs["batch"] = ds_id
        adata.obs["leiden_original"] = adata.obs["leiden"].astype(str)
        adatas.append(adata)

    # Find common feature names.
    common_vars = set(adatas[0].var_names)
    for adata in adatas[1:]:
        common_vars &= set(adata.var_names)

    if not common_vars:
        raise ValueError(
            "Datasets share no common features — cannot run MNN matching. "
            "Ensure all datasets were processed with the same unified PAS BED."
        )

    # Subset and concatenate.
    common_list = sorted(common_vars)
    adatas_sub = [adata[:, common_list].copy() for adata in adatas]
    combined = ad.concat(adatas_sub, label="batch", keys=dataset_ids, merge="same")
    return combined


def _build_shared_embedding(
    combined: ad.AnnData,
    n_components: int = 30,
    random_state: int = 42,
) -> np.ndarray:
    """Compute a shared LSI embedding on the concatenated matrix.

    Applies log-normalization first if the data looks raw (integer-valued).
    Uses TruncatedSVD (LSI) and L2-normalises the resulting embedding.

    Args:
        combined: Concatenated AnnData from all datasets.
        n_components: Number of SVD components (default: 30).
        random_state: Random seed for TruncatedSVD.

    Returns:
        Embedding array of shape (n_cells_total, n_components).
    """
    X = combined.X
    if issparse(X):
        sample = X[:min(10, X.shape[0])].toarray()
    else:
        sample = np.asarray(X[:min(10, X.shape[0])])

    looks_raw = sample.max() > 20 or np.all(sample == sample.astype(int))
    if looks_raw:
        sc.pp.normalize_total(combined, target_sum=1e4)
        sc.pp.log1p(combined)

    X_mat = combined.X
    if issparse(X_mat):
        X_mat = X_mat.toarray()

    n_comp = min(n_components, X_mat.shape[0] - 1, X_mat.shape[1] - 1)
    n_comp = max(n_comp, 2)

    svd = TruncatedSVD(
        n_components=n_comp,
        algorithm="randomized",
        random_state=random_state,
    )
    embedding = svd.fit_transform(X_mat)
    embedding = normalize(embedding, norm="l2")
    return embedding


def _find_mnn_pairs(
    emb_a: np.ndarray,
    emb_b: np.ndarray,
    k: int = 10,
) -> list[tuple[int, int]]:
    """Find mutual nearest neighbour pairs between two embedding subsets.

    A pair (i, j) is mutual if:
    - cell i (from A) is among the k-NN of cell j within A's embedding
      subspace (searched in B's perspective), AND
    - cell j (from B) is among the k-NN of cell i within B's embedding.

    For simplicity, we search across the combined embedding:
    - For each cell in B, find k nearest in A.
    - For each cell in A, find k nearest in B.
    - Intersect.

    Args:
        emb_a: Embedding for dataset A cells, shape (n_a, n_dims).
        emb_b: Embedding for dataset B cells, shape (n_b, n_dims).
        k: Number of neighbours to search.

    Returns:
        List of (index_in_A, index_in_B) mutual nearest neighbour pairs.
    """
    k_eff = min(k, emb_a.shape[0], emb_b.shape[0])

    # A → B: for each cell in A, find k nearest in B.
    nn_ab = NearestNeighbors(n_neighbors=k_eff, algorithm="auto", metric="cosine")
    nn_ab.fit(emb_b)
    # shape (n_a, k_eff)
    neighbours_ab = nn_ab.kneighbors(emb_a, return_distance=False)

    # B → A: for each cell in B, find k nearest in A.
    nn_ba = NearestNeighbors(n_neighbors=k_eff, algorithm="auto", metric="cosine")
    nn_ba.fit(emb_a)
    # shape (n_b, k_eff)
    neighbours_ba = nn_ba.kneighbors(emb_b, return_distance=False)

    # Build set of A→B links.
    ab_set: set[tuple[int, int]] = set()
    for i, neighs in enumerate(neighbours_ab):
        for j in neighs:
            ab_set.add((i, int(j)))

    # Intersect with B→A links (mutuality check).
    mnn_pairs: list[tuple[int, int]] = []
    for j, neighs in enumerate(neighbours_ba):
        for i in neighs:
            if (int(i), j) in ab_set:
                mnn_pairs.append((int(i), j))

    return mnn_pairs


def _vote_similarity_matrix(
    mnn_pairs: list[tuple[int, int]],
    labels_a: np.ndarray,
    labels_b: np.ndarray,
    clusters_a: list[str],
    clusters_b: list[str],
) -> np.ndarray:
    """Convert MNN pairs to a cluster-level similarity matrix.

    Similarity sim(A_i, B_j) = fraction of cluster A_i cells whose MNN
    partner is in cluster B_j.

    Args:
        mnn_pairs: List of (idx_in_A, idx_in_B) mutual pairs.
        labels_a: Cluster label per cell in A (array of strings).
        labels_b: Cluster label per cell in B (array of strings).
        clusters_a: Sorted unique cluster labels for A.
        clusters_b: Sorted unique cluster labels for B.

    Returns:
        Similarity matrix of shape (len(clusters_a), len(clusters_b)).
    """
    cl_to_i = {cl: i for i, cl in enumerate(clusters_a)}
    cl_to_j = {cl: j for j, cl in enumerate(clusters_b)}

    # vote_matrix[i, j] = number of MNN pairs from cluster_a_i to cluster_b_j
    vote_matrix = np.zeros((len(clusters_a), len(clusters_b)), dtype=float)
    for idx_a, idx_b in mnn_pairs:
        i = cl_to_i.get(labels_a[idx_a])
        j = cl_to_j.get(labels_b[idx_b])
        if i is not None and j is not None:
            vote_matrix[i, j] += 1.0

    # Normalise by cluster size in A.
    cluster_sizes_a = np.array(
        [(labels_a == cl).sum() for cl in clusters_a], dtype=float
    )
    cluster_sizes_a[cluster_sizes_a == 0] = 1.0
    sim = vote_matrix / cluster_sizes_a[:, np.newaxis]
    return sim


class MNNStrategy(ClusterMatchStrategy):
    """Cluster matching via mutual nearest neighbours in shared LSI space.

    Re-embeds cells from all datasets jointly using TruncatedSVD, then
    identifies mutual nearest neighbour links across dataset pairs.  Each
    cluster's MNN vote distribution determines which cluster in the other
    dataset it best corresponds to.

    This approach is more sensitive to batch effects than ``marker_overlap``
    but can detect correspondences when marker gene overlap is low (e.g.
    different capture technologies or very different library sizes).

    Parameters:
        n_components: Number of LSI components for the shared embedding.
            More components capture finer structure but increase cost.
        k_neighbors: Number of neighbours for MNN search.
        random_state: Seed for TruncatedSVD reproducibility.
    """

    name = "mnn"

    def __init__(
        self,
        n_components: int = 30,
        k_neighbors: int = 10,
        random_state: int = 42,
    ) -> None:
        self.n_components = n_components
        self.k_neighbors = k_neighbors
        self.random_state = random_state

    def match(
        self,
        h5ad_paths: list[Path],
        dataset_ids: list[str],
        n_top_markers: int = 50,
        n_jobs: int = -1,
    ) -> pd.DataFrame:
        """Match clusters via MNN in shared LSI embedding.

        Note: ``n_top_markers`` is unused by this strategy (MNN does not
        compute marker genes).  It is accepted for API compatibility.

        Args:
            h5ad_paths: Paths to h5ad files (one per dataset).
            dataset_ids: Human-readable dataset identifiers.
            n_top_markers: Unused; accepted for interface compatibility.
            n_jobs: Unused; MNN search is not parallelised at the strategy
                level (sklearn NearestNeighbors uses its own threading).

        Returns:
            DataFrame with columns:
            dataset_id, original_cluster, canonical_cluster,
            match_confidence, matched_to.
        """
        h5ad_paths = [Path(p) for p in h5ad_paths]

        if len(h5ad_paths) == 1:
            return _identity_mapping(h5ad_paths[0], dataset_ids[0])

        # Load and concatenate.
        combined = _concat_adatas(h5ad_paths, dataset_ids)

        # Shared LSI embedding.
        embedding = _build_shared_embedding(
            combined,
            n_components=self.n_components,
            random_state=self.random_state,
        )

        # Attach embedding back to combined for easy slicing.
        batch_col = combined.obs["batch"].values
        leiden_col = combined.obs["leiden_original"].values

        # Build all (dataset_id, original_cluster) nodes.
        all_nodes: list[tuple[str, str]] = []
        ds_cluster_info: dict[str, tuple[np.ndarray, np.ndarray, list[str]]] = {}
        for ds_id in dataset_ids:
            mask = batch_col == ds_id
            emb_ds = embedding[mask]
            labels_ds = leiden_col[mask]
            clusters_ds = sorted(set(labels_ds))
            ds_cluster_info[ds_id] = (emb_ds, labels_ds, clusters_ds)
            for cl in clusters_ds:
                all_nodes.append((ds_id, cl))

        # Match every dataset pair.
        all_edges: list[tuple[tuple[str, str], tuple[str, str], float]] = []
        for ds_a, ds_b in combinations(dataset_ids, 2):
            emb_a, labels_a, clusters_a = ds_cluster_info[ds_a]
            emb_b, labels_b, clusters_b = ds_cluster_info[ds_b]

            mnn_pairs = _find_mnn_pairs(emb_a, emb_b, k=self.k_neighbors)

            if not mnn_pairs:
                continue

            sim = _vote_similarity_matrix(
                mnn_pairs, labels_a, labels_b, clusters_a, clusters_b
            )
            raw_matches = _hungarian_matches(sim, clusters_a, clusters_b)
            for cl_a, cl_b, score in raw_matches:
                all_edges.append(((ds_a, cl_a), (ds_b, cl_b), score))

        components = _transitive_closure(all_nodes, all_edges)
        return _components_to_dataframe(components)


register_match_strategy(MNNStrategy)
