"""Marker-overlap cross-dataset cluster matching strategy.

Algorithm
---------
1. For each dataset, load the h5ad file and compute per-cluster marker PAS via
   Wilcoxon rank-sum (scanpy.tl.rank_genes_groups).  The top-N marker PAS by
   score constitute the cluster's "fingerprint".

2. For every pair of datasets (datasetA, datasetB) compute the Jaccard
   similarity between each cluster_i in A and each cluster_j in B:
       J(A_i, B_j) = |markers_A_i ∩ markers_B_j| / |markers_A_i ∪ markers_B_j|

3. Solve the maximum-weight 1-to-1 assignment between A and B using the
   Hungarian algorithm (scipy.optimize.linear_sum_assignment).  Only accept
   pairs whose Jaccard ≥ a small threshold (≥ 1 shared marker).

4. Build a graph where nodes are (dataset_id, original_cluster) tuples and
   edges connect matched pairs.  Extract connected components via transitive
   closure — if A.0 ~ B.2 and B.2 ~ C.5, all three share one canonical ID.

5. Assign canonical integer IDs (1, 2, …) to each connected component.
   Unmatched nodes (no edges) get unique IDs.  The match_confidence for a
   component is the mean Jaccard of all edges within it.

Parallelism
-----------
Marker computation is the expensive step.  When n_clusters > 10, joblib
parallelises across datasets — each dataset's rank_genes_groups call runs in
its own worker process, exploiting that each h5ad is independent.

Tradeoffs
---------
- Fast: marker computation is the only heavy step; similarity matrix is small
  (n_clusters × n_clusters per pair).
- Requires at least a handful of shared marker PAS to work well.  Datasets
  from very different protocols may share few markers even for the same cell
  type — use ``mnn`` in that case.
- Memory efficient: only the top-N marker names are retained, not the full
  expression matrix.
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed
from scipy.optimize import linear_sum_assignment
from scipy.sparse import issparse

from ema.clustering.cross_dataset.base import ClusterMatchStrategy
from ema.clustering.cross_dataset import register_match_strategy


def _load_and_compute_markers(
    h5ad_path: Path,
    dataset_id: str,
    n_top_markers: int,
) -> dict[str, set[str]]:
    """Load one h5ad and compute top-N marker PAS per cluster.

    Args:
        h5ad_path: Path to the h5ad file.
        dataset_id: Dataset identifier (used only for error messages).
        n_top_markers: Number of top marker PAS to retain per cluster.

    Returns:
        Mapping from original_cluster label → set of top marker PAS names.
    """
    adata = ad.read_h5ad(h5ad_path)

    if "leiden" not in adata.obs.columns:
        raise ValueError(
            f"Dataset '{dataset_id}' ({h5ad_path}) has no 'leiden' column in obs."
        )

    clusters = sorted(adata.obs["leiden"].unique().tolist())

    # rank_genes_groups requires normalised, log-transformed data.
    # Check if data looks raw (integers) and normalise accordingly.
    X = adata.X
    if issparse(X):
        sample = X[:10].toarray()
    else:
        sample = np.asarray(X[:10])

    # Heuristic: if max value >> 10 or values are integers, apply log-norm.
    looks_raw = sample.max() > 20 or np.all(sample == sample.astype(int))
    if looks_raw:
        adata = adata.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # n_genes capped at n_vars to avoid scanpy errors on tiny datasets.
    n_genes = min(n_top_markers, adata.n_vars)

    # rank_genes_groups needs at least 2 groups.
    markers: dict[str, set[str]] = {}
    if len(clusters) < 2:
        # Only one cluster — no meaningful markers; return all features.
        all_features = set(adata.var_names.tolist())
        for cl in clusters:
            markers[str(cl)] = all_features
        return markers

    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        n_genes=n_genes,
        key_added="_cdm_markers",
    )

    for cl in clusters:
        df = sc.get.rank_genes_groups_df(adata, group=str(cl), key="_cdm_markers")
        top_genes = df["names"].head(n_genes).tolist()
        markers[str(cl)] = set(top_genes)

    return markers


def _jaccard(a: set[str], b: set[str]) -> float:
    """Compute Jaccard similarity between two sets.

    Args:
        a: First set of marker names.
        b: Second set of marker names.

    Returns:
        Jaccard index in [0.0, 1.0].  Returns 0.0 if both sets are empty.
    """
    if not a and not b:
        return 0.0
    intersection = len(a & b)
    union = len(a | b)
    return intersection / union if union > 0 else 0.0


def _pairwise_similarity_matrix(
    markers_a: dict[str, set[str]],
    markers_b: dict[str, set[str]],
) -> tuple[np.ndarray, list[str], list[str]]:
    """Build Jaccard similarity matrix between two datasets' clusters.

    Args:
        markers_a: cluster_label → marker set for dataset A.
        markers_b: cluster_label → marker set for dataset B.

    Returns:
        Tuple of (similarity_matrix, cluster_labels_a, cluster_labels_b).
        similarity_matrix shape: (len(markers_a), len(markers_b)).
    """
    labels_a = sorted(markers_a.keys())
    labels_b = sorted(markers_b.keys())
    sim = np.zeros((len(labels_a), len(labels_b)), dtype=float)
    for i, la in enumerate(labels_a):
        for j, lb in enumerate(labels_b):
            sim[i, j] = _jaccard(markers_a[la], markers_b[lb])
    return sim, labels_a, labels_b


def _hungarian_matches(
    sim: np.ndarray,
    labels_a: list[str],
    labels_b: list[str],
    min_confidence: float = 0.0,
) -> list[tuple[str, str, float]]:
    """Run Hungarian assignment on a similarity matrix.

    Args:
        sim: Similarity matrix (n_a × n_b); higher = better.
        labels_a: Cluster labels for rows.
        labels_b: Cluster labels for columns.
        min_confidence: Minimum Jaccard to accept a match (pairs below this
            threshold are treated as unmatched).

    Returns:
        List of (label_a, label_b, jaccard_score) for accepted matches.
    """
    # linear_sum_assignment minimises; negate to maximise.
    row_ind, col_ind = linear_sum_assignment(-sim)
    matches: list[tuple[str, str, float]] = []
    for r, c in zip(row_ind, col_ind):
        score = float(sim[r, c])
        if score > min_confidence:
            matches.append((labels_a[r], labels_b[c], score))
    return matches


def _transitive_closure(
    nodes: list[tuple[str, str]],
    edges: list[tuple[tuple[str, str], tuple[str, str], float]],
) -> list[tuple[list[tuple[str, str]], float]]:
    """Find connected components via union-find.

    Args:
        nodes: All (dataset_id, original_cluster) pairs.
        edges: Matched pairs as ((ds_a, cl_a), (ds_b, cl_b), confidence).

    Returns:
        List of (component_members, mean_confidence) where component_members
        is a list of (dataset_id, original_cluster) and mean_confidence is
        the average Jaccard across all edges within the component.
    """
    parent: dict[tuple[str, str], tuple[str, str]] = {n: n for n in nodes}
    rank: dict[tuple[str, str], int] = {n: 0 for n in nodes}
    edge_weights: dict[frozenset[tuple[str, str]], float] = {}

    def find(x: tuple[str, str]) -> tuple[str, str]:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: tuple[str, str], y: tuple[str, str]) -> None:
        rx, ry = find(x), find(y)
        if rx == ry:
            return
        if rank[rx] < rank[ry]:
            rx, ry = ry, rx
        parent[ry] = rx
        if rank[rx] == rank[ry]:
            rank[rx] += 1

    for (a, b, score) in edges:
        union(a, b)
        key = frozenset({a, b})
        edge_weights[key] = score

    # Group nodes by root.
    component_map: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for node in nodes:
        root = find(node)
        component_map.setdefault(root, []).append(node)

    # Compute mean confidence per component (0.0 for singletons).
    components: list[tuple[list[tuple[str, str]], float]] = []
    for members in component_map.values():
        member_set = set(members)
        relevant_weights = [
            w
            for key, w in edge_weights.items()
            if key <= member_set  # all nodes in key are in this component
        ]
        conf = float(np.mean(relevant_weights)) if relevant_weights else 0.0
        components.append((members, conf))

    return components


class MarkerOverlapStrategy(ClusterMatchStrategy):
    """Cluster matching via top-N marker PAS Jaccard overlap.

    For each dataset, computes per-cluster marker PAS using Wilcoxon
    rank-sum differential expression (via scanpy).  Matches clusters
    across dataset pairs using the Hungarian algorithm on a Jaccard
    similarity matrix.  Builds transitive closure to assign canonical IDs
    that span more than two datasets.

    Suitable for datasets with overlapping PAS features.  If datasets share
    very few PAS (different genomic regions captured), use ``mnn`` instead.
    """

    name = "marker_overlap"

    def match(
        self,
        h5ad_paths: list[Path],
        dataset_ids: list[str],
        n_top_markers: int = 50,
        n_jobs: int = -1,
    ) -> pd.DataFrame:
        """Match clusters across datasets using marker PAS Jaccard overlap.

        Args:
            h5ad_paths: Paths to h5ad files (one per dataset).
            dataset_ids: Human-readable dataset identifiers.
            n_top_markers: Top-N marker PAS per cluster to use.
            n_jobs: Parallel workers for marker computation.

        Returns:
            DataFrame with columns:
            dataset_id, original_cluster, canonical_cluster,
            match_confidence, matched_to.
        """
        h5ad_paths = [Path(p) for p in h5ad_paths]

        # Edge case: single dataset → identity mapping.
        if len(h5ad_paths) == 1:
            return _identity_mapping(h5ad_paths[0], dataset_ids[0])

        # Determine parallelism: use joblib when there are multiple datasets.
        # Each dataset's marker computation is independent.
        # When n_jobs=-1, delegate to ResourceManager so the user's --threads
        # ceiling is respected.  Explicit n_jobs=N is capped at n_datasets.
        if n_jobs == -1:
            from ema.utils import get_resource_manager
            _effective_n_jobs = get_resource_manager().get_n_jobs(
                per_worker_mb=200, stage="marker_overlap"
            )
        else:
            _effective_n_jobs = min(n_jobs, len(h5ad_paths))
        markers_list: list[dict[str, set[str]]] = Parallel(
            n_jobs=_effective_n_jobs,
            backend="loky",
        )(
            delayed(_load_and_compute_markers)(path, ds_id, n_top_markers)
            for path, ds_id in zip(h5ad_paths, dataset_ids)
        )

        # Build all (dataset_id, original_cluster) nodes.
        all_nodes: list[tuple[str, str]] = []
        markers_by_ds: dict[str, dict[str, set[str]]] = {}
        for ds_id, markers in zip(dataset_ids, markers_list):
            markers_by_ds[ds_id] = markers
            for cl in markers:
                all_nodes.append((ds_id, cl))

        # Match every pair of datasets.
        all_edges: list[tuple[tuple[str, str], tuple[str, str], float]] = []
        for (ds_a, ds_b) in combinations(dataset_ids, 2):
            sim, labels_a, labels_b = _pairwise_similarity_matrix(
                markers_by_ds[ds_a], markers_by_ds[ds_b]
            )
            raw_matches = _hungarian_matches(sim, labels_a, labels_b)
            for cl_a, cl_b, score in raw_matches:
                all_edges.append(((ds_a, cl_a), (ds_b, cl_b), score))

        # Transitive closure → canonical components.
        components = _transitive_closure(all_nodes, all_edges)

        return _components_to_dataframe(components)


register_match_strategy(MarkerOverlapStrategy)


def _identity_mapping(h5ad_path: Path, dataset_id: str) -> pd.DataFrame:
    """Return an identity mapping for a single dataset.

    Each cluster maps to its own canonical ID with confidence 1.0 and no
    matched_to partners.

    Args:
        h5ad_path: Path to the h5ad file.
        dataset_id: Dataset identifier.

    Returns:
        DataFrame with columns:
        dataset_id, original_cluster, canonical_cluster,
        match_confidence, matched_to.
    """
    adata = ad.read_h5ad(h5ad_path)
    clusters = sorted(adata.obs["leiden"].unique().tolist())
    rows = []
    for canonical_id, cl in enumerate(clusters, start=1):
        rows.append(
            {
                "dataset_id": dataset_id,
                "original_cluster": str(cl),
                "canonical_cluster": canonical_id,
                "match_confidence": 1.0,
                "matched_to": "[]",
            }
        )
    return pd.DataFrame(rows)


def _components_to_dataframe(
    components: list[tuple[list[tuple[str, str]], float]],
) -> pd.DataFrame:
    """Convert connected components to the canonical mapping DataFrame.

    Args:
        components: List of (members, mean_confidence) from transitive
            closure.  Each member is a (dataset_id, original_cluster) tuple.

    Returns:
        DataFrame with columns:
        dataset_id, original_cluster, canonical_cluster,
        match_confidence, matched_to.
    """
    rows: list[dict] = []
    canonical_id = 1
    for members, confidence in components:
        for ds_id, cl in members:
            # matched_to: all other members in this component.
            matched = [[other_ds, other_cl] for other_ds, other_cl in members
                       if (other_ds, other_cl) != (ds_id, cl)]
            rows.append(
                {
                    "dataset_id": ds_id,
                    "original_cluster": str(cl),
                    "canonical_cluster": canonical_id,
                    "match_confidence": round(confidence, 6),
                    "matched_to": json.dumps(matched),
                }
            )
        canonical_id += 1

    df = pd.DataFrame(rows)
    df = df.sort_values(["canonical_cluster", "dataset_id", "original_cluster"])
    df = df.reset_index(drop=True)
    return df
