"""Abstract base class for cross-dataset cluster matching strategies.

Each strategy maps cluster labels from multiple h5ad files (each with their own
independently-assigned leiden cluster labels) to a shared canonical ID space.
Different datasets may label the same cell population with different integers
(cluster_0 in sampleA may correspond to cluster_5 in sampleB). Strategies find
these correspondences and return a unified mapping DataFrame.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd


class ClusterMatchStrategy(ABC):
    """Base class for cross-dataset cluster matching strategies.

    Subclasses implement ``match`` to find cluster correspondences across
    multiple h5ad files and assign canonical integer cluster IDs.

    Class Attributes:
        name: Registry key used to look up this strategy.
    """

    name: str

    @abstractmethod
    def match(
        self,
        h5ad_paths: list[Path],
        dataset_ids: list[str],
        n_top_markers: int = 50,
        n_jobs: int = -1,
    ) -> pd.DataFrame:
        """Find cluster correspondences across datasets.

        Each dataset has independently-assigned leiden cluster labels
        (``adata.obs['leiden']``). This method identifies which clusters
        across datasets represent the same underlying cell population and
        assigns them a shared canonical integer ID.

        Args:
            h5ad_paths: Absolute paths to h5ad files, one per dataset.
                Each file must have ``adata.obs['leiden']`` populated.
            dataset_ids: Human-readable identifiers for each dataset,
                parallel to ``h5ad_paths``.
            n_top_markers: Number of top marker features (PAS) to use per
                cluster when computing similarity. Strategies that do not use
                marker genes (e.g. ``jaccard``) may ignore this parameter.
            n_jobs: Number of parallel workers for expensive computations.
                ``-1`` uses all available CPUs (passed to joblib).

        Returns:
            DataFrame with the following columns:

            - ``dataset_id`` (str): Dataset identifier from ``dataset_ids``.
            - ``original_cluster`` (str): Original leiden label in that dataset
              (e.g. ``"0"``, ``"1"``).
            - ``canonical_cluster`` (int): Shared integer ID (1, 2, 3, ...)
              assigned to groups of corresponding clusters across datasets.
              Unmatched clusters receive a unique ID not shared with any other
              row.
            - ``match_confidence`` (float): Confidence in [0.0, 1.0] for the
              match. Groups of corresponding clusters report the average
              pairwise similarity. Unmatched clusters report 0.0.
            - ``matched_to`` (str): JSON-encoded list of
              ``[dataset_id, original_cluster]`` pairs that share this
              canonical ID (excluding self). Empty list ``"[]"`` for unmatched
              clusters.

        Notes:
            - Single-dataset input returns an identity mapping: each cluster
              maps to its own canonical ID (starting at 1) with confidence 1.0
              and ``matched_to="[]"``.
            - Canonical IDs are dense integers starting at 1. The mapping is
              consistent within one call but not comparable across separate
              calls.
        """
        ...
