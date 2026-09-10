"""Cell-set Jaccard cross-dataset cluster matching strategy.

Algorithm
---------
For each pair of datasets (A, B), this strategy directly computes the Jaccard
similarity between cluster cell barcodes:

    J(A_i, B_j) = |CB(A_i) ∩ CB(B_j)| / |CB(A_i) ∪ CB(B_j)|

where CB(X_k) is the set of cell barcodes in cluster k of dataset X, with
any dataset prefix stripped so that ``sampleA#ACGT`` and ``sampleB#ACGT``
match on the suffix ``ACGT``.

This is only meaningful when the same physical cells appear in both datasets
(e.g. two runs of the same BAM, or two analysis passes of the same experiment).
For genuinely different samples, use ``marker_overlap`` or ``mnn`` instead.

Hungarian assignment on the Jaccard matrix produces 1-to-1 matches; transitive
closure propagates to three or more datasets.

Tradeoffs
---------
- Very fast: no embedding, no marker computation — just set operations.
- Requires cell barcode overlap; useless for distinct biological samples.
- Perfect for regression testing: running the same sample twice should yield
  confidence ≈ 1.0 if clustering is deterministic.
"""

from __future__ import annotations

import json
import re
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# Real biological barcodes are DNA letters (with optional Cell Ranger
# "-N" suffix). Numeric or arbitrary-string suffixes (e.g. "_0", "_idx_42")
# are NOT biological barcodes and the leading "sampleX_" is part of the
# cell name itself, not a dataset prefix to strip.
_DNA_BARCODE_RE = re.compile(r"^[ACGTN]{6,}(?:-\d+)?$", re.IGNORECASE)

from ema.clustering.cross_dataset.base import ClusterMatchStrategy
from ema.clustering.cross_dataset import register_match_strategy
from ema.clustering.cross_dataset.marker_overlap import (
    _components_to_dataframe,
    _hungarian_matches,
    _identity_mapping,
    _transitive_closure,
)


def _strip_prefix(barcode: str) -> str:
    """Strip dataset prefix from a cell barcode.

    Two conventions are recognised:

    * ``<sample>#<barcode>`` -- Cell Ranger convention.  Always stripped
      because ``#`` never appears inside a real biological barcode.
    * ``<sample>_<barcode>`` -- PeakATail convention.  Only stripped when
      the suffix LOOKS like a real biological barcode (DNA letters of
      length >= 6, optionally followed by ``-<digits>`` Cell Ranger tail).
      This guards against synthetic/test names like ``cellA_0`` where the
      underscore is part of the identifier itself, not a dataset prefix.

    Args:
        barcode: Cell barcode string, possibly prefixed.

    Returns:
        Barcode suffix (the biological barcode without dataset prefix),
        or the original ``barcode`` when no recognisable prefix is found.

    Examples:
        >>> _strip_prefix("sampleA#ACGT-1")
        'ACGT-1'
        >>> _strip_prefix("sampleA_ACGTGCATAGCT")
        'ACGTGCATAGCT'
        >>> _strip_prefix("sampleA_ACGTGCATAGCT-1")
        'ACGTGCATAGCT-1'
        >>> _strip_prefix("ACGT-1")
        'ACGT-1'
        >>> _strip_prefix("cellA_0")           # not a DNA suffix -> kept
        'cellA_0'
    """
    if "#" in barcode:
        return barcode.split("#")[-1]
    if "_" in barcode:
        # rsplit: the sample-id prefix may itself contain "_" (see
        # ema.countmatrix.indexing.split_cb); the barcode half never does.
        candidate = barcode.rsplit("_", 1)[-1]
        if _DNA_BARCODE_RE.match(candidate):
            return candidate
    return barcode


def _cluster_cb_sets(adata: ad.AnnData) -> dict[str, set[str]]:
    """Build a mapping from cluster label to stripped barcode set.

    Args:
        adata: AnnData object with ``obs['leiden']`` populated.

    Returns:
        Mapping from original_cluster label → set of stripped barcodes.

    Raises:
        ValueError: If ``leiden`` column is absent.
    """
    if "leiden" not in adata.obs.columns:
        raise ValueError("AnnData has no 'leiden' column in obs.")

    result: dict[str, set[str]] = {}
    for cl, group in adata.obs.groupby("leiden", observed=True):
        stripped = {_strip_prefix(cb) for cb in group.index}
        result[str(cl)] = stripped
    return result


def _jaccard_sets(a: set[str], b: set[str]) -> float:
    """Compute Jaccard similarity between two barcode sets.

    Args:
        a: First barcode set.
        b: Second barcode set.

    Returns:
        Jaccard index in [0.0, 1.0].  Returns 0.0 if both sets are empty.
    """
    if not a and not b:
        return 0.0
    union = len(a | b)
    return len(a & b) / union if union > 0 else 0.0


def _barcode_similarity_matrix(
    cb_a: dict[str, set[str]],
    cb_b: dict[str, set[str]],
) -> tuple[np.ndarray, list[str], list[str]]:
    """Build Jaccard similarity matrix from two datasets' cluster CB sets.

    Args:
        cb_a: cluster_label → barcode set for dataset A.
        cb_b: cluster_label → barcode set for dataset B.

    Returns:
        Tuple of (similarity_matrix, labels_a, labels_b).
    """
    labels_a = sorted(cb_a.keys())
    labels_b = sorted(cb_b.keys())
    sim = np.zeros((len(labels_a), len(labels_b)), dtype=float)
    for i, la in enumerate(labels_a):
        for j, lb in enumerate(labels_b):
            sim[i, j] = _jaccard_sets(cb_a[la], cb_b[lb])
    return sim, labels_a, labels_b


class JaccardCellSetStrategy(ClusterMatchStrategy):
    """Cluster matching via direct Jaccard overlap of cell barcode sets.

    Compares the sets of cell barcodes (with dataset-prefix stripped) that
    belong to each cluster.  Useful when the same physical cells appear in
    multiple datasets (e.g. same BAM processed twice, or synthetic regression
    tests).

    For datasets with distinct cell populations (different biological samples),
    this strategy will produce no matches and should not be used — use
    ``marker_overlap`` or ``mnn`` instead.

    Tunable hyperparameters:
        None. This strategy performs pure set operations; there are no
        algorithm parameters to expose. ``n_top_markers`` and ``n_jobs`` are
        accepted for interface compatibility but have no effect.
    """

    name = "jaccard"

    def match(
        self,
        h5ad_paths: list[Path],
        dataset_ids: list[str],
        n_top_markers: int = 50,
        n_jobs: int = -1,
    ) -> pd.DataFrame:
        """Match clusters via Jaccard overlap of cell barcode sets.

        Note: ``n_top_markers`` and ``n_jobs`` are unused by this strategy.
        They are accepted for interface compatibility only.

        Args:
            h5ad_paths: Paths to h5ad files (one per dataset).
            dataset_ids: Human-readable dataset identifiers.
            n_top_markers: Unused; accepted for interface compatibility.
            n_jobs: Unused; accepted for interface compatibility.

        Returns:
            DataFrame with columns:
            dataset_id, original_cluster, canonical_cluster,
            match_confidence, matched_to.
        """
        h5ad_paths = [Path(p) for p in h5ad_paths]

        if len(h5ad_paths) == 1:
            return _identity_mapping(h5ad_paths[0], dataset_ids[0])

        # Load CB sets per cluster per dataset.
        cb_by_ds: dict[str, dict[str, set[str]]] = {}
        for path, ds_id in zip(h5ad_paths, dataset_ids):
            adata = ad.read_h5ad(path)
            cb_by_ds[ds_id] = _cluster_cb_sets(adata)

        # Build node list.
        all_nodes: list[tuple[str, str]] = [
            (ds_id, cl)
            for ds_id, cbs in cb_by_ds.items()
            for cl in cbs
        ]

        # Match every pair of datasets.
        all_edges: list[tuple[tuple[str, str], tuple[str, str], float]] = []
        for ds_a, ds_b in combinations(dataset_ids, 2):
            sim, labels_a, labels_b = _barcode_similarity_matrix(
                cb_by_ds[ds_a], cb_by_ds[ds_b]
            )
            raw_matches = _hungarian_matches(sim, labels_a, labels_b)
            for cl_a, cl_b, score in raw_matches:
                all_edges.append(((ds_a, cl_a), (ds_b, cl_b), score))

        components = _transitive_closure(all_nodes, all_edges)
        return _components_to_dataframe(components)


register_match_strategy(JaccardCellSetStrategy)
