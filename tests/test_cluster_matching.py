"""Tests for cross-dataset cluster matching strategies.

All tests use synthetic AnnData objects only — no real BAMs or GTF files.

Test scenarios
--------------
1. Identical composition: 2 datasets with the same 5 clusters and same marker
   PAS patterns.  All strategies should find 5 canonical clusters with high
   confidence.

2. Disjoint composition: 2 datasets with completely non-overlapping cells and
   distinct marker patterns.  No matches → 10 canonical clusters, all
   confidence 0.

3. Single-dataset edge case: identity mapping, confidence 1.0.

4. ``n_top_markers`` sensitivity: marker_overlap should still find correct
   matches at n_top_markers in {5, 50, 200}.

5. JaccardCellSet strategy with identical CBs: confidence ≈ 1.0 per cluster.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Synthetic data helpers
# ---------------------------------------------------------------------------

def _make_synthetic_adata(
    n_clusters: int = 5,
    n_cells_per_cluster: int = 50,
    n_features: int = 200,
    marker_strength: float = 10.0,
    cell_prefix: str = "cell",
    feature_prefix: str = "PAS",
    rng: np.random.Generator | None = None,
    distinct_cells: bool = False,
    distinct_cell_offset: int = 0,
) -> ad.AnnData:
    """Build a synthetic AnnData with structured marker patterns per cluster.

    Each cluster has ``n_features // n_clusters`` "signature" PAS features
    that are expressed at ``marker_strength`` times the background rate.
    All other PAS features are expressed at rate 1.0 (Poisson noise).

    Args:
        n_clusters: Number of clusters to generate.
        n_cells_per_cluster: Cells per cluster.
        n_features: Total number of PAS features (columns).
        marker_strength: Expression multiplier for marker PAS in their cluster.
        cell_prefix: Prefix for cell barcode names.
        feature_prefix: Prefix for PAS feature names.
        rng: Random number generator for reproducibility.
        distinct_cells: If True, use numeric offsets to ensure barcodes are
            unique from any previously generated dataset.
        distinct_cell_offset: Integer offset added to cell indices when
            ``distinct_cells=True``.

    Returns:
        AnnData with:
        - ``X``: int32 count matrix (cells × PAS).
        - ``obs['leiden']``: cluster labels (str of int, 0-indexed).
        - Unique cell barcodes in ``obs_names``.
        - Unique PAS names in ``var_names``.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    n_cells = n_clusters * n_cells_per_cluster
    features_per_cluster = max(1, n_features // n_clusters)

    X = rng.poisson(1.0, size=(n_cells, n_features)).astype(np.float32)

    leiden_labels: list[str] = []
    for cl in range(n_clusters):
        start_cell = cl * n_cells_per_cluster
        end_cell = start_cell + n_cells_per_cluster
        start_feat = cl * features_per_cluster
        end_feat = min(start_feat + features_per_cluster, n_features)
        # Boost marker PAS for this cluster.
        X[start_cell:end_cell, start_feat:end_feat] += rng.poisson(
            marker_strength, size=(n_cells_per_cluster, end_feat - start_feat)
        ).astype(np.float32)
        leiden_labels.extend([str(cl)] * n_cells_per_cluster)

    if distinct_cells:
        cb_names = [
            f"{cell_prefix}_{i + distinct_cell_offset}" for i in range(n_cells)
        ]
    else:
        cb_names = [f"{cell_prefix}_{i}" for i in range(n_cells)]

    var_names = [f"{feature_prefix}_{i}" for i in range(n_features)]

    adata = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame({"leiden": leiden_labels}, index=cb_names),
        var=pd.DataFrame(index=var_names),
    )
    return adata


def _save_adatas(adatas: list[ad.AnnData], tmp_dir: Path) -> list[Path]:
    """Save a list of AnnData objects to temporary h5ad files.

    Args:
        adatas: AnnData objects to save.
        tmp_dir: Directory to write files into.

    Returns:
        List of paths to the saved files.
    """
    paths: list[Path] = []
    for i, adata in enumerate(adatas):
        path = tmp_dir / f"dataset_{i}.h5ad"
        adata.write_h5ad(path)
        paths.append(path)
    return paths


# ---------------------------------------------------------------------------
# Tests: marker_overlap strategy
# ---------------------------------------------------------------------------

class TestMarkerOverlapStrategy:
    """Tests for the marker_overlap cluster matching strategy."""

    def test_identical_datasets_find_five_canonical_clusters(
        self, tmp_path: Path
    ) -> None:
        """Two datasets with identical marker structure → 5 canonical clusters.

        Both datasets contain the same 5 clusters with the same marker PAS
        patterns (and the same cell barcodes). marker_overlap should pair each
        cluster in dataset A with exactly one cluster in dataset B, producing
        5 canonical IDs each shared by exactly 2 rows.
        """
        rng = np.random.default_rng(0)
        adata_a = _make_synthetic_adata(n_clusters=5, n_cells_per_cluster=50, rng=rng)

        # Re-use the same RNG seed to create an identical composition.
        rng2 = np.random.default_rng(0)
        adata_b = _make_synthetic_adata(n_clusters=5, n_cells_per_cluster=50, rng=rng2)

        paths = _save_adatas([adata_a, adata_b], tmp_path)
        ds_ids = ["sampleA", "sampleB"]

        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("marker_overlap")
        result = strategy.match(paths, ds_ids, n_top_markers=50, n_jobs=1)

        # Exactly 10 rows (5 clusters × 2 datasets).
        assert len(result) == 10, f"Expected 10 rows, got {len(result)}"

        # Exactly 5 distinct canonical IDs.
        n_canonical = result["canonical_cluster"].nunique()
        assert n_canonical == 5, (
            f"Expected 5 canonical clusters, got {n_canonical}"
        )

        # Each canonical ID should appear exactly twice (one per dataset).
        counts = result.groupby("canonical_cluster").size()
        assert (counts == 2).all(), (
            f"Some canonical IDs do not have exactly 2 members: {counts}"
        )

        # Each canonical group should contain one row from each dataset.
        for canonical_id, grp in result.groupby("canonical_cluster"):
            ds_in_group = set(grp["dataset_id"])
            assert ds_in_group == {"sampleA", "sampleB"}, (
                f"Canonical cluster {canonical_id} has unexpected datasets: {ds_in_group}"
            )

        # All confidences should be high (marker pattern is strong).
        assert (result["match_confidence"] > 0.8).all(), (
            f"Expected all confidences > 0.8, got:\n{result[['canonical_cluster','match_confidence']]}"
        )

    def test_disjoint_datasets_produce_ten_canonical_clusters(
        self, tmp_path: Path
    ) -> None:
        """Datasets with completely different cells → 10 canonical clusters.

        Each dataset has 5 clusters; cells and marker patterns do not overlap
        at all (distinct_cells=True, distinct feature names for each dataset
        → Jaccard = 0 for all pairs → all confidence = 0 and 10 unique IDs).
        """
        rng_a = np.random.default_rng(1)
        rng_b = np.random.default_rng(999)

        n_features = 200
        # Dataset A: PAS_0..PAS_199; Dataset B: PAS_200..PAS_399 (different PAS)
        adata_a = _make_synthetic_adata(
            n_clusters=5,
            n_cells_per_cluster=50,
            n_features=n_features,
            feature_prefix="PAS_A",
            cell_prefix="cellA",
            rng=rng_a,
        )
        adata_b = _make_synthetic_adata(
            n_clusters=5,
            n_cells_per_cluster=50,
            n_features=n_features,
            feature_prefix="PAS_B",
            cell_prefix="cellB",
            rng=rng_b,
        )

        paths = _save_adatas([adata_a, adata_b], tmp_path)
        ds_ids = ["sampleA", "sampleB"]

        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("marker_overlap")
        result = strategy.match(paths, ds_ids, n_top_markers=50, n_jobs=1)

        # 10 rows total (5 per dataset).
        assert len(result) == 10

        # 10 distinct canonical IDs (no matches → each cluster is its own group).
        n_canonical = result["canonical_cluster"].nunique()
        assert n_canonical == 10, (
            f"Expected 10 canonical clusters (no matches), got {n_canonical}"
        )

        # All confidences should be 0.
        assert (result["match_confidence"] == 0.0).all(), (
            f"Expected all confidences 0.0, got:\n{result[['original_cluster','match_confidence']]}"
        )

        # All matched_to should be empty lists.
        for _, row in result.iterrows():
            assert json.loads(row["matched_to"]) == [], (
                f"Expected empty matched_to, got {row['matched_to']}"
            )

    def test_single_dataset_identity_mapping(self, tmp_path: Path) -> None:
        """Single dataset → identity mapping with confidence 1.0."""
        rng = np.random.default_rng(7)
        adata = _make_synthetic_adata(n_clusters=3, n_cells_per_cluster=30, rng=rng)
        paths = _save_adatas([adata], tmp_path)

        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("marker_overlap")
        result = strategy.match(paths, ["onlyDataset"], n_top_markers=50, n_jobs=1)

        assert len(result) == 3
        assert result["canonical_cluster"].nunique() == 3
        assert (result["match_confidence"] == 1.0).all()
        for _, row in result.iterrows():
            assert json.loads(row["matched_to"]) == []
        assert set(result["dataset_id"]) == {"onlyDataset"}

    @pytest.mark.parametrize("n_top_markers", [5, 50, 200])
    def test_n_top_markers_sensitivity(
        self, tmp_path: Path, n_top_markers: int
    ) -> None:
        """marker_overlap finds 5 canonical clusters across n_top_markers values.

        Uses a very strong marker signal so that even n_top_markers=5 is
        sufficient to distinguish all clusters.
        """
        rng = np.random.default_rng(13)
        adata_a = _make_synthetic_adata(
            n_clusters=5, n_cells_per_cluster=40, n_features=200,
            marker_strength=20.0, rng=rng,
        )
        rng2 = np.random.default_rng(13)
        adata_b = _make_synthetic_adata(
            n_clusters=5, n_cells_per_cluster=40, n_features=200,
            marker_strength=20.0, rng=rng2,
        )

        # n_top_markers > n_features/n_clusters is clipped internally.
        paths = _save_adatas([adata_a, adata_b], tmp_path)
        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("marker_overlap")
        result = strategy.match(paths, ["A", "B"], n_top_markers=n_top_markers, n_jobs=1)

        n_canonical = result["canonical_cluster"].nunique()
        assert n_canonical == 5, (
            f"n_top_markers={n_top_markers}: expected 5 canonical clusters, "
            f"got {n_canonical}"
        )


# ---------------------------------------------------------------------------
# Tests: jaccard strategy
# ---------------------------------------------------------------------------

class TestJaccardCellSetStrategy:
    """Tests for the jaccard cell-set cluster matching strategy."""

    def test_identical_cells_confidence_near_one(self, tmp_path: Path) -> None:
        """Identical cell compositions → confidence ≈ 1.0 and 5 canonical clusters."""
        rng = np.random.default_rng(42)
        adata_a = _make_synthetic_adata(
            n_clusters=5, n_cells_per_cluster=50, rng=rng
        )
        # Identical barcodes and cluster labels (same seed).
        rng2 = np.random.default_rng(42)
        adata_b = _make_synthetic_adata(
            n_clusters=5, n_cells_per_cluster=50, rng=rng2
        )

        paths = _save_adatas([adata_a, adata_b], tmp_path)
        ds_ids = ["dsA", "dsB"]

        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("jaccard")
        result = strategy.match(paths, ds_ids, n_jobs=1)

        assert len(result) == 10
        assert result["canonical_cluster"].nunique() == 5

        # Each canonical ID has exactly one member per dataset.
        counts = result.groupby("canonical_cluster").size()
        assert (counts == 2).all()

        # Since barcodes are identical, Jaccard per cluster = 1.0.
        assert (result["match_confidence"] >= 0.99).all(), (
            f"Expected confidence ≈ 1.0, got:\n{result[['canonical_cluster','match_confidence']]}"
        )

    def test_disjoint_cells_no_matches(self, tmp_path: Path) -> None:
        """Datasets with completely different barcodes → 10 canonical clusters, confidence 0."""
        rng_a = np.random.default_rng(1)
        rng_b = np.random.default_rng(2)

        adata_a = _make_synthetic_adata(
            n_clusters=5, n_cells_per_cluster=30,
            cell_prefix="cellA", rng=rng_a,
        )
        adata_b = _make_synthetic_adata(
            n_clusters=5, n_cells_per_cluster=30,
            cell_prefix="cellB", rng=rng_b,
        )

        paths = _save_adatas([adata_a, adata_b], tmp_path)
        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("jaccard")
        result = strategy.match(paths, ["A", "B"], n_jobs=1)

        assert len(result) == 10
        n_canonical = result["canonical_cluster"].nunique()
        assert n_canonical == 10, (
            f"Expected 10 canonical clusters (no cell overlap), got {n_canonical}"
        )
        assert (result["match_confidence"] == 0.0).all()

    def test_single_dataset_identity_mapping(self, tmp_path: Path) -> None:
        """Single dataset → identity mapping, confidence 1.0."""
        rng = np.random.default_rng(5)
        adata = _make_synthetic_adata(n_clusters=4, n_cells_per_cluster=20, rng=rng)
        paths = _save_adatas([adata], tmp_path)

        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("jaccard")
        result = strategy.match(paths, ["solo"], n_jobs=1)

        assert len(result) == 4
        assert result["canonical_cluster"].nunique() == 4
        assert (result["match_confidence"] == 1.0).all()


# ---------------------------------------------------------------------------
# Tests: MNN strategy
# ---------------------------------------------------------------------------

class TestMNNStrategy:
    """Tests for the MNN cluster matching strategy."""

    def test_identical_datasets_find_canonical_clusters(
        self, tmp_path: Path
    ) -> None:
        """Two datasets with identical composition → clusters are matched."""
        rng = np.random.default_rng(7)
        adata_a = _make_synthetic_adata(
            n_clusters=4, n_cells_per_cluster=40, n_features=100,
            marker_strength=15.0, rng=rng,
        )
        rng2 = np.random.default_rng(7)
        adata_b = _make_synthetic_adata(
            n_clusters=4, n_cells_per_cluster=40, n_features=100,
            marker_strength=15.0, rng=rng2,
        )

        paths = _save_adatas([adata_a, adata_b], tmp_path)
        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("mnn")
        result = strategy.match(paths, ["A", "B"], n_jobs=1)

        assert len(result) == 8
        # At minimum, some clusters should be matched (not 8 canonical IDs).
        n_canonical = result["canonical_cluster"].nunique()
        assert n_canonical <= 8, f"MNN should find some matches; got {n_canonical} IDs"
        # Identical data: expect ≤ 4 canonical clusters.
        assert n_canonical <= 4, (
            f"Expected ≤ 4 canonical clusters for identical data, got {n_canonical}"
        )

    def test_single_dataset_identity_mapping(self, tmp_path: Path) -> None:
        """Single dataset → identity mapping, confidence 1.0."""
        rng = np.random.default_rng(3)
        adata = _make_synthetic_adata(n_clusters=3, n_cells_per_cluster=30, rng=rng)
        paths = _save_adatas([adata], tmp_path)

        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("mnn")
        result = strategy.match(paths, ["single"], n_jobs=1)

        assert len(result) == 3
        assert result["canonical_cluster"].nunique() == 3
        assert (result["match_confidence"] == 1.0).all()

    def test_disjoint_features_raises_no_common_features(
        self, tmp_path: Path
    ) -> None:
        """Datasets with no shared features raise ValueError."""
        rng = np.random.default_rng(11)
        adata_a = _make_synthetic_adata(
            n_clusters=2, n_cells_per_cluster=20, feature_prefix="PAS_A", rng=rng
        )
        adata_b = _make_synthetic_adata(
            n_clusters=2, n_cells_per_cluster=20, feature_prefix="PAS_B",
            rng=np.random.default_rng(12),
        )
        paths = _save_adatas([adata_a, adata_b], tmp_path)

        from ema.clustering.cross_dataset import get_match_strategy

        strategy = get_match_strategy("mnn")
        with pytest.raises(ValueError, match="no common features"):
            strategy.match(paths, ["A", "B"])


# ---------------------------------------------------------------------------
# Tests: registry / __init__.py
# ---------------------------------------------------------------------------

class TestRegistry:
    """Tests for the cross_dataset strategy registry."""

    def test_list_match_strategies_returns_all_three(self) -> None:
        """Registry should contain exactly the three built-in strategies."""
        from ema.clustering.cross_dataset import list_match_strategies

        strategies = list_match_strategies()
        assert set(strategies) == {"marker_overlap", "mnn", "jaccard"}, (
            f"Unexpected registered strategies: {strategies}"
        )

    def test_get_match_strategy_returns_correct_instance(self) -> None:
        """get_match_strategy returns the correct type for each name."""
        from ema.clustering.cross_dataset import get_match_strategy
        from ema.clustering.cross_dataset.marker_overlap import MarkerOverlapStrategy
        from ema.clustering.cross_dataset.mnn import MNNStrategy
        from ema.clustering.cross_dataset.jaccard import JaccardCellSetStrategy

        assert isinstance(get_match_strategy("marker_overlap"), MarkerOverlapStrategy)
        assert isinstance(get_match_strategy("mnn"), MNNStrategy)
        assert isinstance(get_match_strategy("jaccard"), JaccardCellSetStrategy)

    def test_get_match_strategy_raises_on_unknown_name(self) -> None:
        """get_match_strategy raises KeyError for unknown strategy names."""
        from ema.clustering.cross_dataset import get_match_strategy

        with pytest.raises(KeyError, match="Unknown cluster match strategy"):
            get_match_strategy("nonexistent_strategy")

    def test_strategy_names_match_registry_keys(self) -> None:
        """Each strategy's .name attribute matches its registry key."""
        from ema.clustering.cross_dataset import (
            _REGISTRY,
            list_match_strategies,
        )
        for key in list_match_strategies():
            cls = _REGISTRY[key]
            assert cls.name == key, (
                f"Strategy class name attribute '{cls.name}' != registry key '{key}'"
            )


# ---------------------------------------------------------------------------
# Tests: DataFrame schema validation
# ---------------------------------------------------------------------------

class TestDataFrameSchema:
    """Tests that the output DataFrame always has the required schema."""

    def test_output_columns_present(self, tmp_path: Path) -> None:
        """Result DataFrame has all required columns."""
        rng = np.random.default_rng(99)
        adata = _make_synthetic_adata(n_clusters=2, n_cells_per_cluster=20, rng=rng)
        paths = _save_adatas([adata], tmp_path)

        from ema.clustering.cross_dataset import get_match_strategy

        required_cols = {
            "dataset_id",
            "original_cluster",
            "canonical_cluster",
            "match_confidence",
            "matched_to",
        }
        for name in ["marker_overlap", "jaccard", "mnn"]:
            strategy = get_match_strategy(name)
            result = strategy.match(paths, ["ds0"], n_jobs=1)
            missing = required_cols - set(result.columns)
            assert not missing, (
                f"Strategy '{name}' result missing columns: {missing}"
            )

    def test_canonical_cluster_is_positive_integer(self, tmp_path: Path) -> None:
        """canonical_cluster values are positive integers starting at 1."""
        rng = np.random.default_rng(55)
        adata = _make_synthetic_adata(n_clusters=3, n_cells_per_cluster=25, rng=rng)
        paths = _save_adatas([adata], tmp_path)

        from ema.clustering.cross_dataset import get_match_strategy

        for name in ["marker_overlap", "jaccard", "mnn"]:
            result = get_match_strategy(name).match(paths, ["ds"], n_jobs=1)
            assert (result["canonical_cluster"] >= 1).all(), (
                f"Strategy '{name}': canonical_cluster contains values < 1"
            )
            assert result["canonical_cluster"].dtype in (
                np.dtype("int64"), np.dtype("int32"), object
            ) or pd.api.types.is_integer_dtype(result["canonical_cluster"]), (
                f"Strategy '{name}': canonical_cluster should be integer dtype"
            )

    def test_match_confidence_in_unit_interval(self, tmp_path: Path) -> None:
        """match_confidence values are in [0.0, 1.0]."""
        rng = np.random.default_rng(77)
        adata_a = _make_synthetic_adata(n_clusters=3, n_cells_per_cluster=30, rng=rng)
        rng2 = np.random.default_rng(77)
        adata_b = _make_synthetic_adata(n_clusters=3, n_cells_per_cluster=30, rng=rng2)
        paths = _save_adatas([adata_a, adata_b], tmp_path)

        from ema.clustering.cross_dataset import get_match_strategy

        for name in ["marker_overlap", "jaccard"]:
            result = get_match_strategy(name).match(paths, ["A", "B"], n_jobs=1)
            assert (result["match_confidence"] >= 0.0).all(), (
                f"Strategy '{name}': confidence below 0.0"
            )
            assert (result["match_confidence"] <= 1.0).all(), (
                f"Strategy '{name}': confidence above 1.0"
            )

    def test_matched_to_is_valid_json(self, tmp_path: Path) -> None:
        """matched_to column contains valid JSON in every row."""
        rng = np.random.default_rng(33)
        adata_a = _make_synthetic_adata(n_clusters=2, n_cells_per_cluster=25, rng=rng)
        rng2 = np.random.default_rng(33)
        adata_b = _make_synthetic_adata(n_clusters=2, n_cells_per_cluster=25, rng=rng2)
        paths = _save_adatas([adata_a, adata_b], tmp_path)

        from ema.clustering.cross_dataset import get_match_strategy

        for name in ["marker_overlap", "jaccard"]:
            result = get_match_strategy(name).match(paths, ["A", "B"], n_jobs=1)
            for _, row in result.iterrows():
                try:
                    parsed = json.loads(row["matched_to"])
                    assert isinstance(parsed, list), (
                        f"Strategy '{name}': matched_to should decode to a list"
                    )
                except json.JSONDecodeError as exc:
                    pytest.fail(
                        f"Strategy '{name}': matched_to is not valid JSON: {exc}"
                    )
