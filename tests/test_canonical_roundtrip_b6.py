"""B6 regression: canonical_cluster is round-tripped into obs.

Bug B6: canonical_cluster_map.tsv was written but never read back into
clusters.h5ad, so cross-sample comparisons silently used non-comparable
per-sample leiden labels.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")

from ema.clustering.cross_dataset.roundtrip import (
    build_canonical_lookup,
    write_canonical_clusters,
)


def _make_h5ad(path: Path, leiden_labels: list[str]) -> None:
    n = len(leiden_labels)
    X = np.zeros((n, 3), dtype=np.float32)
    obs = pd.DataFrame({"leiden": pd.Categorical(leiden_labels)})
    obs.index = [f"cell{i}" for i in range(n)]
    ad.AnnData(X=X, obs=obs).write(path)


def test_build_canonical_lookup_stringifies_keys() -> None:
    df = pd.DataFrame(
        {
            "dataset_id": ["A", "A", "B"],
            "original_cluster": [0, 1, 0],  # int labels
            "canonical_cluster": [1, 2, 1],
        }
    )
    lut = build_canonical_lookup(df)
    assert lut[("A", "0")] == 1
    assert lut[("A", "1")] == 2
    assert lut[("B", "0")] == 1


def test_build_canonical_lookup_missing_columns_raises() -> None:
    df = pd.DataFrame({"dataset_id": ["A"], "original_cluster": ["0"]})
    with pytest.raises(KeyError):
        build_canonical_lookup(df)


def test_write_canonical_clusters_annotates_obs(tmp_path: Path) -> None:
    a = tmp_path / "A.h5ad"
    b = tmp_path / "B.h5ad"
    _make_h5ad(a, ["0", "0", "1"])
    _make_h5ad(b, ["0", "1", "1"])

    # A.leiden 0 and B.leiden 1 are the SAME canonical population (2).
    match_df = pd.DataFrame(
        {
            "dataset_id": ["A", "A", "B", "B"],
            "original_cluster": ["0", "1", "0", "1"],
            "canonical_cluster": [2, 3, 5, 2],
        }
    )
    summary = write_canonical_clusters(match_df, [(a, "A"), (b, "B")])
    assert summary == {"A": 3, "B": 3}

    aa = ad.read_h5ad(a)
    bb = ad.read_h5ad(b)
    assert list(aa.obs["canonical_cluster"]) == [2, 2, 3]
    assert list(bb.obs["canonical_cluster"]) == [5, 2, 2]
    # Cross-sample: A cluster 0 and B cluster 1 now share canonical id 2.
    assert aa.obs["canonical_cluster"].iloc[0] == bb.obs["canonical_cluster"].iloc[1]


def test_unmatched_leiden_gets_sentinel(tmp_path: Path) -> None:
    a = tmp_path / "A.h5ad"
    _make_h5ad(a, ["0", "7"])  # label "7" not in match_df
    match_df = pd.DataFrame(
        {"dataset_id": ["A"], "original_cluster": ["0"], "canonical_cluster": [1]}
    )
    write_canonical_clusters(match_df, [(a, "A")])
    aa = ad.read_h5ad(a)
    assert list(aa.obs["canonical_cluster"]) == [1, -1]
