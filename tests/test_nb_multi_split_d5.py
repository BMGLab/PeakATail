"""D5 regression: nb_multi sample-split separates selection from inference.

Bug D5: the omnibus "which PAS are testable" filter is decided on the same cells
the LRT then tests (double-dip → ~100% significant). sample_split=True selects
PAS on half A and infers on the disjoint half B. Default is unchanged behaviour.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ema.switch_test.strategies.nb_multi import NbMultiStrategy


def _null_matrix(n_per_cluster: int, n_pas: int, seed: int = 0):
    """Two clusters, counts drawn from the SAME distribution for both clusters
    (no true difference). Under a valid test, few PAS should be significant.
    """
    rng = np.random.RandomState(seed)
    rows, idx, labs = [], [], []
    for grp in ("c1", "c2"):
        for i in range(n_per_cluster):
            rows.append(rng.poisson(5, size=n_pas))
            idx.append(f"{grp}_{i}")
            labs.append(grp)
    cm = pd.DataFrame(rows, index=idx, columns=[str(j) for j in range(n_pas)])
    return cm, pd.Series(labs, index=idx)


def test_sample_split_runs_and_returns_expected_columns():
    cm, labs = _null_matrix(40, 6)
    df = NbMultiStrategy().test(
        cm, labs, min_cells_per_group=5, n_jobs=1, sample_split=True, split_seed=0,
    )
    assert {"pvalue", "qvalue", "test_stat", "df", "dispersion", "n_cells"} <= set(df.columns)
    # Inference ran on ~half the cells → n_cells reflects the held-out split.
    if len(df):
        assert (df["n_cells"] <= 40).all()


def test_default_is_not_sample_split():
    """Default path must be unchanged (regression-safe)."""
    cm, labs = _null_matrix(30, 4)
    df_default = NbMultiStrategy().test(cm, labs, min_cells_per_group=5, n_jobs=1)
    # n_cells on the full data equals total cells for tested PAS.
    if len(df_default):
        assert (df_default["n_cells"] == 60).all()


def test_sample_split_is_deterministic():
    cm, labs = _null_matrix(40, 5)
    a = NbMultiStrategy().test(cm, labs, min_cells_per_group=5, n_jobs=1,
                               sample_split=True, split_seed=7)
    b = NbMultiStrategy().test(cm, labs, min_cells_per_group=5, n_jobs=1,
                               sample_split=True, split_seed=7)
    assert list(a.index) == list(b.index)
    if len(a):
        assert np.allclose(a["pvalue"].values, b["pvalue"].values, equal_nan=True)


def test_sample_split_too_few_cells_raises():
    cm = pd.DataFrame([[1, 2], [3, 4]], index=["a", "b"], columns=["1", "2"])
    labs = pd.Series(["c1", "c2"], index=["a", "b"])
    with pytest.raises(ValueError):
        NbMultiStrategy().test(cm, labs, sample_split=True)
