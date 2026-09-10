"""Regression: select_marker_pas must hand scanpy a CATEGORICAL cluster key.

Bug (caught on a real run): scanpy's rank_genes_groups reads
``adata.obs[groupby].cat.categories`` internally, so a groupby column that
isn't categorical dtype raises "Can only use .cat accessor with a 'category'
dtype". Real clusters.h5ad stores ``canonical_cluster`` as int/object, so
``peakatail switch diff --cluster-key canonical_cluster`` crashed. select_marker_pas
now coerces the column to categorical before calling rank_genes_groups.

These tests assert exactly that coercion by intercepting rank_genes_groups and
checking the dtype it receives — deliberately NOT exercising scanpy's full
wilcoxon internals, which are version-fragile on tiny synthetic data and are
not what this fix is about (the real end-to-end path is covered by the server
run of `peakatail switch diff`).
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import ema.quantification.marker_selector as ms
from ema.quantification.marker_selector import select_marker_pas


def _adata(labels, dtype):
    n = len(labels)
    X = sp.csr_matrix(np.random.default_rng(0).poisson(1.0, (n, 6)).astype("float32"))
    obs = pd.DataFrame(
        {"canonical_cluster": pd.Series(labels, dtype=dtype)},
        index=[f"cell{i}" for i in range(n)],
    )
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[str(i) for i in range(1, 7)]))


class _Stop(Exception):
    pass


def _run_capture(monkeypatch, adata):
    """Call select_marker_pas but intercept rank_genes_groups; return the dtype
    scanpy would have seen for the groupby column."""
    captured = {}

    def fake_rgg(a, groupby, **kw):
        captured["is_categorical"] = isinstance(a.obs[groupby].dtype, pd.CategoricalDtype)
        raise _Stop  # short-circuit before scanpy's real internals

    monkeypatch.setattr(ms.sc.tl, "rank_genes_groups", fake_rgg)
    with pytest.raises(_Stop):
        select_marker_pas(adata, cluster_key="canonical_cluster")
    return captured["is_categorical"]


def test_int_cluster_key_coerced_to_categorical(monkeypatch):
    a = _adata([0, 1] * 8, dtype="int64")
    assert not isinstance(a.obs["canonical_cluster"].dtype, pd.CategoricalDtype)
    assert _run_capture(monkeypatch, a) is True


def test_object_cluster_key_coerced_to_categorical(monkeypatch):
    a = _adata(["cl_A", "cl_B"] * 8, dtype="object")
    assert _run_capture(monkeypatch, a) is True


def test_already_categorical_is_passed_through(monkeypatch):
    a = _adata(pd.Categorical(["cl_A", "cl_B"] * 8), dtype="category")
    assert _run_capture(monkeypatch, a) is True
