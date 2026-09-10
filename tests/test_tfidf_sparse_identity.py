"""The sparse Signac-Method-1 TF-IDF must be bit-identical to the dense
formula it replaced (ema/clustering/strategies/leiden_tfidf.py).

The dense implementation held four dense float64 copies of the cells x PAS
matrix at once (291 GB on the PBMC 10k final run); the sparse one applies
the same scalar sequence to the stored entries only.  These tests pin the
equivalence (values, sparsity pattern, dtype) on integer count matrices,
including the edge cases the dense path handled implicitly: explicit stored
zeros, unsorted / duplicate indices, all-zero rows and columns, and dense
ndarray input.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix, csc_matrix, coo_matrix, issparse

from ema.clustering.strategies.leiden_tfidf import _tfidf_signac_method1


def _tfidf_dense_reference(X, scale_factor=1e4):
    """Verbatim copy of the pre-change dense implementation."""
    if issparse(X):
        X = X.toarray()
    X = X.astype(np.float64)
    cell_totals = X.sum(axis=1, keepdims=True)
    cell_totals[cell_totals == 0] = 1
    tf = X / cell_totals
    n_cells = X.shape[0]
    cells_with_peak = (X > 0).sum(axis=0)
    cells_with_peak[cells_with_peak == 0] = 1
    idf = n_cells / cells_with_peak
    tfidf = np.log1p(tf * idf * scale_factor)
    return csr_matrix(tfidf)


def _random_counts(seed, n_cells=60, n_peaks=200, density=0.08, max_count=120):
    rng = np.random.default_rng(seed)
    dense = np.zeros((n_cells, n_peaks), dtype=np.int64)
    mask = rng.random((n_cells, n_peaks)) < density
    dense[mask] = rng.integers(1, max_count, size=int(mask.sum()))
    dense[3, :] = 0            # an all-zero cell (cell_totals == 0 branch)
    dense[:, 7] = 0            # an all-zero PAS (cells_with_peak == 0 branch)
    dense[5, 11] = 10_000      # PAS-count dynamic range like real data
    return dense


def _assert_identical(a, b):
    a = csr_matrix(a); b = csr_matrix(b)
    a.sort_indices(); b.sort_indices()
    assert a.shape == b.shape
    assert a.dtype == b.dtype == np.float64
    assert np.array_equal(a.indptr, b.indptr)
    assert np.array_equal(a.indices, b.indices)
    assert np.array_equal(a.data, b.data), "TF-IDF values differ bitwise"


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_sparse_matches_dense_on_csr_counts(seed):
    dense = _random_counts(seed)
    _assert_identical(
        _tfidf_signac_method1(csr_matrix(dense)),
        _tfidf_dense_reference(csr_matrix(dense)),
    )


@pytest.mark.parametrize("scale", [1.0, 1e4, 12_345.678])
def test_scale_factor_is_applied_in_the_same_order(scale):
    dense = _random_counts(4)
    _assert_identical(
        _tfidf_signac_method1(csr_matrix(dense), scale_factor=scale),
        _tfidf_dense_reference(csr_matrix(dense), scale_factor=scale),
    )


def test_dense_ndarray_and_other_sparse_formats_accepted():
    dense = _random_counts(5)
    ref = _tfidf_dense_reference(dense)
    for X in (dense, csc_matrix(dense), coo_matrix(dense), csr_matrix(dense, dtype=np.float32)):
        _assert_identical(_tfidf_signac_method1(X), ref)


def test_explicit_zeros_and_duplicate_indices_are_canonicalised():
    dense = _random_counts(6)
    ref = _tfidf_dense_reference(dense)
    # explicit stored zeros
    X = csr_matrix(dense)
    X.data[:5] = 0  # turn 5 stored entries into explicit zeros
    dense2 = X.toarray()
    _assert_identical(_tfidf_signac_method1(X), _tfidf_dense_reference(dense2))
    # duplicate / unsorted COO entries summing to the same matrix
    coo = coo_matrix(dense)
    rows = np.concatenate([coo.row, coo.row[:10]])
    cols = np.concatenate([coo.col, coo.col[:10]])
    vals = np.concatenate([coo.data, np.zeros(10, dtype=coo.data.dtype)])
    perm = np.random.default_rng(0).permutation(len(vals))
    dup = coo_matrix((vals[perm], (rows[perm], cols[perm])), shape=dense.shape)
    _assert_identical(_tfidf_signac_method1(dup), ref)


def test_output_is_sparse_and_never_densifies(monkeypatch):
    """The whole point of the change: no ``toarray`` on the input."""
    dense = _random_counts(7)
    X = csr_matrix(dense)

    def _boom(*a, **k):  # pragma: no cover - only hit on regression
        raise AssertionError("TF-IDF densified its input")

    monkeypatch.setattr(csr_matrix, "toarray", _boom, raising=True)
    monkeypatch.setattr(csr_matrix, "todense", _boom, raising=True)
    out = _tfidf_signac_method1(X)
    assert issparse(out) and out.format == "csr"
    assert out.nnz == csr_matrix(dense).nnz
    assert out.dtype == np.float64
