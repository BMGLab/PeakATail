"""Tests for ``build_count_dfs(adata, which=...)`` — the memory-saving switch
that lets ``run_length`` build only ``pdui_df`` and ``run_diff`` only
``diff_df`` instead of always materializing both dense frames.

The contract under test: whichever frame(s) a ``which`` value builds are
numerically identical to the original ``which="both"`` path, and the frame the
caller did not ask for is ``None`` (so no second full dense copy is made).
"""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import anndata as ad

from ema.switch_test.runner import build_count_dfs


def _make_adata(sparse: bool = False, int_names: bool = True) -> ad.AnnData:
    # 4 cells x 3 PAS, deliberately asymmetric so a stray transpose is caught.
    X = np.array(
        [[1.0, 0.0, 2.0],
         [0.0, 3.0, 0.0],
         [4.0, 5.0, 0.0],
         [0.0, 0.0, 6.0]],
        dtype=float,
    )
    var_idx = ["1", "2", "3"] if int_names else ["pasA", "pasB", "pasC"]
    var = pd.DataFrame(index=var_idx)
    obs = pd.DataFrame(index=["c0", "c1", "c2", "c3"])
    return ad.AnnData(X=sp.csr_matrix(X) if sparse else X, obs=obs, var=var)


@pytest.mark.parametrize("sparse", [False, True])
def test_both_matches_original_shapes_and_values(sparse):
    a = _make_adata(sparse=sparse)
    pdui, diff, cells, pas = build_count_dfs(a, which="both")
    # diff_df is cells x PAS, pdui_df is PAS x cells (the transpose).
    assert diff.shape == (4, 3)
    assert pdui.shape == (3, 4)
    assert list(diff.index) == cells == ["c0", "c1", "c2", "c3"]
    assert list(diff.columns) == pas == [1, 2, 3]
    assert list(pdui.index) == [1, 2, 3]
    assert list(pdui.columns) == cells
    np.testing.assert_array_equal(pdui.values, diff.values.T)


@pytest.mark.parametrize("sparse", [False, True])
def test_which_diff_builds_only_diff(sparse):
    a = _make_adata(sparse=sparse)
    both_pdui, both_diff, _, _ = build_count_dfs(a, which="both")
    pdui, diff, cells, pas = build_count_dfs(a, which="diff")
    assert pdui is None                       # the copy we skip
    assert diff is not None
    pd.testing.assert_frame_equal(diff, both_diff)   # identical to original path
    assert cells == ["c0", "c1", "c2", "c3"]
    assert pas == [1, 2, 3]


@pytest.mark.parametrize("sparse", [False, True])
def test_which_pdui_builds_only_pdui(sparse):
    a = _make_adata(sparse=sparse)
    both_pdui, both_diff, _, _ = build_count_dfs(a, which="both")
    pdui, diff, cells, pas = build_count_dfs(a, which="pdui")
    assert diff is None                       # the copy we skip
    assert pdui is not None
    # numerically identical to the both-path pdui_df (direct-transpose densify)
    pd.testing.assert_frame_equal(pdui, both_pdui)
    np.testing.assert_array_equal(pdui.values, both_diff.values.T)


def test_non_integer_var_names_fall_back_to_str_index():
    a = _make_adata(int_names=False)
    pdui, diff, cells, pas = build_count_dfs(a, which="both")
    assert pas == ["pasA", "pasB", "pasC"]
    assert list(diff.columns) == ["pasA", "pasB", "pasC"]


def test_invalid_which_raises():
    a = _make_adata()
    with pytest.raises(ValueError, match="which must be"):
        build_count_dfs(a, which="all")


def test_default_is_both():
    a = _make_adata()
    pdui, diff, _, _ = build_count_dfs(a)   # no which= => original behaviour
    assert pdui is not None and diff is not None
