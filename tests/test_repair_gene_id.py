"""Regression tests for ``_repair_gene_id`` categorical handling.

``adata.var['gene_id']`` frequently arrives as a pandas Categorical.  Repairing
a NaN entry with a value that is not already a category (a gene id absent from
this cell type, or the ``"-"`` unknown/strand placeholder) raised
``TypeError: Cannot setitem on a Categorical with a new category`` under
pandas 2.x, which crashed every ``peakatail switch diff`` run whose var carried a
categorical ``gene_id`` (observed: fisher/mwu produced 0 output).
"""
import pandas as pd

from ema.switch_test.runner import _repair_gene_id


def test_categorical_gene_id_repair_with_dash():
    # exact reproducer of the crash: replacement "-" is not an existing category
    s = pd.Series(pd.Categorical(["ENSG1", None, "ENSG2"]), index=["10", "11", "12"])
    assert isinstance(s.dtype, pd.CategoricalDtype)
    out, n = _repair_gene_id(s, {"11": "-"})
    assert n == 1
    assert out.loc["11"] == "-"
    assert out.loc["10"] == "ENSG1"
    assert out.loc["12"] == "ENSG2"


def test_categorical_repair_new_gene_category():
    s = pd.Series(pd.Categorical(["ENSG1", None]), index=["1", "2"])
    out, n = _repair_gene_id(s, {"2": "ENSG_NEW"})
    assert n == 1
    assert out.loc["2"] == "ENSG_NEW"


def test_noncategorical_object_series_still_works():
    s = pd.Series(["ENSG1", None], index=["1", "2"], dtype=object)
    out, n = _repair_gene_id(s, {"2": "-"})
    assert n == 1
    assert out.loc["2"] == "-"


def test_empty_pas_gene_map_is_noop():
    s = pd.Series(pd.Categorical(["ENSG1", None]), index=["1", "2"])
    out, n = _repair_gene_id(s, {})
    assert n == 0
    # unchanged (still has the NaN)
    assert out.isna().sum() == 1
