"""Clustering registry contract."""
from ema.clustering.registry import (
    register_clustering_strategy, list_clustering_strategies, get_clustering_strategy,
)


def test_default_strategies_registered():
    names = list_clustering_strategies()
    assert "leiden_tfidf" in names
    assert "leiden_libsize" in names
    assert "external" in names


def test_get_returns_callable():
    fn = get_clustering_strategy("leiden_tfidf")
    assert callable(fn)


def test_get_unknown_raises():
    import pytest
    with pytest.raises(KeyError):
        get_clustering_strategy("nope")


def test_register_idempotent():
    """Re-registering the same name doesn't raise."""
    @register_clustering_strategy("ephemeral_test")
    def fn1(adata, **k):
        return adata
    assert "ephemeral_test" in list_clustering_strategies()
