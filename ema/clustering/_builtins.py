"""Register the 3 built-in clustering strategies onto the registry.

Algorithm code itself stays in ema/clustering/clustering.py. This file
just thin-wraps the existing if/elif branches into registry entries.
"""
from __future__ import annotations

from ema.clustering.registry import register_clustering_strategy


@register_clustering_strategy("leiden_tfidf")
def _leiden_tfidf(adata, **kwargs):
    from ema.clustering.clustering import _do_leiden_tfidf
    return _do_leiden_tfidf(adata, **kwargs)


@register_clustering_strategy("leiden_libsize")
def _leiden_libsize(adata, **kwargs):
    from ema.clustering.clustering import _do_leiden_libsize
    return _do_leiden_libsize(adata, **kwargs)


@register_clustering_strategy("external")
def _external(adata, **kwargs):
    from ema.clustering.clustering import _do_external
    return _do_external(adata, **kwargs)
