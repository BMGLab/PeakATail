"""Clustering strategy registry.

Follows the same pattern as ema.strategies for peak finding.
Each strategy handles: normalization -> dimensionality reduction -> clustering.
"""

_REGISTRY = {}


def register(name):
    """Decorator to register a ClusteringStrategy by name."""
    def decorator(cls):
        _REGISTRY[name] = cls
        return cls
    return decorator


def get_strategy(name, **kwargs):
    """Get a clustering strategy instance by name."""
    if name not in _REGISTRY:
        raise ValueError(
            f"Unknown clustering strategy '{name}'. "
            f"Available: {list(_REGISTRY.keys())}"
        )
    return _REGISTRY[name](**kwargs)


def list_strategies():
    """Return list of registered clustering strategy names."""
    return list(_REGISTRY.keys())


# Import strategies to trigger registration
from ema.clustering.strategies.leiden_libsize import LeidenLibsizeStrategy  # noqa: E402
from ema.clustering.strategies.leiden_tfidf import LeidenTfidfStrategy  # noqa: E402
from ema.clustering.strategies.external import ExternalStrategy  # noqa: E402
