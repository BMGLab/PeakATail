"""Registry pattern for clustering strategies (mirrors switch_test/strategies)."""
from __future__ import annotations

from typing import Callable, Dict

_REGISTRY: Dict[str, Callable] = {}


def register_clustering_strategy(name: str):
    """Decorator: register a clustering function under ``name``."""
    def decorator(fn: Callable) -> Callable:
        _REGISTRY[name] = fn
        return fn
    return decorator


def get_clustering_strategy(name: str) -> Callable:
    if name not in _REGISTRY:
        raise KeyError(
            f"Unknown clustering strategy {name!r}. Available: "
            f"{sorted(_REGISTRY.keys())}"
        )
    return _REGISTRY[name]


def list_clustering_strategies() -> list[str]:
    return sorted(_REGISTRY.keys())


# Trigger registration of built-ins. Imported here so any caller that
# imports the registry sees the defaults populated.
from ema.clustering import _builtins  # noqa: E402,F401
