"""Cross-dataset cluster matching strategy registry.

This module provides a pluggable registry for strategies that find cluster
correspondences across multiple h5ad files and assign canonical cluster IDs.

Each strategy implements ``ClusterMatchStrategy.match`` and is registered
under a short name string.  New strategies can be added without modifying the
orchestration layer — just create a new module and call
``register_match_strategy``.

Usage
-----
    from ema.clustering.cross_dataset import (
        get_match_strategy,
        list_match_strategies,
    )

    strategy = get_match_strategy("marker_overlap")
    mapping = strategy.match(h5ad_paths, dataset_ids, n_top_markers=50)

Registry contract
-----------------
The registry stores ``type[ClusterMatchStrategy]`` (the class, not an
instance).  ``get_match_strategy`` instantiates the class with no arguments
— strategies must be usable with default parameters.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ema.clustering.cross_dataset.base import ClusterMatchStrategy

_REGISTRY: dict[str, type["ClusterMatchStrategy"]] = {}


def register_match_strategy(cls: type["ClusterMatchStrategy"]) -> type["ClusterMatchStrategy"]:
    """Register a ClusterMatchStrategy class under its ``name`` attribute.

    This function is intended to be called at module level in each strategy
    module, after the class definition.

    Args:
        cls: A concrete subclass of ``ClusterMatchStrategy`` with a ``name``
            class attribute.

    Returns:
        The class unchanged (allows use as a decorator if desired).

    Raises:
        AttributeError: If ``cls`` does not have a ``name`` attribute.
    """
    _REGISTRY[cls.name] = cls
    return cls


def get_match_strategy(
    name: str,
    mnn_components: int = 30,
    mnn_k_neighbors: int = 10,
) -> "ClusterMatchStrategy":
    """Instantiate and return a registered cluster match strategy by name.

    MNN-specific kwargs (``mnn_components``, ``mnn_k_neighbors``) are forwarded
    only to :class:`MNNStrategy`; they are silently ignored for other strategies
    so callers can pass them unconditionally.

    Args:
        name: Strategy name as registered (e.g. ``"marker_overlap"``).
        mnn_components: Number of LSI components for the MNN shared embedding.
            Only used when ``name="mnn"``. Default 30.
        mnn_k_neighbors: Number of nearest neighbours for MNN search.
            Only used when ``name="mnn"``. Default 10.

    Returns:
        A fresh instance of the requested strategy.

    Raises:
        KeyError: If ``name`` is not a registered strategy.
    """
    if name not in _REGISTRY:
        available = list(_REGISTRY.keys())
        raise KeyError(
            f"Unknown cluster match strategy '{name}'. "
            f"Available strategies: {available}"
        )
    cls = _REGISTRY[name]
    if name == "mnn":
        return cls(n_components=mnn_components, k_neighbors=mnn_k_neighbors)
    return cls()


def list_match_strategies() -> list[str]:
    """Return a sorted list of all registered cluster match strategy names.

    Returns:
        List of strategy name strings in alphabetical order.
    """
    return sorted(_REGISTRY.keys())


# Import strategy modules to trigger registration.
# Order does not matter for correctness; alphabetical for readability.
from ema.clustering.cross_dataset.jaccard import JaccardCellSetStrategy  # noqa: E402, F401
from ema.clustering.cross_dataset.marker_overlap import MarkerOverlapStrategy  # noqa: E402, F401
from ema.clustering.cross_dataset.mnn import MNNStrategy  # noqa: E402, F401
