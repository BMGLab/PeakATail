"""PDUI strategy registry for PeakATail.

Usage::

    from ema.quantification.strategies import (
        get_pdui_strategy,
        list_pdui_strategies,
        register_pdui_strategy,
    )

    strategy = get_pdui_strategy("shannon")
    result_df = strategy.compute(count_matrix, pas_isoform_map)

    print(list_pdui_strategies())   # ['classic', 'proportion', 'shannon']
"""

from __future__ import annotations

from ema.quantification.strategies.base import PDUIStrategy

_REGISTRY: dict[str, type[PDUIStrategy]] = {}


def register_pdui_strategy(cls: type[PDUIStrategy]) -> type[PDUIStrategy]:
    """Class decorator that registers a PDUIStrategy subclass by its ``name``.

    Args:
        cls: A concrete subclass of :class:`PDUIStrategy` with a non-empty
            ``name`` class attribute.

    Returns:
        The class unchanged (decorator pattern).

    Raises:
        AttributeError: If ``cls.name`` is missing or empty.
    """
    if not getattr(cls, "name", None):
        raise AttributeError(
            f"PDUIStrategy subclass {cls.__qualname__} must define a non-empty "
            f"class attribute 'name'."
        )
    _REGISTRY[cls.name] = cls
    return cls


def get_pdui_strategy(name: str) -> PDUIStrategy:
    """Return a fresh instance of the named strategy.

    Args:
        name: Registry key (e.g. ``"classic"``, ``"proportion"``,
            ``"shannon"``).

    Returns:
        An instance of the corresponding :class:`PDUIStrategy` subclass.

    Raises:
        ValueError: If *name* is not registered.
    """
    if name not in _REGISTRY:
        available = list(_REGISTRY.keys())
        raise ValueError(
            f"Unknown PDUI strategy '{name}'. Available: {available}"
        )
    return _REGISTRY[name]()


def list_pdui_strategies() -> list[str]:
    """Return a sorted list of all registered strategy names."""
    return sorted(_REGISTRY.keys())


# Trigger registration by importing strategy modules
from ema.quantification.strategies.classic import ClassicPDUIStrategy  # noqa: E402, F401
from ema.quantification.strategies.proportion import ProportionPDUIStrategy  # noqa: E402, F401
from ema.quantification.strategies.shannon import ShannonPDUIStrategy  # noqa: E402, F401
