"""Registry for differential APA testing strategies.

Usage
-----
    from ema.switch_test.strategies import get_diff_strategy, list_diff_strategies

    strategy = get_diff_strategy("nb_pairwise")
    results  = strategy.test(count_matrix, cluster_labels,
                             cluster1="0", cluster2="1")

Registration
------------
Strategy classes self-register by calling
:func:`register_diff_strategy` **after** their class definition
(or via the ``@register_diff_strategy`` decorator style).

The imports at the bottom of this module trigger all built-in registrations.
"""

from __future__ import annotations

from typing import Type

from ema.switch_test.strategies.base import DiffAPAStrategy

_REGISTRY: dict[str, Type[DiffAPAStrategy]] = {}


def register_diff_strategy(cls: Type[DiffAPAStrategy]) -> Type[DiffAPAStrategy]:
    """Register a :class:`DiffAPAStrategy` subclass in the global registry.

    Intended for use as a class decorator::

        @register_diff_strategy
        class MyStrategy(DiffAPAStrategy):
            name = "my_strategy"
            supports_multi_condition = False
            ...

    Parameters
    ----------
    cls:
        Strategy class with a ``name`` class attribute.

    Returns
    -------
    The same class (decorator passthrough).
    """
    if not hasattr(cls, "name") or not isinstance(cls.name, str):
        raise TypeError(
            f"{cls.__qualname__} must define a string 'name' class attribute "
            "before registration."
        )
    _REGISTRY[cls.name] = cls
    return cls


def get_diff_strategy(name: str) -> DiffAPAStrategy:
    """Instantiate and return the named strategy.

    Parameters
    ----------
    name:
        Registry key (e.g. ``"fisher"``, ``"nb_pairwise"``, ``"nb_multi"``).

    Returns
    -------
    A fresh instance of the requested strategy.

    Raises
    ------
    KeyError
        If ``name`` is not in the registry.
    """
    if name not in _REGISTRY:
        available = list(_REGISTRY.keys())
        raise KeyError(
            f"Unknown diff-APA strategy '{name}'. "
            f"Available strategies: {available}"
        )
    return _REGISTRY[name]()


def list_diff_strategies() -> list[str]:
    """Return sorted list of registered strategy names."""
    return sorted(_REGISTRY.keys())


# ---------------------------------------------------------------------------
# Trigger built-in strategy registrations by importing the modules.
# Keep these at the bottom so the registry helpers above are defined first.
# ---------------------------------------------------------------------------
from ema.switch_test.strategies.fisher import FisherStrategy  # noqa: E402, F401
from ema.switch_test.strategies.nb_pairwise import NbPairwiseStrategy  # noqa: E402, F401
from ema.switch_test.strategies.nb_multi import NbMultiStrategy  # noqa: E402, F401
