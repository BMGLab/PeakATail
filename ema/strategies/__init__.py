_REGISTRY = {}


def register(name):
    """Decorator to register a PeakFinderStrategy by name."""
    def decorator(cls):
        _REGISTRY[name] = cls
        return cls
    return decorator


def get_strategy(name, **kwargs):
    """Get a strategy instance by name."""
    if name not in _REGISTRY:
        raise ValueError(f"Unknown strategy '{name}'. Available: {list(_REGISTRY.keys())}")
    return _REGISTRY[name](**kwargs)


def list_strategies():
    """Return list of registered strategy names."""
    return list(_REGISTRY.keys())


# Import strategies to trigger registration
from ema.strategies.original import OriginalStrategy  # noqa: E402
from ema.strategies.lambda_poisson import LambdaPoissonStrategy  # noqa: E402
from ema.strategies.sierra_iterative import SierraIterativeStrategy  # noqa: E402
from ema.strategies.lambda_gradient import LambdaGradientStrategy  # noqa: E402
from ema.strategies.clip_seeded import ClipSeededStrategy  # noqa: E402

# Alias so all 5 registry helpers share the `list_*_strategies()` shape.
list_peak_strategies = list_strategies
