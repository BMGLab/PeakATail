"""Visualization strategies registry (mirrors ema/strategies, ema/clustering/registry, ...).

Adding a new plot:
    @register_viz_strategy
    class MyPlot(VizStrategy):
        name = "myplot_matplotlib"
        plot_type = "myplot"
        engine = "matplotlib"
        def render(self, data, output_basepath):
            # write matplotlib figure(s) under output_basepath.{png,svg}
            return [output_basepath.with_suffix(".png"), output_basepath.with_suffix(".svg")]

Then in your pipeline:
    from ema.viz import render_all
    render_all(plot_type="myplot", data=df, output_basepath=Path("out/figures/myplot"),
               engines=["matplotlib", "plotly"])
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from ema.viz.base import VizStrategy

_REGISTRY: dict[str, VizStrategy] = {}


def register_viz_strategy(cls: type[VizStrategy]) -> type[VizStrategy]:
    """Class decorator: instantiate and register under cls.name."""
    instance = cls()
    _REGISTRY[instance.name] = instance
    return cls


def get_viz_strategy(name: str) -> VizStrategy:
    if name not in _REGISTRY:
        raise KeyError(
            f"Unknown viz strategy {name!r}. Available: {sorted(_REGISTRY.keys())}"
        )
    return _REGISTRY[name]


def list_viz_strategies(
    plot_type: Optional[str] = None,
    engine: Optional[str] = None,
) -> list[str]:
    """Return sorted strategy names. Optional filters by plot_type and/or engine."""
    out = []
    for name, strat in _REGISTRY.items():
        if plot_type is not None and strat.plot_type != plot_type:
            continue
        if engine is not None and strat.engine != engine:
            continue
        out.append(name)
    return sorted(out)


def render_all(
    plot_type: str,
    data: Any,
    output_basepath: Path,
    engines: list[str],
) -> list[Path]:
    """Render every registered strategy for this plot_type whose engine is in engines.

    Returns the union of file paths written.
    """
    written: list[Path] = []
    for name in list_viz_strategies(plot_type=plot_type):
        strat = _REGISTRY[name]
        if strat.engine not in engines:
            continue
        try:
            paths = strat.render(data, output_basepath)
            written.extend(paths)
        except Exception as e:
            import logging
            logging.getLogger(__name__).warning(
                "viz strategy %r failed: %s", name, e
            )
    return written


# Trigger registration of built-ins (auto-discovery via importing submodules).
# V1 will list its files here; V2/V3/V4 will append theirs in their merges.
def _autoimport():
    import importlib
    import pkgutil
    pkg = importlib.import_module(__name__)
    for _, modname, _ in pkgutil.iter_modules(pkg.__path__):
        if modname.startswith("_") or modname in ("base",):
            continue
        try:
            importlib.import_module(f"{__name__}.{modname}")
        except ImportError:
            pass


_autoimport()
