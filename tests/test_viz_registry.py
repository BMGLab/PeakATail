"""Registry contract tests for ema.viz."""
from ema.viz import (
    register_viz_strategy, list_viz_strategies,
    get_viz_strategy, render_all,
)
from ema.viz.base import VizStrategy
from pathlib import Path
import pytest


def test_list_returns_sorted_list():
    names = list_viz_strategies()
    assert isinstance(names, list)
    assert names == sorted(names)


def test_filter_by_plot_type():
    names = list_viz_strategies(plot_type="umap")
    assert all("umap" in n for n in names)


def test_filter_by_engine():
    names = list_viz_strategies(engine="matplotlib")
    for n in names:
        s = get_viz_strategy(n)
        assert s.engine == "matplotlib"


def test_get_unknown_raises():
    with pytest.raises(KeyError):
        get_viz_strategy("nonexistent_xyz_123")


def test_render_all_dispatch(tmp_path):
    """render_all should iterate every registered strategy matching plot_type+engines."""
    # placeholder — no strategies yet; just verify the function signature
    paths = render_all(plot_type="umap", data=None, output_basepath=tmp_path / "umap",
                       engines=["matplotlib"])
    assert isinstance(paths, list)


def test_register_idempotent():
    """Re-registering the same name is a no-op."""
    @register_viz_strategy
    class _Tmp(VizStrategy):
        name = "test_tmp_strategy"
        plot_type = "test"
        engine = "matplotlib"
        def render(self, data, output_basepath):
            return []
    @register_viz_strategy
    class _Tmp2(VizStrategy):
        name = "test_tmp_strategy"
        plot_type = "test"
        engine = "matplotlib"
        def render(self, data, output_basepath):
            return []
    assert "test_tmp_strategy" in list_viz_strategies()
