"""PAS overlap diagram — matplotlib backend.

Uses ``upsetplot`` for all dataset counts (handles 2–N sets uniformly).
For exactly 2 datasets a simple 3-bar Venn-style bar chart is also produced
alongside the UpSet plot since upsetplot's layout can look sparse for 2 sets.

Data contract::

    data: dict[str, set[str]]
        Keys are dataset IDs; values are sets of PAS identifiers.
        At least 2 datasets required to produce a plot.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib

log = logging.getLogger(__name__)


def _build_upset_data(pas_sets: dict[str, set[str]]) -> "upsetplot.UpSet":
    """Build an UpSet plot object from a dict of sets."""
    from upsetplot import from_memberships, UpSet

    # Build list of (membership_tuple, element) for upsetplot
    all_pas = set().union(*pas_sets.values())
    ds_names = list(pas_sets.keys())
    memberships: list[tuple[str, ...]] = []
    for pas in all_pas:
        member = tuple(ds for ds in ds_names if pas in pas_sets[ds])
        if member:
            memberships.append(member)

    series = from_memberships(memberships)
    return UpSet(series, subset_size="count", show_counts=True)


def _venn2_bar(ax_venn: plt.Axes, sets: dict[str, set[str]]) -> None:
    """Simple 3-bar Venn-equivalent for exactly 2 datasets."""
    names = list(sets.keys())
    a, b = sets[names[0]], sets[names[1]]
    only_a = len(a - b)
    only_b = len(b - a)
    both = len(a & b)

    categories = [f"Only {names[0]}", "Both", f"Only {names[1]}"]
    values = [only_a, both, only_b]
    colors = ["#4c72b0", "#55a868", "#dd8452"]

    bars = ax_venn.bar(categories, values, color=colors, edgecolor="white")
    for bar, val in zip(bars, values):
        ax_venn.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.5,
            str(val),
            ha="center", va="bottom", fontsize=9,
        )
    ax_venn.set_ylabel("# PAS")
    ax_venn.set_title("PAS overlap (2 datasets)")
    ax_venn.spines["top"].set_visible(False)
    ax_venn.spines["right"].set_visible(False)


@register_viz_strategy
class PasOverlapMatplotlib(VizStrategy):
    """PAS set overlap — upsetplot (matplotlib) for any number of datasets."""

    name = "pas_overlap_matplotlib"
    plot_type = "pas_overlap"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render PAS overlap.

        Args:
            data: ``dict[str, set[str]]`` mapping dataset_id -> set of PAS IDs.
            output_basepath: Base path (no extension).

        Returns:
            List of written file paths (PNG + SVG).
        """
        pas_sets: dict[str, set[str]] = data
        n_ds = len(pas_sets)

        if n_ds < 2:
            log.warning(
                "pas_overlap_matplotlib: need ≥2 datasets, got %d — skipping", n_ds
            )
            return []

        paths: list[Path] = []

        if n_ds == 2:
            # 2-dataset: produce simple Venn-bar
            fig, ax = plt.subplots(figsize=(6, 4))
            _venn2_bar(ax, pas_sets)
            plt.tight_layout()
            paths.extend(save_matplotlib(fig, output_basepath))
        else:
            # ≥3 datasets: UpSet plot
            try:
                upset = _build_upset_data(pas_sets)
                upset_fig = upset.plot()
                # upsetplot returns a dict of axes; grab the figure from any axis
                _ax_any = next(iter(upset_fig.values())) if isinstance(upset_fig, dict) else None
                if _ax_any is not None:
                    fig = _ax_any.get_figure()
                else:
                    fig = plt.gcf()
                fig.suptitle("PAS overlap across datasets", y=1.02, fontsize=11)
                paths.extend(save_matplotlib(fig, output_basepath))
            except Exception as exc:
                log.warning("upsetplot failed (%s) — falling back to bar chart", exc)
                fig, ax = plt.subplots(figsize=(8, 5))
                _intersection_bar_fallback(ax, pas_sets)
                plt.tight_layout()
                paths.extend(save_matplotlib(fig, output_basepath))

        return paths


def _intersection_bar_fallback(ax: plt.Axes, pas_sets: dict[str, set[str]]) -> None:
    """Fallback: per-dataset unique + total bar when upsetplot is unavailable."""
    names = list(pas_sets.keys())
    union_all = set().union(*pas_sets.values())
    intersection_all = set.intersection(*[s for s in pas_sets.values()])

    totals = [len(pas_sets[n]) for n in names]
    ax.bar(names, totals, color="#4c72b0", label="total PAS")
    ax.axhline(len(intersection_all), color="#c44e52", linestyle="--",
               label=f"shared by all ({len(intersection_all)})")
    ax.set_ylabel("# PAS")
    ax.set_title("PAS per dataset (overlap bar)")
    ax.legend(fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
