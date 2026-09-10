"""PAS overlap diagram — plotly backend.

Plotly has no native UpSet widget, so we render a horizontal bar chart of
intersection sizes grouped by membership pattern (UpSet-style decomposition).
Each bar represents one unique membership combination; bars are coloured by
the number of datasets in the intersection.

Data contract::

    data: dict[str, set[str]]
        Keys are dataset IDs; values are sets of PAS identifiers.
        At least 2 datasets required to produce a plot.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly
from ema.viz._meta import write_figure_meta

log = logging.getLogger(__name__)

# Colour palette for membership size (1..N datasets)
_PALETTE = [
    "#aec6cf", "#779ecb", "#4c72b0", "#2d4e8c",
    "#1b2f5c", "#0d1b3e", "#060d1f",
]


def _build_intersection_table(
    pas_sets: dict[str, set[str]],
) -> list[tuple[str, int, int]]:
    """Decompose pas_sets into (membership_label, count, n_members).

    Returns list sorted descending by count.
    """
    ds_names = list(pas_sets.keys())
    all_pas = set().union(*pas_sets.values())

    # Map each PAS to its membership frozenset
    membership_counts: dict[frozenset[str], int] = {}
    for pas in all_pas:
        member = frozenset(ds for ds in ds_names if pas in pas_sets[ds])
        membership_counts[member] = membership_counts.get(member, 0) + 1

    rows: list[tuple[str, int, int]] = []
    for fset, count in membership_counts.items():
        label = " & ".join(sorted(fset)) if fset else "(none)"
        rows.append((label, count, len(fset)))

    rows.sort(key=lambda r: -r[1])
    return rows


@register_viz_strategy
class PasOverlapPlotly(VizStrategy):
    """PAS set overlap — interactive bar of intersection sizes (plotly)."""

    name = "pas_overlap_plotly"
    plot_type = "pas_overlap"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render PAS overlap intersection bar chart.

        Args:
            data: ``dict[str, set[str]]`` mapping dataset_id -> set of PAS IDs.
            output_basepath: Base path (no extension).

        Returns:
            List of written file paths (HTML + SVG).
        """
        pas_sets: dict[str, set[str]] = data
        n_ds = len(pas_sets)

        if n_ds < 2:
            log.warning(
                "pas_overlap_plotly: need ≥2 datasets, got %d — skipping", n_ds
            )
            return []

        rows = _build_intersection_table(pas_sets)
        labels = [r[0] for r in rows]
        counts = [r[1] for r in rows]
        n_members = [r[2] for r in rows]

        # Colour by number of members in intersection
        colors = [_PALETTE[min(n - 1, len(_PALETTE) - 1)] for n in n_members]

        fig = go.Figure(
            go.Bar(
                x=counts,
                y=labels,
                orientation="h",
                marker_color=colors,
                text=[str(c) for c in counts],
                textposition="outside",
                hovertemplate=(
                    "<b>%{y}</b><br>PAS count: %{x}<extra></extra>"
                ),
            )
        )
        fig.update_layout(
            title="PAS overlap — intersection sizes",
            xaxis_title="# PAS",
            yaxis_title="Membership",
            template="plotly_white",
            height=max(300, 40 * len(rows) + 100),
            margin=dict(l=200, r=40, t=60, b=40),
        )
        paths = save_plotly(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "n_datasets": n_ds,
            "dataset_ids": list(pas_sets.keys()),
            "dataset_sizes": {k: len(v) for k, v in pas_sets.items()},
        })
        return paths
