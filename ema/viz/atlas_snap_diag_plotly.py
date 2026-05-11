"""Atlas-snap diagnostic diagram — plotly backend.

Two-panel interactive figure:
  Left  — bar chart of snapped vs. unsnapped peak counts.
  Right — histogram of snap distances (bp) for the snapped peaks.

Data contract::

    data: dict with keys:
        "snapped"        : int  — number of peaks that matched within distance
        "unsnapped"      : int  — number of peaks beyond the threshold
        "snap_distances" : list[int]  — distance in bp for each snapped peak
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

log = logging.getLogger(__name__)


@register_viz_strategy
class AtlasSnapDiagPlotly(VizStrategy):
    """Atlas-snap diagnostic — 2-panel plotly figure."""

    name = "atlas_snap_diag_plotly"
    plot_type = "atlas_snap_diag"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render the 2-panel atlas-snap diagnostic.

        Args:
            data: Dict with keys ``snapped``, ``unsnapped``, ``snap_distances``.
            output_basepath: Base path (no extension); HTML + SVG are written.

        Returns:
            List of written file paths.
        """
        snapped: int = int(data.get("snapped", 0))
        unsnapped: int = int(data.get("unsnapped", 0))
        distances: list[int] = list(data.get("snap_distances", []))

        total = snapped + unsnapped
        pct_label = f"{100 * snapped / total:.1f}%" if total > 0 else "N/A"

        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=(
                f"Atlas snap result ({pct_label} snapped)",
                "Snap-distance distribution",
            ),
            column_widths=[0.35, 0.65],
        )

        # ---- Left: snapped / unsnapped bar ----
        fig.add_trace(
            go.Bar(
                x=["Snapped", "Unsnapped"],
                y=[snapped, unsnapped],
                marker_color=["#55a868", "#c44e52"],
                text=[str(snapped), str(unsnapped)],
                textposition="outside",
                hovertemplate="<b>%{x}</b>: %{y}<extra></extra>",
                showlegend=False,
            ),
            row=1, col=1,
        )

        # ---- Right: distance histogram ----
        if distances:
            import statistics
            median_dist = statistics.median(distances)
            fig.add_trace(
                go.Histogram(
                    x=distances,
                    nbinsx=min(50, max(10, int(len(distances) ** 0.5))),
                    marker_color="#4c72b0",
                    name="snap distance",
                    hovertemplate="Distance %{x} bp: %{y} peaks<extra></extra>",
                    showlegend=False,
                ),
                row=1, col=2,
            )
            # Median vline via shape
            fig.add_vline(
                x=median_dist,
                line_dash="dash",
                line_color="#c44e52",
                annotation_text=f"median {median_dist:.0f} bp",
                annotation_position="top right",
                row=1, col=2,  # type: ignore[call-arg]
            )
        else:
            fig.add_annotation(
                text="No snapped peaks",
                x=0.75, y=0.5, xref="paper", yref="paper",
                showarrow=False, font=dict(size=14, color="gray"),
            )

        fig.update_layout(
            template="plotly_white",
            height=400,
            bargap=0.3,
        )
        fig.update_yaxes(title_text="# peaks", row=1, col=1)
        fig.update_xaxes(title_text="Distance to atlas (bp)", row=1, col=2)
        fig.update_yaxes(title_text="# peaks", row=1, col=2)

        return save_plotly(fig, output_basepath)
