"""Resource timeline (RAM + CPU over time) — plotly backend.

Input ``data``: ``dict`` with two keys::

    {
        "samples": list[dict],       # periodic resource samples
        "annotations": list[dict],   # stage-transition markers (optional)
    }

Each sample record::

    {"elapsed_s": float, "rss_gb": float, "cpu_pct": float}

Each annotation record::

    {"label": str, "elapsed_s": float}

Renders an interactive HTML chart with dual y-axes (RAM GB left, CPU% right)
and vertical dashed lines for stage transitions with hover labels.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from ema.viz import register_viz_strategy
from ema.viz._io import save_plotly
from ema.viz.base import VizStrategy

log = logging.getLogger(__name__)


@register_viz_strategy
class ResourceTimelinePlotly(VizStrategy):
    """Interactive dual-axis chart: RAM (GB) and CPU% over elapsed time."""

    name = "resource_timeline_plotly"
    plot_type = "resource_timeline"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render the resource timeline as an interactive plotly figure.

        Args:
            data: Dict with ``"samples"`` and ``"annotations"`` keys.
            output_basepath: Base path (without extension) for output files.

        Returns:
            List of paths written (HTML and SVG).  Returns ``[]`` if samples
            list is empty.
        """
        import plotly.graph_objects as go

        samples: list[dict[str, Any]] = data.get("samples", [])
        annotations: list[dict[str, Any]] = data.get("annotations", [])

        if not samples:
            log.warning("resource_timeline_plotly: empty samples, skipping render")
            return []

        elapsed = [s["elapsed_s"] for s in samples]
        rss_gb = [s["rss_gb"] for s in samples]
        cpu_pct = [s["cpu_pct"] for s in samples]

        fig = go.Figure()

        # RAM trace — primary y-axis
        fig.add_trace(
            go.Scatter(
                x=elapsed,
                y=rss_gb,
                name="RAM (GB)",
                mode="lines",
                line=dict(color="#4c72b0", width=2),
                fill="tozeroy",
                fillcolor="rgba(76, 114, 176, 0.1)",
                yaxis="y1",
                hovertemplate="t=%{x:.1f}s<br>RAM=%{y:.3f} GB<extra></extra>",
            )
        )

        # CPU trace — secondary y-axis
        fig.add_trace(
            go.Scatter(
                x=elapsed,
                y=cpu_pct,
                name="CPU (%)",
                mode="lines",
                line=dict(color="#dd8452", width=1.5, dash="dash"),
                yaxis="y2",
                hovertemplate="t=%{x:.1f}s<br>CPU=%{y:.1f}%<extra></extra>",
            )
        )

        # Stage-transition vertical lines via shapes + annotations
        shapes = []
        ann_labels = []
        for ann in annotations:
            x_val = ann.get("elapsed_s")
            label = ann.get("label", "")
            if x_val is None:
                continue
            shapes.append(
                dict(
                    type="line",
                    x0=x_val, x1=x_val,
                    y0=0, y1=1,
                    xref="x", yref="paper",
                    line=dict(color="#888888", width=1, dash="dot"),
                )
            )
            ann_labels.append(
                dict(
                    x=x_val,
                    y=1.0,
                    xref="x",
                    yref="paper",
                    text=label,
                    showarrow=False,
                    font=dict(size=9, color="#666666"),
                    textangle=-90,
                    xanchor="left",
                    yanchor="top",
                )
            )

        _cpu_max = max((s["cpu_pct"] for s in samples), default=100.0)
        fig.update_layout(
            template="plotly_white",
            title="Resource usage over time",
            xaxis=dict(title="Elapsed time (s)"),
            yaxis=dict(
                title="RAM (GB)",
                title_font=dict(color="#4c72b0"),
                tickfont=dict(color="#4c72b0"),
                rangemode="tozero",
            ),
            yaxis2=dict(
                title="CPU (%)",
                title_font=dict(color="#dd8452"),
                tickfont=dict(color="#dd8452"),
                overlaying="y",
                side="right",
                range=[0, max(105, _cpu_max * 1.1)],
            ),
            legend=dict(x=0.01, y=0.99, xanchor="left", yanchor="top"),
            shapes=shapes,
            annotations=ann_labels,
            height=480,
        )

        return save_plotly(fig, output_basepath)
