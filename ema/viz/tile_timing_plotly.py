"""Tile timing bar chart — plotly backend.

Input ``data``: ``list[dict]`` — per-tile timing records with schema::

    {
        "dataset_id": str,
        "chrom": str,
        "tile_idx": int,
        "tile_start": int,
        "tile_end": int,
        "wall_seconds": float,
    }

Renders interactive HTML with hoverable tooltips showing chrom, tile range, and
wall seconds.  Outliers (>2σ) are coloured red.  Multi-dataset runs produce
vertically-stacked subplots.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

from ema.viz import register_viz_strategy
from ema.viz._io import save_plotly
from ema.viz.base import VizStrategy

log = logging.getLogger(__name__)

_NORMAL_COLOR = "#4c72b0"
_OUTLIER_COLOR = "#c44e52"


def _label(rec: dict[str, Any]) -> str:
    return f"{rec['chrom']}:{rec['tile_idx']}"


def _outlier_mask(values: np.ndarray) -> np.ndarray:
    """Return boolean mask of values > mean + 2σ."""
    if len(values) < 2:
        return np.zeros(len(values), dtype=bool)
    threshold = values.mean() + 2.0 * values.std(ddof=0)
    return values > threshold


@register_viz_strategy
class TileTimingPlotly(VizStrategy):
    """Interactive bar chart of tile wall-seconds with outlier colouring."""

    name = "tile_timing_plotly"
    plot_type = "tile_timing"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render tile timing interactive chart.

        Args:
            data: List of per-tile timing dicts.  Empty list returns ``[]``.
            output_basepath: Base path (without extension) for output files.

        Returns:
            List of paths written (HTML and SVG).
        """
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots

        if not data:
            log.warning("tile_timing_plotly: empty data, skipping render")
            return []

        # Group by dataset_id (preserving insertion order)
        grouped: dict[str, list[dict[str, Any]]] = {}
        for rec in data:
            ds = rec.get("dataset_id", "default")
            grouped.setdefault(ds, []).append(rec)

        n_datasets = len(grouped)
        subplot_titles = list(grouped.keys())
        fig = make_subplots(
            rows=n_datasets,
            cols=1,
            subplot_titles=[f"Tile timing — {t}" for t in subplot_titles],
            vertical_spacing=0.15 / max(n_datasets, 1),
        )

        for row_idx, (ds_id, recs) in enumerate(grouped.items(), start=1):
            if not recs:
                continue
            sorted_recs = sorted(recs, key=lambda r: r["wall_seconds"], reverse=True)
            labels = [_label(r) for r in sorted_recs]
            values = np.array([r["wall_seconds"] for r in sorted_recs], dtype=float)
            is_outlier = _outlier_mask(values)
            colors = [_OUTLIER_COLOR if o else _NORMAL_COLOR for o in is_outlier]

            hover_texts = [
                (
                    f"<b>{r['chrom']}</b> tile {r['tile_idx']}<br>"
                    f"range: {r.get('tile_start', '?'):,}–{r.get('tile_end', '?'):,}<br>"
                    f"wall: <b>{r['wall_seconds']:.2f}s</b>"
                    + ("<br><b>OUTLIER (>2σ)</b>" if o else "")
                )
                for r, o in zip(sorted_recs, is_outlier)
            ]

            fig.add_trace(
                go.Bar(
                    x=labels,
                    y=values,
                    marker_color=colors,
                    text=hover_texts,
                    hovertemplate="%{text}<extra></extra>",
                    name=ds_id,
                    showlegend=False,
                ),
                row=row_idx,
                col=1,
            )

        fig.update_layout(
            template="plotly_white",
            height=max(400, 350 * n_datasets),
            title_text="Tile wall-second timing",
        )
        fig.update_xaxes(tickangle=45, tickfont_size=9)
        fig.update_yaxes(title_text="wall seconds")

        return save_plotly(fig, output_basepath)
