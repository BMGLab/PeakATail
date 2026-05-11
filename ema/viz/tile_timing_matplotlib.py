"""Tile timing bar chart — matplotlib backend.

Input ``data``: ``list[dict]`` — per-tile timing records with schema::

    {
        "dataset_id": str,
        "chrom": str,
        "tile_idx": int,
        "tile_start": int,
        "tile_end": int,
        "wall_seconds": float,
    }

For a single-dataset run a single bar chart is produced, sorted descending by
``wall_seconds`` and with tiles more than 2 standard deviations above the mean
highlighted in red.

For multi-dataset runs one subplot per dataset is produced (shared y-axis).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ema.viz import register_viz_strategy
from ema.viz._io import save_matplotlib
from ema.viz._meta import write_figure_meta
from ema.viz.base import VizStrategy

log = logging.getLogger(__name__)

# Colour palette
_NORMAL_COLOR = "#4c72b0"
_OUTLIER_COLOR = "#c44e52"


def _label(rec: dict[str, Any]) -> str:
    """Short bar label: chrom:tile_idx."""
    return f"{rec['chrom']}:{rec['tile_idx']}"


def _plot_dataset(
    ax: Any,
    records: list[dict[str, Any]],
    ds_id: str,
) -> None:
    """Render one bar chart for a single dataset into *ax*."""
    if not records:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(ds_id)
        return

    # Sort descending by wall_seconds so slowest tiles are on the left
    sorted_recs = sorted(records, key=lambda r: r["wall_seconds"], reverse=True)
    labels = [_label(r) for r in sorted_recs]
    values = np.array([r["wall_seconds"] for r in sorted_recs], dtype=float)

    # Outlier detection: > mean + 2σ
    mean_v = values.mean()
    std_v = values.std(ddof=0) if len(values) > 1 else 0.0
    threshold = mean_v + 2.0 * std_v
    colors = [_OUTLIER_COLOR if v > threshold else _NORMAL_COLOR for v in values]

    x_pos = np.arange(len(labels))
    ax.bar(x_pos, values, color=colors)

    # X-axis labels — rotate and thin if crowded
    step = max(1, len(labels) // 30)
    ax.set_xticks(x_pos[::step])
    ax.set_xticklabels(labels[::step], rotation=45, ha="right", fontsize=7)

    ax.set_ylabel("wall seconds")
    ax.set_title(f"Tile timing — {ds_id}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend patch for outliers
    if any(c == _OUTLIER_COLOR for c in colors):
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor=_OUTLIER_COLOR, label=">2σ outlier"),
            Patch(facecolor=_NORMAL_COLOR, label="normal"),
        ]
        ax.legend(handles=legend_elements, fontsize=8, frameon=False, loc="upper right")


@register_viz_strategy
class TileTimingMatplotlib(VizStrategy):
    """Bar chart of tile wall-seconds, outliers highlighted, per dataset."""

    name = "tile_timing_matplotlib"
    plot_type = "tile_timing"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render tile timing bars.

        Args:
            data: List of per-tile timing dicts.  Empty list returns ``[]``.
            output_basepath: Base path (without extension) for output files.

        Returns:
            List of paths written (PNG and SVG).
        """
        if not data:
            log.warning("tile_timing_matplotlib: empty data, skipping render")
            return []

        # Group by dataset_id
        grouped: dict[str, list[dict[str, Any]]] = {}
        for rec in data:
            ds = rec.get("dataset_id", "default")
            grouped.setdefault(ds, []).append(rec)

        n_datasets = len(grouped)
        fig_height = max(4, 4 * n_datasets)
        fig, axes = plt.subplots(
            n_datasets, 1, figsize=(14, fig_height), squeeze=False
        )

        for ax, (ds_id, recs) in zip(axes.flatten(), grouped.items()):
            _plot_dataset(ax, recs, ds_id)

        plt.tight_layout()
        paths = save_matplotlib(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "n_datasets": len(grouped),
            "dataset_ids": list(grouped.keys()),
            "n_tiles_total": len(data),
        })
        return paths
