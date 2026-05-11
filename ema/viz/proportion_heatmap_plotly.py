"""Per-PAS proportion heatmap (PAS x cluster) — plotly backend."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly
from ema.viz._meta import write_figure_meta
from ema.viz.proportion_heatmap_matplotlib import _build_pas_cluster_matrix

log = logging.getLogger(__name__)


@register_viz_strategy
class ProportionHeatmapPlotly(VizStrategy):
    """PAS x cluster mean-proportion heatmap (plotly interactive)."""

    name = "proportion_heatmap_plotly"
    plot_type = "proportion_heatmap"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        if isinstance(data, dict):
            pdui_df = data["pdui_df"]
            adata = data["adata"]
            cluster_key = data.get("cluster_key", "leiden")
            top_n = int(data.get("top_n", 50))
        else:
            pdui_df, adata = data
            cluster_key, top_n = "leiden", 50

        pivot, clusters = _build_pas_cluster_matrix(
            pdui_df, adata, cluster_key, top_n,
        )
        if pivot is None or pivot.empty:
            log.warning(
                "proportion_heatmap_plotly: no rows to plot (cluster_key=%r)",
                cluster_key,
            )
            return []

        fig = px.imshow(
            pivot.fillna(0.0),
            x=clusters,
            y=[str(p) for p in pivot.index],
            color_continuous_scale="Viridis",
            zmin=0.0, zmax=1.0,
            labels={"x": "cluster", "y": "PAS id", "color": "mean proportion"},
            title="Mean PAS proportion per cluster",
            aspect="auto",
        )
        fig.update_layout(
            template="plotly_white",
            width=max(500, pivot.shape[1] * 60 + 200),
            height=max(400, min(1200, 18 * pivot.shape[0])),
        )

        paths = save_plotly(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "cluster_key": cluster_key,
            "n_pas_shown": int(pivot.shape[0]),
            "n_clusters": int(pivot.shape[1]),
            "top_n_selected_by": "cross-cluster variance",
        })
        return paths
