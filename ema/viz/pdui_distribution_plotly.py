"""PDUI distribution per cluster — plotly backend.

Renders interactive per-cluster violin plots of a PDUI score.

Data shape: ``(adata, score_key)`` where ``adata.obs[score_key]`` contains
per-cell PDUI scores (float in [0, 1]) and ``adata.obs["leiden"]`` contains
cluster assignments.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly

log = logging.getLogger(__name__)


@register_viz_strategy
class PduiDistributionPlotly(VizStrategy):
    """Per-cluster PDUI violin plot (plotly interactive).

    Data shape:
        ``(adata, score_key)`` tuple.
        ``adata`` — AnnData with ``obs["leiden"]`` and ``obs[score_key]``.
        ``score_key`` — column name in ``adata.obs`` holding PDUI scores.
    """

    name = "pdui_distribution_plotly"
    plot_type = "pdui_distribution"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render interactive per-cluster PDUI violin plot.

        Args:
            data: ``(adata, score_key)`` tuple.
            output_basepath: Path stem.

        Returns:
            List of paths written (HTML + SVG).
        """
        adata, score_key = data
        if score_key not in adata.obs.columns:
            log.warning(
                "pdui_distribution_plotly: score key %r not in obs columns", score_key
            )
            return []

        if "leiden" not in adata.obs.columns:
            log.warning("pdui_distribution_plotly: no 'leiden' column in obs")
            return []

        df = pd.DataFrame(
            {
                "cluster": adata.obs["leiden"].astype(str).values,
                score_key: adata.obs[score_key].astype(float).values,
            }
        )

        # Sort clusters naturally
        cluster_order = sorted(
            df["cluster"].unique(),
            key=lambda x: int(x) if x.isdigit() else x,
        )

        fig = px.violin(
            df,
            x="cluster",
            y=score_key,
            color="cluster",
            box=True,
            points="outliers",
            category_orders={"cluster": cluster_order},
            title=f"PDUI distribution per cluster ({score_key})",
            labels={"cluster": "cluster", score_key: score_key},
        )
        fig.update_layout(
            template="plotly_white",
            showlegend=False,
            width=max(500, len(cluster_order) * 80 + 200),
            height=500,
        )

        return save_plotly(fig, output_basepath)
