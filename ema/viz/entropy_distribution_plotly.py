"""Shannon entropy distribution per cluster — plotly backend."""
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

log = logging.getLogger(__name__)


@register_viz_strategy
class EntropyDistributionPlotly(VizStrategy):
    """Per-cluster Shannon-entropy violin plot (plotly interactive)."""

    name = "entropy_distribution_plotly"
    plot_type = "entropy_distribution"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        adata, score_key = data
        if score_key not in adata.obs.columns:
            log.warning(
                "entropy_distribution_plotly: score key %r not in obs columns",
                score_key,
            )
            return []
        if "leiden" not in adata.obs.columns:
            log.warning("entropy_distribution_plotly: no 'leiden' column in obs")
            return []

        df = pd.DataFrame({
            "cluster": adata.obs["leiden"].astype(str).values,
            score_key: adata.obs[score_key].astype(float).values,
        }).dropna(subset=[score_key])

        cluster_order = sorted(
            df["cluster"].unique(),
            key=lambda x: int(x) if x.isdigit() else x,
        )

        fig = px.violin(
            df, x="cluster", y=score_key, color="cluster",
            box=True, points="outliers",
            category_orders={"cluster": cluster_order},
            title=f"Shannon entropy per cluster ({score_key})",
            labels={"cluster": "cluster", score_key: score_key},
        )
        fig.update_layout(
            template="plotly_white",
            showlegend=False,
            width=max(500, len(cluster_order) * 80 + 200),
            height=500,
        )

        paths = save_plotly(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "cluster_key": "leiden",
            "score_key": score_key,
            "n_clusters": len(cluster_order),
            "n_observations": int(len(df)),
        })
        return paths
