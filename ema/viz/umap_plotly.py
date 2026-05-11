"""UMAP scatter — interactive plotly."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import plotly.express as px

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly


@register_viz_strategy
class UmapPlotly(VizStrategy):
    name = "umap_plotly"
    plot_type = "umap"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        adata, ds_id = data
        coords = adata.obsm["X_umap"]
        df = {
            "UMAP 1": coords[:, 0],
            "UMAP 2": coords[:, 1],
            "cell": list(adata.obs_names),
        }
        if "leiden" in adata.obs:
            df["cluster"] = list(adata.obs["leiden"].astype(str))
        fig = px.scatter(df, x="UMAP 1", y="UMAP 2",
                         color="cluster" if "cluster" in df else None,
                         hover_data=["cell"],
                         title=f"UMAP — {ds_id}")
        fig.update_traces(marker=dict(size=5, opacity=0.85))
        fig.update_layout(template="plotly_white", legend_title_text="cluster")
        return save_plotly(fig, output_basepath)
