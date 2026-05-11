"""UMAP scatter — interactive plotly."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import plotly.express as px

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly
from ema.viz._meta import write_figure_meta


@register_viz_strategy
class UmapPlotly(VizStrategy):
    name = "umap_plotly"
    plot_type = "umap"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        adata, ds_id = data
        coords = adata.obsm["X_umap"]
        color_key = "leiden" if "leiden" in adata.obs else None
        df = {
            "UMAP 1": coords[:, 0],
            "UMAP 2": coords[:, 1],
            "cell": list(adata.obs_names),
        }
        if color_key is not None:
            df["cluster"] = list(adata.obs["leiden"].astype(str))
        fig = px.scatter(df, x="UMAP 1", y="UMAP 2",
                         color="cluster" if "cluster" in df else None,
                         hover_data=["cell"],
                         title=f"UMAP — {ds_id}")
        fig.update_traces(marker=dict(size=5, opacity=0.85))
        fig.update_layout(template="plotly_white", legend_title_text="cluster")
        paths = save_plotly(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "dataset_id": ds_id,
            "n_cells": int(adata.n_obs),
            "n_genes_used": int(adata.n_vars),
            "color_key": color_key,
        })
        return paths
