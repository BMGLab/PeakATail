"""UMAP scatter plot — matplotlib backend."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib
from ema.viz._meta import write_figure_meta


@register_viz_strategy
class UmapMatplotlib(VizStrategy):
    name = "umap_matplotlib"
    plot_type = "umap"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        adata, ds_id = data
        coords = adata.obsm["X_umap"]
        labels = adata.obs.get("leiden")
        color_key = "leiden" if labels is not None else None
        fig, ax = plt.subplots(figsize=(7, 6))
        if labels is not None:
            cats = sorted(set(labels))
            cmap = plt.cm.tab20(np.linspace(0, 1, len(cats)))
            for i, c in enumerate(cats):
                mask = labels == c
                ax.scatter(coords[mask, 0], coords[mask, 1], s=8, c=[cmap[i]],
                           label=f"cluster {c}", alpha=0.8)
            ax.legend(loc="best", fontsize=8, frameon=False)
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=8, alpha=0.8)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_title(f"UMAP — {ds_id}")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        paths = save_matplotlib(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "dataset_id": ds_id,
            "n_cells": int(adata.n_obs),
            "n_genes_used": int(adata.n_vars),
            "color_key": color_key,
        })
        return paths
