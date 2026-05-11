"""UMAP via scanpy.pl.umap — matplotlib backend with scanpy's defaults."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib


@register_viz_strategy
class UmapScanpy(VizStrategy):
    name = "umap_scanpy"
    plot_type = "umap"
    engine = "scanpy"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        adata, ds_id = data
        # sc.pl.umap returns nothing useful; use show=False and grab current figure
        sc.pl.umap(adata, color="leiden" if "leiden" in adata.obs else None,
                   show=False, return_fig=False, title=f"UMAP — {ds_id}")
        fig = plt.gcf()
        return save_matplotlib(fig, output_basepath)
