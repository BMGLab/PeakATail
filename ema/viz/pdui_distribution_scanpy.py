"""PDUI distribution per cluster — scanpy violin wrapper.

Uses ``scanpy.pl.violin`` which produces a publication-quality grouped violin
plot.  Output is captured from the current matplotlib figure and saved via
``save_matplotlib``.

Data shape: ``(adata, score_key)`` — same as the matplotlib variant.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib
from ema.viz._meta import write_figure_meta

log = logging.getLogger(__name__)


@register_viz_strategy
class PduiDistributionScanpy(VizStrategy):
    """Per-cluster PDUI violin plot via scanpy.pl.violin.

    Data shape:
        ``(adata, score_key)`` tuple.
        ``adata`` — AnnData with ``obs["leiden"]`` and ``obs[score_key]``.
        ``score_key`` — column name in ``adata.obs`` holding PDUI scores.
    """

    name = "pdui_distribution_scanpy"
    plot_type = "pdui_distribution"
    engine = "scanpy"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render per-cluster PDUI violin using scanpy.

        Args:
            data: ``(adata, score_key)`` tuple.
            output_basepath: Path stem.

        Returns:
            List of paths written (PNG + SVG).
        """
        try:
            import scanpy as sc
        except ImportError:
            log.warning("pdui_distribution_scanpy: scanpy not available, skipping")
            return []

        adata, score_key = data
        if score_key not in adata.obs.columns:
            log.warning(
                "pdui_distribution_scanpy: score key %r not in obs columns", score_key
            )
            return []
        if "leiden" not in adata.obs.columns:
            log.warning("pdui_distribution_scanpy: no 'leiden' column in obs")
            return []

        plt.close("all")
        try:
            sc.pl.violin(
                adata,
                keys=score_key,
                groupby="leiden",
                show=False,
                rotation=45,
            )
        except Exception as exc:
            log.warning("pdui_distribution_scanpy: sc.pl.violin failed: %s", exc)
            return []

        fig = plt.gcf()
        fig.suptitle(f"PDUI distribution per cluster ({score_key})", y=1.01)
        paths = save_matplotlib(fig, output_basepath)
        n_clusters = int(adata.obs["leiden"].nunique()) if "leiden" in adata.obs else 0
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "cluster_key": "leiden",
            "score_key": score_key,
            "n_clusters": n_clusters,
            "n_observations": int(adata.n_obs),
        })
        return paths
