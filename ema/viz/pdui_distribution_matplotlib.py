"""PDUI distribution per cluster — matplotlib backend.

Renders per-cluster violin/box plots of a PDUI score column from an AnnData.

Data shape: ``(adata, score_key)`` where ``adata.obs[score_key]`` contains
per-cell PDUI scores (float in [0, 1]) and ``adata.obs["leiden"]`` contains
cluster assignments.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib
from ema.viz._meta import write_figure_meta

log = logging.getLogger(__name__)


@register_viz_strategy
class PduiDistributionMatplotlib(VizStrategy):
    """Per-cluster PDUI violin plot (matplotlib).

    Data shape:
        ``(adata, score_key)`` tuple.
        ``adata`` — AnnData with ``obs["leiden"]`` and ``obs[score_key]``.
        ``score_key`` — column name in ``adata.obs`` holding PDUI scores.
    """

    name = "pdui_distribution_matplotlib"
    plot_type = "pdui_distribution"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render per-cluster PDUI violin plot.

        Args:
            data: ``(adata, score_key)`` tuple.
            output_basepath: Path stem.

        Returns:
            List of paths written (PNG + SVG).
        """
        adata, score_key = data
        if score_key not in adata.obs.columns:
            log.warning(
                "pdui_distribution_matplotlib: score key %r not in obs columns", score_key
            )
            return []

        scores: pd.Series = adata.obs[score_key].astype(float)
        clusters = adata.obs.get("leiden")
        if clusters is None:
            log.warning("pdui_distribution_matplotlib: no 'leiden' column in obs")
            return []

        labels = sorted(clusters.unique(), key=lambda x: int(x) if str(x).isdigit() else x)
        groups = [scores[clusters == lbl].values for lbl in labels]
        # Filter out empty groups
        valid = [(lbl, grp) for lbl, grp in zip(labels, groups) if len(grp) > 0]
        if not valid:
            log.warning("pdui_distribution_matplotlib: no valid groups")
            return []
        labels, groups = zip(*valid)  # type: ignore[assignment]

        fig, ax = plt.subplots(figsize=(max(5, len(labels) * 0.9), 5))
        parts = ax.violinplot(
            list(groups),
            positions=range(len(labels)),
            showmedians=True,
            showextrema=True,
        )
        # Colour the violin bodies
        cmap = plt.cm.tab10(np.linspace(0, 1, len(labels)))
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(cmap[i])
            body.set_alpha(0.75)

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_xlabel("cluster")
        ax.set_ylabel(score_key)
        ax.set_title(f"PDUI distribution per cluster ({score_key})")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        paths = save_matplotlib(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "cluster_key": "leiden",
            "score_key": score_key,
            "n_clusters": len(labels),
            "n_observations": int(adata.n_obs),
        })
        return paths
