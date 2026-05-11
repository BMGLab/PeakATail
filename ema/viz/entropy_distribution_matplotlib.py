"""Shannon entropy distribution per cluster — matplotlib backend.

Mirrors :mod:`ema.viz.pdui_distribution_matplotlib` but visualises the
Shannon-entropy score column produced by
:class:`~ema.quantification.strategies.shannon.ShannonPDUIStrategy`.

Data shape: ``(adata, score_key)`` where ``adata.obs[score_key]`` holds
per-cell mean entropy and ``adata.obs["leiden"]`` holds cluster labels.
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
class EntropyDistributionMatplotlib(VizStrategy):
    """Per-cluster Shannon-entropy violin plot (matplotlib)."""

    name = "entropy_distribution_matplotlib"
    plot_type = "entropy_distribution"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        adata, score_key = data
        if score_key not in adata.obs.columns:
            log.warning(
                "entropy_distribution_matplotlib: score key %r not in obs columns",
                score_key,
            )
            return []

        scores: pd.Series = adata.obs[score_key].astype(float)
        clusters = adata.obs.get("leiden")
        if clusters is None:
            log.warning("entropy_distribution_matplotlib: no 'leiden' column in obs")
            return []

        labels = sorted(
            clusters.unique(), key=lambda x: int(x) if str(x).isdigit() else x
        )
        groups = [scores[clusters == lbl].dropna().values for lbl in labels]
        valid = [(lbl, grp) for lbl, grp in zip(labels, groups) if len(grp) > 0]
        if not valid:
            log.warning("entropy_distribution_matplotlib: no valid groups")
            return []
        labels, groups = zip(*valid)  # type: ignore[assignment]

        fig, ax = plt.subplots(figsize=(max(5, len(labels) * 0.9), 5))
        parts = ax.violinplot(
            list(groups),
            positions=range(len(labels)),
            showmedians=True,
            showextrema=True,
        )
        cmap = plt.cm.viridis(np.linspace(0, 1, len(labels)))
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(cmap[i])
            body.set_alpha(0.75)

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_xlabel("cluster")
        ax.set_ylabel(score_key)
        ax.set_title(f"Shannon entropy per cluster ({score_key})")
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
