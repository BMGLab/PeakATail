"""Differential agreement heatmap — matplotlib backend.

Renders a pairwise overlap heatmap showing how many significant PAS are
shared between each pair of differential strategies.

Receives ``dict[strategy_name, set[sig_pas_ids]]``.
Only renders if more than one strategy is present.
"""
from __future__ import annotations

import logging
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

log = logging.getLogger(__name__)


@register_viz_strategy
class DiffAgreementMatplotlib(VizStrategy):
    """Pairwise overlap heatmap across differential strategies (matplotlib).

    Data shape:
        ``dict[str, set[str]]`` — mapping strategy name → set of significant
        PAS IDs declared significant by that strategy.
    """

    name = "diff_agreement_matplotlib"
    plot_type = "diff_agreement"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render pairwise Jaccard overlap heatmap.

        Args:
            data: ``dict[strategy_name, set[sig_pas_ids]]``.
            output_basepath: Path stem.

        Returns:
            List of paths written (PNG + SVG), empty list if <2 strategies.
        """
        if not isinstance(data, dict) or len(data) < 2:
            log.info("diff_agreement_matplotlib: need >= 2 strategies, skipping")
            return []

        names = sorted(data.keys())
        n = len(names)
        matrix = np.zeros((n, n), dtype=float)

        for i, a in enumerate(names):
            for j, b in enumerate(names):
                sa, sb = data[a], data[b]
                union = len(sa | sb)
                matrix[i, j] = len(sa & sb) / union if union > 0 else 0.0

        fig, ax = plt.subplots(figsize=(max(4, n + 1), max(4, n + 1)))
        im = ax.imshow(matrix, vmin=0, vmax=1, cmap="Blues")
        fig.colorbar(im, ax=ax, label="Jaccard overlap")

        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(names, rotation=45, ha="right")
        ax.set_yticklabels(names)

        for i in range(n):
            for j in range(n):
                ax.text(
                    j, i, f"{matrix[i, j]:.2f}",
                    ha="center", va="center",
                    fontsize=9,
                    color="white" if matrix[i, j] > 0.6 else "black",
                )

        ax.set_title("Strategy agreement (Jaccard, significant PAS)")
        plt.tight_layout()
        paths = save_matplotlib(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "sig_set_names": names,
            "sig_set_sizes": {name: len(data[name]) for name in names},
        })
        return paths
