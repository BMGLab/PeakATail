"""3'UTR length shift heatmap — matplotlib backend.

Receives a pd.DataFrame indexed by ``gene_id`` with columns for each
cluster pair (e.g. ``"0_vs_1"``, ``"0_vs_2"``).  Each cell is a PDUI
delta (positive = lengthening, negative = shortening).
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

log = logging.getLogger(__name__)

_MAX_GENES = 50  # cap rows to keep the figure readable


@register_viz_strategy
class LengthShiftsMatplotlib(VizStrategy):
    """Heatmap of PDUI delta across genes × cluster pairs (matplotlib).

    Data shape:
        pd.DataFrame indexed by ``gene_id`` with one column per cluster pair.
        Values are PDUI delta (float, typically in [-1, 1]).
    """

    name = "length_shifts_matplotlib"
    plot_type = "length_shifts"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render length-shift heatmap.

        Args:
            data: DataFrame (genes × cluster-pairs) of PDUI deltas.
            output_basepath: Path stem.

        Returns:
            List of paths written (PNG + SVG).
        """
        df: pd.DataFrame = data
        if df is None or df.empty:
            log.warning("length_shifts_matplotlib: empty DataFrame, skipping")
            return []

        # Subset to the top genes by absolute magnitude
        if len(df) > _MAX_GENES:
            top_idx = df.abs().max(axis=1).nlargest(_MAX_GENES).index
            df = df.loc[top_idx]

        n_genes, n_pairs = df.shape
        cell_h = max(0.25, min(0.6, 12 / n_genes))
        cell_w = max(0.8, min(2.0, 10 / n_pairs))
        fig_h = max(4, n_genes * cell_h + 2)
        fig_w = max(4, n_pairs * cell_w + 3)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        vmax = float(df.abs().values.max() or 1.0)
        im = ax.imshow(df.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        fig.colorbar(im, ax=ax, label="ΔPDUI")

        ax.set_xticks(range(n_pairs))
        ax.set_xticklabels(list(df.columns), rotation=45, ha="right")
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(list(df.index), fontsize=max(4, min(9, 200 // n_genes)))
        ax.set_title("3'UTR length shifts (ΔPDUI per gene per cluster pair)")

        plt.tight_layout()
        return save_matplotlib(fig, output_basepath)
