"""Volcano plot — matplotlib backend.

Receives a pd.DataFrame with columns ["log2fc", "qvalue", "pas_id"].
Significant PAS (qvalue < threshold) are highlighted in a contrasting colour.
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

_FDR = 0.05
_LOG2FC_THRESH = 1.0


@register_viz_strategy
class VolcanoMatplotlib(VizStrategy):
    """Volcano plot rendered with matplotlib (paper-quality static PNG + SVG).

    Data shape:
        pd.DataFrame with columns ``["log2fc", "qvalue", "pas_id"]``.
        Rows = PAS tested in one cluster pair.
    """

    name = "volcano_matplotlib"
    plot_type = "volcano"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render volcano plot.

        Args:
            data: DataFrame with ``log2fc``, ``qvalue``, ``pas_id`` columns.
            output_basepath: Path stem (suffix will be added by save_matplotlib).

        Returns:
            List of paths written (PNG + SVG).
        """
        df: pd.DataFrame = data
        if df is None or df.empty:
            log.warning("volcano_matplotlib: empty DataFrame, skipping")
            return []

        neg_log10_q = -np.log10(df["qvalue"].clip(lower=1e-300))

        sig_mask = (df["qvalue"] < _FDR) & (df["log2fc"].abs() >= _LOG2FC_THRESH)
        up = sig_mask & (df["log2fc"] >= _LOG2FC_THRESH)
        down = sig_mask & (df["log2fc"] <= -_LOG2FC_THRESH)
        ns = ~sig_mask

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(
            df.loc[ns, "log2fc"], neg_log10_q[ns],
            s=10, color="#aaaaaa", alpha=0.5, linewidths=0, label="NS",
        )
        ax.scatter(
            df.loc[up, "log2fc"], neg_log10_q[up],
            s=14, color="#c44e52", alpha=0.8, linewidths=0, label=f"Up (n={up.sum()})",
        )
        ax.scatter(
            df.loc[down, "log2fc"], neg_log10_q[down],
            s=14, color="#4c72b0", alpha=0.8, linewidths=0, label=f"Down (n={down.sum()})",
        )
        # Reference lines
        ax.axhline(-np.log10(_FDR), color="grey", lw=0.8, ls="--")
        ax.axvline(_LOG2FC_THRESH, color="grey", lw=0.8, ls="--")
        ax.axvline(-_LOG2FC_THRESH, color="grey", lw=0.8, ls="--")

        ax.set_xlabel("log₂ fold change")
        ax.set_ylabel("-log₁₀(q-value)")
        ax.set_title("Differential APA — volcano")
        ax.legend(loc="upper left", fontsize=8, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        return save_matplotlib(fig, output_basepath)
