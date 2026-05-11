"""Volcano plot — plotly backend with hover gene names.

Receives a pd.DataFrame with columns ["log2fc", "qvalue", "pas_id"].
Interactive HTML + SVG export via kaleido.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly

log = logging.getLogger(__name__)

_FDR = 0.05
_LOG2FC_THRESH = 1.0


@register_viz_strategy
class VolcanoPlotly(VizStrategy):
    """Volcano plot rendered with plotly (interactive HTML + SVG).

    Data shape:
        pd.DataFrame with columns ``["log2fc", "qvalue", "pas_id"]``.
        Rows = PAS tested in one cluster pair.
    """

    name = "volcano_plotly"
    plot_type = "volcano"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render interactive volcano plot.

        Args:
            data: DataFrame with ``log2fc``, ``qvalue``, ``pas_id`` columns.
            output_basepath: Path stem (suffix will be added by save_plotly).

        Returns:
            List of paths written (HTML + SVG).
        """
        df: pd.DataFrame = data
        if df is None or df.empty:
            log.warning("volcano_plotly: empty DataFrame, skipping")
            return []

        neg_log10_q = -np.log10(df["qvalue"].clip(lower=1e-300))
        sig_mask = (df["qvalue"] < _FDR) & (df["log2fc"].abs() >= _LOG2FC_THRESH)

        def _color(row_sig: bool, row_lfc: float) -> str:
            if not row_sig:
                return "#aaaaaa"
            return "#c44e52" if row_lfc > 0 else "#4c72b0"

        colors = [
            _color(sig_mask.iloc[i], df["log2fc"].iloc[i])
            for i in range(len(df))
        ]

        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=df["log2fc"].tolist(),
                y=neg_log10_q.tolist(),
                mode="markers",
                marker=dict(color=colors, size=6, opacity=0.75),
                text=df["pas_id"].tolist() if "pas_id" in df.columns else None,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "log₂FC: %{x:.3f}<br>"
                    "-log₁₀(q): %{y:.2f}<extra></extra>"
                ),
                name="PAS",
            )
        )
        # Reference lines
        fig.add_hline(y=-np.log10(_FDR), line_dash="dash", line_color="grey", line_width=1)
        fig.add_vline(x=_LOG2FC_THRESH, line_dash="dash", line_color="grey", line_width=1)
        fig.add_vline(x=-_LOG2FC_THRESH, line_dash="dash", line_color="grey", line_width=1)

        n_up = int((sig_mask & (df["log2fc"] >= _LOG2FC_THRESH)).sum())
        n_down = int((sig_mask & (df["log2fc"] <= -_LOG2FC_THRESH)).sum())
        fig.update_layout(
            template="plotly_white",
            title=f"Differential APA — volcano (up={n_up}, down={n_down})",
            xaxis_title="log₂ fold change",
            yaxis_title="-log₁₀(q-value)",
            showlegend=False,
            width=700,
            height=580,
        )

        return save_plotly(fig, output_basepath)
