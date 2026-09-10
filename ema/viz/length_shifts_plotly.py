"""3'UTR length shift heatmap — plotly backend.

Interactive heatmap with hover details showing gene name, cluster pair,
and PDUI delta value.

Receives a pd.DataFrame indexed by ``gene_id`` with columns per cluster pair.
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
from ema.viz._meta import write_figure_meta

log = logging.getLogger(__name__)

_MAX_GENES = 50


@register_viz_strategy
class LengthShiftsPlotly(VizStrategy):
    """Heatmap of PDUI delta across genes × cluster pairs (plotly interactive).

    Data shape:
        pd.DataFrame indexed by ``gene_id`` with one column per cluster pair.
        Values are PDUI delta (float, typically in [-1, 1]).
    """

    name = "length_shifts_plotly"
    plot_type = "length_shifts"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render interactive length-shift heatmap.

        Args:
            data: DataFrame (genes × cluster-pairs) of PDUI deltas.
            output_basepath: Path stem.

        Returns:
            List of paths written (HTML + SVG).
        """
        df: pd.DataFrame = data
        if df is None or df.empty:
            log.warning("length_shifts_plotly: empty DataFrame, skipping")
            return []

        if len(df) > _MAX_GENES:
            top_idx = df.abs().max(axis=1).nlargest(_MAX_GENES).index
            df = df.loc[top_idx]

        vmax = float(df.abs().values.max() or 1.0)
        pairs = list(df.columns)
        genes = list(df.index)

        # Build custom hover text
        hover = [
            [
                f"Gene: {gene}<br>Pair: {pair}<br>ΔPDUI: {df.loc[gene, pair]:.4f}"
                for pair in pairs
            ]
            for gene in genes
        ]

        fig = go.Figure(
            go.Heatmap(
                z=df.values.tolist(),
                x=pairs,
                y=genes,
                colorscale="RdBu_r",
                zmid=0,
                zmin=-vmax,
                zmax=vmax,
                text=hover,
                hovertemplate="%{text}<extra></extra>",
                colorbar=dict(title="ΔPDUI"),
            )
        )
        fig.update_layout(
            template="plotly_white",
            title="3'UTR length shifts (ΔPDUI per gene per cluster pair)",
            xaxis=dict(title="Cluster pair", tickangle=45),
            yaxis=dict(title="Gene", autorange="reversed"),
            width=max(500, len(pairs) * 100 + 250),
            height=max(400, len(genes) * 14 + 150),
        )

        paths = save_plotly(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "cluster_key": "leiden",
            "n_genes_shown": len(genes),
            "n_cluster_pairs": len(pairs),
            "cluster_pairs": pairs,
        })
        return paths
