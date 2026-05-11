"""Match confidence heatmap — plotly backend.

Interactive heatmap with hover details: cluster source, canonical target,
and match confidence value.

Data shape: same DataFrame as cluster_match_sankey (from ClusterMatchStrategy.match()).
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


@register_viz_strategy
class MatchConfidencePlotly(VizStrategy):
    """Interactive heatmap of cluster match confidence scores (plotly).

    Data shape:
        pd.DataFrame with columns:
        ``dataset_id``, ``original_cluster``, ``canonical_cluster``,
        ``match_confidence``.
    """

    name = "match_confidence_plotly"
    plot_type = "match_confidence"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render interactive match confidence heatmap.

        Args:
            data: DataFrame from ClusterMatchStrategy.match().
            output_basepath: Path stem.

        Returns:
            List of paths written (HTML + SVG).
        """
        df: pd.DataFrame = data
        if df is None or df.empty:
            log.warning("match_confidence_plotly: empty DataFrame, skipping")
            return []

        required = {"dataset_id", "original_cluster", "canonical_cluster", "match_confidence"}
        missing = required - set(df.columns)
        if missing:
            log.warning(
                "match_confidence_plotly: missing columns %s, skipping", missing
            )
            return []

        df = df.copy()
        df["row_label"] = df["dataset_id"].astype(str) + ":" + df["original_cluster"].astype(str)

        canonical_ids = sorted(df["canonical_cluster"].unique())
        row_labels = list(df["row_label"].unique())
        n_rows = len(row_labels)
        n_cols = len(canonical_ids)

        matrix = np.zeros((n_rows, n_cols))
        hover_text = [[" " for _ in range(n_cols)] for _ in range(n_rows)]

        row_idx = {r: i for i, r in enumerate(row_labels)}
        col_idx = {c: j for j, c in enumerate(canonical_ids)}

        for _, row in df.iterrows():
            i = row_idx[row["row_label"]]
            j = col_idx[row["canonical_cluster"]]
            conf = float(row["match_confidence"])
            matrix[i, j] = conf
            hover_text[i][j] = (
                f"Source: {row['row_label']}<br>"
                f"Canonical: {row['canonical_cluster']}<br>"
                f"Confidence: {conf:.4f}"
            )

        col_labels = [f"can_{c}" for c in canonical_ids]

        fig = go.Figure(
            go.Heatmap(
                z=matrix.tolist(),
                x=col_labels,
                y=row_labels,
                colorscale="YlOrRd",
                zmin=0,
                zmax=1,
                text=hover_text,
                hovertemplate="%{text}<extra></extra>",
                colorbar=dict(title="confidence"),
            )
        )
        fig.update_layout(
            template="plotly_white",
            title="Cluster match confidence",
            xaxis=dict(title="Canonical cluster", tickangle=45),
            yaxis=dict(title="Dataset:cluster", autorange="reversed"),
            width=max(450, n_cols * 80 + 250),
            height=max(350, n_rows * 28 + 150),
        )

        return save_plotly(fig, output_basepath)
