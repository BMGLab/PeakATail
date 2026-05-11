"""Match confidence heatmap — matplotlib backend.

Renders a heatmap of match_confidence values as a matrix:
rows = ``(dataset_id, original_cluster)`` combined label,
columns = canonical cluster IDs.

Data shape: same DataFrame as cluster_match_sankey (from ClusterMatchStrategy.match()).
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
class MatchConfidenceMatplotlib(VizStrategy):
    """Heatmap of cluster match confidence scores (matplotlib).

    Data shape:
        pd.DataFrame with columns:
        ``dataset_id``, ``original_cluster``, ``canonical_cluster``,
        ``match_confidence``.
    """

    name = "match_confidence_matplotlib"
    plot_type = "match_confidence"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render match confidence heatmap.

        Rows = ``{dataset}:{cluster}`` labels.
        Columns = canonical cluster IDs.
        Cell value = match_confidence (0 when no match).

        Args:
            data: DataFrame from ClusterMatchStrategy.match().
            output_basepath: Path stem.

        Returns:
            List of paths written (PNG + SVG).
        """
        df: pd.DataFrame = data
        if df is None or df.empty:
            log.warning("match_confidence_matplotlib: empty DataFrame, skipping")
            return []

        required = {"dataset_id", "original_cluster", "canonical_cluster", "match_confidence"}
        missing = required - set(df.columns)
        if missing:
            log.warning(
                "match_confidence_matplotlib: missing columns %s, skipping", missing
            )
            return []

        df = df.copy()
        df["row_label"] = df["dataset_id"].astype(str) + ":" + df["original_cluster"].astype(str)

        canonical_ids = sorted(df["canonical_cluster"].unique())
        row_labels = list(df["row_label"].unique())
        n_rows = len(row_labels)
        n_cols = len(canonical_ids)

        # Build matrix: confidence for each (row, canonical) pair; 0 if no match
        matrix = np.zeros((n_rows, n_cols))
        row_idx = {r: i for i, r in enumerate(row_labels)}
        col_idx = {c: j for j, c in enumerate(canonical_ids)}
        for _, row in df.iterrows():
            i = row_idx[row["row_label"]]
            j = col_idx[row["canonical_cluster"]]
            matrix[i, j] = float(row["match_confidence"])

        cell_h = max(0.3, min(0.7, 10 / n_rows))
        fig_h = max(4, n_rows * cell_h + 2)
        fig_w = max(4, n_cols * 1.1 + 2)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(matrix, vmin=0, vmax=1, cmap="YlOrRd", aspect="auto")
        fig.colorbar(im, ax=ax, label="match confidence")

        ax.set_xticks(range(n_cols))
        ax.set_xticklabels([f"can_{c}" for c in canonical_ids], rotation=45, ha="right")
        ax.set_yticks(range(n_rows))
        ax.set_yticklabels(
            row_labels,
            fontsize=max(5, min(9, 200 // n_rows)),
        )

        for i in range(n_rows):
            for j in range(n_cols):
                val = matrix[i, j]
                if val > 0:
                    ax.text(
                        j, i, f"{val:.2f}",
                        ha="center", va="center",
                        fontsize=max(5, min(8, 180 // max(n_rows, n_cols))),
                        color="white" if val > 0.7 else "black",
                    )

        ax.set_title("Cluster match confidence")
        plt.tight_layout()
        paths = save_matplotlib(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "n_source_rows": n_rows,
            "n_canonical_clusters": n_cols,
            "datasets": sorted(df["dataset_id"].unique().tolist()),
        })
        return paths
