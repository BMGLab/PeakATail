"""Per-PAS proportion heatmap (PAS x cluster) — matplotlib backend.

The proportion strategy produces a long-format DataFrame::

    gene_id  transcript_id  pas_id  rank  cell  proportion

This viz collapses cells -> clusters by averaging, picks the top-N most
variable PAS (by cluster-level variance), and renders a row-normalised
PAS x cluster heatmap so a researcher can see which PAS shift between
cell populations.

Data shape:
    ``(pdui_df, adata)`` tuple — the long DataFrame plus the AnnData
    that carries ``obs["leiden"]`` (or whichever cluster column).
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


def _build_pas_cluster_matrix(
    pdui_df: pd.DataFrame,
    adata,
    cluster_key: str,
    top_n: int,
) -> tuple[pd.DataFrame, list[str]] | tuple[None, None]:
    """Return (pas x cluster mean proportion table, sorted cluster labels).

    Returns ``(None, None)`` if the inputs lack the columns we need.
    """
    if "proportion" not in pdui_df.columns or "cell" not in pdui_df.columns:
        return None, None
    if "pas_id" not in pdui_df.columns:
        return None, None
    if cluster_key not in adata.obs.columns:
        return None, None

    cell_to_cluster = adata.obs[cluster_key].astype(str).to_dict()
    df = pdui_df.dropna(subset=["proportion"]).copy()
    df["cluster"] = df["cell"].map(cell_to_cluster)
    df = df.dropna(subset=["cluster"])
    if df.empty:
        return None, None

    # Mean proportion per (pas_id, cluster)
    pivot = (
        df.groupby(["pas_id", "cluster"], observed=True)["proportion"]
          .mean()
          .unstack("cluster")
    )
    # Pick the top-N PAS by cross-cluster variance (most discriminating).
    variances = pivot.var(axis=1, skipna=True).fillna(0.0)
    top_idx = variances.sort_values(ascending=False).head(top_n).index
    pivot = pivot.loc[top_idx]
    cluster_order = sorted(
        pivot.columns.tolist(),
        key=lambda x: int(x) if str(x).isdigit() else x,
    )
    pivot = pivot[cluster_order]
    return pivot, cluster_order


@register_viz_strategy
class ProportionHeatmapMatplotlib(VizStrategy):
    """PAS x cluster mean-proportion heatmap (matplotlib).

    Tunable knobs are passed via a data dict::

        render_all("proportion_heatmap",
                   {"pdui_df": df, "adata": adata,
                    "cluster_key": "leiden", "top_n": 50}, ...)

    Falls back to ``(pdui_df, adata)`` tuple form with sensible defaults.
    """

    name = "proportion_heatmap_matplotlib"
    plot_type = "proportion_heatmap"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        if isinstance(data, dict):
            pdui_df = data["pdui_df"]
            adata = data["adata"]
            cluster_key = data.get("cluster_key", "leiden")
            top_n = int(data.get("top_n", 50))
        else:
            pdui_df, adata = data
            cluster_key, top_n = "leiden", 50

        pivot, clusters = _build_pas_cluster_matrix(
            pdui_df, adata, cluster_key, top_n,
        )
        if pivot is None or pivot.empty:
            log.warning(
                "proportion_heatmap_matplotlib: no rows to plot (cluster_key=%r)",
                cluster_key,
            )
            return []

        # Display values in [0, 1]; clip strictly to avoid color flicker on NaN.
        mat = pivot.fillna(0.0).values
        fig_h = max(5, 0.18 * pivot.shape[0])
        fig_w = max(5, 0.8 * pivot.shape[1] + 2)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)

        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels(clusters)
        ax.set_xlabel("cluster")
        # Only annotate row labels when readable (<= 50 rows).
        if pivot.shape[0] <= 50:
            ax.set_yticks(range(pivot.shape[0]))
            ax.set_yticklabels([str(p) for p in pivot.index])
        else:
            ax.set_yticks([])
        ax.set_ylabel(f"PAS (top-{pivot.shape[0]} by variance)")
        ax.set_title("Mean PAS proportion per cluster")
        fig.colorbar(im, ax=ax, label="mean proportion")

        paths = save_matplotlib(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "cluster_key": cluster_key,
            "n_pas_shown": int(pivot.shape[0]),
            "n_clusters": int(pivot.shape[1]),
            "top_n_selected_by": "cross-cluster variance",
        })
        return paths
