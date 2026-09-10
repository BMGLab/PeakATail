"""Cluster match Sankey / chord diagram — matplotlib backend.

Since matplotlib does not have a native Sankey for categorical flows,
this implementation renders a stacked bar "alluvial"-style chart that
visually shows how original clusters in each dataset map to canonical
clusters.

Data shape: pd.DataFrame from ``ClusterMatchStrategy.match()`` with columns
``["dataset_id", "original_cluster", "canonical_cluster", "match_confidence",
"matched_to"]``.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib
from ema.viz._meta import write_figure_meta

log = logging.getLogger(__name__)


@register_viz_strategy
class ClusterMatchSankeyMatplotlib(VizStrategy):
    """Alluvial-style cluster correspondence chart (matplotlib chord fallback).

    Data shape:
        pd.DataFrame with columns:
        ``dataset_id``, ``original_cluster``, ``canonical_cluster``,
        ``match_confidence``, ``matched_to``.
    """

    name = "cluster_match_sankey_matplotlib"
    plot_type = "cluster_match_sankey"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render alluvial cluster-match chart.

        Args:
            data: DataFrame from ClusterMatchStrategy.match().
            output_basepath: Path stem.

        Returns:
            List of paths written (PNG + SVG).
        """
        df: pd.DataFrame = data
        if df is None or df.empty:
            log.warning("cluster_match_sankey_matplotlib: empty DataFrame, skipping")
            return []

        required = {"dataset_id", "original_cluster", "canonical_cluster", "match_confidence"}
        missing = required - set(df.columns)
        if missing:
            log.warning(
                "cluster_match_sankey_matplotlib: missing columns %s, skipping", missing
            )
            return []

        datasets = sorted(df["dataset_id"].unique())
        canonical_ids = sorted(df["canonical_cluster"].unique())
        n_datasets = len(datasets)
        n_canonical = len(canonical_ids)

        # Build a colour palette per canonical cluster
        cmap = plt.cm.tab20(np.linspace(0, 1, max(n_canonical, 1)))
        can_color = {cid: cmap[i] for i, cid in enumerate(canonical_ids)}

        fig, axes = plt.subplots(
            1, n_datasets,
            figsize=(max(6, n_datasets * 3), max(4, n_canonical * 0.8 + 2)),
            sharey=False,
        )
        if n_datasets == 1:
            axes = [axes]

        for ax, ds in zip(axes, datasets):
            sub = df[df["dataset_id"] == ds].sort_values("canonical_cluster")
            bottoms = np.zeros(1)
            for _, row in sub.iterrows():
                color = can_color.get(row["canonical_cluster"], "#cccccc")
                ax.bar(
                    0,
                    1,
                    bottom=bottoms[0],
                    color=color,
                    alpha=float(row["match_confidence"]) * 0.6 + 0.4,
                    edgecolor="white",
                    linewidth=0.8,
                    width=0.6,
                )
                ax.text(
                    0,
                    bottoms[0] + 0.5,
                    f"c{row['original_cluster']}→{row['canonical_cluster']}",
                    ha="center", va="center", fontsize=7,
                )
                bottoms[0] += 1

            ax.set_title(ds, fontsize=9)
            ax.set_xlim(-0.5, 0.5)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_visible(False)

        # Legend
        patches = [
            mpatches.Patch(color=can_color[cid], label=f"canonical {cid}")
            for cid in canonical_ids
        ]
        fig.legend(
            handles=patches,
            loc="lower center",
            ncol=min(6, n_canonical),
            fontsize=7,
            frameon=False,
        )
        fig.suptitle("Cross-dataset cluster correspondence", y=1.02)
        plt.tight_layout()

        paths = save_matplotlib(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "n_datasets": n_datasets,
            "dataset_ids": datasets,
            "n_canonical_clusters": n_canonical,
            "n_match_rows": len(df),
        })
        return paths
