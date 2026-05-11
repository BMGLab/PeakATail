"""Cluster match Sankey diagram — plotly backend.

Renders a proper Sankey diagram showing how original per-dataset clusters
flow into canonical shared clusters.

Data shape: pd.DataFrame from ``ClusterMatchStrategy.match()`` with columns
``["dataset_id", "original_cluster", "canonical_cluster", "match_confidence",
"matched_to"]``.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import plotly.graph_objects as go
import pandas as pd

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly

log = logging.getLogger(__name__)


def _hex_to_rgba(hex_color: str, alpha: float = 0.7) -> str:
    """Convert hex color to rgba string for plotly."""
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha:.2f})"


_PALETTE = [
    "#4c72b0", "#dd8452", "#55a868", "#c44e52", "#8172b2",
    "#937860", "#da8bc3", "#8c8c8c", "#ccb974", "#64b5cd",
]


@register_viz_strategy
class ClusterMatchSankeyPlotly(VizStrategy):
    """Proper Sankey diagram of cluster correspondences (plotly).

    Data shape:
        pd.DataFrame with columns:
        ``dataset_id``, ``original_cluster``, ``canonical_cluster``,
        ``match_confidence``, ``matched_to``.
    """

    name = "cluster_match_sankey_plotly"
    plot_type = "cluster_match_sankey"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render Sankey flow diagram for cluster matches.

        Each source node is ``{dataset_id}:cluster_{original_cluster}`` and
        each target node is ``canonical_{canonical_cluster}``.  Link width is
        proportional to match_confidence.

        Args:
            data: DataFrame from ClusterMatchStrategy.match().
            output_basepath: Path stem.

        Returns:
            List of paths written (HTML + SVG).
        """
        df: pd.DataFrame = data
        if df is None or df.empty:
            log.warning("cluster_match_sankey_plotly: empty DataFrame, skipping")
            return []

        required = {"dataset_id", "original_cluster", "canonical_cluster", "match_confidence"}
        missing = required - set(df.columns)
        if missing:
            log.warning(
                "cluster_match_sankey_plotly: missing columns %s, skipping", missing
            )
            return []

        # Build node index
        source_nodes = [
            f"{row['dataset_id']}:c{row['original_cluster']}"
            for _, row in df.iterrows()
        ]
        target_nodes = [
            f"canonical_{row['canonical_cluster']}"
            for _, row in df.iterrows()
        ]
        all_nodes = list(dict.fromkeys(source_nodes + target_nodes))
        node_idx = {n: i for i, n in enumerate(all_nodes)}

        canonical_ids = sorted(df["canonical_cluster"].unique())
        can_color = {
            cid: _PALETTE[i % len(_PALETTE)]
            for i, cid in enumerate(canonical_ids)
        }

        sources, targets, values, link_colors = [], [], [], []
        for (_, row), src_label, tgt_label in zip(
            df.iterrows(), source_nodes, target_nodes
        ):
            sources.append(node_idx[src_label])
            targets.append(node_idx[tgt_label])
            values.append(max(0.01, float(row["match_confidence"])))
            hex_c = can_color.get(row["canonical_cluster"], "#aaaaaa")
            link_colors.append(_hex_to_rgba(hex_c, alpha=0.55))

        node_colors = []
        for node in all_nodes:
            if node.startswith("canonical_"):
                cid_str = node.replace("canonical_", "")
                try:
                    cid = type(list(canonical_ids)[0])(cid_str)
                except Exception:
                    cid = cid_str
                node_colors.append(can_color.get(cid, "#aaaaaa"))
            else:
                node_colors.append("#b0b0b0")

        fig = go.Figure(
            go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=20,
                    thickness=20,
                    line=dict(color="white", width=0.5),
                    label=all_nodes,
                    color=node_colors,
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=link_colors,
                    hovertemplate=(
                        "From %{source.label}<br>"
                        "To %{target.label}<br>"
                        "Confidence: %{value:.3f}<extra></extra>"
                    ),
                ),
            )
        )
        fig.update_layout(
            title="Cross-dataset cluster correspondence (Sankey)",
            template="plotly_white",
            width=900,
            height=max(400, len(all_nodes) * 20 + 200),
            font=dict(size=11),
        )

        return save_plotly(fig, output_basepath)
