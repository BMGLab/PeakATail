"""Differential agreement heatmap — plotly backend.

Renders an interactive pairwise Jaccard overlap heatmap across multiple
differential APA strategies.

Receives ``dict[strategy_name, set[sig_pas_ids]]``.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import plotly.graph_objects as go

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly
from ema.viz._meta import write_figure_meta

log = logging.getLogger(__name__)


@register_viz_strategy
class DiffAgreementPlotly(VizStrategy):
    """Interactive pairwise overlap heatmap across differential strategies (plotly).

    Data shape:
        ``dict[str, set[str]]`` — mapping strategy name → set of significant
        PAS IDs declared significant by that strategy.
    """

    name = "diff_agreement_plotly"
    plot_type = "diff_agreement"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render interactive Jaccard overlap heatmap.

        Args:
            data: ``dict[strategy_name, set[sig_pas_ids]]``.
            output_basepath: Path stem.

        Returns:
            List of paths written (HTML + SVG), empty if <2 strategies.
        """
        if not isinstance(data, dict) or len(data) < 2:
            log.info("diff_agreement_plotly: need >= 2 strategies, skipping")
            return []

        names = sorted(data.keys())
        n = len(names)
        matrix = np.zeros((n, n), dtype=float)
        text = [[" " for _ in range(n)] for _ in range(n)]

        for i, a in enumerate(names):
            for j, b in enumerate(names):
                sa, sb = data[a], data[b]
                union = len(sa | sb)
                val = len(sa & sb) / union if union > 0 else 0.0
                matrix[i, j] = val
                intersection = len(sa & sb)
                text[i][j] = (
                    f"Jaccard: {val:.3f}<br>"
                    f"Intersection: {intersection}<br>"
                    f"|{a}|={len(sa)}, |{b}|={len(sb)}"
                )

        fig = go.Figure(
            go.Heatmap(
                z=matrix.tolist(),
                x=names,
                y=names,
                colorscale="Blues",
                zmin=0,
                zmax=1,
                text=text,
                hovertemplate="%{text}<extra></extra>",
                colorbar=dict(title="Jaccard"),
            )
        )
        fig.update_layout(
            template="plotly_white",
            title="Strategy agreement (Jaccard, significant PAS)",
            width=max(400, n * 100 + 150),
            height=max(400, n * 100 + 150),
            xaxis=dict(title="Strategy"),
            yaxis=dict(title="Strategy", autorange="reversed"),
        )

        paths = save_plotly(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "sig_set_names": names,
            "sig_set_sizes": {name: len(data[name]) for name in names},
        })
        return paths
