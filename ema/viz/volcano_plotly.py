"""Volcano plot — plotly backend with hover gene names.

Receives a pd.DataFrame with columns ["log2fc", "qvalue", "pas_id"].
Interactive HTML + SVG export via kaleido.

FDR threading
-------------
The ``render`` method accepts ``fdr`` and ``log2fc_thresh`` via the ``data``
dict (keys ``"fdr"`` and ``"log2fc_thresh"``).  This is the least-invasive
approach: the call site in ``pipeline_hooks.py::render_switch_diff_outputs``
only needs to include those keys when calling ``render_all``.  The module-level
defaults are preserved as fallbacks so existing calls with a bare DataFrame
continue to work.

TODO (pipeline_hooks.py owner — Bug 2 wire-up):
    In ``render_switch_diff_outputs``, change::

        render_all("volcano", df, ...)

    to::

        render_all("volcano", {"df": df, "fdr": fdr, "log2fc_thresh": log2fc_thresh}, ...)

    where ``fdr`` and ``log2fc_thresh`` come from the user's CLI options.
    Until that change lands, the strategy falls back to 0.05 / 1.0.
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

# Module-level fallback defaults (used when data dict does not supply overrides).
_FDR_DEFAULT: float = 0.05
_LOG2FC_THRESH_DEFAULT: float = 1.0


@register_viz_strategy
class VolcanoPlotly(VizStrategy):
    """Volcano plot rendered with plotly (interactive HTML + SVG).

    Data shape:
        Either a ``pd.DataFrame`` with columns ``["log2fc", "qvalue", "pas_id"]``
        (legacy), or a dict with key ``"df"`` carrying that DataFrame plus
        optional keys ``"fdr"`` (float) and ``"log2fc_thresh"`` (float) to
        override the significance thresholds drawn on the figure.

        When ``gene_id`` is present as a column in the DataFrame, hover text
        shows ``gene_id : pas_id``; otherwise raw ``pas_id`` is shown.

    Rows = PAS tested in one cluster pair.
    """

    name = "volcano_plotly"
    plot_type = "volcano"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render interactive volcano plot.

        Args:
            data: Either a DataFrame with ``log2fc``, ``qvalue``, ``pas_id``
                columns (legacy), or a dict with key ``"df"`` plus optional
                ``"fdr"`` and ``"log2fc_thresh"`` overrides.
            output_basepath: Path stem (suffix will be added by save_plotly).

        Returns:
            List of paths written (HTML + SVG).
        """
        # --- unpack data dict or bare DataFrame ---
        if isinstance(data, dict):
            df: pd.DataFrame = data["df"]
            fdr: float = float(data.get("fdr", _FDR_DEFAULT))
            log2fc_thresh: float = float(data.get("log2fc_thresh", _LOG2FC_THRESH_DEFAULT))
        else:
            df = data
            fdr = _FDR_DEFAULT
            log2fc_thresh = _LOG2FC_THRESH_DEFAULT

        if df is None or df.empty:
            log.warning("volcano_plotly: empty DataFrame, skipping")
            return []

        neg_log10_q = -np.log10(df["qvalue"].clip(lower=1e-300))
        sig_mask = (df["qvalue"] < fdr) & (df["log2fc"].abs() >= log2fc_thresh)

        def _color(row_sig: bool, row_lfc: float) -> str:
            if not row_sig:
                return "#aaaaaa"
            return "#c44e52" if row_lfc > 0 else "#4c72b0"

        colors = [
            _color(sig_mask.iloc[i], df["log2fc"].iloc[i])
            for i in range(len(df))
        ]

        # Build hover text: prefer "gene_id : pas_id" when gene_id is available.
        has_gene = "gene_id" in df.columns
        hover_texts: list[str] = []
        for i in range(len(df)):
            pas = str(df["pas_id"].iloc[i]) if "pas_id" in df.columns else str(df.index[i])
            if has_gene:
                gene = df["gene_id"].iloc[i]
                label = f"{gene} : {pas}" if pd.notna(gene) and str(gene).strip() else pas
            else:
                label = pas
            hover_texts.append(label)

        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=df["log2fc"].tolist(),
                y=neg_log10_q.tolist(),
                mode="markers",
                marker=dict(color=colors, size=6, opacity=0.75),
                text=hover_texts,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "log₂FC: %{x:.3f}<br>"
                    "-log₁₀(q): %{y:.2f}<extra></extra>"
                ),
                name="PAS",
            )
        )
        # Reference lines using the actual user-supplied thresholds.
        fig.add_hline(y=-np.log10(fdr), line_dash="dash", line_color="grey", line_width=1)
        fig.add_vline(x=log2fc_thresh, line_dash="dash", line_color="grey", line_width=1)
        fig.add_vline(x=-log2fc_thresh, line_dash="dash", line_color="grey", line_width=1)

        n_up = int((sig_mask & (df["log2fc"] >= log2fc_thresh)).sum())
        n_down = int((sig_mask & (df["log2fc"] <= -log2fc_thresh)).sum())
        fig.update_layout(
            template="plotly_white",
            title=f"Differential APA — volcano (up={n_up}, down={n_down})",
            xaxis_title="log₂ fold change",
            yaxis_title="-log₁₀(q-value)",
            showlegend=False,
            width=700,
            height=580,
        )

        paths = save_plotly(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "fdr": fdr,
            "log2fc_thresh": log2fc_thresh,
            "n_tested": len(df),
            "n_significant": int(sig_mask.sum()),
            "n_up": n_up,
            "n_down": n_down,
        })
        return paths
