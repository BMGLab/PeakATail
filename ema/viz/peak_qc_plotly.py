from __future__ import annotations

from pathlib import Path
from typing import Any

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly


@register_viz_strategy
class PeakQCPlotly(VizStrategy):
    name = "peak_qc_plotly"
    plot_type = "peak_qc"
    engine = "plotly"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        ppc = data.get("peaks_per_chrom", {})
        per_cell_pas = data.get("per_cell_pas", [])
        peak_widths = data.get("peak_widths", [])
        per_cell_reads = data.get("per_cell_reads", [])
        fig = make_subplots(rows=2, cols=2, subplot_titles=(
            "Peaks per chromosome", "PAS per cell",
            "Peak widths", "Reads per cell"))
        if ppc:
            fig.add_trace(go.Bar(x=list(ppc.keys()), y=list(ppc.values())), row=1, col=1)
        if len(per_cell_pas) > 0:
            fig.add_trace(go.Histogram(x=list(per_cell_pas), nbinsx=40), row=1, col=2)
        if len(peak_widths) > 0:
            fig.add_trace(go.Histogram(x=list(peak_widths), nbinsx=40), row=2, col=1)
        if len(per_cell_reads) > 0:
            fig.add_trace(go.Histogram(x=list(per_cell_reads), nbinsx=40), row=2, col=2)
        fig.update_layout(template="plotly_white", showlegend=False, height=700)
        return save_plotly(fig, output_basepath)
