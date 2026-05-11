"""4-panel grid: peaks/chrom, cells/PAS dist, peak width dist, reads/cell dist."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib


@register_viz_strategy
class PeakQCMatplotlib(VizStrategy):
    name = "peak_qc_matplotlib"
    plot_type = "peak_qc"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        ppc = data.get("peaks_per_chrom", {})
        per_cell_pas = data.get("per_cell_pas", [])
        peak_widths = data.get("peak_widths", [])
        per_cell_reads = data.get("per_cell_reads", [])
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        # peaks per chrom
        chroms = list(ppc.keys())
        axes[0, 0].bar(chroms, [ppc[c] for c in chroms], color="#4c72b0")
        axes[0, 0].set_title("Peaks per chromosome")
        axes[0, 0].tick_params(axis="x", rotation=45)
        # per-cell PAS dist
        if len(per_cell_pas) > 0:
            axes[0, 1].hist(per_cell_pas, bins=40, color="#dd8452", edgecolor="white")
            axes[0, 1].set_title("PAS detected per cell")
            axes[0, 1].set_xlabel("# PAS")
        # peak widths
        if len(peak_widths) > 0:
            axes[1, 0].hist(peak_widths, bins=40, color="#55a868", edgecolor="white")
            axes[1, 0].set_title("Peak widths (bp)")
            axes[1, 0].set_xlabel("width")
        # per-cell reads
        if len(per_cell_reads) > 0:
            axes[1, 1].hist(per_cell_reads, bins=40, color="#c44e52", edgecolor="white")
            axes[1, 1].set_title("Reads per cell")
            axes[1, 1].set_xlabel("# reads")
        for ax in axes.flat:
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        plt.tight_layout()
        return save_matplotlib(fig, output_basepath)
