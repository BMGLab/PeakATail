"""Atlas-snap diagnostic diagram — matplotlib backend.

Two-panel figure:
  Left  — bar chart of snapped vs. unsnapped peak counts.
  Right — histogram of snap distances (bp) for the snapped peaks.

Data contract::

    data: dict with keys:
        "snapped"        : int  — number of peaks that matched within distance
        "unsnapped"      : int  — number of peaks beyond the threshold
        "snap_distances" : list[int]  — distance in bp for each snapped peak
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib
from ema.viz._meta import write_figure_meta

log = logging.getLogger(__name__)


@register_viz_strategy
class AtlasSnapDiagMatplotlib(VizStrategy):
    """Atlas-snap diagnostic — 2-panel matplotlib figure."""

    name = "atlas_snap_diag_matplotlib"
    plot_type = "atlas_snap_diag"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render the 2-panel atlas-snap diagnostic.

        Args:
            data: Dict with keys ``snapped``, ``unsnapped``, ``snap_distances``.
            output_basepath: Base path (no extension); PNG + SVG are written.

        Returns:
            List of written file paths.
        """
        snapped: int = int(data.get("snapped", 0))
        unsnapped: int = int(data.get("unsnapped", 0))
        distances: list[int] = list(data.get("snap_distances", []))

        fig, (ax_bar, ax_hist) = plt.subplots(1, 2, figsize=(10, 4))

        # ---- Left panel: snapped vs unsnapped bar ----
        categories = ["Snapped", "Unsnapped"]
        values = [snapped, unsnapped]
        colors = ["#55a868", "#c44e52"]
        bars = ax_bar.bar(categories, values, color=colors, edgecolor="white", width=0.5)
        for bar, val in zip(bars, values):
            ax_bar.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(values) * 0.01 + 0.5,
                str(val),
                ha="center", va="bottom", fontsize=10, fontweight="bold",
            )
        total = snapped + unsnapped
        if total > 0:
            pct = 100 * snapped / total
            ax_bar.set_title(f"Atlas snap result\n{pct:.1f}% peaks snapped", fontsize=10)
        else:
            ax_bar.set_title("Atlas snap result\n(no peaks)", fontsize=10)
        ax_bar.set_ylabel("# peaks")
        ax_bar.spines["top"].set_visible(False)
        ax_bar.spines["right"].set_visible(False)

        # ---- Right panel: distance histogram ----
        if distances:
            arr = np.array(distances, dtype=int)
            n_bins = min(50, max(10, int(np.sqrt(len(arr)))))
            ax_hist.hist(arr, bins=n_bins, color="#4c72b0", edgecolor="white")
            ax_hist.axvline(
                float(np.median(arr)), color="#c44e52", linestyle="--",
                label=f"median {np.median(arr):.0f} bp",
            )
            ax_hist.legend(fontsize=8, frameon=False)
            ax_hist.set_xlabel("Distance to atlas (bp)")
            ax_hist.set_ylabel("# peaks")
            ax_hist.set_title("Snap-distance distribution", fontsize=10)
        else:
            ax_hist.text(
                0.5, 0.5, "No snapped peaks",
                ha="center", va="center", transform=ax_hist.transAxes,
                fontsize=12, color="gray",
            )
            ax_hist.set_title("Snap-distance distribution", fontsize=10)
        ax_hist.spines["top"].set_visible(False)
        ax_hist.spines["right"].set_visible(False)

        plt.tight_layout()
        paths = save_matplotlib(fig, output_basepath)
        total = snapped + unsnapped
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "n_snapped": snapped,
            "n_unsnapped": unsnapped,
            "snap_rate": round(snapped / total, 4) if total > 0 else None,
            "n_snap_distances": len(distances),
        })
        return paths
