"""Resource timeline (RAM + CPU over time) — matplotlib backend.

Input ``data``: ``dict`` with two keys::

    {
        "samples": list[dict],       # periodic resource samples
        "annotations": list[dict],   # stage-transition markers (optional)
    }

Each sample record::

    {"elapsed_s": float, "rss_gb": float, "cpu_pct": float}

Each annotation record::

    {"label": str, "elapsed_s": float}

The plot has two y-axes sharing the same x-axis:
- Left y-axis: RSS (GB) — blue line
- Right y-axis: CPU % — orange line

Stage-transition annotations are rendered as vertical dashed grey lines with
small rotated labels.
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
from ema.viz._io import save_matplotlib
from ema.viz.base import VizStrategy

log = logging.getLogger(__name__)

_RAM_COLOR = "#4c72b0"   # blue
_CPU_COLOR = "#dd8452"   # orange
_ANN_COLOR = "#888888"   # grey for annotation lines


@register_viz_strategy
class ResourceTimelineMatplotlib(VizStrategy):
    """Dual-axis line chart: RAM (GB) and CPU% over elapsed time."""

    name = "resource_timeline_matplotlib"
    plot_type = "resource_timeline"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render the resource timeline.

        Args:
            data: Dict with ``"samples"`` and ``"annotations"`` keys.
            output_basepath: Base path (without extension) for output files.

        Returns:
            List of paths written (PNG and SVG).  Returns ``[]`` if samples
            list is empty.
        """
        samples: list[dict[str, Any]] = data.get("samples", [])
        annotations: list[dict[str, Any]] = data.get("annotations", [])

        if not samples:
            log.warning("resource_timeline_matplotlib: empty samples, skipping render")
            return []

        elapsed = np.array([s["elapsed_s"] for s in samples], dtype=float)
        rss_gb = np.array([s["rss_gb"] for s in samples], dtype=float)
        cpu_pct = np.array([s["cpu_pct"] for s in samples], dtype=float)

        fig, ax_ram = plt.subplots(figsize=(12, 5))
        ax_cpu = ax_ram.twinx()

        line_ram = ax_ram.plot(
            elapsed, rss_gb,
            color=_RAM_COLOR, linewidth=1.8, label="RAM (GB)",
        )
        ax_ram.fill_between(elapsed, rss_gb, alpha=0.12, color=_RAM_COLOR)
        ax_ram.set_ylabel("RAM (GB)", color=_RAM_COLOR, fontsize=10)
        ax_ram.tick_params(axis="y", labelcolor=_RAM_COLOR)
        ax_ram.set_ylim(bottom=0)

        line_cpu = ax_cpu.plot(
            elapsed, cpu_pct,
            color=_CPU_COLOR, linewidth=1.5, linestyle="--", label="CPU %",
        )
        ax_cpu.set_ylabel("CPU (%)", color=_CPU_COLOR, fontsize=10)
        ax_cpu.tick_params(axis="y", labelcolor=_CPU_COLOR)
        ax_cpu.set_ylim(0, max(105, cpu_pct.max() * 1.1))

        ax_ram.set_xlabel("Elapsed time (s)", fontsize=10)
        ax_ram.set_title("Resource usage over time")

        # Stage-transition annotations
        for ann in annotations:
            x_val = ann.get("elapsed_s")
            label = ann.get("label", "")
            if x_val is None:
                continue
            ax_ram.axvline(x=x_val, color=_ANN_COLOR, linestyle=":", linewidth=1.0, alpha=0.8)
            ax_ram.text(
                x_val, ax_ram.get_ylim()[1],
                f" {label}",
                rotation=90, va="top", ha="left",
                fontsize=7, color=_ANN_COLOR, alpha=0.9,
            )

        # Combined legend
        lines = line_ram + line_cpu
        labels = [line.get_label() for line in lines]
        ax_ram.legend(lines, labels, loc="upper left", fontsize=9, frameon=False)

        ax_ram.spines["top"].set_visible(False)
        ax_cpu.spines["top"].set_visible(False)

        plt.tight_layout()
        return save_matplotlib(fig, output_basepath)
