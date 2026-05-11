"""Volcano plot — matplotlib backend.

Receives a pd.DataFrame with columns ["log2fc", "qvalue", "pas_id"].
Significant PAS (qvalue < threshold) are highlighted in a contrasting colour.

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

# Module-level fallback defaults (used when data dict does not supply overrides).
_FDR_DEFAULT: float = 0.05
_LOG2FC_THRESH_DEFAULT: float = 1.0


@register_viz_strategy
class VolcanoMatplotlib(VizStrategy):
    """Volcano plot rendered with matplotlib (paper-quality static PNG + SVG).

    Data shape:
        Either a ``pd.DataFrame`` with columns ``["log2fc", "qvalue", "pas_id"]``
        (legacy), or a dict with key ``"df"`` carrying that DataFrame plus
        optional keys ``"fdr"`` (float) and ``"log2fc_thresh"`` (float) to
        override the significance thresholds drawn on the figure.

        When ``gene_id`` is present as a column in the DataFrame the top-N
        most significant points are annotated with their gene symbol.  N is
        controlled by the ``n_label`` key in the data dict (default 10).

    Rows = PAS tested in one cluster pair.
    """

    name = "volcano_matplotlib"
    plot_type = "volcano"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render volcano plot.

        Args:
            data: Either a DataFrame with ``log2fc``, ``qvalue``, ``pas_id``
                columns (legacy), or a dict with key ``"df"`` plus optional
                ``"fdr"``, ``"log2fc_thresh"``, and ``"n_label"`` overrides.
            output_basepath: Path stem (suffix will be added by save_matplotlib).

        Returns:
            List of paths written (PNG + SVG).
        """
        # --- unpack data dict or bare DataFrame ---
        if isinstance(data, dict):
            df: pd.DataFrame = data["df"]
            fdr: float = float(data.get("fdr", _FDR_DEFAULT))
            log2fc_thresh: float = float(data.get("log2fc_thresh", _LOG2FC_THRESH_DEFAULT))
            n_label: int = int(data.get("n_label", 10))
        else:
            df = data
            fdr = _FDR_DEFAULT
            log2fc_thresh = _LOG2FC_THRESH_DEFAULT
            n_label = 10

        if df is None or df.empty:
            log.warning("volcano_matplotlib: empty DataFrame, skipping")
            return []

        neg_log10_q = -np.log10(df["qvalue"].clip(lower=1e-300))

        sig_mask = (df["qvalue"] < fdr) & (df["log2fc"].abs() >= log2fc_thresh)
        up = sig_mask & (df["log2fc"] >= log2fc_thresh)
        down = sig_mask & (df["log2fc"] <= -log2fc_thresh)
        ns = ~sig_mask

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(
            df.loc[ns, "log2fc"], neg_log10_q[ns],
            s=10, color="#aaaaaa", alpha=0.5, linewidths=0, label="NS",
        )
        ax.scatter(
            df.loc[up, "log2fc"], neg_log10_q[up],
            s=14, color="#c44e52", alpha=0.8, linewidths=0, label=f"Up (n={up.sum()})",
        )
        ax.scatter(
            df.loc[down, "log2fc"], neg_log10_q[down],
            s=14, color="#4c72b0", alpha=0.8, linewidths=0, label=f"Down (n={down.sum()})",
        )
        # Reference lines using the actual user-supplied thresholds.
        ax.axhline(-np.log10(fdr), color="grey", lw=0.8, ls="--")
        ax.axvline(log2fc_thresh, color="grey", lw=0.8, ls="--")
        ax.axvline(-log2fc_thresh, color="grey", lw=0.8, ls="--")

        ax.set_xlabel("log₂ fold change")
        ax.set_ylabel("-log₁₀(q-value)")
        ax.set_title("Differential APA — volcano")
        ax.legend(loc="upper left", fontsize=8, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # --- top-N gene labels (Bug 4) ---
        if n_label > 0 and "gene_id" in df.columns:
            # Score = -log10(q) * |log2fc|; rank descending.
            score = neg_log10_q * df["log2fc"].abs()
            top_idx = score.nlargest(n_label).index
            for idx in top_idx:
                gene = df.loc[idx, "gene_id"]
                if pd.isna(gene) or str(gene).strip() == "":
                    continue
                x_val = float(df.loc[idx, "log2fc"])
                y_val = float(neg_log10_q[idx])
                ax.annotate(
                    str(gene),
                    xy=(x_val, y_val),
                    xytext=(4, 2),
                    textcoords="offset points",
                    fontsize=6,
                    color="#333333",
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
                    ha="left",
                    va="bottom",
                    clip_on=True,
                )

        paths = save_matplotlib(fig, output_basepath)
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "fdr": fdr,
            "log2fc_thresh": log2fc_thresh,
            "n_tested": len(df),
            "n_significant": int(sig_mask.sum()),
            "n_up": int(up.sum()),
            "n_down": int(down.sum()),
        })
        return paths
