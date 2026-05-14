"""Gene-track view — matplotlib backend.

Renders a stacked-subplot figure for one gene showing:

  1. Gene structure track (exon bars per isoform, when isoform data is present)
  2. Per-cluster PAS coverage (reads/cell normalised) as vertical bar charts

Each cluster row colourises bars by within-gene proportion using the viridis
sequential colourmap.  Bars are annotated with the proportion percent (skipped
for NaN rows).

Data shape:
    Either a :class:`~ema.viz._gene_track_helpers.GenePanel` directly, or a
    dict ``{"panel": GenePanel, "top_n_clusters": int}`` to cap how many
    clusters are rendered.
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

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib
from ema.viz._meta import write_figure_meta

log = logging.getLogger(__name__)

# Maximum clusters rendered even when panel has more.
_MAX_CLUSTERS = 12
# Height (inches) per isoform row in the gene structure track.
_ISOFORM_ROW_HEIGHT = 0.35
# Height (inches) per cluster coverage row.
_CLUSTER_ROW_HEIGHT = 1.0
# Minimum total figure height (inches).
_MIN_FIG_HEIGHT = 2.5


@register_viz_strategy
class GeneTrackMatplotlib(VizStrategy):
    """Gene-track: PAS positions, per-cluster reads/cell, within-gene proportions.

    Data shape:
        A :class:`~ema.viz._gene_track_helpers.GenePanel` instance, or a dict
        with key ``"panel"`` (GenePanel) and optional ``"top_n_clusters"``
        (int, default 12) capping how many clusters are rendered.
    """

    name = "gene_track_matplotlib"
    plot_type = "gene_track"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render the gene-track figure.

        Args:
            data: GenePanel or dict with ``"panel"`` + optional ``"top_n_clusters"``.
            output_basepath: Path stem (PNG + SVG will be written next to it).

        Returns:
            List of paths written (PNG + SVG).
        """
        from ema.viz._gene_track_helpers import GenePanel

        # --- unpack data ---
        top_n_clusters = _MAX_CLUSTERS
        if isinstance(data, dict):
            panel: GenePanel = data["panel"]
            top_n_clusters = int(data.get("top_n_clusters", _MAX_CLUSTERS))
        elif isinstance(data, GenePanel):
            panel = data
        else:
            log.warning(
                "gene_track_matplotlib: unsupported data type %s; skipping",
                type(data).__name__,
            )
            return []

        if not panel.pas_ids:
            log.warning(
                "gene_track_matplotlib: gene %r has no PAS; skipping",
                panel.gene_id,
            )
            return []

        # --- determine cluster subset to render ---
        cluster_indices = _select_top_clusters(panel, top_n_clusters)
        n_clusters_rendered = len(cluster_indices)

        if n_clusters_rendered == 0:
            log.warning(
                "gene_track_matplotlib: gene %r has no clusters with reads; "
                "skipping",
                panel.gene_id,
            )
            return []

        # --- compute layout ---
        has_structure = bool(panel.isoforms)
        n_isoforms = len(panel.isoforms) if has_structure else 0

        # Shared x range with 5% padding on each side.
        gene_span = max(panel.end - panel.start, 1)
        pad = int(gene_span * 0.05) + 1
        x_min = panel.start - pad
        x_max = panel.end + pad

        # Minimum visible bar width = 0.5% of gene span (so 1-bp PAS don't
        # disappear when the gene span is large).  Actual PAS widths from
        # ``pas_starts``/``pas_ends`` are used otherwise so a 200-bp merged
        # PAS visibly differs from a 50-bp singleton.
        min_visible_w = max(int(gene_span * 0.005), 1)
        # Per-PAS widths, with sensible fallbacks if the helper didn't
        # populate pas_starts/pas_ends (legacy panels).
        if panel.pas_starts and panel.pas_ends and len(panel.pas_starts) == len(panel.pas_positions):
            pas_widths = [
                max(int(e - s), min_visible_w)
                for s, e in zip(panel.pas_starts, panel.pas_ends)
            ]
            pas_left_edges = list(panel.pas_starts)
        else:
            pas_widths = [min_visible_w] * len(panel.pas_positions)
            pas_left_edges = [int(p) - min_visible_w // 2 for p in panel.pas_positions]

        # Subplot heights in order: [gene_structure (optional), cluster0, cluster1, ...]
        subplot_heights: list[float] = []
        if has_structure:
            subplot_heights.append(_ISOFORM_ROW_HEIGHT * n_isoforms)
        subplot_heights.extend([_CLUSTER_ROW_HEIGHT] * n_clusters_rendered)

        total_height = max(sum(subplot_heights), _MIN_FIG_HEIGHT)
        fig_width = max(9, min(14, gene_span / 1000 + 8))

        # Create figure with height_ratios so the gene structure track is compact.
        fig, axes = plt.subplots(
            nrows=len(subplot_heights),
            ncols=1,
            figsize=(fig_width, total_height),
            gridspec_kw={"height_ratios": subplot_heights, "hspace": 0.35},
            sharex=True,
        )
        if len(subplot_heights) == 1:
            axes = [axes]  # type: ignore[list-item]

        ax_idx = 0

        # --- 1. gene structure track ---
        if has_structure:
            ax_struct = axes[ax_idx]
            ax_idx += 1
            _draw_gene_structure(ax_struct, panel, x_min, x_max, n_isoforms)
            _annotate_pas_positions(
                ax_struct, panel, pas_left_edges, pas_widths, n_isoforms,
            )

        # --- 2. per-cluster coverage rows ---
        # Cap the reads_per_cell y-axis at 95th-pctile × 1.1 across all rendered
        # clusters so an outlier PAS doesn't visually crush the other bars.
        all_rpc_values = []
        for ci in cluster_indices:
            all_rpc_values.extend(panel.reads_per_cell[ci].tolist())
        finite_rpc = [v for v in all_rpc_values if np.isfinite(v) and v > 0]
        if finite_rpc:
            y_cap = float(np.percentile(finite_rpc, 95)) * 1.1
        else:
            y_cap = 1.0
        y_cap = max(y_cap, 1e-6)

        cmap = plt.cm.viridis  # type: ignore[attr-defined]

        for rank, ci in enumerate(cluster_indices):
            ax = axes[ax_idx + rank]
            cluster_label = panel.clusters[ci]
            n_cells = panel.n_cells_per_cluster[ci]
            rpc_row = panel.reads_per_cell[ci]          # shape: (n_pas,)
            prop_row = panel.proportions[ci]             # shape: (n_pas,)

            for j, (pos, left_edge, width, rpc, prop) in enumerate(
                zip(
                    panel.pas_positions,
                    pas_left_edges,
                    pas_widths,
                    rpc_row,
                    prop_row,
                )
            ):
                # colour by within-gene proportion (NaN → grey)
                if np.isfinite(prop):
                    colour = cmap(float(np.clip(prop, 0.0, 1.0)))
                else:
                    colour = "#aaaaaa"

                bar_height = float(rpc) if np.isfinite(rpc) and rpc > 0 else 0.0
                # Draw bar at the PAS's actual genomic span (left edge +
                # width) so wide merged PAS show as wide bars.
                ax.bar(
                    left_edge,
                    bar_height,
                    width=width,
                    align="edge",
                    color=colour,
                    linewidth=0,
                    zorder=2,
                )

                # Proportion annotation centred above the bar (skip NaN).
                if np.isfinite(prop) and bar_height > 0:
                    pct_str = f"{prop * 100:.0f}%"
                    ax.annotate(
                        pct_str,
                        xy=(pos, min(bar_height, y_cap)),
                        xytext=(0, 2),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=5,
                        color="#333333",
                        clip_on=True,
                    )

            # Y-axis styling.
            ax.set_ylim(0, y_cap)
            ax.set_ylabel(
                f"cluster {cluster_label}\n(n={n_cells})",
                fontsize=7,
                labelpad=4,
            )
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(axis="y", labelsize=6)
            ax.set_xlim(x_min, x_max)
            ax.yaxis.set_major_formatter(
                matplotlib.ticker.FormatStrFormatter("%.2g")  # type: ignore[attr-defined]
            )

        # --- shared x-axis ---
        bottom_ax = axes[-1]
        bottom_ax.set_xlabel("Genomic position (bp)", fontsize=9)
        bottom_ax.tick_params(axis="x", labelsize=7)

        # --- figure title ---
        n_pas = len(panel.pas_ids)
        title = (
            f"Gene {panel.gene_id} — "
            f"{panel.chrom}:{panel.start}-{panel.end} "
            f"({'−' if panel.strand == '-' else '+'} strand) — "
            f"{n_pas} PAS"
        )
        fig.suptitle(title, fontsize=9, y=1.01)

        # --- colorbar legend (proportion) ---
        sm = plt.cm.ScalarMappable(  # type: ignore[attr-defined]
            cmap=cmap,
            norm=matplotlib.colors.Normalize(vmin=0, vmax=1),  # type: ignore[attr-defined]
        )
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, fraction=0.015, pad=0.04, shrink=0.5)
        cbar.set_label("Within-gene proportion", fontsize=7)
        cbar.ax.tick_params(labelsize=6)

        paths = save_matplotlib(fig, output_basepath)

        # --- meta sidecar ---
        top_prop_per_cluster: dict[str, float] = {}
        for ci in cluster_indices:
            prop_row = panel.proportions[ci]
            finite = prop_row[np.isfinite(prop_row)]
            top_prop_per_cluster[panel.clusters[ci]] = (
                float(finite.max()) if len(finite) > 0 else float("nan")
            )

        write_figure_meta(
            output_basepath,
            {
                "viz_strategy": self.name,
                "description": (
                    "Gene track: PAS positions, per-cluster reads/cell, "
                    "within-gene proportions."
                ),
                "gene_id": panel.gene_id,
                "chrom": panel.chrom,
                "start": panel.start,
                "end": panel.end,
                "strand": panel.strand,
                "n_pas": n_pas,
                "n_clusters_rendered": n_clusters_rendered,
                "n_isoforms": n_isoforms,
                "top_proportion_per_cluster": top_prop_per_cluster,
            },
        )
        return paths


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _select_top_clusters(panel: Any, top_n: int) -> list[int]:
    """Return indices into panel.clusters for up to top_n non-empty clusters.

    A cluster is "non-empty" for this gene if it has at least one read across
    any PAS.  When there are more non-empty clusters than ``top_n`` the top
    ones by total reads at the gene are kept.

    Args:
        panel: GenePanel.
        top_n: Maximum clusters to render.

    Returns:
        Ordered list of integer indices into ``panel.clusters``.
    """
    n_clusters = len(panel.clusters)
    # reads row sums per cluster.
    row_totals = panel.reads.sum(axis=1)  # shape: (n_clusters,)
    non_empty = [i for i in range(n_clusters) if row_totals[i] > 0]
    if len(non_empty) <= top_n:
        return non_empty
    # Sort non-empty by descending total reads and take top_n.
    non_empty_sorted = sorted(non_empty, key=lambda i: row_totals[i], reverse=True)
    return non_empty_sorted[:top_n]


def _draw_gene_structure(
    ax: "plt.Axes",  # type: ignore[name-defined]
    panel: Any,
    x_min: float,
    x_max: float,
    n_isoforms: int,
) -> None:
    """Draw exon bars and intron lines for each isoform.

    Args:
        ax: Matplotlib Axes to draw into.
        panel: GenePanel with isoforms field.
        x_min, x_max: X display range.
        n_isoforms: Number of isoforms (used for y positioning).
    """
    exon_colour = "#4477aa"
    intron_colour = "#888888"
    exon_height = 0.6  # fraction of the row height

    for row_idx, (tid, exons) in enumerate(panel.isoforms):
        y_centre = row_idx + 0.5
        if not exons:
            continue
        # Draw intron backbone (thin line spanning all exons).
        all_starts = [e[0] for e in exons]
        all_ends = [e[1] for e in exons]
        backbone_start = min(all_starts)
        backbone_end = max(all_ends)
        ax.plot(
            [backbone_start, backbone_end],
            [y_centre, y_centre],
            color=intron_colour,
            lw=0.8,
            zorder=1,
        )
        # Draw strand arrows along the backbone.
        arrow_step = max(int((backbone_end - backbone_start) / 8), 1)
        arrow_dir = 1 if panel.strand == "+" else -1
        for pos in range(backbone_start, backbone_end, arrow_step):
            ax.annotate(
                "",
                xy=(pos + arrow_dir * arrow_step * 0.3, y_centre),
                xytext=(pos, y_centre),
                arrowprops=dict(
                    arrowstyle="->",
                    color=intron_colour,
                    lw=0.5,
                ),
            )

        # Draw exon rectangles.
        for exon_start, exon_end in exons:
            rect = mpatches.FancyBboxPatch(
                (exon_start, y_centre - exon_height / 2),
                exon_end - exon_start,
                exon_height,
                boxstyle="square,pad=0",
                facecolor=exon_colour,
                edgecolor="none",
                zorder=3,
            )
            ax.add_patch(rect)
        # Label with transcript id (abbreviated).
        short_tid = tid if len(tid) <= 18 else tid[:15] + "..."
        ax.text(
            x_min + (x_max - x_min) * 0.005,
            y_centre,
            short_tid,
            va="center",
            ha="left",
            fontsize=5,
            color="#333333",
        )

    ax.set_ylim(0, n_isoforms)
    ax.set_yticks([])
    ax.set_xlim(x_min, x_max)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.set_ylabel("Isoforms", fontsize=7, labelpad=4)
    ax.tick_params(axis="x", bottom=False)


def _annotate_pas_positions(
    ax: "plt.Axes",  # type: ignore[name-defined]
    panel: Any,
    pas_left_edges: list[int],
    pas_widths: list[int],
    n_isoforms: int,
) -> None:
    """Overlay each PAS region as a coloured stripe + pas_id tick on the structure axis.

    Makes the actual genomic position (and width, post-merger) of every PAS
    visible against the isoform structure --- the user can immediately see
    which exon / 3'UTR each PAS sits in, and whether two PAS shown side-by-
    side in the coverage rows are genuinely close on the genome or just
    appear close because the rendered bars are clipped to a small axis.

    Args:
        ax: The gene-structure subplot axes.
        panel: GenePanel with ``pas_ids`` aligned to ``pas_left_edges``.
        pas_left_edges: Genomic left-edge of each PAS region.
        pas_widths: Width (bp) of each PAS region.
        n_isoforms: Isoform row count (used to size annotations).
    """
    pas_colour = "#C44E52"     # warm red — high-contrast against blue exons
    label_y = n_isoforms + 0.05

    for pas_id, left, width in zip(panel.pas_ids, pas_left_edges, pas_widths):
        # Shaded stripe spanning the PAS's actual genomic extent across
        # every isoform row.  Alpha is intentionally low so the underlying
        # exon/intron structure remains visible.
        ax.axvspan(
            left, left + width,
            ymin=0.0, ymax=1.0,
            color=pas_colour, alpha=0.18, lw=0,
            zorder=4,
        )
        # Tick mark + pas_id label above the structure track.
        centre = left + width / 2.0
        ax.text(
            centre, label_y,
            f"{pas_id}",
            ha="center", va="bottom",
            fontsize=5, color=pas_colour,
            rotation=0, clip_on=False,
        )
