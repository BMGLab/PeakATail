"""Gene track visualisation — plotly interactive backend.

Renders, for one gene:

  1. **Gene structure track** (optional) — one horizontal segment trace per
     isoform when ``panel.isoforms`` is non-empty.
  2. **Per-cluster PAS coverage bars** — one subplot row per cluster; bar
     height = reads / cell, bar colour = within-gene proportion (Viridis).
  3. **Shared genomic x-axis** labelled with chromosome and strand.

All subplots share the same x-axis so zooming in the browser on one panel
immediately zooms all others — important when comparing cluster-specific usage
across the same PAS set.

Registration
------------
Decorated with ``@register_viz_strategy``.  Auto-discovery in
``ema/viz/__init__.py`` picks this up without any extra hookup.

Data shapes
-----------
Accepts either a :class:`~ema.viz._gene_track_helpers.GenePanel` directly, or
a ``dict`` with the following structure::

    {
        "panel": GenePanel,
        "top_n_clusters": int,  # default 12 — cap on rendered cluster rows
    }
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_plotly
from ema.viz._meta import write_figure_meta
from ema.viz._gene_track_helpers import GenePanel

log = logging.getLogger(__name__)

# Maximum clusters rendered by default (overridable via data dict).
_DEFAULT_TOP_N_CLUSTERS: int = 12

# Distinct palette for isoform tracks — qualitative, visually separable.
_ISOFORM_COLORS: list[str] = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
]


@register_viz_strategy
class GeneTrackPlotly(VizStrategy):
    """Interactive gene track: isoforms + per-cluster PAS coverage (plotly).

    Produces a vertically stacked multi-panel figure where all subplots share
    the genomic x-axis.  Hovering any bar shows PAS id, genomic position,
    cluster, cell count, raw reads, reads/cell, and within-gene proportion.
    """

    name = "gene_track_plotly"
    plot_type = "gene_track"
    engine = "plotly"

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render interactive gene track figure.

        Args:
            data: A :class:`~ema.viz._gene_track_helpers.GenePanel` or a dict
                with keys ``"panel"`` (GenePanel) and optional
                ``"top_n_clusters"`` (int, default 12).
            output_basepath: Stem path; ``.html`` / ``.svg`` / ``.meta.json``
                suffixes are appended by the IO helpers.

        Returns:
            List of paths written (HTML, SVG, meta JSON).
        """
        panel, top_n = self._unpack(data)
        if panel is None:
            log.warning("gene_track_plotly: received None panel, skipping")
            return []

        cluster_indices, capped = self._select_clusters(panel, top_n)
        n_clusters_rendered = len(cluster_indices)
        n_isoforms = len(panel.isoforms)
        has_isoforms = n_isoforms > 0

        if n_clusters_rendered == 0:
            log.warning(
                "gene_track_plotly: gene %s has no clusters with reads, skipping",
                panel.gene_id,
            )
            return []

        fig = self._build_figure(panel, cluster_indices, has_isoforms)
        paths = save_plotly(fig, output_basepath)

        # Compute top_proportion_per_cluster for meta sidecar.
        top_prop: dict[str, float] = {}
        for i in cluster_indices:
            row = panel.proportions[i, :]
            finite = row[np.isfinite(row)]
            top_prop[panel.clusters[i]] = float(finite.max()) if finite.size > 0 else 0.0

        meta: dict[str, Any] = {
            "viz_strategy": self.name,
            "description": (
                "Gene track (interactive): PAS positions, per-cluster reads/cell, "
                "within-gene proportions. Hover for exact counts."
            ),
            "gene_id": panel.gene_id,
            "chrom": panel.chrom,
            "start": panel.start,
            "end": panel.end,
            "strand": panel.strand,
            "n_pas": len(panel.pas_ids),
            "n_clusters_rendered": n_clusters_rendered,
            "n_isoforms": n_isoforms,
            "top_proportion_per_cluster": top_prop,
        }
        if capped:
            meta["cluster_cap_applied"] = True
            meta["cluster_cap_limit"] = top_n
            meta["total_clusters_available"] = len(panel.clusters)

        write_figure_meta(output_basepath, meta)
        return paths

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _unpack(data: Any) -> tuple[GenePanel | None, int]:
        """Return (panel, top_n_clusters) from raw data argument."""
        if isinstance(data, dict):
            panel: GenePanel | None = data.get("panel")
            top_n: int = int(data.get("top_n_clusters", _DEFAULT_TOP_N_CLUSTERS))
        elif isinstance(data, GenePanel):
            panel = data
            top_n = _DEFAULT_TOP_N_CLUSTERS
        else:
            log.error(
                "gene_track_plotly: unexpected data type %s", type(data).__name__
            )
            return None, _DEFAULT_TOP_N_CLUSTERS
        return panel, top_n

    @staticmethod
    def _select_clusters(
        panel: GenePanel, top_n: int
    ) -> tuple[list[int], bool]:
        """Pick up to *top_n* cluster indices, ranked by total reads for the gene.

        Returns:
            (selected_indices_in_display_order, was_capped)
        """
        # Total reads per cluster for this gene.
        total_reads = panel.reads.sum(axis=1)  # shape (n_clusters,)
        ranked = sorted(
            range(len(panel.clusters)),
            key=lambda i: total_reads[i],
            reverse=True,
        )
        capped = len(ranked) > top_n
        selected = sorted(ranked[:top_n])  # restore display order
        return selected, capped

    def _build_figure(
        self,
        panel: GenePanel,
        cluster_indices: list[int],
        has_isoforms: bool,
    ) -> go.Figure:
        """Assemble the multi-panel plotly Figure."""
        n_isoforms = len(panel.isoforms)
        n_rows = len(cluster_indices) + (1 if has_isoforms else 0)

        # Row heights: isoform track shorter than coverage bars.
        if has_isoforms:
            isoform_weight = max(0.5, n_isoforms * 0.15)
            coverage_weight = 1.0
            total_w = isoform_weight + coverage_weight * len(cluster_indices)
            row_heights = (
                [isoform_weight / total_w]
                + [coverage_weight / total_w] * len(cluster_indices)
            )
        else:
            row_heights = [1.0 / len(cluster_indices)] * len(cluster_indices)

        subplot_titles = self._build_subplot_titles(panel, cluster_indices, has_isoforms)

        fig = make_subplots(
            rows=n_rows,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.03,
            subplot_titles=subplot_titles,
            row_heights=row_heights,
        )

        current_row = 1

        # --- Panel 1: gene structure track (optional) ---
        if has_isoforms:
            self._add_isoform_traces(fig, panel, row=current_row)
            current_row += 1

        # --- Panels 2..N: per-cluster PAS coverage ---
        for list_pos, ci in enumerate(cluster_indices):
            self._add_cluster_bar_traces(fig, panel, cluster_idx=ci, row=current_row)
            current_row += 1

        # --- Layout ---
        height = max(
            400,
            70 * len(cluster_indices) + 60 * n_isoforms + 80,
        )
        fig.update_layout(
            template="plotly_white",
            title=dict(
                text=(
                    f"Gene {panel.gene_id} — "
                    f"{panel.chrom}:{panel.start}-{panel.end} "
                    f"({panel.strand} strand) — "
                    f"{len(panel.pas_ids)} PAS"
                ),
                font=dict(size=13),
            ),
            height=height,
            width=900,
            margin=dict(l=80, r=30, t=70, b=50),
            showlegend=False,
            # Colorbar shared reference axis for the coverage bars.
            coloraxis=dict(
                colorscale="Viridis",
                cmin=0.0,
                cmax=1.0,
                colorbar=dict(
                    title="proportion",
                    thickness=12,
                    len=0.4,
                    y=0.2,
                ),
            ),
        )

        # Label the bottom-most x-axis.
        bottom_xaxis = f"xaxis{n_rows}" if n_rows > 1 else "xaxis"
        fig.update_layout(
            **{
                bottom_xaxis: dict(
                    title=dict(
                        text=f"genomic position ({panel.chrom}, {panel.strand} strand)",
                        font=dict(size=11),
                    )
                )
            }
        )

        return fig

    @staticmethod
    def _build_subplot_titles(
        panel: GenePanel,
        cluster_indices: list[int],
        has_isoforms: bool,
    ) -> list[str]:
        titles: list[str] = []
        if has_isoforms:
            titles.append("gene structure")
        for ci in cluster_indices:
            n_cells = panel.n_cells_per_cluster[ci]
            titles.append(f"cluster {panel.clusters[ci]} ({n_cells} cells)")
        return titles

    @staticmethod
    def _add_isoform_traces(fig: go.Figure, panel: GenePanel, row: int) -> None:
        """Add one horizontal-segment trace per isoform to *row*."""
        for iso_idx, (transcript_id, exons) in enumerate(panel.isoforms):
            color = _ISOFORM_COLORS[iso_idx % len(_ISOFORM_COLORS)]
            # Draw exon blocks as filled rectangles via scatter with fill.
            for exon_start, exon_end in exons:
                y_val = iso_idx
                fig.add_trace(
                    go.Scatter(
                        x=[exon_start, exon_end, exon_end, exon_start, exon_start],
                        y=[y_val - 0.3, y_val - 0.3, y_val + 0.3, y_val + 0.3, y_val - 0.3],
                        mode="lines",
                        fill="toself",
                        fillcolor=color,
                        line=dict(color=color, width=0),
                        hoverinfo="text",
                        hovertext=f"{transcript_id}: {exon_start}-{exon_end}",
                        showlegend=False,
                        name=transcript_id,
                    ),
                    row=row,
                    col=1,
                )
            # Intron backbone — thin horizontal line spanning the full isoform range.
            if exons:
                iso_min = min(s for s, _ in exons)
                iso_max = max(e for _, e in exons)
                fig.add_trace(
                    go.Scatter(
                        x=[iso_min, iso_max],
                        y=[iso_idx, iso_idx],
                        mode="lines",
                        line=dict(color=color, width=1, dash="dot"),
                        hoverinfo="text",
                        hovertext=transcript_id,
                        showlegend=False,
                        name=f"{transcript_id}_intron",
                    ),
                    row=row,
                    col=1,
                )
        # Annotate transcript IDs on the left margin.
        for iso_idx, (transcript_id, _) in enumerate(panel.isoforms):
            fig.add_annotation(
                xref="paper", x=0,
                yref=f"y{row}",
                y=iso_idx,
                text=f"<b>{transcript_id[:20]}</b>",
                showarrow=False,
                xanchor="right",
                font=dict(size=8),
                row=row, col=1,
            )

    @staticmethod
    def _add_cluster_bar_traces(
        fig: go.Figure,
        panel: GenePanel,
        cluster_idx: int,
        row: int,
    ) -> None:
        """Add a bar trace for cluster *cluster_idx* to *row*.

        Bar height = reads_per_cell; bar colour mapped to proportion via Viridis.
        """
        ci = cluster_idx
        n_pas = len(panel.pas_ids)
        n_cells = panel.n_cells_per_cluster[ci]
        cluster_label = panel.clusters[ci]

        x_vals: list[int] = panel.pas_positions
        y_vals: list[float] = panel.reads_per_cell[ci, :].tolist()
        proportions_row = panel.proportions[ci, :]

        # Build custom_data columns aligned with each bar.
        custom_data: list[list[Any]] = []
        for j in range(n_pas):
            prop = float(proportions_row[j]) if np.isfinite(proportions_row[j]) else float("nan")
            custom_data.append([
                panel.pas_ids[j],          # 0: pas_id
                panel.pas_positions[j],    # 1: position
                cluster_label,             # 2: cluster
                n_cells,                   # 3: n_cells
                float(panel.reads[ci, j]), # 4: reads
                float(panel.reads_per_cell[ci, j]),  # 5: reads_per_cell
                prop,                      # 6: proportion
            ])

        # Colour each bar individually by proportion (mapped onto Viridis).
        # We pass marker.color as a numeric array and reference coloraxis.
        prop_vals: list[float] = [
            float(proportions_row[j]) if np.isfinite(proportions_row[j]) else 0.0
            for j in range(n_pas)
        ]

        fig.add_trace(
            go.Bar(
                x=x_vals,
                y=y_vals,
                customdata=custom_data,
                hovertemplate=(
                    "<b>PAS %{customdata[0]}</b><br>"
                    "position: %{customdata[1]}<br>"
                    "cluster: %{customdata[2]}<br>"
                    "n_cells: %{customdata[3]}<br>"
                    "reads: %{customdata[4]:.1f}<br>"
                    "reads / cell: %{customdata[5]:.4f}<br>"
                    "proportion within gene: %{customdata[6]:.1%}"
                    "<extra></extra>"
                ),
                marker=dict(
                    color=prop_vals,
                    coloraxis="coloraxis",
                    line=dict(width=0),
                ),
                name=f"cluster {cluster_label}",
                showlegend=False,
            ),
            row=row,
            col=1,
        )
