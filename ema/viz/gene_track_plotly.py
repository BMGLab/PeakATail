"""Gene track visualisation — plotly interactive backend.

This is the INTERACTIVE counterpart to ``gene_track_matplotlib`` and renders
the same GenePanel. Where the matplotlib backend is the publication figure
(static, print-ready), this one is the exploration figure: everything the
static view shows, plus controls a reader can drive in the browser.

Rendered elements (mirroring ``gene_track_matplotlib``)
------------------------------------------------------
  1. **Gene structure track** — one row per isoform: a thin intron backbone
     with strand-direction arrows and the transcript id, overlaid with the
     typed gene model when CDS/UTR annotation is available (tall CDS boxes,
     short 5'UTR boxes, a distinctly-coloured/outlined + labelled 3'UTR —
     the APA-relevant region). Falls back to uniform exon blocks when only
     plain ``exon`` features were loaded.
  2. **Per-cluster PAS tracks** — one subplot row per cluster; bar colour is
     the within-gene proportion (Viridis, shared 0..1 colour axis), matching
     the static figure.
  3. **PAS bands** — a vertical band per PAS spanning every row, with the PAS
     id annotated at the top, so a position lines up across all clusters.
  4. **PAS distance table** — the same ``pas_distance_table`` the matplotlib
     backend draws, as a real plotly table row.
  5. **Title** — ``Gene SYMBOL (ENSG) — chrN:start-end (strand, span) — N PAS``,
     identical in wording to the static figure.

Interactive controls (what this backend adds)
---------------------------------------------
  * **Metric selector** (dropdown) — switch the bar height between
    ``within-gene proportion``, ``reads / cell`` and ``raw reads`` without
    re-rendering. Bar colour stays proportion-mapped so the colour axis keeps
    one meaning across metrics.
  * **Proportion filter** (slider) — dim every PAS whose within-gene
    proportion falls below the threshold, so minor sites drop out of the eye's
    path while the axis/layout stays fixed. Implemented as an opacity mask, so
    it composes with the metric selector rather than fighting it.
  * **Full metadata on hover** — every bar carries its PAS record (id,
    ``chrom:start-end``, width, summit, 5'->3' rank, gap and summit-distance to
    the next PAS, assigned UTR transcripts) *and* the cluster/cell record
    (cluster label, cells in cluster, condition, reads, reads/cell,
    proportion).
  * Shared x-axis: zooming/panning one track zooms all of them.

Registration
------------
Decorated with ``@register_viz_strategy``; auto-discovery in
``ema/viz/__init__.py`` picks this up. ``render_all("gene_track", panel, base,
engines=[...])`` renders this and/or the matplotlib backend from the SAME
panel, which is what lets a UI offer the two as a user-selectable toggle.

Data shapes
-----------
Accepts either a :class:`~ema.viz._gene_track_helpers.GenePanel` directly, or
a ``dict`` with::

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
from ema.viz._gene_track_helpers import (
    GenePanel, TranscriptRegions, format_count, pas_distance_table,
)

log = logging.getLogger(__name__)

# Maximum clusters rendered by default (overridable via data dict).
_DEFAULT_TOP_N_CLUSTERS: int = 12

# Distinct palette for isoform tracks — qualitative, visually separable.
_ISOFORM_COLORS: list[str] = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
]

# PAS band colour — the same salmon wash the static figure uses to tie a PAS
# position across every cluster row.
_PAS_BAND_COLOR = "rgba(214, 39, 40, 0.10)"
_PAS_BAND_LINE = "rgba(214, 39, 40, 0.55)"

# Gene-model region colours (matched to gene_track_matplotlib's palette).
# CDS: tall, neutral slate-blue -- the coding backbone of the transcript.
_CDS_COLOR = "#3a6b91"
# 5'UTR: short, muted -- present but not the focus of an APA view.
_UTR5_COLOR = "#a9c4d8"
# 3'UTR: short but distinctly coloured + outlined -- this is the region APA
# acts on, so it needs to read as "different" at a glance, not just "smaller".
_UTR3_COLOR = "#f2b56b"
_UTR3_OUTLINE = "#c8781f"

# Bar opacity for PAS filtered OUT by the proportion slider. Not zero: a
# filtered PAS is de-emphasised, not hidden — hiding it would silently change
# what the reader thinks the gene contains.
_DIMMED_OPACITY = 0.12

# Thresholds offered by the proportion filter slider.
_PROP_STEPS: list[float] = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50]

# Metrics offered by the dropdown: (label, GenePanel attribute, y-axis title).
_METRICS: list[tuple[str, str, str]] = [
    ("within-gene proportion", "proportions", "proportion"),
    ("reads / cell", "reads_per_cell", "reads / cell"),
    ("raw reads", "reads", "reads"),
]


def _fmt_span(bp: int) -> str:
    """Match the matplotlib backend's span wording exactly."""
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:,.2f} Mb"
    if bp >= 1_000:
        return f"{bp / 1_000:,.1f} kb"
    return f"{bp:,} bp"


@register_viz_strategy
class GeneTrackPlotly(VizStrategy):
    """Interactive gene track: isoforms + per-cluster PAS usage (plotly).

    Vertically stacked multi-panel figure sharing the genomic x-axis, with a
    metric dropdown, a proportion filter slider, PAS bands spanning all rows,
    the PAS distance table, and full PAS/cell metadata on hover.
    """

    name = "gene_track_plotly"
    plot_type = "gene_track"
    engine = "plotly"

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render the interactive gene track figure.

        Args:
            data: A :class:`GenePanel` or a dict with ``"panel"`` and optional
                ``"top_n_clusters"``.
            output_basepath: Stem path; ``.html`` / ``.svg`` / ``.meta.json``
                suffixes are appended by the IO helpers.

        Returns:
            List of paths written.
        """
        panel, top_n = self._unpack(data)
        if panel is None:
            log.warning("gene_track_plotly: received None panel, skipping")
            return []

        cluster_indices, capped = self._select_clusters(panel, top_n)
        if not cluster_indices:
            log.warning(
                "gene_track_plotly: gene %s has no clusters with reads, skipping",
                panel.gene_id,
            )
            return []

        dist_df = self._distance_table(panel)
        fig = self._build_figure(panel, cluster_indices, dist_df)

        # HTML is the artifact that matters for this backend; a missing kaleido
        # (static export) must not cost us the interactive file.
        paths = save_plotly(fig, output_basepath, formats=("html",))
        try:
            paths += save_plotly(fig, output_basepath, formats=("svg",))
        except Exception as exc:  # pragma: no cover - depends on kaleido
            log.info(
                "gene_track_plotly: static SVG export unavailable (%s); "
                "interactive HTML written regardless",
                exc,
            )

        top_prop: dict[str, float] = {}
        for i in cluster_indices:
            row = panel.proportions[i, :]
            finite = row[np.isfinite(row)]
            top_prop[panel.clusters[i]] = float(finite.max()) if finite.size > 0 else 0.0

        meta: dict[str, Any] = {
            "viz_strategy": self.name,
            "description": (
                "Gene track (interactive): PAS positions, per-cluster usage with a "
                "switchable metric (proportion / reads per cell / raw reads), a "
                "proportion filter, PAS bands, the PAS distance table, and full "
                "PAS + cell metadata on hover."
            ),
            "gene_id": panel.gene_id,
            "gene_name": panel.gene_name,
            "chrom": panel.chrom,
            "start": panel.start,
            "end": panel.end,
            "strand": panel.strand,
            "n_pas": len(panel.pas_ids),
            "n_clusters_rendered": len(cluster_indices),
            "n_isoforms": len(panel.isoforms),
            "gene_model_regions_available": bool(panel.isoform_regions),
            "cluster_key": panel.cluster_key,
            "metrics_available": [m[0] for m in _METRICS],
            "proportion_filter_steps": _PROP_STEPS,
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
        """Return (panel, top_n_clusters) from the raw data argument."""
        if isinstance(data, dict):
            panel: GenePanel | None = data.get("panel")
            top_n: int = int(data.get("top_n_clusters", _DEFAULT_TOP_N_CLUSTERS))
        elif isinstance(data, GenePanel):
            panel = data
            top_n = _DEFAULT_TOP_N_CLUSTERS
        else:
            log.error("gene_track_plotly: unexpected data type %s", type(data).__name__)
            return None, _DEFAULT_TOP_N_CLUSTERS
        return panel, top_n

    @staticmethod
    def _select_clusters(panel: GenePanel, top_n: int) -> tuple[list[int], bool]:
        """Pick up to *top_n* non-empty cluster indices by total reads.

        Aligned with GeneTrackMatplotlib._select_top_clusters: clusters whose
        row of ``reads`` sums to zero are excluded before capping, so both
        backends render the SAME cluster set for a given panel.
        """
        total_reads = panel.reads.sum(axis=1)
        non_empty = [i for i in range(len(panel.clusters)) if total_reads[i] > 0]
        ranked = sorted(non_empty, key=lambda i: total_reads[i], reverse=True)
        capped = len(ranked) > top_n
        return sorted(ranked[:top_n]), capped

    @staticmethod
    def _distance_table(panel: GenePanel) -> Any:
        """The same table the matplotlib backend draws; None if unavailable."""
        try:
            return pas_distance_table(panel)
        except Exception as exc:  # pragma: no cover - defensive
            log.warning("gene_track_plotly: distance table unavailable: %s", exc)
            return None

    # -- per-PAS metadata ------------------------------------------------

    @staticmethod
    def _pas_metadata(panel: GenePanel, dist_df: Any) -> list[dict[str, Any]]:
        """One metadata record per PAS, in panel (coordinate) order.

        Pulls rank / gap / summit-distance / assigned UTR transcripts out of
        the distance table so the hover carries exactly what the static
        figure's table reports -- one source of truth, not a re-derivation.
        """
        n = len(panel.pas_ids)
        starts = list(panel.pas_starts) if len(panel.pas_starts) == n else [p for p in panel.pas_positions]
        ends = list(panel.pas_ends) if len(panel.pas_ends) == n else [p + 1 for p in panel.pas_positions]

        by_pas_id: dict[int, dict[str, Any]] = {}
        if dist_df is not None and len(dist_df) > 0:
            for rec in dist_df.to_dict("records"):
                by_pas_id[int(rec["pas_id"])] = rec

        out: list[dict[str, Any]] = []
        for j in range(n):
            pid = int(panel.pas_ids[j])
            rec = by_pas_id.get(pid, {})

            def _num(key: str) -> Any:
                v = rec.get(key)
                try:
                    if v is None or (isinstance(v, float) and not np.isfinite(v)):
                        return "—"
                except Exception:
                    pass
                return v if v is not None and str(v) != "<NA>" else "—"

            out.append(
                {
                    "pas_id": pid,
                    "start": int(starts[j]),
                    "end": int(ends[j]),
                    "summit": int(panel.pas_positions[j]),
                    "width": int(ends[j]) - int(starts[j]),
                    "rank": _num("rank"),
                    "gap_to_next": _num("gap_to_next_bp"),
                    "summit_dist_to_next": _num("summit_dist_to_next_bp"),
                    "utr": rec.get("utr_transcripts") or ";".join(panel.pas_isoforms.get(pid, [])) or "—",
                }
            )
        return out

    # -- figure assembly -------------------------------------------------

    def _build_figure(
        self,
        panel: GenePanel,
        cluster_indices: list[int],
        dist_df: Any,
    ) -> go.Figure:
        """Assemble the multi-panel interactive figure."""
        has_isoforms = len(panel.isoforms) > 0
        has_table = dist_df is not None and len(dist_df) > 0
        n_clusters = len(cluster_indices)

        n_rows = n_clusters + (1 if has_isoforms else 0) + (1 if has_table else 0)

        # Row heights: isoform track shorter than the coverage rows; the table
        # gets a share proportional to its row count (bounded).
        weights: list[float] = []
        specs: list[list[dict[str, Any]]] = []
        if has_isoforms:
            weights.append(max(0.55, len(panel.isoforms) * 0.22))
            specs.append([{"type": "xy"}])
        weights += [1.0] * n_clusters
        specs += [[{"type": "xy"}] for _ in range(n_clusters)]
        if has_table:
            weights.append(min(2.5, 0.5 + 0.16 * len(dist_df)))
            specs.append([{"type": "table"}])
        total_w = sum(weights)
        row_heights = [w / total_w for w in weights]

        titles = self._build_subplot_titles(panel, cluster_indices, has_isoforms, has_table)

        fig = make_subplots(
            rows=n_rows,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.038,
            subplot_titles=titles,
            row_heights=row_heights,
            specs=specs,
        )

        pas_meta = self._pas_metadata(panel, dist_df)
        row = 1
        if has_isoforms:
            self._add_isoform_traces(fig, panel, row=row)
            row += 1

        first_cluster_row = row
        bar_trace_idx: list[int] = []
        for ci in cluster_indices:
            bar_trace_idx.append(len(fig.data))
            self._add_cluster_bar_traces(fig, panel, pas_meta, cluster_idx=ci, row=row)
            row += 1
        last_cluster_row = row - 1

        if has_table:
            self._add_distance_table(fig, dist_df, row=row)

        # PAS bands span every xy row (the table row has no x-axis).
        del first_cluster_row  # rows 1..last_cluster_row are exactly the xy rows
        self._add_pas_bands(fig, panel, pas_meta, last_cluster_row, has_isoforms)

        self._apply_layout(fig, panel, cluster_indices, bar_trace_idx, n_rows, has_table)
        return fig

    @staticmethod
    def _build_subplot_titles(
        panel: GenePanel,
        cluster_indices: list[int],
        has_isoforms: bool,
        has_table: bool,
    ) -> list[str]:
        """Row titles carrying the cell metadata for each cluster."""
        titles: list[str] = []
        if has_isoforms:
            titles.append("gene structure")
        for ci in cluster_indices:
            n_cells = panel.n_cells_per_cluster[ci]
            label = panel.clusters[ci]
            # Numeric labels are Leiden ids and read badly bare ("3"); a
            # descriptive label (stage / cell type) is already self-describing.
            shown = f"{panel.cluster_key} {label}" if str(label).isdigit() else str(label)
            cond = ""
            if panel.group_conditions and ci < len(panel.group_conditions):
                if panel.group_conditions[ci]:
                    cond = f" · {panel.group_conditions[ci]}"
            titles.append(f"{shown} ({n_cells} cells){cond}")
        if has_table:
            titles.append("PAS distances (5'→3')")
        return titles

    @staticmethod
    def _add_isoform_traces(fig: go.Figure, panel: GenePanel, row: int) -> None:
        """Gene model per isoform: intron backbone + strand arrows + region blocks.

        When ``panel.isoform_regions`` is populated (GTF had CDS/UTR feature
        lines) draws the standard gene-model shape per transcript: thin UTR
        boxes, tall CDS boxes, and a distinctly-coloured/outlined 3'UTR (the
        APA-relevant region). Falls back to the legacy uniform exon block
        when regions aren't available for a transcript (e.g. non-coding
        transcripts with exon-only annotation), or for the whole gene when
        ``isoform_regions`` wasn't loaded at all.
        """
        regions_by_tid = {r.transcript_id: r for r in panel.isoform_regions}

        for iso_idx, (transcript_id, exons) in enumerate(panel.isoforms):
            color = _ISOFORM_COLORS[iso_idx % len(_ISOFORM_COLORS)]
            y_val = iso_idx
            regions = regions_by_tid.get(transcript_id)

            if exons:
                iso_min = min(s for s, _ in exons)
                iso_max = max(e for _, e in exons)
                # Intron backbone.
                fig.add_trace(
                    go.Scatter(
                        x=[iso_min, iso_max],
                        y=[y_val, y_val],
                        mode="lines",
                        line=dict(color=color, width=1),
                        hoverinfo="text",
                        hovertext=f"{transcript_id} ({iso_min:,}-{iso_max:,})",
                        showlegend=False,
                    ),
                    row=row, col=1,
                )
                # Strand-direction arrows along the backbone.
                step = max((iso_max - iso_min) / 9.0, 1.0)
                marker = "triangle-right" if panel.strand != "-" else "triangle-left"
                arrow_x = list(np.arange(iso_min + step / 2, iso_max, step))
                if arrow_x:
                    fig.add_trace(
                        go.Scatter(
                            x=arrow_x,
                            y=[y_val] * len(arrow_x),
                            mode="markers",
                            marker=dict(symbol=marker, size=6, color=color),
                            hoverinfo="skip",
                            showlegend=False,
                        ),
                        row=row, col=1,
                    )

            if regions is not None and regions.has_typed_regions:
                GeneTrackPlotly._add_region_blocks(fig, regions, transcript_id, y_val, row)
            else:
                # Legacy fallback: uniform exon blocks.
                for exon_start, exon_end in exons:
                    fig.add_trace(
                        go.Scatter(
                            x=[exon_start, exon_end, exon_end, exon_start, exon_start],
                            y=[y_val - 0.32, y_val - 0.32, y_val + 0.32, y_val + 0.32, y_val - 0.32],
                            mode="lines",
                            fill="toself",
                            fillcolor=color,
                            line=dict(color=color, width=0),
                            hoverinfo="text",
                            hovertext=(
                                f"<b>{transcript_id}</b><br>exon {exon_start:,}-{exon_end:,}"
                                f"<br>width: {exon_end - exon_start:,} bp"
                            ),
                            showlegend=False,
                        ),
                        row=row, col=1,
                    )

        fig.update_yaxes(
            row=row, col=1,
            tickmode="array",
            tickvals=list(range(len(panel.isoforms))),
            ticktext=[t[:22] for t, _ in panel.isoforms],
            tickfont=dict(size=9),
            range=[-0.7, len(panel.isoforms) - 0.3],
            showgrid=False, zeroline=False,
        )

    @staticmethod
    def _add_region_blocks(
        fig: go.Figure,
        regions: TranscriptRegions,
        transcript_id: str,
        y_val: float,
        row: int,
    ) -> None:
        """Draw CDS / 5'UTR / 3'UTR boxes for one transcript's gene model.

        CDS is drawn tall (y +/- 0.34), UTRs short (y +/- 0.18) -- the
        standard gene-model convention -- and 3'UTR gets its own colour plus
        a thicker outline and a text label since it's the region APA acts on.
        """
        specs = [
            ("cds", regions.cds, "CDS", _CDS_COLOR, _CDS_COLOR, 0.34, 0, False),
            ("utr5", regions.utr5, "5'UTR", _UTR5_COLOR, _UTR5_COLOR, 0.18, 0, False),
            ("utr3", regions.utr3, "3'UTR", _UTR3_COLOR, _UTR3_OUTLINE, 0.18, 1.4, True),
        ]
        for kind, intervals, label, fill, outline, half_h, outline_w, annotate in specs:
            for seg_start, seg_end in intervals:
                fig.add_trace(
                    go.Scatter(
                        x=[seg_start, seg_end, seg_end, seg_start, seg_start],
                        y=[y_val - half_h, y_val - half_h, y_val + half_h, y_val + half_h, y_val - half_h],
                        mode="lines",
                        fill="toself",
                        fillcolor=fill,
                        line=dict(color=outline, width=outline_w),
                        hoverinfo="text",
                        hovertext=(
                            f"<b>{transcript_id}</b><br>{label} {seg_start:,}-{seg_end:,}"
                            f"<br>width: {seg_end - seg_start:,} bp"
                        ),
                        showlegend=False,
                    ),
                    row=row, col=1,
                )
                if annotate:
                    fig.add_annotation(
                        x=(seg_start + seg_end) / 2, y=y_val + half_h,
                        text=label, showarrow=False, yshift=7,
                        font=dict(size=7, color=_UTR3_OUTLINE),
                        row=row, col=1,
                    )

    @staticmethod
    def _add_cluster_bar_traces(
        fig: go.Figure,
        panel: GenePanel,
        pas_meta: list[dict[str, Any]],
        cluster_idx: int,
        row: int,
    ) -> None:
        """One bar trace for a cluster, carrying full PAS + cell metadata."""
        ci = cluster_idx
        n_pas = len(panel.pas_ids)
        n_cells = panel.n_cells_per_cluster[ci]
        cluster_label = str(panel.clusters[ci])
        condition = ""
        if panel.group_conditions and ci < len(panel.group_conditions):
            condition = str(panel.group_conditions[ci] or "")
        cluster_total_reads = float(panel.reads[ci, :].sum())

        proportions_row = panel.proportions[ci, :]
        prop_vals = [
            float(proportions_row[j]) if np.isfinite(proportions_row[j]) else 0.0
            for j in range(n_pas)
        ]

        # Real PAS widths so a merged/wide peak renders wide (matches static).
        gene_span = max(panel.end - panel.start, 1)
        min_w = max(int(gene_span * 0.005), 1)
        bar_widths = [max(int(m["end"] - m["start"]), min_w) for m in pas_meta]

        # Read-count fields go through format_count() into customdata as
        # already-formatted strings ("12.3k", "1.23M") -- raw magnitudes like
        # "235.8" or "12345.0" force the reader to count digits/commas
        # themselves. Proportion stays a plotly-formatted percentage (:.1%)
        # since 0-100% is unambiguous on its own.
        customdata: list[list[Any]] = []
        for j in range(n_pas):
            m = pas_meta[j]
            customdata.append([
                m["pas_id"], panel.chrom, m["start"], m["end"], m["width"], m["summit"],
                m["rank"], m["gap_to_next"], m["summit_dist_to_next"], m["utr"],
                cluster_label, n_cells, condition or "—",
                format_count(panel.reads[ci, j]), format_count(panel.reads_per_cell[ci, j]),
                prop_vals[j], format_count(cluster_total_reads),
            ])

        hovertemplate = (
            "<b>PAS %{customdata[0]}</b> — %{customdata[1]}:%{customdata[2]:,}-%{customdata[3]:,}<br>"
            "width: %{customdata[4]:,} bp · summit: %{customdata[5]:,}<br>"
            "rank 5'→3': %{customdata[6]} · gap to next: %{customdata[7]} bp"
            " · summit dist: %{customdata[8]} bp<br>"
            "UTR: %{customdata[9]}<br>"
            "<b>cluster %{customdata[10]}</b> · cells: %{customdata[11]}"
            " · condition: %{customdata[12]}<br>"
            "reads: %{customdata[13]} · reads/cell: %{customdata[14]}<br>"
            "within-gene proportion (0–100%): %{customdata[15]:.1%}"
            " (cluster total reads: %{customdata[16]})"
            "<extra></extra>"
        )

        fig.add_trace(
            go.Bar(
                x=list(panel.pas_positions),
                y=prop_vals,  # default metric == proportion (matches static)
                width=bar_widths,
                customdata=customdata,
                hovertemplate=hovertemplate,
                marker=dict(color=prop_vals, coloraxis="coloraxis", line=dict(width=0)),
                opacity=1.0,
                name=f"cluster {cluster_label}",
                showlegend=False,
            ),
            row=row, col=1,
        )

    @staticmethod
    def _add_pas_bands(
        fig: go.Figure,
        panel: GenePanel,
        pas_meta: list[dict[str, Any]],
        xy_rows: int,
        has_isoforms: bool,
    ) -> None:
        """Vertical band per PAS spanning every track row, id annotated on top.

        This is what makes a position readable ACROSS clusters -- the static
        figure's defining feature.

        Shapes are appended to ``layout.shapes`` with EXPLICIT axis references
        rather than via ``add_vrect``/``add_vline``. Those helpers run plotly's
        subplot-emptiness check, which reads ``trace.xaxis`` off every trace in
        the figure -- and a ``go.Table`` trace (the PAS distance table) has no
        ``xaxis``, so the helper raises no matter which row is targeted.
        """
        def _ax(r: int) -> tuple[str, str]:
            return ("x" if r == 1 else f"x{r}", "y" if r == 1 else f"y{r}")

        shapes = list(fig.layout.shapes or [])
        for m in pas_meta:
            x0, x1 = m["start"], m["end"]
            if x1 <= x0:
                x1 = x0 + 1
            for r in range(1, xy_rows + 1):
                xref, yref = _ax(r)
                shapes.append(dict(
                    type="rect", xref=xref, yref=f"{yref} domain",
                    x0=x0, x1=x1, y0=0, y1=1,
                    fillcolor=_PAS_BAND_COLOR, line=dict(width=0), layer="below",
                ))
                shapes.append(dict(
                    type="line", xref=xref, yref=f"{yref} domain",
                    x0=m["summit"], x1=m["summit"], y0=0, y1=1,
                    line=dict(color=_PAS_BAND_LINE, width=0.8, dash="dot"),
                    layer="below",
                ))
        fig.update_layout(shapes=shapes)

        # PAS id labels above the top track (explicit refs, same reason).
        annotations = list(fig.layout.annotations or [])
        xref, yref = _ax(1)
        for m in pas_meta:
            annotations.append(dict(
                x=m["summit"], y=1.02, xref=xref, yref=f"{yref} domain",
                text=str(m["pas_id"]), showarrow=False, yanchor="bottom",
                font=dict(size=8, color="#b23"),
            ))
        fig.update_layout(annotations=annotations)

    @staticmethod
    def _add_distance_table(fig: go.Figure, dist_df: Any, row: int) -> None:
        """The PAS distance table as a real (scrollable) plotly table row."""
        df = dist_df.copy()
        for c in df.columns:
            df[c] = df[c].astype(object).where(df[c].notna(), "—")
        fig.add_trace(
            go.Table(
                header=dict(
                    values=[f"<b>{c}</b>" for c in df.columns],
                    fill_color="#eef2f5",
                    font=dict(size=10),
                    align="right",
                ),
                cells=dict(
                    values=[df[c].tolist() for c in df.columns],
                    fill_color="white",
                    font=dict(size=9),
                    align="right",
                    height=18,
                ),
            ),
            row=row, col=1,
        )

    # -- layout + interactive controls ------------------------------------

    def _apply_layout(
        self,
        fig: go.Figure,
        panel: GenePanel,
        cluster_indices: list[int],
        bar_trace_idx: list[int],
        n_rows: int,
        has_table: bool,
    ) -> None:
        """Title, colour axis, and the metric dropdown + proportion slider."""
        gene_span = max(panel.end - panel.start, 1)
        gene_label = (
            f"{panel.gene_name} ({panel.gene_id})" if panel.gene_name else panel.gene_id
        )
        title = (
            f"Gene {gene_label} — "
            f"chr{panel.chrom}:{panel.start:,}-{panel.end:,} "
            f"({'−' if panel.strand == '-' else '+'} strand, {_fmt_span(gene_span)}) — "
            f"{len(panel.pas_ids)} PAS"
        )
        if panel.subtitle:
            title = f"{title}<br><span style='font-size:11px'>{panel.subtitle}</span>"

        n_cluster_rows = len(cluster_indices)
        height = max(
            520,
            150 * n_cluster_rows + 26 * len(panel.isoforms) + (26 * min(len(panel.pas_ids), 14) if has_table else 0) + 200,
        )

        fig.update_layout(
            template="plotly_white",
            # Title is pinned to the CONTAINER top; the metric dropdown and the
            # proportion slider live just above the plotting area (paper
            # y~1.01). Without that split they render on top of each other.
            title=dict(
                text=title, font=dict(size=15),
                x=0.02, xanchor="left", y=0.985, yanchor="top", yref="container",
            ),
            height=height,
            autosize=True,
            margin=dict(l=90, r=40, t=200, b=60),
            showlegend=False,
            bargap=0.0,
            hovermode="closest",
            coloraxis=dict(
                colorscale="Viridis", cmin=0.0, cmax=1.0,
                colorbar=dict(
                    title=dict(text="Within-gene<br>proportion", side="right"),
                    thickness=12, len=0.45, y=0.5,
                ),
            ),
            updatemenus=self._metric_menu(panel, cluster_indices, bar_trace_idx),
            sliders=self._proportion_slider(panel, cluster_indices, bar_trace_idx),
        )

        # Only the bottom-most x-axis carries the label; a table row (when
        # present) has no x-axis, so the last xy row is the anchor.
        xy_rows = n_rows - (1 if has_table else 0)
        bottom = f"xaxis{xy_rows}" if xy_rows > 1 else "xaxis"
        # The axis title is dropped when the distance table is present: the
        # table row starts immediately under the last track, so the title
        # lands on top of its header. The figure title already states the
        # chromosome and strand, and the table repeats both as columns.
        if not has_table:
            fig.update_layout(**{
                bottom: dict(
                    title=dict(
                        text=f"genomic position — chr{panel.chrom} ({'−' if panel.strand == '-' else '+'} strand)",
                        font=dict(size=11),
                    ),
                )
            })
        for i in range(1, xy_rows + 1):
            fig.update_yaxes(row=i, col=1, tickfont=dict(size=8))
        # Subplot titles are annotations; keep them small and left-aligned.
        for ann in fig.layout.annotations:
            if ann.text in ("gene structure", "PAS distances (5'→3')") or "cells)" in (ann.text or ""):
                ann.font = dict(size=10)
                ann.xanchor = "left"
                ann.x = 0.0

    @staticmethod
    def _metric_menu(
        panel: GenePanel,
        cluster_indices: list[int],
        bar_trace_idx: list[int],
    ) -> list[dict[str, Any]]:
        """Dropdown switching bar height between proportion / reads-per-cell / reads.

        Bar COLOUR stays proportion-mapped in every mode, so the shared colour
        axis keeps exactly one meaning no matter which height is displayed.
        """
        buttons: list[dict[str, Any]] = []
        for label, attr, axis_title in _METRICS:
            mat = getattr(panel, attr)
            ys = []
            for ci in cluster_indices:
                row = mat[ci, :]
                ys.append([float(v) if np.isfinite(v) else 0.0 for v in row])
            relayout: dict[str, Any] = {}
            for pos in range(len(cluster_indices)):
                axis = "yaxis" if pos == 0 and not panel.isoforms else f"yaxis{pos + 1 + (1 if panel.isoforms else 0)}"
                relayout[f"{axis}.title.text"] = axis_title if pos == len(cluster_indices) // 2 else ""
            buttons.append(dict(
                label=label,
                method="update",
                args=[{"y": ys}, relayout, bar_trace_idx],
            ))
        return [dict(
            type="dropdown",
            buttons=buttons,
            direction="down",
            showactive=True,
            x=0.0, xanchor="left", y=1.045, yanchor="bottom",
            pad=dict(t=4, b=4),
            bgcolor="#ffffff",
            bordercolor="#c8ccd2",
            font=dict(size=11),
        )]

    @staticmethod
    def _proportion_slider(
        panel: GenePanel,
        cluster_indices: list[int],
        bar_trace_idx: list[int],
    ) -> list[dict[str, Any]]:
        """Slider dimming PAS below a within-gene-proportion threshold.

        Masks OPACITY rather than the y values, so it composes with the metric
        dropdown (which owns y) instead of overwriting it -- and a filtered PAS
        stays visible-but-muted rather than vanishing, so the reader never
        silently loses a site from the gene's picture.
        """
        steps: list[dict[str, Any]] = []
        for thr in _PROP_STEPS:
            opacities = []
            for ci in cluster_indices:
                row = panel.proportions[ci, :]
                opacities.append([
                    1.0 if (np.isfinite(v) and float(v) >= thr) else _DIMMED_OPACITY
                    for v in row
                ])
            steps.append(dict(
                label=f"≥{thr:.0%}" if thr > 0 else "all",
                method="restyle",
                args=[{"marker.opacity": opacities}, bar_trace_idx],
            ))
        return [dict(
            active=0,
            steps=steps,
            x=0.34, xanchor="left", y=1.045, yanchor="bottom",
            len=0.42,
            pad=dict(t=4, b=4),
            currentvalue=dict(prefix="show PAS with proportion ", font=dict(size=11)),
            font=dict(size=10),
        )]
