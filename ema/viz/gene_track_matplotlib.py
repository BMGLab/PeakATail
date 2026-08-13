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
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

from ema.viz import register_viz_strategy
from ema.viz.base import VizStrategy
from ema.viz._io import save_matplotlib
from ema.viz._meta import write_figure_meta
from ema.viz._gene_track_helpers import format_count, pas_distance_table

log = logging.getLogger(__name__)

# Maximum clusters rendered even when panel has more.
_MAX_CLUSTERS = 12
# Above this many PAS the in-figure distance table degrades to a summary plus
# the largest gaps; the complete table is always written to the CSV sidecar.
_MAX_TABLE_ROWS = 12
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

        # Shared x range: include BOTH the PAS coords (pasbed-derived) and
        # the GTF isoform structure span so a gene whose exons extend
        # outside the detected-PAS window still renders in full.  Without
        # this, ``panel.start``/``panel.end`` came from pasbed only and
        # exons could fall off-screen on the left or right side.
        x_lo, x_hi = panel.start, panel.end
        if panel.isoforms:
            for _, exons in panel.isoforms:
                if not exons:
                    continue
                ex_lo = min(s for s, _ in exons)
                ex_hi = max(e for _, e in exons)
                x_lo = min(x_lo, ex_lo)
                x_hi = max(x_hi, ex_hi)
        gene_span = max(x_hi - x_lo, 1)
        pad = int(gene_span * 0.05) + 1
        x_min = x_lo - pad
        x_max = x_hi + pad

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

        # Subplot heights in order: [gene_structure (optional), cluster0, ...,
        # distance_table (optional)]
        subplot_heights: list[float] = []
        if has_structure:
            subplot_heights.append(_ISOFORM_ROW_HEIGHT * n_isoforms)
        subplot_heights.extend([_CLUSTER_ROW_HEIGHT] * n_clusters_rendered)
        n_track_axes = len(subplot_heights)

        dist_df = pas_distance_table(panel) if panel.show_distance_table else None
        if dist_df is not None and not dist_df.empty:
            n_tbl_rows = min(len(dist_df), _MAX_TABLE_ROWS) + 1  # +1 header
            # Extra 0.55in of headroom: the shared x-axis label lives here too.
            subplot_heights.append(0.80 + 0.22 * n_tbl_rows)

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

        # One colour per distinct condition (e.g. healthy / primary tumour /
        # metastasis), assigned in first-appearance order so the palette follows
        # the track order rather than the alphabet.
        cond_colours: dict[str, str] = {}
        if panel.group_conditions:
            _palette = [
                "#1b7837", "#2166ac", "#b2182b", "#762a83",
                "#e08214", "#4d4d4d", "#01665e", "#8c510a",
            ]
            _seen: list[str] = []
            for c in panel.group_conditions:
                if c and c not in _seen:
                    _seen.append(c)
            cond_colours = {c: _palette[i % len(_palette)] for i, c in enumerate(_seen)}

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
            # Only prefix numeric labels (Leiden ids). A descriptive label --
            # a stage name, a cell type -- speaks for itself, and writing
            # "cluster Normal" for --cluster-key stage was actively wrong.
            if cluster_label.isdigit():
                track_label = f"{panel.cluster_key} {cluster_label}"
            else:
                track_label = cluster_label
            ax.set_ylabel(
                f"{track_label}\n(n={n_cells})",
                fontsize=7,
                labelpad=4,
                color=cond_colours.get(
                    panel.group_conditions[ci] if panel.group_conditions else "", "black"
                ),
            )
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(axis="y", labelsize=6)
            # Encode condition as a faint wash behind the whole track. A
            # coloured left spine (the previous encoding) is a vertical rule at
            # x=0 with the same width and orientation as a PAS bar, and readers
            # mistook it for data. A background band cannot be confused with a
            # peak because it spans the axes rather than rising from the floor.
            if panel.group_conditions:
                cond = panel.group_conditions[ci]
                if cond in cond_colours:
                    r, g, b = mcolors.to_rgb(cond_colours[cond])
                    ax.set_facecolor((r, g, b, 0.07))
            ax.set_xlim(x_min, x_max)
            # Same human-readable count formatting as the plotly hover
            # ("12.3k" / "1.23M") rather than raw/scientific digits, so the
            # two backends read consistently.
            ax.yaxis.set_major_formatter(
                matplotlib.ticker.FuncFormatter(lambda v, _: format_count(v))  # type: ignore[attr-defined]
            )

        # --- shared x-axis ---
        # Show distance from gene start in the unit that best fits the
        # viewing window:
        #   <  5 kb span   → bp        (tight zoom: a single peak / small gene)
        #   5 kb – 1 Mb    → kb        (normal gene scale)
        #   >= 1 Mb        → Mb        (very long genes; keeps ticks short)
        # Anchoring the axis at gene_start (so it reads 0 → gene_length)
        # avoids the cognitive trap of seeing absolute coordinates like
        # ``155,276 kb`` and mistaking the view for a chromosome-scale
        # window.  Absolute coordinates remain in the figure title for IGV
        # / UCSC cross-reference.
        # The last *track* axis, not axes[-1]: a distance-table axis may sit
        # below it. sharex=True suppresses tick labels on every axis but the
        # bottom one, so re-enable them here.
        bottom_ax = axes[n_track_axes - 1]
        bottom_ax.tick_params(axis="x", labelbottom=True)
        gene_anchor = x_lo
        if gene_span >= 1_000_000:
            unit_label, unit_div, fmt_str = "Mb", 1_000_000, "{val:,.2f}"
        elif gene_span >= 5_000:
            unit_label, unit_div, fmt_str = "kb", 1_000, "{val:,.1f}"
        else:
            unit_label, unit_div, fmt_str = "bp", 1, "{val:,.0f}"

        def _format_tick(x: float, _: int) -> str:
            return fmt_str.format(val=(x - gene_anchor) / unit_div)

        bottom_ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(_format_tick)  # type: ignore[attr-defined]
        )
        bottom_ax.set_xlabel(
            f"Distance from gene start ({unit_label})  "
            f"— anchor chr{panel.chrom}:{gene_anchor:,}",
            fontsize=9,
        )
        bottom_ax.tick_params(axis="x", labelsize=7)

        # --- PAS distance table ---
        if dist_df is not None and not dist_df.empty:
            _draw_distance_table(axes[-1], dist_df, panel)

        # --- figure title ---
        # Make the gene-scope explicit:
        #   1. lead with the gene symbol (e.g. "CLIC2") when available, so the
        #      reader immediately sees this is a one-gene view, not a region
        #      or chromosome view;
        #   2. include the Ensembl ID in parens for unambiguous lookup;
        #   3. show absolute genomic range AND gene length in a matching
        #      unit so the axis below (0..gene_length) is anchored to a
        #      concrete coordinate.
        n_pas = len(panel.pas_ids)
        if gene_span >= 1_000_000:
            length_str = f"{gene_span / 1_000_000:,.2f} Mb"
        elif gene_span >= 1_000:
            length_str = f"{gene_span / 1_000:,.1f} kb"
        else:
            length_str = f"{gene_span:,} bp"
        if panel.gene_name:
            gene_label = f"{panel.gene_name} ({panel.gene_id})"
        else:
            gene_label = panel.gene_id
        title = (
            f"Gene {gene_label} — "
            f"chr{panel.chrom}:{panel.start:,}-{panel.end:,} "
            f"({'−' if panel.strand == '-' else '+'} strand, "
            f"{length_str}) — {n_pas} PAS"
        )
        # A panel is usually restricted to one cell type; name it once in the
        # header rather than repeating it on every track label.
        if panel.subtitle:
            title = f"{title}\n{panel.subtitle}"
        fig.suptitle(title, fontsize=9, y=1.01)

        # Condition legend (healthy / primary tumour / metastasis, ...).
        if cond_colours:
            # Swatches match the track wash: a filled band with a saturated
            # edge, not a line (a line reads as a bar).
            handles = [
                mpatches.Patch(
                    facecolor=(*mcolors.to_rgb(col), 0.25),
                    edgecolor=col,
                    linewidth=0.8,
                    label=cond,
                )
                for cond, col in cond_colours.items()
            ]
            fig.legend(
                handles=handles,
                loc="upper right",
                bbox_to_anchor=(0.995, 1.0),
                fontsize=6,
                frameon=False,
                title="Condition",
                title_fontsize=6,
            )

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
                "gene_name": panel.gene_name,
                "chrom": panel.chrom,
                "start": panel.start,
                "end": panel.end,
                "strand": panel.strand,
                "view_start": int(x_lo),
                "view_end": int(x_hi),
                "view_unit": unit_label,
                "n_pas": n_pas,
                "n_clusters_rendered": n_clusters_rendered,
                "n_isoforms": n_isoforms,
                "gene_model_regions_available": bool(panel.isoform_regions),
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


def _fmt_bp(v: Any) -> str:
    """Format a base-pair count, rendering NA as an em dash."""
    if v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NA:
        return "—"
    try:
        return f"{int(v):,}"
    except (TypeError, ValueError):
        return "—"


def _draw_distance_table(
    ax: "plt.Axes",  # type: ignore[name-defined]
    dist_df: "pd.DataFrame",  # type: ignore[name-defined]
    panel: Any,
) -> None:
    """Render the PAS distance table beneath the tracks.

    Genes with more than ``_MAX_TABLE_ROWS`` PAS get a summary plus the largest
    gaps instead of every row -- a 38-row table is unreadable at figure scale.
    The truncation is stated in the caption line, never silent, and the full
    table is written alongside the figure as a CSV.
    """
    ax.axis("off")
    n = len(dist_df)
    truncated = n > _MAX_TABLE_ROWS

    if truncated:
        # Keep the rows a reader would ask about: the widest gaps.
        show = (
            dist_df.dropna(subset=["gap_to_next_bp"])
            .sort_values("summit_dist_to_next_bp", ascending=False)
            .head(_MAX_TABLE_ROWS)
            .sort_values("rank")
        )
    else:
        show = dist_df

    # Single-line headers: matplotlib table cells do not grow to fit a second
    # line, so "\n" in a header is silently clipped.
    has_utr = "utr_transcripts" in show.columns
    col_labels = ["#", "PAS id", "start", "end", "width bp",
                  "summit", "gap→next bp", "summit→next bp"]
    if has_utr:
        col_labels += ["3'UTR (transcript)"]

    def _utr_cell(v: str) -> str:
        """One transcript per PAS is the common case; summarise when several."""
        tx = [t for t in str(v).split(";") if t]
        if not tx:
            return "—"
        return tx[0] if len(tx) == 1 else f"{tx[0]} +{len(tx) - 1}"

    cells = []
    for _, r in show.iterrows():
        row = [
            str(r["rank"]), str(r["pas_id"]), f"{int(r['start']):,}",
            f"{int(r['end']):,}", f"{int(r['width_bp']):,}",
            f"{int(r['summit_pos']):,}",
            _fmt_bp(r["gap_to_next_bp"]), _fmt_bp(r["summit_dist_to_next_bp"]),
        ]
        if has_utr:
            row.append(_utr_cell(r["utr_transcripts"]))
        cells.append(row)

    # Explicit bbox rather than loc=: the axis above this one owns the shared
    # x-axis label, which is drawn *below* its axes and would collide with a
    # table anchored to the top of ours.
    tbl = ax.table(cellText=cells, colLabels=col_labels, cellLoc="right",
                   bbox=[0.02, 0.16, 0.96, 0.68])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(6)
    for (row, _col), cell in tbl.get_celld().items():
        cell.set_linewidth(0.3)
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#f0f0f0")

    gaps = dist_df["summit_dist_to_next_bp"].dropna()
    if len(gaps):
        stats = (f"n_PAS={n};  adjacent-summit distance "
                 f"min={int(gaps.min()):,}  median={int(gaps.median()):,}  "
                 f"max={int(gaps.max()):,} bp")
    else:
        stats = f"n_PAS={n}; single PAS, no adjacent pair"
    note = (f"PAS ordered 5'→3' ({panel.strand} strand). {stats}.")
    if truncated:
        note += (f"  Showing the {_MAX_TABLE_ROWS} largest of {n - 1} adjacent "
                 f"gaps; full table in the accompanying _pas_distances.csv.")
    if "utr_transcripts" in dist_df.columns:
        note += ("  3'UTR column: transcript(s) in which this PAS is the proximal or "
                 "distal site (per-isoform length output); blank = neither.")
    ax.text(0.5, 0.0, note, transform=ax.transAxes, ha="center", va="bottom",
            fontsize=6, color="#444444")


def _draw_gene_structure(
    ax: "plt.Axes",  # type: ignore[name-defined]
    panel: Any,
    x_min: float,
    x_max: float,
    n_isoforms: int,
) -> None:
    """Draw exon bars and intron lines for each isoform.

    The exon boxes are drawn at the GTF coordinates exactly, BUT clamped to
    a minimum visible width of ~0.4% of the view span so that small exons
    (commonly 80–200 bp on a 50–900 kb gene = 0.1–0.4% of axis width)
    don't render as 1-pixel slivers and disappear visually.  The intron
    backbone behind them is at true coordinates, so the visual lengthening
    only widens the exon box around its true left edge and has no impact
    on PAS positions or anything else.

    Args:
        ax: Matplotlib Axes to draw into.
        panel: GenePanel with isoforms field.
        x_min, x_max: X display range.
        n_isoforms: Number of isoforms (used for y positioning).
    """
    exon_colour = "#4477aa"
    intron_colour = "#888888"
    exon_height = 0.6  # fraction of the row height

    # Gene-model region colours -- matched to gene_track_plotly's palette so
    # the two backends render the same picture.
    cds_colour = "#3a6b91"
    utr5_colour = "#a9c4d8"
    utr3_colour = "#f2b56b"
    utr3_outline = "#c8781f"
    cds_height = 0.75    # CDS tall
    utr_height = 0.40    # UTR short

    view_span = max(x_max - x_min, 1)
    min_visible_exon_w = max(view_span * 0.004, 1)

    regions_by_tid = {r.transcript_id: r for r in getattr(panel, "isoform_regions", [])}

    def _draw_block(start: int, end: int, y_centre: float, height: float,
                     facecolour: str, edgecolour: str = "none", lw: float = 0.0) -> None:
        visual_width = max(end - start, min_visible_exon_w)
        rect = mpatches.FancyBboxPatch(
            (start, y_centre - height / 2),
            visual_width,
            height,
            boxstyle="square,pad=0",
            facecolor=facecolour,
            edgecolor=edgecolour,
            linewidth=lw,
            zorder=3,
        )
        ax.add_patch(rect)

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

        regions = regions_by_tid.get(tid)
        if regions is not None and regions.has_typed_regions:
            # Full gene model: CDS tall, UTRs short, 3'UTR distinctly
            # coloured + outlined + labelled (the APA-relevant region).
            for seg_start, seg_end in regions.cds:
                _draw_block(seg_start, seg_end, y_centre, cds_height, cds_colour)
            for seg_start, seg_end in regions.utr5:
                _draw_block(seg_start, seg_end, y_centre, utr_height, utr5_colour)
            for seg_start, seg_end in regions.utr3:
                _draw_block(seg_start, seg_end, y_centre, utr_height, utr3_colour,
                             edgecolour=utr3_outline, lw=0.9)
                centre = (seg_start + seg_end) / 2.0
                ax.text(
                    centre, y_centre + utr_height / 2 + 0.03, "3'UTR",
                    ha="center", va="bottom", fontsize=4.5, color=utr3_outline,
                    clip_on=False,
                )
        else:
            # Legacy fallback: uniform exon rectangles.  Visual width is
            # clamped up to a minimum so tiny exons remain visible; actual
            # coordinates are unaltered.
            for exon_start, exon_end in exons:
                _draw_block(exon_start, exon_end, y_centre, exon_height, exon_colour)

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
        # every isoform row.  Low alpha so the underlying exon/intron
        # structure stays visible (especially small 80–200 bp exons that
        # the GTF parser correctly draws but which can otherwise get
        # washed out under the PAS tint).
        ax.axvspan(
            left, left + width,
            ymin=0.0, ymax=1.0,
            color=pas_colour, alpha=0.10, lw=0,
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
