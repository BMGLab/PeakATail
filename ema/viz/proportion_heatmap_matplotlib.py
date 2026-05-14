"""Per-PAS proportion heatmap (PAS x cluster) — matplotlib backend.

The proportion strategy produces a long-format DataFrame::

    gene_id  transcript_id  pas_id  rank  cell  proportion

This viz collapses cells -> clusters by averaging, picks the top-N most
variable PAS (by cluster-level variance), and renders a row-normalised
PAS x cluster heatmap so a researcher can see which PAS shift between
cell populations.

Data shape:
    ``(pdui_df, adata)`` tuple — the long DataFrame plus the AnnData
    that carries ``obs["leiden"]`` (or whichever cluster column).
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


def _build_pas_cluster_matrix(
    pdui_df: pd.DataFrame,
    adata,
    cluster_key: str,
    top_n: int,
) -> tuple[pd.DataFrame, list[str], pd.DataFrame] | tuple[None, None, None]:
    """Return (pas x cluster mean proportion table, sorted cluster labels, pos meta).

    Returns ``(None, None, None)`` if the inputs lack the columns we need.

    The returned position metadata has one row per PAS surviving the
    top-N variance cut, with columns ``pas_id, gene_id, chrom, start, end,
    strand``.  Used by the caller to label rows with genomic coordinates
    and to insert gene separator lines on the heatmap.
    """
    if "proportion" not in pdui_df.columns or "cell" not in pdui_df.columns:
        return None, None, None
    if "pas_id" not in pdui_df.columns:
        return None, None, None
    if cluster_key not in adata.obs.columns:
        return None, None, None

    cell_to_cluster = adata.obs[cluster_key].astype(str).to_dict()
    df = pdui_df.dropna(subset=["proportion"]).copy()
    df["cluster"] = df["cell"].map(cell_to_cluster)
    df = df.dropna(subset=["cluster"])
    if df.empty:
        return None, None, None

    # Mean proportion per (pas_id, cluster)
    pivot = (
        df.groupby(["pas_id", "cluster"], observed=True)["proportion"]
          .mean()
          .unstack("cluster")
    )
    # Pick the top-N PAS by cross-cluster variance (most discriminating).
    variances = pivot.var(axis=1, skipna=True).fillna(0.0)
    top_idx = variances.sort_values(ascending=False).head(top_n).index
    pivot = pivot.loc[top_idx]

    cluster_order = sorted(
        pivot.columns.tolist(),
        key=lambda x: int(x) if str(x).isdigit() else x,
    )
    pivot = pivot[cluster_order]

    # Position metadata for each surviving pas_id (one row per pas_id).
    pos_cols = [c for c in ("gene_id", "chrom", "start", "end", "strand")
                if c in pdui_df.columns]
    if pos_cols:
        pos_meta = (
            pdui_df[["pas_id"] + pos_cols]
            .drop_duplicates(subset=["pas_id"])
            .set_index("pas_id")
            .reindex(pivot.index)
            .reset_index()
        )
    else:
        pos_meta = pd.DataFrame({"pas_id": pivot.index})

    # Re-order rows: group by gene_id, ascending start position within each
    # gene.  Genes themselves are ordered by their first-occurring PAS in
    # the original variance ranking — preserves "most-interesting-first"
    # while keeping same-gene PAS contiguous.
    if "gene_id" in pos_meta.columns and "start" in pos_meta.columns:
        # gene order = order of first appearance in the (already variance-
        # sorted) pivot index
        pos_meta["_gene_rank"] = (
            pos_meta["gene_id"].map(
                {g: i for i, g in enumerate(pos_meta["gene_id"].drop_duplicates())}
            )
        )
        pos_meta = pos_meta.sort_values(["_gene_rank", "start"]).drop(columns="_gene_rank")
        pivot = pivot.loc[pos_meta["pas_id"].tolist()]

    return pivot, cluster_order, pos_meta


@register_viz_strategy
class ProportionHeatmapMatplotlib(VizStrategy):
    """PAS x cluster mean-proportion heatmap (matplotlib).

    Tunable knobs are passed via a data dict::

        render_all("proportion_heatmap",
                   {"pdui_df": df, "adata": adata,
                    "cluster_key": "leiden", "top_n": 50}, ...)

    Falls back to ``(pdui_df, adata)`` tuple form with sensible defaults.
    """

    name = "proportion_heatmap_matplotlib"
    plot_type = "proportion_heatmap"
    engine = "matplotlib"

    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        if isinstance(data, dict):
            pdui_df = data["pdui_df"]
            adata = data["adata"]
            cluster_key = data.get("cluster_key", "leiden")
            top_n = int(data.get("top_n", 50))
        else:
            pdui_df, adata = data
            cluster_key, top_n = "leiden", 50

        pivot, clusters, pos_meta = _build_pas_cluster_matrix(
            pdui_df, adata, cluster_key, top_n,
        )
        if pivot is None or pivot.empty:
            log.warning(
                "proportion_heatmap_matplotlib: no rows to plot (cluster_key=%r)",
                cluster_key,
            )
            return []

        # Display values in [0, 1]; clip strictly to avoid color flicker on NaN.
        mat = pivot.fillna(0.0).values
        # Slightly more vertical room per row so the longer (pas | gene |
        # chr:start-end) labels stay readable.
        fig_h = max(5, 0.24 * pivot.shape[0])
        fig_w = max(7, 0.8 * pivot.shape[1] + 4)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)

        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels(clusters)
        ax.set_xlabel("cluster")

        # Build informative y-labels: "pas_id  GENE  chrN:start-end (strand)".
        # When position metadata is missing we fall back to bare pas_id, which
        # reproduces the legacy behaviour.
        def _format_row(pas_id, row) -> str:
            parts = [str(pas_id)]
            if "gene_id" in row and pd.notna(row.get("gene_id")):
                parts.append(str(row["gene_id"]))
            if (
                "chrom" in row and "start" in row and "end" in row
                and pd.notna(row.get("chrom"))
            ):
                strand = row.get("strand", "")
                strand_str = f" ({strand})" if isinstance(strand, str) and strand else ""
                parts.append(
                    f"chr{row['chrom']}:{int(row['start']):,}-{int(row['end']):,}"
                    f"{strand_str}"
                )
            return "  ".join(parts)

        if pivot.shape[0] <= 50 and pos_meta is not None and not pos_meta.empty:
            ax.set_yticks(range(pivot.shape[0]))
            pos_by_pas = pos_meta.set_index("pas_id")
            labels = [
                _format_row(p, pos_by_pas.loc[p] if p in pos_by_pas.index else {})
                for p in pivot.index
            ]
            ax.set_yticklabels(labels, fontsize=7, family="monospace")
        elif pivot.shape[0] <= 50:
            ax.set_yticks(range(pivot.shape[0]))
            ax.set_yticklabels([str(p) for p in pivot.index], fontsize=7)
        else:
            ax.set_yticks([])

        # Insert separator lines between gene groups so same-gene rows
        # group visually.  Skips when no gene metadata is available.
        if (
            pos_meta is not None
            and "gene_id" in pos_meta.columns
            and pivot.shape[0] >= 2
        ):
            ordered_genes = pos_meta["gene_id"].tolist()
            for i in range(1, len(ordered_genes)):
                if ordered_genes[i] != ordered_genes[i - 1]:
                    ax.axhline(
                        y=i - 0.5, color="#FFFFFF", lw=0.8, alpha=0.85,
                    )

        ax.set_ylabel(
            f"PAS (top-{pivot.shape[0]} by variance, grouped by gene & "
            f"sorted by genomic position)"
        )
        ax.set_title("Mean PAS proportion per cluster")
        fig.colorbar(im, ax=ax, label="mean proportion")

        paths = save_matplotlib(fig, output_basepath)
        n_genes = (
            int(pos_meta["gene_id"].nunique())
            if pos_meta is not None and "gene_id" in pos_meta.columns
            else None
        )
        write_figure_meta(output_basepath, {
            "viz_strategy": self.name,
            "cluster_key": cluster_key,
            "n_pas_shown": int(pivot.shape[0]),
            "n_clusters": int(pivot.shape[1]),
            "n_genes": n_genes,
            "top_n_selected_by": "cross-cluster variance",
            "row_order": (
                "grouped by gene_id (variance-rank first PAS), "
                "within-gene by genomic start"
            ),
            "y_label_format": "pas_id  gene_id  chrN:start-end (strand)",
        })
        return paths
