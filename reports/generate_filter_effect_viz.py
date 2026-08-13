#!/usr/bin/env python
"""Figures 21–22 — UMAP and Sankey per filter-effect scenario.

Same conventions as the sweep's clustering figures (15–19):

* unlabelled cells are drawn **grey and kept**, never dropped — discarding them
  makes the cluster -> cell-type mapping look far cleaner than it is;
* each panel is a **single named dataset**, because ``leiden`` ids are assigned
  per dataset and pooling merges unrelated clusters into one node;
* the cell-type-labelled fraction is stated on every panel.

    python reports/generate_filter_effect_viz.py \
        --in reports/cumulative_analysis/filter_effect \
        --out reports/figures/corrected
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402

logger = logging.getLogger(__name__)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.2,
    "axes.spines.top": False, "axes.spines.right": False,
})

ORDER = ["baseline", "atlas_filter", "annot_filter_3utr", "ip_filter"]
TITLE = {
    "baseline": "baseline (keep all)",
    "atlas_filter": "atlas non-match (8.7% excluded)",
    "annot_filter_3utr": "3′UTR-only (73.1% excluded)",
    "ip_filter": "internal priming (8.3% excluded)",
}
UNLABELLED = "unlabelled"
GREY = "#d0d0d0"


def short_ct(name: str) -> str:
    return (str(name).replace("CELL_TYPES_WANSLEEBEN_HOGAN_2013:", "")
            .replace("CELL_TYPES_WANSLEEBEN_HOGAN_2013_", "")
            .replace("lung_epithelial_lineage_signatures:", "LE:")
            .replace("lung_epithelial_lineage_signatures_", "LE:")
            .replace("stem_cell_signatures_merged:", "SC:")
            .replace("stem_cell_signatures_merged_", "SC:"))


def celltype_colors(cts) -> dict:
    """Stable colours across every panel, so a type keeps its colour."""
    named = sorted(c for c in cts if c != UNLABELLED)
    cmap = plt.get_cmap("tab20")
    out = {c: cmap(i % 20) for i, c in enumerate(named)}
    out[UNLABELLED] = GREY
    return out


# ---------------------------------------------------------------------------
# 21 — UMAP grid
# ---------------------------------------------------------------------------


def fig_umap(u: pd.DataFrame, out_dir: Path) -> None:
    ds = u["dataset"].iloc[0]
    ct_colors = celltype_colors(u["celltype"].unique())
    scen = [s for s in ORDER if s in set(u["scenario"])]
    fig, axes = plt.subplots(2, len(scen), figsize=(3.3 * len(scen), 6.6))
    axes = np.atleast_2d(axes)

    for j, s in enumerate(scen):
        d = u[u["scenario"] == s]
        # top: by PAS-leiden cluster
        ax = axes[0, j]
        clusters = sorted(d["leiden"].astype(str).unique(), key=lambda v: int(v) if v.isdigit() else 1e9)
        cmap = plt.get_cmap("tab20")
        for i, cl in enumerate(clusters):
            m = d["leiden"].astype(str) == cl
            ax.scatter(d.loc[m, "umap_x"], d.loc[m, "umap_y"], s=1.6,
                       color=cmap(i % 20), linewidths=0, rasterized=True)
        ax.set_title(f"{TITLE[s]}\n{len(d):,} cells · {len(clusters)} clusters",
                     fontsize=8)
        if j == 0:
            ax.set_ylabel("by PAS cluster", fontsize=9)

        # bottom: by GEX cell type, unlabelled grey and drawn first
        ax = axes[1, j]
        un = d["celltype"] == UNLABELLED
        ax.scatter(d.loc[un, "umap_x"], d.loc[un, "umap_y"], s=1.6, color=GREY,
                   linewidths=0, rasterized=True)
        for ct, g in d.loc[~un].groupby("celltype"):
            ax.scatter(g["umap_x"], g["umap_y"], s=1.6, color=ct_colors[ct],
                       linewidths=0, rasterized=True)
        frac = float((~un).mean())
        ax.set_title(f"{frac:.0%} cell-type-labelled", fontsize=8)
        if j == 0:
            ax.set_ylabel("by GEX cell type", fontsize=9)

    for ax in axes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)

    handles = [mpatches.Patch(color=GREY, label="unlabelled (kept, not dropped)")]
    fig.legend(handles=handles, loc="lower center", ncol=1, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(
        f"Filter-effect scenarios on {ds}: same cells, same pipeline, one PAS filter changed.\n"
        "Each embedding is computed independently, so layout differences BETWEEN panels are not "
        "evidence of a clustering difference — see figure 20 for that. Read each panel on its own.",
        fontsize=9.5, y=1.005,
    )
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"21_filter_effect_umap.{ext}")
        logger.info("wrote 21_filter_effect_umap.%s", ext)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 22 — Sankey grid
# ---------------------------------------------------------------------------


def _ribbon(ax, y0a, y1a, y0b, y1b, x0, x1, color, alpha=0.55):
    verts = [(x0, y0a), (x0 + (x1 - x0) * .5, y0a), (x0 + (x1 - x0) * .5, y0b), (x1, y0b),
             (x1, y1b), (x0 + (x1 - x0) * .5, y1b), (x0 + (x1 - x0) * .5, y1a), (x0, y1a)]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
             MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
    ax.add_patch(mpatches.PathPatch(MPath(verts, codes), facecolor=color,
                                    edgecolor="none", alpha=alpha))


def fig_sankey(counts: pd.DataFrame, out_dir: Path, dataset: str) -> None:
    d = counts[counts["dataset"] == dataset]
    scen = [s for s in ORDER if s in set(d["scenario"])]
    ct_colors = celltype_colors(d["celltype"].unique())
    fig, axes = plt.subplots(1, len(scen), figsize=(3.6 * len(scen), 7.2))
    axes = np.atleast_1d(axes)
    gap = 0.012

    for ax, s in zip(axes, scen):
        sub = d[d["scenario"] == s]
        total = sub["n_cells"].sum()
        clusters = sorted(sub["leiden"].astype(str).unique(),
                          key=lambda v: int(v) if v.isdigit() else 1e9)
        cts = ([c for c in sorted(sub["celltype"].unique()) if c != UNLABELLED]
               + ([UNLABELLED] if (sub["celltype"] == UNLABELLED).any() else []))

        lsize = sub.groupby("leiden")["n_cells"].sum()
        rsize = sub.groupby("celltype")["n_cells"].sum()
        ly, ry = {}, {}
        y = 0.0
        for cl in clusters:
            h = lsize.get(cl, 0) / total * (1 - gap * len(clusters))
            ly[cl] = [y, y + h]; y += h + gap
        y = 0.0
        for ct in cts:
            h = rsize.get(ct, 0) / total * (1 - gap * len(cts))
            ry[ct] = [y, y + h]; y += h + gap

        lcur = {c: ly[c][0] for c in clusters}
        rcur = {c: ry[c][0] for c in cts}
        for cl in clusters:
            for ct in cts:
                n = sub[(sub["leiden"].astype(str) == cl) & (sub["celltype"] == ct)]["n_cells"].sum()
                if not n:
                    continue
                ha = n / total * (1 - gap * len(clusters))
                hb = n / total * (1 - gap * len(cts))
                _ribbon(ax, lcur[cl], lcur[cl] + ha, rcur[ct], rcur[ct] + hb,
                        0.12, 0.88, ct_colors[ct],
                        alpha=0.25 if ct == UNLABELLED else 0.6)
                lcur[cl] += ha; rcur[ct] += hb

        for cl in clusters:
            ax.add_patch(mpatches.Rectangle((0.07, ly[cl][0]), 0.05,
                                            ly[cl][1] - ly[cl][0],
                                            color="#555555"))
            if ly[cl][1] - ly[cl][0] > 0.025:
                ax.text(0.06, np.mean(ly[cl]), cl, ha="right", va="center", fontsize=6)
        for ct in cts:
            ax.add_patch(mpatches.Rectangle((0.88, ry[ct][0]), 0.05,
                                            ry[ct][1] - ry[ct][0],
                                            color=ct_colors[ct]))
            if ry[ct][1] - ry[ct][0] > 0.022:
                ax.text(0.94, np.mean(ry[ct]), short_ct(ct)[:22], ha="left",
                        va="center", fontsize=5.6)

        # purity: share of each cluster going to its dominant *labelled* type
        lab = sub[sub["celltype"] != UNLABELLED]
        if not lab.empty:
            top = lab.groupby(["leiden", "celltype"])["n_cells"].sum().reset_index()
            best = top.loc[top.groupby("leiden")["n_cells"].idxmax()]
            tot = lab.groupby("leiden")["n_cells"].sum()
            purity = float((best.set_index("leiden")["n_cells"] / tot).median())
        else:
            purity = float("nan")
        ax.set_title(f"{TITLE[s]}\n{len(clusters)} clusters · median purity {purity:.2f}",
                     fontsize=8)
        ax.set_xlim(0, 1.25)
        ax.set_ylim(-0.02, 1.02)
        ax.axis("off")

    fig.suptitle(
        f"PAS cluster → GEX cell type, per filter scenario ({dataset}).\n"
        "Ribbon colour is the destination cell type; pale grey ribbons are unlabelled cells, kept "
        "rather than dropped. Purity rises mechanically as clusters shrink — it does NOT rank the "
        "scenarios.",
        fontsize=9.5, y=1.0,
    )
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"22_filter_effect_sankey.{ext}")
        logger.info("wrote 22_filter_effect_sankey.%s", ext)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_dir", required=True, type=Path)
    ap.add_argument("--out", dest="out_dir", required=True, type=Path)
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    a.out_dir.mkdir(parents=True, exist_ok=True)

    u = pd.read_csv(a.in_dir / "fe_umap_cells.tsv.gz", sep="\t")
    fig_umap(u, a.out_dir)
    c = pd.read_csv(a.in_dir / "fe_cluster_celltype_counts.tsv", sep="\t")
    fig_sankey(c, a.out_dir, dataset=u["dataset"].iloc[0])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
