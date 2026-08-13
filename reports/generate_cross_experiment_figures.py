#!/usr/bin/env python
"""Figures for the experiment-driven Laughney report.

Consumes the TSVs written by ``ema.benchmark.cross_experiment.run_all`` and
emits the figures ``reports/REPORT_ANALYSIS_PLAN.md`` calls for — the ones that
actually discriminate between arms, rather than restating per-run numbers.

    python scripts/analyst_figures.py --in <cross_experiment out dir> \
        --out reports/figures/cross_experiment

Every figure degrades gracefully: a missing input TSV is logged and skipped, so
this runs against a partially complete sweep (e.g. before phase 3 lands).
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

logger = logging.getLogger(__name__)

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

OBS = "#1f4e79"
NULL_SHIFT = "#e07b39"
NULL_UNIF = "#bbbbbb"


def _load(in_dir: Path, name: str) -> Optional[pd.DataFrame]:
    path = in_dir / name
    if not path.exists():
        logger.warning("skip: %s not found", name)
        return None
    # index_col=False: a stray trailing delimiter would otherwise make pandas
    # promote column 0 to the index and silently shift every column left.
    df = pd.read_csv(path, sep="\t", index_col=False)
    return df if not df.empty else None


def _save(fig, out_dir: Path, name: str) -> None:
    """Write both raster and vector, matching the manifest's convention.

    ``name`` is the numbered stem (e.g. ``24_atlas_null_control``); figures
    24-36 are the cross-experiment set. Slots 20-23 are reserved for the
    pending filter-effect figures and must not be reused.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        path = out_dir / f"{name}.{ext}"
        fig.savefig(path)
        logger.info("wrote %s", path)
    plt.close(fig)


# ---------------------------------------------------------------------------
# F1 — the atlas benchmark is mostly free
# ---------------------------------------------------------------------------


def fig_atlas_null(df: pd.DataFrame, out_dir: Path) -> None:
    """Observed vs null precision per cutoff — why precision ~1.0 means little."""
    runs = sorted(df["run"].unique())
    cutoffs = sorted(df["cutoff_bp"].unique())
    fig, axes = plt.subplots(1, len(runs), figsize=(2.5 * len(runs), 3.2), sharey=True)
    axes = np.atleast_1d(axes)
    x = np.arange(len(cutoffs))
    for ax, run in zip(axes, runs):
        sub = df[df["run"] == run].set_index("cutoff_bp").loc[cutoffs]
        ax.bar(x - 0.27, sub["precision_observed"], 0.27, color=OBS, label="observed")
        ax.bar(x, sub["precision_null_shift"], 0.27, color=NULL_SHIFT, label="null: shifted 5–50 kb")
        ax.bar(x + 0.27, sub["precision_null_uniform"], 0.27, color=NULL_UNIF, label="null: uniform")
        ax.set_xticks(x)
        ax.set_xticklabels([str(c) for c in cutoffs])
        ax.set_title(run, fontsize=8)
        ax.set_xlabel("cutoff (bp)")
        ax.set_ylim(0, 1.05)
    axes[0].set_ylabel("precision vs PolyASite v3.0")
    axes[0].legend(fontsize=7, loc="lower right", framealpha=0.9)
    fig.suptitle(
        "Atlas precision is ~70% free at ±50 bp: a randomly displaced PAS scores almost as well",
        fontsize=9.5,
    )
    _save(fig, out_dir, "24_atlas_null_control")

    # excess-over-null only — the part that carries information
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    for run in runs:
        sub = df[df["run"] == run].set_index("cutoff_bp").loc[cutoffs]
        ax.plot(cutoffs, sub["excess_over_shift"], marker="o", label=run, lw=1.6)
    ax.set_xscale("log")
    ax.set_xticks(cutoffs)
    ax.set_xticklabels([str(c) for c in cutoffs])
    ax.set_xlabel("cutoff (bp)")
    ax.set_ylabel("excess precision over shifted null")
    ax.set_title("Excess over null: strategies are indistinguishable at every cutoff")
    ax.legend(fontsize=7)
    _save(fig, out_dir, "25_atlas_excess_over_null")


def fig_pas_overlap(df: pd.DataFrame, out_dir: Path) -> None:
    """Pairwise PAS-set Jaccard — same score, different sites."""
    runs = sorted(set(df["run_a"]) | set(df["run_b"]))
    mat = pd.DataFrame(np.eye(len(runs)), index=runs, columns=runs)
    for _, r in df.iterrows():
        mat.loc[r["run_a"], r["run_b"]] = r["jaccard"]
        mat.loc[r["run_b"], r["run_a"]] = r["jaccard"]
    fig, ax = plt.subplots(figsize=(4.6, 3.9))
    im = ax.imshow(mat.values, cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(len(runs)))
    ax.set_xticklabels(runs, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(runs)))
    ax.set_yticklabels(runs, fontsize=7)
    for i in range(len(runs)):
        for j in range(len(runs)):
            v = mat.values[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if v < 0.6 else "black")
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="Jaccard (±50 bp)")
    ax.set_title("PAS-set agreement between peak strategies", fontsize=9.5)
    _save(fig, out_dir, "26_pas_set_jaccard")


def fig_yield_vs_biology(
    yield_df: pd.DataFrame, cluster_df: pd.DataFrame, out_dir: Path
) -> None:
    """PAS yield against cell-type recovery — the flat line."""
    merged = yield_df.merge(
        cluster_df, left_on="run", right_on="branch", how="inner"
    )
    if merged.empty or "ARI_vs_celltype" not in merged:
        return
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.scatter(merged["n_pas"], merged["ARI_vs_celltype"], s=70, color=OBS, zorder=3)
    for _, r in merged.iterrows():
        ax.annotate(r["run"], (r["n_pas"], r["ARI_vs_celltype"]),
                    textcoords="offset points", xytext=(6, 4), fontsize=7)
    lo, hi = merged["ARI_vs_celltype"].min(), merged["ARI_vs_celltype"].max()
    pad = max(0.05, (hi - lo) * 3)
    ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel("PAS called")
    ax.set_ylabel("ARI vs GEX cell type")
    ax.set_title("Doubling the PAS catalogue does not improve cell-type recovery", fontsize=9.5)
    _save(fig, out_dir, "27_yield_vs_celltype_ari")


# ---------------------------------------------------------------------------
# F2/F3 — which knob actually matters
# ---------------------------------------------------------------------------

# label -> every branch on that axis (the swept arms *including* the default, so
# the reported delta is the full range max-min, not an endpoint difference; this
# matters because resolution is non-monotone and peaks at the default).
KNOBS = [
    ("cluster method\ntfidf → libsize", ["A2_trim_default", "A3_libsize"]),
    ("resolution\n0.5 / 1.0 / 2.0", ["A3_res0.5", "A2_trim_default", "A3_res2.0"]),
    ("include_extended\nfalse → true", ["A2_trim_default", "A2_trim_d5000_ext"]),
    ("n_neighbors\n15 / 30 / 50", ["A3_nn15", "A2_trim_default", "A3_nn50"]),
    (
        "utr_multiplier\n1.5 / 2.0 / 3.0",
        ["A2_trim_mult1.5", "A2_trim_default", "A2_trim_mult3.0"],
    ),
    (
        "max_gene_distance\n1000 → 10000",
        [
            "A2_trim_d1000_ext",
            "A2_trim_d2000_ext",
            "A2_trim_d3000_ext",
            "A2_trim_d5000_ext",
            "A2_trim_d10000_ext",
        ],
    ),
]


def fig_knob_ranking(
    branch_df: pd.DataFrame, strategy_df: Optional[pd.DataFrame], out_dir: Path
) -> None:
    """The headline: every swept knob ranked by its effect on the biology."""
    b = branch_df.set_index("branch")
    rows = []
    for label, arms in KNOBS:
        vals = [
            b.loc[a, "ARI_vs_celltype"]
            for a in arms
            if a in b.index and pd.notna(b.loc[a, "ARI_vs_celltype"])
        ]
        if len(vals) >= 2:
            rows.append({"knob": label, "delta_ari": float(max(vals) - min(vals))})
    if strategy_df is not None and "ARI_vs_celltype" in strategy_df:
        s = strategy_df["ARI_vs_celltype"]
        rows.append({"knob": "peak strategy\nlg / lp / si", "delta_ari": float(s.max() - s.min())})
    if not rows:
        return
    df = pd.DataFrame(rows).sort_values("delta_ari")
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    colors = ["#c62828" if v >= 0.05 else "#8c8c8c" for v in df["delta_ari"]]
    ax.barh(df["knob"], df["delta_ari"], color=colors)
    for y, v in enumerate(df["delta_ari"]):
        ax.text(v + 0.002, y, f"{v:.3f}", va="center", fontsize=8)
    ax.set_xlabel("|Δ ARI vs GEX cell type| across the swept range")
    ax.set_title("Only the normalisation method moves the biology", fontsize=10)
    ax.set_xlim(0, max(df["delta_ari"]) * 1.25)
    _save(fig, out_dir, "31_knob_ranking")


def fig_resolution_curve(branch_df: pd.DataFrame, out_dir: Path) -> None:
    """ARI and silhouette disagree about resolution — show both."""
    mapping = {"A3_res0.5": 0.5, "A2_trim_default": 1.0, "A3_res2.0": 2.0}
    sub = branch_df[branch_df["branch"].isin(mapping)].copy()
    if len(sub) < 3:
        return
    sub["resolution"] = sub["branch"].map(mapping)
    sub = sub.sort_values("resolution")
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.plot(sub["resolution"], sub["ARI_vs_celltype"], marker="o", color=OBS,
            label="ARI vs GEX cell type", lw=1.8)
    ax.set_xlabel("leiden resolution")
    ax.set_ylabel("ARI vs GEX cell type", color=OBS)
    ax.tick_params(axis="y", labelcolor=OBS)
    ax2 = ax.twinx()
    ax2.plot(sub["resolution"], sub["silhouette_lsi"], marker="s", color=NULL_SHIFT,
             label="silhouette (LSI)", lw=1.8, ls="--")
    ax2.set_ylabel("silhouette in LSI space", color=NULL_SHIFT)
    ax2.tick_params(axis="y", labelcolor=NULL_SHIFT)
    ax2.grid(False)
    for _, r in sub.iterrows():
        ax.annotate(f"{r['n_clusters']:.0f} clusters",
                    (r["resolution"], r["ARI_vs_celltype"]),
                    textcoords="offset points", xytext=(0, -14), fontsize=7, ha="center")
    ax.set_title("Cell-type agreement peaks at res=1.0;\nsilhouette prefers ever-coarser clusters",
                 fontsize=9.5)
    _save(fig, out_dir, "32_resolution_ari_silhouette")


def fig_celltype_shift(df: pd.DataFrame, out_dir: Path) -> None:
    """% of cells whose cell-type call changes, per branch."""
    agg = (
        df.groupby("branch")["frac_cells_reassigned"]
        .agg(["mean", "std"])
        .sort_values("mean")
    )
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.barh(agg.index, agg["mean"] * 100, xerr=(agg["std"].fillna(0) * 100),
            color=OBS, error_kw={"lw": 0.8, "ecolor": "#999999"})
    ax.set_xlabel("% of cells whose cell-type call changes (vs base cohort)")
    ax.set_title("Blast radius of each parameter on the cell-type call", fontsize=9.5)
    ax.tick_params(axis="y", labelsize=7)
    _save(fig, out_dir, "29_celltype_reassignment")


# max_gene_distance encoded by each trim branch (all with include_extended=true,
# so the series varies one factor).
TRIM_DISTANCE = {
    "A2_trim_d1000_ext": 1000,
    "A2_trim_d2000_ext": 2000,
    "A2_trim_d3000_ext": 3000,
    "A2_trim_d5000_ext": 5000,
    "A2_trim_d10000_ext": 10000,
}


def fig_trim_sensitivity(
    branch_df: pd.DataFrame, yield_df: Optional[pd.DataFrame], out_dir: Path
) -> None:
    """Trim window against cluster count, cell-type ARI and PAS yield.

    The point of the figure is that the first two lines are flat: a 10x change
    in the annotation window does not move the clustering or its agreement with
    cell identity.
    """
    sub = branch_df[branch_df["branch"].isin(TRIM_DISTANCE)].copy()
    if len(sub) < 3:
        return
    sub["distance"] = sub["branch"].map(TRIM_DISTANCE)
    sub = sub.sort_values("distance")
    ncols = 3 if yield_df is not None else 2
    fig, axes = plt.subplots(1, ncols, figsize=(3.2 * ncols, 3.2))

    axes[0].plot(sub["distance"], sub["ARI_vs_celltype"], marker="o", color=OBS, lw=1.8)
    axes[0].set_ylabel("ARI vs GEX cell type")
    # Anchor the scale to the knob that *does* matter, or a flat line looks
    # dramatic on an auto-scaled axis.
    ref = branch_df.set_index("branch")
    if {"A2_trim_default", "A3_libsize"} <= set(ref.index):
        lo = float(ref.loc["A3_libsize", "ARI_vs_celltype"])
        hi = float(ref.loc["A2_trim_default", "ARI_vs_celltype"])
        axes[0].axhspan(lo, hi, color="#c62828", alpha=0.10)
        axes[0].annotate(
            "range spanned by\ntfidf → libsize",
            xy=(sub["distance"].iloc[1], (lo + hi) / 2),
            fontsize=6.5, color="#c62828", ha="center",
        )
        axes[0].set_ylim(lo - 0.02, hi + 0.02)
    axes[0].set_title("Cell-type agreement", fontsize=9)

    axes[1].plot(sub["distance"], sub["n_clusters"], marker="s", color=NULL_SHIFT, lw=1.8)
    axes[1].set_ylabel("mean clusters per dataset")
    axes[1].set_ylim(sub["n_clusters"].mean() - 3, sub["n_clusters"].mean() + 3)
    axes[1].set_title("Cluster count", fontsize=9)

    if yield_df is not None:
        y = yield_df[yield_df["run"].isin(TRIM_DISTANCE)].copy()
        if not y.empty:
            y["distance"] = y["run"].map(TRIM_DISTANCE)
            y = y.sort_values("distance")
            axes[2].plot(y["distance"], y["n_pas"], marker="^", color="#2e7d32", lw=1.8)
            axes[2].set_ylabel("PAS called")
            axes[2].set_title("PAS yield\n(see reproducibility caveat)", fontsize=9)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks(sorted(TRIM_DISTANCE.values()))
        ax.set_xticklabels([str(v) for v in sorted(TRIM_DISTANCE.values())], fontsize=7)
        ax.set_xlabel("max_gene_distance (bp)")
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.suptitle(
        "A 10× change in the annotation window moves neither the clustering nor the biology",
        fontsize=9.5, y=0.99,
    )
    _save(fig, out_dir, "28_trim_sensitivity")


def fig_reannotate_reproducibility(df: pd.DataFrame, out_dir: Path) -> None:
    """Branches with identical trim parameters that produced different output.

    Rendered as a table-figure because the evidence *is* the table: same
    declared parameters, three distinct md5s, grouped by write time.
    """
    need = {"branch", "pas_gene_md5", "pas_gene_rows", "n_genes", "mtime"}
    if not need <= set(df.columns):
        return
    df = df.copy()
    order = df.groupby("pas_gene_md5")["mtime"].min().sort_values().index.tolist()
    df["_g"] = df["pas_gene_md5"].map({m: i for i, m in enumerate(order)})
    df = df.sort_values(["_g", "branch"])
    palette = ["#1f4e79", "#e07b39", "#2e7d32", "#8e24aa"]

    fig, ax = plt.subplots(figsize=(9.2, 0.42 * len(df) + 1.9))
    ax.axis("off")
    cols = ["branch", "trim params", "clustering", "pas_gene md5",
            "rows", "genes", "written"]
    cells, colours = [], []
    for _, r in df.iterrows():
        c = palette[int(r["_g"]) % len(palette)]
        cells.append([
            r["branch"],
            f"{r.get('declared_max_gene_distance','')}/"
            f"{r.get('declared_utr_multiplier','')}/"
            f"{r.get('declared_include_extended','')}",
            f"{str(r.get('cluster_method','')).replace('leiden_','')}"
            f" r{r.get('resolution','')} nn{r.get('n_neighbors','')}",
            r["pas_gene_md5"],
            f"{int(r['pas_gene_rows']):,}",
            f"{int(r['n_genes']):,}",
            r["mtime"],
        ])
        colours.append([c] * len(cols))
    table = ax.table(cellText=cells, colLabels=cols, loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    table.scale(1, 1.5)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#f0f0f0")
        elif col in (3, 4, 5):
            cell.set_text_props(color=colours[row - 1][col], weight="bold")
    n_groups = df["pas_gene_md5"].nunique()
    ax.set_title(
        "Identical declared trim parameters, "
        f"{n_groups} different PAS→gene assignment tables\n"
        "(colour = md5 group; groups track write time, not parameters — "
        "a concurrency race, not a parameter effect)",
        fontsize=10, pad=14,
    )
    _save(fig, out_dir, "30_reannotate_reproducibility")


# ---------------------------------------------------------------------------
# F5 — switch biology honesty panels
# ---------------------------------------------------------------------------


def fig_trend_confound(df: pd.DataFrame, out_dir: Path) -> None:
    """Slope with vs without the detection-rate confound, and stage counts."""
    if "slope" not in df or "slope_informative" not in df:
        return
    sub = df.dropna(subset=["slope", "slope_informative"])
    if sub.empty:
        return
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.6))

    ok = sub["interpretable"].astype(bool) if "interpretable" in sub else np.ones(len(sub), bool)
    ax.scatter(sub.loc[ok, "slope"], sub.loc[ok, "slope_informative"], s=45,
               color=OBS, label="≥3 stages", zorder=3)
    ax.scatter(sub.loc[~ok, "slope"], sub.loc[~ok, "slope_informative"], s=45,
               facecolors="none", edgecolors="#c62828",
               label="2 stages (uninterpretable)", zorder=3)
    lim = np.nanmax(np.abs(np.r_[sub["slope"], sub["slope_informative"]])) * 1.15
    ax.plot([-lim, lim], [-lim, lim], color="#999999", lw=0.9, ls="--")
    ax.axhline(0, color="#cccccc", lw=0.8)
    ax.axvline(0, color="#cccccc", lw=0.8)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("PDUI slope (all gene–cell pairs)")
    ax.set_ylabel("PDUI slope (informative pairs only)")
    ax.set_title("Does 3'UTR shortening survive the\ndetection-rate control?", fontsize=9.5)
    ax.legend(fontsize=7)

    counts = df["n_stages"].value_counts().sort_index()
    ax2.bar(counts.index.astype(str), counts.values,
            color=["#c62828" if int(i) < 3 else OBS for i in counts.index])
    ax2.set_xlabel("stages available for the trend fit")
    ax2.set_ylabel("cell types")
    ax2.set_title("Spearman ρ is ±1 by construction\nwherever only 2 stages exist", fontsize=9.5)
    for x, v in zip(range(len(counts)), counts.values):
        ax2.text(x, v + 0.1, str(v), ha="center", fontsize=8)
    _save(fig, out_dir, "34_trend_depth_confound")


def fig_fisher_power(df: pd.DataFrame, out_dir: Path) -> None:
    """n_cells vs n_significant — the pseudoreplication signature."""
    sub = df.dropna(subset=["n_cells_total", "n_significant"])
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    size = 30
    if "median_abs_delta_sig" in sub and sub["median_abs_delta_sig"].notna().any():
        d = sub["median_abs_delta_sig"].fillna(sub["median_abs_delta_sig"].median())
        size = 20 + 220 * (d - d.min()) / max(d.max() - d.min(), 1e-9)
    ax.scatter(sub["n_cells_total"], sub["n_significant"], s=size, alpha=0.75,
               color=OBS, edgecolors="white", lw=0.5)
    ax.set_xscale("log")
    ax.set_yscale("symlog")
    ax.set_xlabel("cells in the contrast (log)")
    ax.set_ylabel("significant PAS at FDR 0.05 (symlog)")
    ax.set_title("Fisher significance tracks cell count, not effect size\n"
                 "(point size = median |Δ proportion| of the hits)", fontsize=9.5)
    _save(fig, out_dir, "35_fisher_pseudoreplication")


def fig_recurrent(df: pd.DataFrame, out_dir: Path, top_n: int = 20) -> None:
    """Recurrent vs private switch genes."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.8),
                                  gridspec_kw={"width_ratios": [1, 1.3]})
    counts = df["n_celltypes"].value_counts().sort_index()
    private = float((df["n_celltypes"] == 1).mean())
    ax.bar(counts.index.astype(str), counts.values, color=OBS)
    ax.set_xlabel("cell types the gene switches in")
    ax.set_ylabel("genes")
    # Real cell-type-specific APA would be concentrated at 1. It is not: the
    # distribution has a long flat tail out to all 24 cell types.
    ax.set_title(
        f"Switch calls are not sparse: only {private:.0%} are private,\n"
        "with a long tail out to all cell types",
        fontsize=9.5,
    )

    top = df.head(top_n).iloc[::-1]
    colors = ["#2e7d32" if c else "#c62828" for c in top["consistent_direction"]]
    ax2.barh(top["gene_id"].astype(str), top["n_celltypes"], color=colors)
    ax2.set_xlabel("cell types with a significant switch")
    ax2.tick_params(axis="y", labelsize=6.5)
    ax2.set_title("Top recurrent switch genes\n(green = consistent direction)", fontsize=9.5)
    _save(fig, out_dir, "36_recurrent_switches")


def fig_pdui_heatmap(df: pd.DataFrame, out_dir: Path) -> None:
    """Cell type x stage PDUI, greying out cell types with <3 stages."""
    if "stage_index" not in df:
        return
    piv = df.pivot_table(index="celltype", columns="stage", values="mean_pdui")
    order = (
        df.sort_values("stage_index")["stage"].drop_duplicates().tolist()
    )
    piv = piv[[c for c in order if c in piv.columns]]
    n_stages = df.groupby("celltype")["stage"].nunique()
    piv = piv.loc[n_stages.sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(4.8, 0.26 * len(piv) + 1.6))
    im = ax.imshow(piv.values, aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels(piv.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels(
        [
            f"{c[:38]}  ({n_stages[c]}st)" + ("  ⚠" if n_stages[c] < 3 else "")
            for c in piv.index
        ],
        fontsize=6.5,
    )
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="mean PDUI")
    ax.set_title("Mean PDUI by cell type × stage\n(⚠ = <3 stages, trend uninterpretable)",
                 fontsize=9.5)
    _save(fig, out_dir, "33_pdui_celltype_stage")


# ---------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_dir", required=True, type=Path)
    ap.add_argument("--out", dest="out_dir", required=True, type=Path)
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    i, o = args.in_dir, args.out_dir

    if (d := _load(i, "strategy_atlas_null_control.tsv")) is not None:
        fig_atlas_null(d, o)
    if (d := _load(i, "strategy_pas_overlap.tsv")) is not None:
        fig_pas_overlap(d, o)

    strat_yield = _load(i, "strategy_pas_yield.tsv")
    strat_clust = _load(i, "strategy_clustering_summary.tsv")
    if strat_yield is not None and strat_clust is not None:
        fig_yield_vs_biology(strat_yield, strat_clust, o)

    branch = _load(i, "branch_clustering_summary.tsv")
    branch_yield = _load(i, "branch_pas_yield.tsv")
    if branch is not None:
        fig_knob_ranking(branch, strat_clust, o)
        fig_resolution_curve(branch, o)
        fig_trim_sensitivity(branch, branch_yield, o)
    if (d := _load(i, "branch_celltype_shift.tsv")) is not None:
        fig_celltype_shift(d, o)
    if (d := _load(i, "reannotate_reproducibility.tsv")) is not None:
        fig_reannotate_reproducibility(d, o)

    if (d := _load(i, "switch_trend_recomputed.tsv")) is not None:
        fig_trend_confound(d, o)
    if (d := _load(i, "switch_depth_by_stage.tsv")) is not None:
        fig_pdui_heatmap(d, o)
    if (d := _load(i, "switch_fisher_power.tsv")) is not None:
        fig_fisher_power(d, o)
    if (d := _load(i, "switch_recurrent_genes.tsv")) is not None:
        fig_recurrent(d, o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
