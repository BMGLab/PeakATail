#!/usr/bin/env python
"""Single-cell clustering figures for the corrected sweep: UMAP, Sankey, sizes.

The corrected figure set (01-14) covered strategy and switch analyses but dropped
the single-cell visualisations entirely, which left the report without a UMAP or
a Sankey -- the two figures a reader of a single-cell paper looks for first.

Everything here is built from real run artifacts, exported by
``_server_scripts/s6_clustering_viz_data.py`` on the analysis host:

  <run>/07_clustering/<dataset>/clusters.h5ad             -- obsm['X_umap'], obs['leiden']
  B1_cohort_full/B2_gex_celltyping/<ds>_pas_labeled.h5ad  -- obs['celltype']

Cell types are the real signature-scored labels (28 of them), cleaned to display
names -- ``CELL_TYPES_WANSLEEBEN_HOGAN_2013:MACROPHAGE_M2`` renders as
``Macrophage M2``.  No fixture, mock or synthetic data is used anywhere.

One honesty constraint runs through every figure: **only ~55% of clustered cells
carry a cell-type label.** B2 labelled a subset, so an unlabelled bucket exists in
every panel. It is drawn explicitly in grey rather than dropped, because silently
discarding 45% of cells would make the PAS-to-cell-type mapping look far cleaner
than it is.

Figures, each covering every experiment rather than the cohort alone:

  15_umap_cohort_stages     cohort UMAP, one dataset per stage, cluster + cell type
  16_umap_clustering_grid   the clustering variants on one shared dataset
  17_sankey_clustering_grid PAS cluster -> cell type flow, per clustering variant
  18_sankey_peak_strategies same flow, per peak-calling strategy
  19_cluster_sizes          cluster size distributions across all 19 experiments

Usage::

    python reports/generate_clustering_figures.py [--extra DIR] [--out DIR]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
EXTRA = HERE / "cumulative_analysis" / "extra"
OUT = HERE / "figures" / "corrected"

UNLABELLED = "unlabelled"
#: clustering variants, in the order the report discusses them
CLUSTER_VARIANTS = [
    ("reannotate/A3_res0.5", "res 0.5"),
    ("reannotate/A2_trim_default", "res 1.0 (default)"),
    ("reannotate/A3_res2.0", "res 2.0"),
    ("reannotate/A3_nn15", "n_neighbors 15"),
    ("reannotate/A3_nn50", "n_neighbors 50"),
    ("reannotate/A3_libsize", "leiden_libsize"),
]
PEAK_RUNS = [
    ("grid/lg_annotate", "lambda_gradient"),
    ("grid/lp_annotate", "lambda_poisson"),
    ("grid/si_annotate", "sierra_iterative"),
    ("grid/lg_ip_filter", "lambda_gradient\n(ip=filter)"),
]

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 200, "savefig.bbox": "tight",
    "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
    "svg.hashsalt": "peakatail-clustering",
})

_written: list[str] = []
_skipped: list[str] = []


def save(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png")
    fig.savefig(OUT / f"{name}.svg", metadata={"Date": None})
    plt.close(fig)
    _written.append(name)
    print(f"  wrote {name}.png/.svg")


def load(name: str, **kw) -> pd.DataFrame | None:
    p = EXTRA / name
    if not p.exists():
        print(f"  SKIP: missing {p}")
        return None
    return pd.read_csv(p, sep="\t", **kw)


def celltype_palette(labels: list[str]) -> dict[str, tuple]:
    """Stable colours across every figure; unlabelled always neutral grey."""
    named = sorted(x for x in labels if x != UNLABELLED)
    cmap = plt.get_cmap("tab20")
    cmap2 = plt.get_cmap("tab20b")
    pal = {}
    for i, lab in enumerate(named):
        pal[lab] = cmap(i % 20) if i < 20 else cmap2((i - 20) % 20)
    pal[UNLABELLED] = (0.82, 0.82, 0.82, 1.0)
    return pal


def scatter_umap(ax, df: pd.DataFrame, colour_by: str, pal: dict, title: str,
                 legend: bool = False, max_legend: int = 14) -> None:
    """Unlabelled cells go down first, faint, so labelled types stay visible.

    Roughly 45% of clustered cells have no cell-type call. Dropping them would
    overstate how cleanly PAS clusters map to cell types, so they are drawn --
    just underneath, and never in the legend's ranked slots.
    """
    unl = df[df[colour_by] == UNLABELLED]
    lab_df = df[df[colour_by] != UNLABELLED]
    if len(unl):
        ax.scatter(unl["umap1"], unl["umap2"], s=1.6, alpha=0.22,
                   c=[pal[UNLABELLED]], linewidths=0, label=None)
    order = lab_df[colour_by].value_counts()
    top = list(order.index)[:max_legend]
    for lab in list(order.index)[::-1]:            # rare types drawn on top
        sub = lab_df[lab_df[colour_by] == lab]
        ax.scatter(sub["umap1"], sub["umap2"], s=3.4, alpha=0.9,
                   c=[pal.get(lab, (0.6, 0.6, 0.6, 1.0))], linewidths=0,
                   label=lab if lab in top else None)
    ax.set_xticks([]); ax.set_yticks([])
    ax.grid(False)
    ax.set_title(title, fontsize=9)
    if legend:
        h, l = ax.get_legend_handles_labels()
        rank = {lab: i for i, lab in enumerate(top)}
        pairs = sorted(zip(l, h), key=lambda t: rank.get(t[0], 99))
        if pairs:
            ax.legend([x[1] for x in pairs], [x[0] for x in pairs],
                      loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=6,
                      markerscale=3.5, frameon=False,
                      title="cell type (grey = unlabelled)", title_fontsize=6)


# ===========================================================================
# 15. cohort UMAP, one dataset per stage
# ===========================================================================
def fig_umap_cohort() -> None:
    u = load("umap_cells.tsv.gz", compression="gzip")
    if u is None:
        _skipped.append("15_umap_cohort_stages"); return
    u = u[u["experiment"] == "B1_cohort_full"]
    if u.empty:
        _skipped.append("15_umap_cohort_stages"); return
    order = ["Normal", "StageIA", "StageIVprimary", "MetBone"]
    datasets = sorted(u["dataset"].unique(),
                      key=lambda d: next((i for i, s in enumerate(order) if s in d), 99))
    pal = celltype_palette(sorted(u["celltype"].unique()))

    fig, axes = plt.subplots(2, len(datasets), figsize=(3.5 * len(datasets), 7.4),
                             squeeze=False, layout="constrained")
    for j, ds in enumerate(datasets):
        d = u[u["dataset"] == ds]
        lp = {k: plt.get_cmap("tab20")(i % 20)
              for i, k in enumerate(sorted(d["leiden"].unique(), key=lambda x: int(x)))}
        scatter_umap(axes[0][j], d, "leiden", lp,
                     f"{ds}\n{len(d):,} cells · {d['leiden'].nunique()} PAS clusters")
        lab = d[d["celltype"] != UNLABELLED]
        scatter_umap(axes[1][j], d, "celltype", pal,
                     f"cell type · {len(lab):,}/{len(d):,} labelled "
                     f"({len(lab)/len(d):.0%})")
    axes[0][0].set_ylabel("coloured by PAS-Leiden cluster", fontsize=9)
    axes[1][0].set_ylabel("coloured by GEX cell type", fontsize=9)

    # One figure-level legend spanning every panel: a per-axes legend would show
    # only the last dataset's types and silently omit the rest.
    from matplotlib.lines import Line2D
    present = (u[u["celltype"] != UNLABELLED]["celltype"].value_counts())
    handles = [Line2D([], [], marker="o", ls="", ms=5, color=pal[k]) for k in present.index]
    handles.append(Line2D([], [], marker="o", ls="", ms=5, color=pal[UNLABELLED]))
    fig.legend(handles, list(present.index) + ["unlabelled"],
               loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=6.5,
               frameon=False, title="GEX cell type", title_fontsize=7)
    fig.suptitle(
        "PAS-space UMAP across disease stage (cohort run) — same cells, coloured by "
        "PAS cluster (top) and by expression-defined cell type (bottom)", fontsize=10)
    save(fig, "15_umap_cohort_stages")


# ===========================================================================
# 16. clustering variants on one shared dataset
# ===========================================================================
def fig_umap_variants() -> None:
    u = load("umap_cells.tsv.gz", compression="gzip")
    if u is None:
        _skipped.append("16_umap_clustering_grid"); return
    avail = [(e, lab) for e, lab in CLUSTER_VARIANTS if e in set(u["experiment"])]
    if not avail:
        _skipped.append("16_umap_clustering_grid"); return
    pal = celltype_palette(sorted(u["celltype"].unique()))

    fig, axes = plt.subplots(2, len(avail), figsize=(3.1 * len(avail), 7.0),
                             squeeze=False, layout="constrained")
    for j, (exp, lab) in enumerate(avail):
        d = u[u["experiment"] == exp]
        ds = d["dataset"].iloc[0]
        lp = {k: plt.get_cmap("tab20")(i % 20)
              for i, k in enumerate(sorted(d["leiden"].unique(), key=lambda x: int(x)))}
        scatter_umap(axes[0][j], d, "leiden", lp,
                     f"{lab}\n{d['leiden'].nunique()} clusters")
        scatter_umap(axes[1][j], d, "celltype", pal, "", legend=(j == len(avail) - 1))
    axes[0][0].set_ylabel("PAS-Leiden cluster", fontsize=9)
    axes[1][0].set_ylabel("GEX cell type", fontsize=9)
    ds = u[u["experiment"] == avail[0][0]]["dataset"].iloc[0]
    fig.suptitle(
        f"Clustering variants on one shared dataset ({ds}) — the embedding is stable; "
        "only how finely it is cut changes", fontsize=10)
    save(fig, "16_umap_clustering_grid")


# ===========================================================================
# Sankey: PAS cluster -> cell type
# ===========================================================================
def _sankey(ax, counts: pd.DataFrame, pal: dict, title: str,
            min_frac: float = 0.01, show_labels: bool = True) -> None:
    """One Sankey panel. Left nodes = PAS clusters, right nodes = cell types."""
    tot = counts["n_cells"].sum()
    if tot == 0:
        ax.set_axis_off(); return
    counts = counts[counts["n_cells"] >= max(1, min_frac * tot / 10)]
    clusters = (counts.groupby("leiden")["n_cells"].sum()
                .sort_values(ascending=False))
    types = (counts.groupby("celltype")["n_cells"].sum()
             .sort_values(ascending=False))

    gap = 0.012
    def layout(series):
        pos, y = {}, 0.0
        span = 1.0 - gap * max(len(series) - 1, 0)
        for k, v in series.items():
            h = span * v / series.sum()
            pos[k] = (y, h)
            y += h + gap
        return pos
    lpos, rpos = layout(clusters), layout(types)
    lcur = {k: v[0] for k, v in lpos.items()}
    rcur = {k: v[0] for k, v in rpos.items()}

    for k, (y0, h) in lpos.items():
        ax.add_patch(Rectangle((0.0, y0), 0.045, h, color="#555555"))
        if show_labels and h > 0.028:
            ax.text(-0.012, y0 + h / 2, k, ha="right", va="center", fontsize=6)
    for k, (y0, h) in rpos.items():
        ax.add_patch(Rectangle((0.955, y0), 0.045, h, color=pal.get(k, "#999999")))
        if show_labels and h > 0.022:
            ax.text(1.012, y0 + h / 2, k, ha="left", va="center", fontsize=6)

    for _, r in counts.sort_values("n_cells", ascending=False).iterrows():
        cl, ct, n = r["leiden"], r["celltype"], r["n_cells"]
        if cl not in lpos or ct not in rpos:
            continue
        hl = lpos[cl][1] * n / clusters[cl]
        hr = rpos[ct][1] * n / types[ct]
        y0, y1 = lcur[cl], rcur[ct]
        lcur[cl] += hl; rcur[ct] += hr
        t = np.linspace(0, 1, 60)
        smooth = t * t * (3 - 2 * t)
        top = (y0 + hl) + ((y1 + hr) - (y0 + hl)) * smooth
        bot = y0 + (y1 - y0) * smooth
        ax.fill_between(0.045 + t * 0.91, bot, top,
                        color=pal.get(ct, "#999999"), alpha=0.55, linewidth=0)
    ax.set_xlim(-0.16, 1.16); ax.set_ylim(-0.02, 1.02)
    ax.set_axis_off()
    ax.set_title(title, fontsize=8.5)


def ref_dataset(sub: pd.DataFrame) -> str:
    """The dataset used for a per-experiment panel: the shared reference if the
    experiment has it, else its largest dataset."""
    have = set(sub["dataset"])
    for cand in ("GSM3516666-Normal", "GSM3516675-Normal"):
        if cand in have:
            return cand
    return sub.groupby("dataset")["n_cells"].sum().idxmax()


def _sankey_panel(spec: list[tuple[str, str]], name: str, suptitle: str) -> None:
    c = load("cluster_celltype_counts.tsv")
    if c is None:
        _skipped.append(name); return
    avail = [(e, lab) for e, lab in spec if e in set(c["experiment"])]
    if not avail:
        _skipped.append(name); return
    pal = celltype_palette(sorted(c["celltype"].unique()))
    ncol = min(3, len(avail))
    nrow = int(np.ceil(len(avail) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 5.2 * nrow),
                             squeeze=False, layout="constrained")
    for i, (exp, lab) in enumerate(avail):
        ax = axes[i // ncol][i % ncol]
        # One dataset, never pooled: `leiden` ids are assigned per dataset, so
        # cluster "0" in one GSM has nothing to do with cluster "0" in another.
        # Pooling would merge unrelated clusters into a single node.
        sub = c[c["experiment"] == exp]
        ds = ref_dataset(sub)
        d = (sub[sub["dataset"] == ds]
             .groupby(["leiden", "celltype"], as_index=False)["n_cells"].sum())
        rate = sub[sub["dataset"] == ds]["celltype_match_rate"].mean()
        _sankey(ax, d, pal, f"{lab}  ·  {ds}\n{d['leiden'].nunique()} PAS clusters → "
                            f"{d[d.celltype != UNLABELLED]['celltype'].nunique()} cell types "
                            f"· {rate:.0%} labelled")
    for j in range(len(avail), nrow * ncol):
        axes[j // ncol][j % ncol].set_axis_off()
    fig.suptitle(suptitle, fontsize=10)
    save(fig, name)


def fig_sankey_variants() -> None:
    _sankey_panel(
        CLUSTER_VARIANTS, "17_sankey_clustering_grid",
        "PAS cluster → GEX cell type, per clustering variant, on one shared dataset — "
        "grey is the unlabelled fraction, kept visible rather than dropped")


def fig_sankey_strategies() -> None:
    _sankey_panel(
        PEAK_RUNS, "18_sankey_peak_strategies",
        "PAS cluster → GEX cell type, per peak-calling strategy, on one shared dataset")


# ===========================================================================
# 19. cluster sizes across every experiment
# ===========================================================================
def fig_cluster_sizes() -> None:
    c = load("cluster_celltype_counts.tsv")
    if c is None:
        _skipped.append("19_cluster_sizes"); return
    sizes = c.groupby(["experiment", "dataset", "leiden"], as_index=False)["n_cells"].sum()
    per = sizes.groupby(["experiment", "dataset"], as_index=False).agg(
        n_clusters=("leiden", "nunique"), n_cells=("n_cells", "sum"),
        largest=("n_cells", "max"))
    per["largest_frac"] = per["largest"] / per["n_cells"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.0), layout="constrained")

    ax = axes[0]
    order = per.groupby("experiment")["n_clusters"].median().sort_values()
    data = [per[per["experiment"] == e]["n_clusters"] for e in order.index]
    bp = ax.boxplot(data, vert=False, tick_labels=[e.split("/")[-1] for e in order.index],
                    showmeans=True, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("#4C72B0"); patch.set_alpha(0.6)
    ax.set_xlabel("clusters per dataset")
    ax.tick_params(axis="y", labelsize=6.5)
    ax.set_title(f"a  Cluster count, all {per['experiment'].nunique()} experiments")

    ax = axes[1]
    for e, lab in CLUSTER_VARIANTS:
        s = sizes[sizes["experiment"] == e]
        if s.empty:
            continue
        v = np.sort(s["n_cells"].values)[::-1]
        ax.plot(np.arange(1, len(v) + 1), v, marker="o", ms=2.5, lw=1, label=lab)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("cluster rank"); ax.set_ylabel("cells in cluster")
    ax.legend(fontsize=6.5)
    ax.set_title("b  Cluster size distribution, clustering variants")

    ax = axes[2]
    ax.scatter(per["n_clusters"], per["largest_frac"], s=14, alpha=0.6, color="#DD8452")
    ax.set_xlabel("clusters per dataset")
    ax.set_ylabel("share of cells in the largest cluster")
    rho = per["n_clusters"].corr(per["largest_frac"], method="spearman")
    ax.set_title(f"c  Fragmentation vs cluster count (ρ={rho:+.2f})")

    fig.suptitle("PAS-space cluster structure across every experiment in the sweep",
                 fontsize=10)
    save(fig, "19_cluster_sizes")


FIGURES = [fig_umap_cohort, fig_umap_variants, fig_sankey_variants,
           fig_sankey_strategies, fig_cluster_sizes]


def main(argv: list[str] | None = None) -> int:
    global EXTRA, OUT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--extra", type=Path, default=EXTRA)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args(argv)
    EXTRA, OUT = a.extra, a.out
    print(f"extra: {EXTRA}\nout  : {OUT}\n")
    for fn in FIGURES:
        print(f"[{fn.__name__}]")
        try:
            fn()
        except Exception as exc:
            print(f"  ERROR {fn.__name__}: {type(exc).__name__}: {exc}")
            _skipped.append(fn.__name__)
    print(f"\nwrote {len(_written)} figures -> {OUT}")
    if _skipped:
        print(f"skipped {len(_skipped)}: {', '.join(_skipped)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
