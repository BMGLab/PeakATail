#!/usr/bin/env python
"""Figures 37–41 — per-cell-type switch deep dive.

Consumes the tables from ``reports/_server_scripts/s7_switch_deepdive.py`` and
renders what the aggregate panels leave out: both differential strategies side
by side per cell type, the distribution of switch effects (volcanoes), the named
strongest switches, both length metrics per stage, and pathway enrichment of the
switch-gene sets.

    python reports/generate_switch_deepdive_figures.py \
        --in reports/cumulative_analysis/deepdive \
        --out reports/figures/corrected

Slots 37–41. As with 24–36, a missing input table is logged and skipped.
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

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
    "axes.spines.top": False, "axes.spines.right": False,
})

FISHER = "#1f4e79"
NB = "#c62828"
GREY = "#9e9e9e"
GREEN = "#2e7d32"


def short(name: str) -> str:
    """Trim the signature prefixes so 24 cell-type labels fit on an axis."""
    return (name.replace("CELL_TYPES_WANSLEEBEN_HOGAN_2013_", "")
                .replace("lung_epithelial_lineage_signatures_", "LE:")
                .replace("stem_cell_signatures_merged_", "SC:"))


def _load(in_dir: Path, name: str) -> Optional[pd.DataFrame]:
    p = in_dir / name
    if not p.exists():
        logger.warning("skip: %s not found", name)
        return None
    df = pd.read_csv(p, sep="\t", index_col=False)
    return df if not df.empty else None


def _save(fig, out_dir: Path, stem: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"{stem}.{ext}")
        logger.info("wrote %s.%s", stem, ext)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 37 — both differential strategies, per cell type
# ---------------------------------------------------------------------------


def fig_diff_strategies(df: pd.DataFrame, out_dir: Path) -> None:
    d = df[df["contrast"] == "Normal_vs_StageI"].copy()
    if d.empty:
        d = df[df["contrast"] == df["contrast"].mode().iloc[0]].copy()
    contrast = d["contrast"].iloc[0]
    d["ct"] = d["celltype"].map(short)
    d = d.sort_values("frac_sig_fisher", ascending=False)
    y = np.arange(len(d))

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 0.30 * len(d) + 2.2),
                             gridspec_kw={"width_ratios": [1.5, 1.1, 1.0]})

    ax = axes[0]
    ax.barh(y + 0.20, d["frac_sig_fisher"], 0.38, color=FISHER, label="fisher (q≤0.05)")
    ax.barh(y + 0.20, d["frac_sig_fisher_large"], 0.38, color="#7fb3d5",
            label="fisher, also |Δp|≥0.30")
    ax.barh(y - 0.20, d["frac_sig_nb_multi"], 0.38, color=NB, label="nb_multi (q≤0.05)")
    ax.set_yticks(y)
    ax.set_yticklabels(d["ct"], fontsize=6.5)
    ax.set_xlabel("fraction of tested PAS called significant")
    ax.set_xlim(0, 1.05)
    ax.axvline(1.0, color="k", lw=0.7, ls=":")
    ax.set_title("Significance rate per strategy", fontsize=9.5)
    ax.legend(fontsize=6.5, loc="lower right", framealpha=0.9)
    ax.invert_yaxis()

    ax = axes[1]
    ax.barh(y + 0.20, d["n_tests_fisher"], 0.38, color=FISHER, label="fisher")
    ax.barh(y - 0.20, d["n_tests_nb_multi"], 0.38, color=NB, label="nb_multi")
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xscale("log")
    ax.set_xlabel("PAS tested (log)")
    ax.set_title("How many PAS each\nstrategy even tests", fontsize=9.5)
    ax.legend(fontsize=6.5, loc="lower right")
    ax.invert_yaxis()

    ax = axes[2]
    ax.barh(y, d["jaccard_sig"], 0.6, color=GREY)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel("Jaccard of significant sets")
    ax.set_xlim(0, max(0.05, float(d["jaccard_sig"].max()) * 1.2))
    ax.set_title("Agreement between\nthe two strategies", fontsize=9.5)
    ax.invert_yaxis()

    fig.suptitle(
        f"Both differential strategies, per cell type ({contrast}). "
        "nb_multi calls ~100% of the ~3% of PAS it tests significant; "
        "the two agree on ~1%.",
        fontsize=9.5, y=1.005,
    )
    _save(fig, out_dir, "37_diff_strategies_by_celltype")


# ---------------------------------------------------------------------------
# 38 — volcano grid
# ---------------------------------------------------------------------------


def fig_volcano_grid(df: pd.DataFrame, out_dir: Path,
                     strategy: Optional[pd.DataFrame] = None,
                     contrast: str = "Normal_vs_StageI",
                     max_panels: int = 12, seed: int = 0) -> None:
    """Volcano per cell type, drawn at a *single* sampling rate.

    The harvest keeps every large-effect significant point but subsamples the
    rest, so plotting the file as-is shows a hole at Δp≈0 that is an artefact of
    that sampling, not a property of the data. Here the retained large-effect
    points are thinned to the same rate as the background, using the true
    per-panel counts, so the rendered density is faithful.
    """
    rng = np.random.default_rng(seed)
    d = df[df["contrast"] == contrast]
    if d.empty:
        contrast = df["contrast"].mode().iloc[0]
        d = df[df["contrast"] == contrast]
    counts = None
    if strategy is not None:
        s = strategy[strategy["contrast"] == contrast]
        counts = s.set_index("celltype")[["n_tests_fisher", "n_sig_fisher_large"]]
    order = d.groupby("celltype").size().sort_values(ascending=False).index[:max_panels]
    ncol = 4
    nrow = int(np.ceil(len(order) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.5 * nrow),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, ct in zip(axes, order):
        s = d[d["celltype"] == ct]
        big = s["is_large_sig"].astype(bool)
        small, large = s[~big], s[big]
        if counts is not None and ct in counts.index:
            n_tests = float(counts.loc[ct, "n_tests_fisher"])
            n_large = float(counts.loc[ct, "n_sig_fisher_large"])
            n_small_true = max(n_tests - n_large, 1.0)
            rate = min(1.0, len(small) / n_small_true)
            keep = max(int(round(n_large * rate)), 1)
            if keep < len(large):
                large = large.iloc[rng.choice(len(large), keep, replace=False)]
        ax.scatter(small["delta_proportion"], small["neglog10q"],
                   s=2, color=GREY, alpha=0.35, linewidths=0, rasterized=True)
        ax.scatter(large["delta_proportion"], large["neglog10q"],
                   s=3, color=FISHER, alpha=0.55, linewidths=0, rasterized=True)
        ax.axhline(-np.log10(0.05), color=NB, lw=0.7, ls="--")
        for v in (-0.30, 0.30):
            ax.axvline(v, color=NB, lw=0.7, ls=":")
        ax.set_title(short(ct)[:26], fontsize=7)
        ax.set_xlim(-1.05, 1.05)
    for ax in axes[len(order):]:
        ax.axis("off")
    fig.supxlabel("Δ proportion  (proximal-up →)", fontsize=9)
    fig.supylabel("−log₁₀ q", fontsize=9)
    fig.suptitle(
        f"Per-cell-type differential APA, {contrast}. Dashed = q 0.05, dotted = |Δp| 0.30. "
        "Blue = significant with |Δp|≥0.30, thinned to the background sampling rate "
        "so the density is faithful.",
        fontsize=9.5, y=1.002,
    )
    fig.tight_layout()
    _save(fig, out_dir, "38_volcano_by_celltype")


# ---------------------------------------------------------------------------
# 39 — both length metrics by stage
# ---------------------------------------------------------------------------

STAGES = ["Normal", "StageI", "IVprimary", "Met"]


def fig_length_metrics(shannon: pd.DataFrame, depth: Optional[pd.DataFrame],
                       out_dir: Path) -> None:
    idx = {s: i for i, s in enumerate(STAGES)}
    sh = shannon.copy()
    sh["si"] = sh["stage"].map(idx)
    sh = sh.dropna(subset=["si"])

    npanel = 3 if depth is not None else 2
    fig, axes = plt.subplots(1, npanel, figsize=(4.1 * npanel, 3.4))

    ax = axes[0]
    for ct, g in sh.groupby("celltype"):
        g = g.sort_values("si")
        ax.plot(g["si"], g["mean_entropy"], color=GREY, lw=0.8, alpha=0.5)
    m = sh.groupby("si")["mean_entropy"].mean()
    ax.plot(m.index, m.values, color=FISHER, lw=2.4, marker="o", label="mean")
    ax.set_xticks(range(len(STAGES)))
    ax.set_xticklabels(STAGES, rotation=30, ha="right")
    ax.set_ylabel("Shannon entropy (as reported)")
    lo, hi = sh["frac_uninformative"].min(), sh["frac_uninformative"].max()
    ax.set_title("shannon — as reported", fontsize=9.5)
    # The reported mean sits just above 1.0 because that is the constant every
    # zero-coverage row is assigned; saying so on the figure face stops the
    # near-flat line being read as "entropy is stable".
    ax.annotate(
        f"{lo:.0%}–{hi:.0%} of rows are\nthe constant 1.0 padding\n"
        "(zero-coverage gene–cell pairs)",
        xy=(0.03, 0.03), xycoords="axes fraction", fontsize=6.8, color=NB, va="bottom",
    )
    ax.legend(fontsize=7, loc="upper left")

    ax = axes[1]
    for ct, g in sh.groupby("celltype"):
        g = g.sort_values("si")
        ax.plot(g["si"], g["mean_entropy_informative"], color=GREY, lw=0.8, alpha=0.5)
    m = sh.groupby("si")["mean_entropy_informative"].mean()
    ax.plot(m.index, m.values, color=GREEN, lw=2.4, marker="s", label="mean")
    ax.set_xticks(range(len(STAGES)))
    ax.set_xticklabels(STAGES, rotation=30, ha="right")
    ax.set_ylabel("Shannon entropy (informative rows)")
    ax.set_title("shannon — conditioned on coverage", fontsize=9.5)
    ax.legend(fontsize=7)

    if depth is not None:
        d = depth.copy()
        d["si"] = d["stage"].map(idx)
        d = d.dropna(subset=["si"])
        ax = axes[2]
        m1 = d.groupby("si")["mean_pdui"].mean()
        m2 = d.groupby("si")["mean_pdui_informative"].mean()
        ax.plot(m1.index, m1.values, color=FISHER, lw=2.2, marker="o",
                label="classic PDUI, as reported")
        ax2 = ax.twinx()
        ax2.plot(m2.index, m2.values, color=GREEN, lw=2.2, marker="s", ls="--",
                 label="classic PDUI, informative")
        ax2.grid(False)
        ax.set_xticks(range(len(STAGES)))
        ax.set_xticklabels(STAGES, rotation=30, ha="right")
        ax.set_ylabel("mean PDUI (as reported)", color=FISHER)
        ax2.set_ylabel("mean PDUI (informative)", color=GREEN)
        ax.set_title("classic — the same split", fontsize=9.5)

    fig.suptitle(
        "Length metrics by stage. Only classic-as-reported trends, and conditioning on "
        "coverage removes it; shannon trends in neither form — but 94–97% of its rows "
        "are constant padding, the same defect that invalidates proportion.",
        fontsize=9.5, y=1.02,
    )
    fig.tight_layout()
    _save(fig, out_dir, "39_length_metrics_by_stage")


# ---------------------------------------------------------------------------
# 40 — named top switches
# ---------------------------------------------------------------------------


def fig_top_switches(df: pd.DataFrame, out_dir: Path, top_n: int = 25) -> None:
    d = df.copy()
    d["ct"] = d["celltype"].map(short)
    rec = (d.groupby("gene_name")
             .agg(n_hits=("celltype", "size"),
                  n_celltypes=("celltype", "nunique"),
                  n_contrasts=("contrast", "nunique"),
                  max_delta=("delta_proportion", lambda s: float(s.abs().max())),
                  frac_proximal=("direction", lambda s: float((s == "proximal-up").mean())))
             .sort_values(["n_celltypes", "n_hits"], ascending=False)
             .head(top_n).iloc[::-1])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.4, 0.30 * len(rec) + 1.8),
                                  gridspec_kw={"width_ratios": [1.25, 1.0]})
    ax.barh(rec.index, rec["n_celltypes"], color=FISHER)
    ax.set_xlabel("cell types with a large-effect switch (|Δp| ≥ 0.30, q ≤ 0.05)")
    ax.set_title(f"Top {top_n} recurrent switch genes", fontsize=9.5)
    ax.tick_params(axis="y", labelsize=7)

    # Direction consistency: 0.5 means the same gene is called both ways.
    colors = ["#2e7d32" if (v >= 0.8 or v <= 0.2) else NB for v in rec["frac_proximal"]]
    ax2.barh(rec.index, rec["frac_proximal"], color=colors)
    ax2.axvline(0.5, color="k", lw=0.8, ls="--")
    ax2.set_xlim(0, 1)
    ax2.set_yticklabels([])
    ax2.set_xlabel("fraction of calls that are proximal-up")
    ax2.set_title("Direction consistency\n(red = called both ways)", fontsize=9.5)
    fig.tight_layout()
    _save(fig, out_dir, "40_top_switch_genes")


# ---------------------------------------------------------------------------
# 41 — Hallmark enrichment
# ---------------------------------------------------------------------------


def fig_enrichment(df: pd.DataFrame, out_dir: Path, top_n: int = 12) -> None:
    if "qvalue" not in df.columns:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 0.34 * top_n + 2.0), sharex=True)
    for ax, group in zip(axes, ["recurrent", "private"]):
        d = df[df["set_group"] == group].nsmallest(top_n, "pvalue").iloc[::-1]
        if d.empty:
            ax.axis("off")
            continue
        sig = d["qvalue"] <= 0.05
        ax.barh(d["gene_set"].str.replace("HALLMARK_", "", regex=False),
                d["odds_ratio"], color=[GREEN if s else GREY for s in sig])
        ax.axvline(1.0, color="k", lw=0.8, ls="--")
        ax.set_xlabel("odds ratio vs tested-gene background")
        n_sig = int(sig.sum())
        ax.set_title(f"{group} switch genes (n={int(d['n_query'].iloc[0])})\n"
                     f"{n_sig}/{top_n} shown at FDR ≤ 0.05", fontsize=9.5)
        ax.tick_params(axis="y", labelsize=7)
    fig.suptitle(
        "Hallmark enrichment of switch-gene sets, against the genes actually tested "
        "(not all genes — that background would manufacture enrichment).",
        fontsize=9.5, y=1.01,
    )
    fig.tight_layout()
    _save(fig, out_dir, "41_hallmark_switch_enrichment")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in", dest="in_dir", required=True, type=Path)
    ap.add_argument("--out", dest="out_dir", required=True, type=Path)
    ap.add_argument("--depth", type=Path, default=None,
                    help="switch_depth_by_stage.tsv, for the classic panel of fig 39")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    i, o = a.in_dir, a.out_dir

    if (d := _load(i, "diff_strategy_by_celltype.tsv")) is not None:
        fig_diff_strategies(d, o)
    strat = _load(i, "diff_strategy_by_celltype.tsv")
    if (d := _load(i, "volcano_points.tsv")) is not None:
        fig_volcano_grid(d, o, strategy=strat)
    sh = _load(i, "shannon_by_stage.tsv")
    if sh is not None:
        depth = None
        if a.depth and a.depth.exists():
            depth = pd.read_csv(a.depth, sep="\t", index_col=False)
        fig_length_metrics(sh, depth, o)
    if (d := _load(i, "top_switch_genes.tsv")) is not None:
        fig_top_switches(d, o)
    if (d := _load(i, "hallmark_switch_enrichment.tsv")) is not None:
        fig_enrichment(d, o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
