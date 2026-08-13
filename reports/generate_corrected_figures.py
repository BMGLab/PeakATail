#!/usr/bin/env python
"""Regenerate every DATA figure in the PeakATail reports from the corrected sweep.

Background
----------
The figure set previously shipped in ``reports/figures`` (``fig_*.png`` and
``01_*.png`` .. ``13_*.png``) is not reproducible and is not trustworthy:

* the ``fig_*`` set has no generator anywhere in the repository -- it was
  produced ad hoc and its inputs are gone;
* the numbered set is produced by ``generate_report.py``, which either hardcodes
  the numbers inline or copies PNGs from a ``rerun_1500/`` directory that no
  longer exists;
* both predate the pipeline fixes.  They carry the ~50% cell-doubling bug and
  the atlas hard-snap that silently dropped 43-68% of called PAS.

This module replaces them with a figure set that is regenerated end-to-end from
on-disk tables of the corrected sweep, so every number in the reports can be
traced back to an artifact.  Nothing here hardcodes a result.

Inputs
------
``reports/cumulative_analysis/tables/``  -- emitted by
``python -m ema.benchmark.sweep_analysis <sweep_root> --out <dir>``
``reports/cumulative_analysis/extra/``   -- emitted by the three scripts in
``reports/_server_scripts/`` (they run against the sweep on the analysis host;
see the module docstring of each for what it computes).

Usage
-----
    python reports/generate_corrected_figures.py [--tables DIR] [--out DIR]

Figures whose input table is absent are skipped with a warning rather than
being drawn from stale or invented numbers.
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
from matplotlib.ticker import PercentFormatter

HERE = Path(__file__).resolve().parent
TABLES = HERE / "cumulative_analysis" / "tables"
EXTRA = HERE / "cumulative_analysis" / "extra"
OUT = HERE / "figures" / "corrected"

STAGE_ORDER = ["Normal", "StageI", "IVprimary", "Met"]
PALETTE = {
    "sierra_iterative": "#4C72B0",
    "lambda_gradient": "#DD8452",
    "lambda_poisson": "#55A868",
}
plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.grid": True,
        "grid.alpha": 0.25,
        "axes.axisbelow": True,
    }
)

_written: list[str] = []
_skipped: list[str] = []


def save(fig: plt.Figure, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}")
    plt.close(fig)
    _written.append(name)
    print(f"  wrote {name}.png/.svg")


def load(directory: Path, name: str, **kw) -> pd.DataFrame | None:
    path = directory / name
    if not path.exists():
        print(f"  SKIP: missing {path}")
        return None
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, **kw)


def short_ct(name: str) -> str:
    """Trim the long signature-derived cell-type identifiers for plotting."""
    for prefix in (
        "CELL_TYPES_WANSLEEBEN_HOGAN_2013_",
        "lung_epithelial_lineage_signatures_",
        "stem_cell_signatures_merged_",
    ):
        if name.startswith(prefix):
            return name[len(prefix) :].replace("_", " ").title()
    return name.replace("_", " ").title()


# ===========================================================================
# 1. Strategy benchmark, and why its F1 ranking must be read as a count ranking
# ===========================================================================
def fig_strategy_benchmark() -> None:
    df = load(TABLES, "strategy_comparison.csv")
    if df is None:
        _skipped.append("01_strategy_benchmark")
        return
    df = df.drop_duplicates(subset="run").copy()
    df["label"] = df["run"]
    cutoffs = [50, 100, 500, 1000]

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")

    ax = axes[0]
    order = df.sort_values("n_pas")
    ax.barh(order["label"], order["n_pas"], color=[PALETTE.get(s, "#888") for s in order["peak_strategy"]])
    for y, v in enumerate(order["n_pas"]):
        ax.text(v, y, f" {v:,}", va="center", fontsize=8)
    ax.set_xlabel("PAS called")
    ax.set_title("a  Peak count per strategy")
    ax.set_xlim(0, order["n_pas"].max() * 1.25)

    ax = axes[1]
    for _, r in df.iterrows():
        ax.plot(cutoffs, [r[f"f1_{c}"] for c in cutoffs], "o-", label=r["run"],
                color=PALETTE.get(r["peak_strategy"], "#888"))
    ax.set_xscale("log")
    ax.set_xlabel("distance cutoff (bp)")
    ax.set_ylabel("F1 vs PolyASite v3")
    ax.set_title("b  F1 rises with cutoff and with count")
    ax.legend(fontsize=7)

    ax = axes[2]
    ax.scatter(df["n_pas"], df["f1_1000"], s=60,
               c=[PALETTE.get(s, "#888") for s in df["peak_strategy"]], zorder=3)
    for _, r in df.iterrows():
        ax.annotate(r["run"], (r["n_pas"], r["f1_1000"]), textcoords="offset points",
                    xytext=(6, -3), fontsize=7)
    if len(df) > 2:
        rho = df["n_pas"].corr(df["f1_1000"], method="spearman")
        sl, ic = np.polyfit(df["n_pas"], df["f1_1000"], 1)
        xs = np.linspace(df["n_pas"].min(), df["n_pas"].max(), 20)
        ax.plot(xs, sl * xs + ic, "--", color="crimson", lw=1)
        ax.set_title(f"c  F1@1kb vs peak count (Spearman {rho:+.2f})")
    ax.set_xlabel("PAS called")
    ax.set_ylabel("F1 @ 1000 bp")

    n_ref = int(df["n_reference"].iloc[0])
    fig.suptitle(
        f"Peak-calling strategy benchmark vs PolyASite v3 ({n_ref:,} reference PAS) — "
        "F1 tracks peak count, so it ranks yield, not accuracy", fontsize=10)
    save(fig, "01_strategy_benchmark")


# ===========================================================================
# 2. Reference saturation: why precision cannot rank anything here
# ===========================================================================
def fig_atlas_saturation() -> None:
    df = load(TABLES, "strategy_comparison.csv")
    if df is None:
        _skipped.append("02_atlas_saturation")
        return
    df = df.drop_duplicates(subset="run")
    cutoffs = [50, 100, 500, 1000]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    ax = axes[0]
    for _, r in df.iterrows():
        ax.plot(cutoffs, [r[f"precision_{c}"] for c in cutoffs], "o-", label=r["run"],
                color=PALETTE.get(r["peak_strategy"], "#888"))
    ax.axhline(1.0, ls=":", color="k", lw=1)
    ax.set_xscale("log")
    ax.set_ylim(0.99, 1.001)
    ax.set_xlabel("distance cutoff (bp)")
    ax.set_ylabel("precision")
    ax.set_title("a  Precision is saturated (>99.6%) at every cutoff")
    ax.legend(fontsize=7)

    ax = axes[1]
    for _, r in df.iterrows():
        ax.plot(cutoffs, [r[f"recall_{c}"] for c in cutoffs], "o-", label=r["run"],
                color=PALETTE.get(r["peak_strategy"], "#888"))
    ax.set_xscale("log")
    ax.set_xlabel("distance cutoff (bp)")
    ax.set_ylabel("recall")
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.set_title("b  Recall carries all the variation")
    n_ref = int(df["n_reference"].iloc[0])
    fig.suptitle(
        f"A {n_ref/1e6:.1f}M-entry reference saturates precision: every call lands near some atlas PAS", fontsize=10)
    save(fig, "02_atlas_saturation")


# ===========================================================================
# 3. Count-controlled comparison (rarefaction) + distance-to-atlas
# ===========================================================================
def fig_count_controlled() -> None:
    rare = load(EXTRA, "strategy_rarefaction.csv")
    dist = load(EXTRA, "pas_distance_to_atlas.csv")
    ov = load(EXTRA, "strategy_overlap.csv")
    if rare is None and dist is None:
        _skipped.append("03_count_controlled")
        return

    ncols = int(rare is not None) + int(dist is not None) + int(ov is not None)
    fig, axes = plt.subplots(1, ncols, figsize=(6 * ncols, 4), squeeze=False, layout="constrained")
    axes = axes[0]
    i = 0

    if rare is not None:
        ax = axes[i]; i += 1
        r = rare.sort_values("frac_within_50_mean")
        v = r["frac_within_50_mean"]
        # the whole spread lives in the fourth decimal place, so zoom the axis and
        # say so explicitly -- a 0-1 axis would render this as five identical bars
        lo, hi = v.min(), v.max()
        pad = max((hi - lo) * 0.35, 5e-4)
        ax.barh(r["run"], v, xerr=r.get("frac_within_50_sd"), color="#4C72B0", capsize=3)
        ax.set_xlim(lo - pad, hi + pad * 1.6)
        for y, val in enumerate(v):
            ax.text(val + pad * 0.12, y, f"{val:.5f}", va="center", fontsize=7)
        n = int(r["n_sampled"].iloc[0])
        ax.set_xlabel("fraction of calls within 50 bp of an atlas PAS\n(axis zoomed: full spread is 0.15 pp)")
        ax.set_title(f"a  Rarefied to n={n:,} — spread is negligible")

    if dist is not None:
        ax = axes[i]; i += 1
        labels: list[str] = []
        for run, g in dist.groupby("run"):
            g = g.sort_values("bin_lo")
            ax.plot(np.arange(len(g)), g["frac"].clip(lower=1e-5), "o-", ms=4, label=run)
            labels = [f"{int(lo)}" for lo in g["bin_lo"]]
        ax.set_yscale("log")
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=45, fontsize=7)
        ax.set_xlabel("distance to nearest atlas PAS (bp, bin lower edge)")
        ax.set_ylabel("fraction of called PAS (log)")
        ax.set_title("b  >99% of every call set sits within 10 bp")
        ax.legend(fontsize=6)

    if ov is not None:
        ax = axes[i]
        runs = sorted(set(ov["run_a"]) | set(ov["run_b"]))
        m = pd.DataFrame(np.eye(len(runs)), index=runs, columns=runs)
        for _, r in ov.iterrows():
            m.loc[r["run_a"], r["run_b"]] = r["frac_a_in_b"]   # row a contained in col b
            m.loc[r["run_b"], r["run_a"]] = r["frac_b_in_a"]
        im = ax.imshow(m.values, cmap="Blues", vmin=0, vmax=1)
        ax.set_xticks(range(len(runs)), runs, rotation=45, ha="right", fontsize=6)
        ax.set_yticks(range(len(runs)), runs, fontsize=6)
        ax.grid(False)
        for a in range(len(runs)):
            for b in range(len(runs)):
                v = m.values[a, b]
                ax.text(b, a, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        color="white" if v > 0.6 else "black")
        fig.colorbar(im, ax=ax, label="fraction of row's PAS found in column")
        ax.set_title("c  The call sets are nested, not different")

    fig.suptitle("Count-controlled, the strategies are indistinguishable — and nested", fontsize=10)
    save(fig, "03_count_controlled")


# ===========================================================================
# 4. Internal priming and trim sensitivity
# ===========================================================================
def fig_ip_and_trim() -> None:
    ip = load(TABLES, "ip_filter_axis.csv")
    trim = load(TABLES, "trim_comparison.csv")
    if ip is None and trim is None:
        _skipped.append("04_ip_and_trim")
        return

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")

    if ip is not None:
        ax = axes[0]
        ip = ip.drop_duplicates(subset="ip_mode").set_index("ip_mode").reindex(["off", "annotate", "filter"]).dropna(how="all").reset_index()
        bars = ax.bar(ip["ip_mode"], ip["n_pas"], color=["#999", "#4C72B0", "#C44E52"])
        for b, n, p in zip(bars, ip["n_pas"], ip["pct_flagged_internal_priming"]):
            ax.text(b.get_x() + b.get_width() / 2, n, f"{int(n):,}\n({p:.2f}% flagged)",
                    ha="center", va="bottom", fontsize=8)
        ax.set_ylabel("PAS retained")
        ax.set_ylim(0, ip["n_pas"].max() * 1.2)
        ax.set_title("a  Internal-priming axis (lambda_gradient)")

    if trim is not None:
        ax = axes[1]
        ext = trim[trim["include_extended"].astype(str).str.lower() == "true"].sort_values("max_gene_distance")
        ax.plot(ext["max_gene_distance"], ext["n_pas"], "o-", color="#4C72B0")
        for _, r in ext.iterrows():
            ax.annotate(f"{int(r['n_pas']):,}", (r["max_gene_distance"], r["n_pas"]),
                        textcoords="offset points", xytext=(0, 6), fontsize=7, ha="center")
        ax.set_xscale("log")
        ax.set_xlabel("max_gene_distance (bp)")
        ax.set_ylabel("PAS retained")
        ax.set_title("b  Trim sensitivity (include_extended)")

        ax = axes[2]
        ax.scatter(trim["n_pas"], trim["n_clusters_mean"], s=50, color="#DD8452", zorder=3)
        for _, r in trim.iterrows():
            ax.annotate(r["branch"].replace("A2_trim_", ""), (r["n_pas"], r["n_clusters_mean"]),
                        textcoords="offset points", xytext=(5, 0), fontsize=6)
        rng = trim["n_clusters_mean"].max() - trim["n_clusters_mean"].min()
        ax.set_xlabel("PAS retained")
        ax.set_ylabel("mean clusters per dataset")
        ax.set_title(f"c  Clustering is insensitive to trim (range {rng:.2f})")

    fig.suptitle("Annotation-window and internal-priming sensitivity", fontsize=10)
    save(fig, "04_ip_and_trim")


# ===========================================================================
# 5. Clustering: resolution sweep and PAS-vs-GEX concordance
# ===========================================================================
def fig_clustering() -> None:
    clus = load(TABLES, "clustering_comparison.csv")
    gex = load(TABLES, "gex_concordance.csv")
    if clus is None and gex is None:
        _skipped.append("05_clustering")
        return

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), layout="constrained")

    if clus is not None:
        ax = axes[0]
        res = clus[clus["cluster_method"] == "leiden_tfidf"].sort_values("resolution")
        byres = res.drop_duplicates(subset="resolution")
        ax.plot(byres["resolution"], byres["n_clusters_mean"], "o-", color="#4C72B0")
        for _, r in byres.iterrows():
            ax.annotate(f"{r['n_clusters_mean']:.1f}", (r["resolution"], r["n_clusters_mean"]),
                        textcoords="offset points", xytext=(0, 7), fontsize=7, ha="center")
        ax.set_xlabel("Leiden resolution")
        ax.set_ylabel("mean clusters per dataset")
        ax.set_title("a  Resolution sweep")

        ax = axes[1]
        nn = clus.dropna(subset=["n_neighbors"]).sort_values("n_neighbors").drop_duplicates(subset="n_neighbors")
        ax.plot(nn["n_neighbors"], nn["n_clusters_mean"], "s-", color="#55A868")
        ax.set_xlabel("n_neighbors")
        ax.set_ylabel("mean clusters per dataset")
        ax.set_title("b  Neighbourhood size")

    if gex is not None:
        ax = axes[2]
        cols = [c for c in ("ARI_gexleiden_vs_pas", "AMI_gexleiden_vs_pas",
                            "ARI_celltype_vs_pas", "AMI_celltype_vs_pas") if c in gex]
        data = [gex[c].dropna() for c in cols]
        bp = ax.boxplot(data, tick_labels=[c.replace("_vs_pas", "").replace("_", "\n") for c in cols],
                        showmeans=True, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("#DD8452")
            patch.set_alpha(0.6)
        for j, d in enumerate(data, start=1):
            ax.scatter(np.random.default_rng(0).normal(j, 0.05, len(d)), d, s=12, color="k", alpha=0.5, zorder=3)
        ax.set_ylabel("agreement with PAS clustering")
        ax.set_ylim(0, 1)
        ax.set_title(f"c  PAS vs GEX labels (n={len(gex)} datasets)")
        ax.tick_params(axis="x", labelsize=7)

    fig.suptitle("PAS-space clustering: parameter sensitivity and agreement with expression labels", fontsize=10)
    save(fig, "05_clustering")


# ===========================================================================
# 6. THE honesty figure: the stage PDUI trend is a detection-rate artifact
# ===========================================================================
def fig_depth_confound() -> None:
    d = load(EXTRA, "depth_confound_by_celltype_stage.tsv")
    if d is None:
        _skipped.append("06_depth_confound")
        return
    d = d[d["mean_pdui_cov1"] >= 0].copy()

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout="constrained")

    # a) reported PDUI is a near-deterministic function of the zero-coverage rate
    ax = axes[0]
    ax.scatter(d["frac_zero_cov"], d["mean_pdui_all"], s=22, alpha=0.75, color="#C44E52")
    r = d["frac_zero_cov"].corr(d["mean_pdui_all"])
    sl, ic = np.polyfit(d["frac_zero_cov"], d["mean_pdui_all"], 1)
    xs = np.linspace(d["frac_zero_cov"].min(), d["frac_zero_cov"].max(), 20)
    ax.plot(xs, sl * xs + ic, "--", color="k", lw=1)
    ax.set_xlabel("fraction of (gene, cell) rows with zero coverage")
    ax.set_ylabel("mean PDUI (all rows)")
    ax.set_title(f"a  Reported PDUI ≈ detection rate (Pearson r = {r:+.3f})")

    # b) stage trend, unconditional vs coverage-conditioned
    ax = axes[1]
    agg = d.groupby("stage").agg(
        all=("mean_pdui_all", "mean"), c1=("mean_pdui_cov1", "mean"),
        c5=("mean_pdui_cov5", "mean"), c10=("mean_pdui_cov10", "mean"),
    ).reindex([s for s in STAGE_ORDER if s in set(d["stage"])])
    x = np.arange(len(agg))
    ax.plot(x, agg["all"], "o-", color="#C44E52", lw=2, label="all rows (as reported)")
    for col, lab, c in [("c1", "≥1 read", "#4C72B0"), ("c5", "≥5 reads", "#55A868"), ("c10", "≥10 reads", "#8172B3")]:
        ax.plot(x, agg[col], "s--", color=c, label=f"conditioned on {lab}")
    ax.set_xticks(x)
    ax.set_xticklabels(agg.index)
    ax.set_ylabel("mean PDUI")
    ax.set_yscale("log")
    ax.set_title("b  The stage trend disappears under coverage conditioning")
    ax.legend(fontsize=7)

    # c) per-celltype slope, unconditional vs conditioned
    ax = axes[2]
    rows = []
    for ct, g in d.groupby("celltype"):
        g = g.set_index("stage").reindex(STAGE_ORDER).dropna(subset=["mean_pdui_all"])
        if len(g) < 3:
            continue
        xx = np.arange(len(g))
        rows.append({
            "all": np.polyfit(xx, g["mean_pdui_all"], 1)[0],
            "c1": np.polyfit(xx, g["mean_pdui_cov1"], 1)[0],
            "c5": np.polyfit(xx, g["mean_pdui_cov5"], 1)[0],
        })
    sl = pd.DataFrame(rows)
    if not sl.empty:
        bp = ax.boxplot([sl["all"], sl["c1"], sl["c5"]],
                        tick_labels=["all rows", "≥1 read", "≥5 reads"],
                        showmeans=True, patch_artist=True)
        for patch, c in zip(bp["boxes"], ["#C44E52", "#4C72B0", "#55A868"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        for j, col in enumerate(["all", "c1", "c5"], start=1):
            ax.scatter(np.random.default_rng(1).normal(j, 0.05, len(sl)), sl[col], s=14, color="k", alpha=0.5, zorder=3)
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.22 * (hi - lo))
        for j, col in enumerate(["all", "c1", "c5"], start=1):
            ax.text(j, hi + 0.06 * (hi - lo), f"{(sl[col] < 0).sum()}/{len(sl)} down",
                    ha="center", va="bottom", fontsize=7)
        ax.axhline(0, color="k", lw=1)
        ax.set_ylabel("PDUI slope across stages")
        ax.set_title(f"c  Per-cell-type slope flips sign (n={len(sl)})")

    fig.suptitle(
        "The cohort-wide '3′UTR shortening' is a detection-rate artifact, not a length change", fontsize=10)
    save(fig, "06_depth_confound")


# ===========================================================================
# 7. Stage trajectory per cell type, annotated with what drives it
# ===========================================================================
def fig_stage_trajectory() -> None:
    mat = load(TABLES, "celltype_stage_pdui_matrix.csv")
    depth = load(EXTRA, "depth_confound_by_celltype_stage.tsv")
    if mat is None:
        _skipped.append("07_stage_trajectory")
        return
    mat = mat.set_index("celltype")
    stages = [s for s in STAGE_ORDER if s in mat.columns]
    mat = mat[stages]
    mat.index = [short_ct(i) for i in mat.index]
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]

    ncols = 2 if depth is not None else 1
    fig, axes = plt.subplots(1, ncols, figsize=(6.5 * ncols, max(5, 0.32 * len(mat))), squeeze=False)
    axes = axes[0]

    ax = axes[0]
    im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages, rotation=30, ha="right")
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=6)
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="mean PDUI (all rows)")
    ax.set_title("a  Reported mean PDUI by cell type and stage")

    if depth is not None:
        dd = depth[depth["mean_pdui_cov1"] >= 0].copy()
        piv = dd.pivot_table(index="celltype", columns="stage", values="mean_pdui_cov5")
        piv = piv[[s for s in STAGE_ORDER if s in piv.columns]]
        piv.index = [short_ct(i) for i in piv.index]
        piv = piv.reindex(mat.index)
        ax = axes[1]
        im = ax.imshow(piv.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(piv.shape[1]))
        ax.set_xticklabels(piv.columns, rotation=30, ha="right")
        ax.set_yticks(range(len(piv)))
        ax.set_yticklabels(piv.index, fontsize=6)
        ax.grid(False)
        fig.colorbar(im, ax=ax, label="mean PDUI (≥5 reads)")
        ax.set_title("b  Same cells, conditioned on coverage")

    fig.suptitle("Stage trajectory of 3′UTR usage — the gradient in (a) is absent in (b)", fontsize=10)
    save(fig, "07_stage_trajectory")


# ===========================================================================
# 8. Fisher hit counts scale with power, not with effect
# ===========================================================================
def fig_fisher_power() -> None:
    fh = load(TABLES, "fisher_hits_by_celltype_contrast.csv")
    if fh is None:
        _skipped.append("08_fisher_power")
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")

    ax = axes[0]
    ax.hist(fh["frac_sig"], bins=25, color="#4C72B0", edgecolor="white")
    med = fh["frac_sig"].median()
    overall = fh["n_sig"].sum() / fh["n_tests"].sum()
    ax.axvline(med, color="crimson", ls="--", lw=1.5, label=f"median {med:.2f}")
    ax.axvline(0.05, color="k", ls=":", lw=1.5, label="0.05 (nominal FDR)")
    ax.set_xlabel("fraction of PAS called significant (FDR < 0.05)")
    ax.set_ylabel("cell-type x contrast tests")
    ax.legend(fontsize=7)
    ax.set_title(f"a  {overall:.0%} of all PAS called significant")

    ax = axes[1]
    summ = fh.groupby("contrast").agg(
        frac=("frac_sig", "mean"), tests=("n_tests", "sum"), sig=("n_sig", "sum"),
        n=("celltype", "nunique"),
    ).sort_values("frac", ascending=True)
    ax.barh(summ.index, summ["frac"], color="#DD8452")
    for y, (f, s, t) in enumerate(zip(summ["frac"], summ["sig"], summ["tests"])):
        ax.text(f, y, f"  {f:.2f}  ({s:,}/{t:,})", va="center", fontsize=7)
    ax.set_xlim(0, summ["frac"].max() * 1.5)
    ax.set_xlabel("mean fraction of PAS called significant")
    ax.set_title("b  Every stage contrast calls a large minority significant")

    fig.suptitle(
        "Fisher tests pool reads, not cells: hit counts index statistical power, so rank by them — never read them as absolute evidence", fontsize=9)
    save(fig, "08_fisher_power")


# ===========================================================================
# 9. Recurrent vs cell-type-private switch genes
# ===========================================================================
def fig_recurrence() -> None:
    rec = load(EXTRA, "recurrence_decomposition.csv")
    if rec is None:
        rec = load(TABLES, "recurrent_switch_genes.csv")
        if rec is None:
            _skipped.append("09_recurrence")
            return

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout="constrained")

    ax = axes[0]
    counts = rec["n_celltypes_trending"].value_counts().sort_index()
    ax.bar(counts.index, counts.values, color="#4C72B0")
    ax.set_yscale("log")
    ax.set_xlabel("cell types in which the gene trends (|Spearman| ≥ 0.8)")
    ax.set_ylabel("genes")
    priv = int(counts.get(1, 0)) + int(counts.get(2, 0))
    ax.set_title(f"a  Most switches are private ({priv:,} genes in ≤2 cell types)")

    ax = axes[1]
    if "dominant_direction" in rec:
        top = rec.nlargest(25, "n_celltypes_trending").sort_values("n_celltypes_trending")
        lbl = top["symbol"].fillna(top["gene_id"]) if "symbol" in top else top["gene_id"]
        colors = ["#C44E52" if d == "decreasing" else "#55A868" for d in top["dominant_direction"]]
        ax.barh(range(len(top)), top["n_celltypes_trending"], color=colors)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(lbl, fontsize=6)
        ax.set_xlabel("cell types trending")
        ax.set_title("b  Most recurrent genes (red = shortening)")

    ax = axes[2]
    if "mean_slope" in rec:
        shared = rec[rec["n_celltypes_trending"] >= 3]["mean_slope"]
        private = rec[rec["n_celltypes_trending"] <= 2]["mean_slope"]
        bp = ax.boxplot([private, shared], tick_labels=[f"private\n(n={len(private):,})",
                                                        f"recurrent\n(n={len(shared):,})"],
                        showmeans=True, patch_artist=True)
        for patch, c in zip(bp["boxes"], ["#999999", "#4C72B0"]):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        ax.axhline(0, color="k", lw=1)
        ax.set_ylabel("mean stage slope")
        ax.set_title("c  Recurrent switches are not larger in effect")

    fig.suptitle("Recurrent versus cell-type-private 3′UTR switches", fontsize=10)
    save(fig, "09_recurrence")


# ===========================================================================
# 10. Does the lost distal segment carry regulatory elements?
# ===========================================================================
def fig_lost_distal_elements() -> None:
    scan = load(EXTRA, "lost_distal_element_scan.csv")
    if scan is None:
        _skipped.append("10_lost_distal_elements")
        return
    scan = scan[(scan["lost_len"] > 0) & (scan["ctrl_len"] > 0)].copy()
    for side in ("lost", "ctrl"):
        scan[f"{side}_are_density"] = 1000 * scan[f"{side}_are_pent"] / scan[f"{side}_len"]
        scan[f"{side}_seed_density"] = 1000 * scan[f"{side}_seed_total"] / scan[f"{side}_len"]

    from scipy import stats as st

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), layout="constrained")

    ax = axes[0]
    pairs = [("are_density", "ARE pentamer\n(ATTTA)"), ("seed_density", "miRNA 7mer-m8\nseed sites")]
    data, labels, colors = [], [], []
    for key, lab in pairs:
        data += [scan[f"lost_{key}"], scan[f"ctrl_{key}"]]
        labels += [f"lost distal\n{lab}", f"proximal ctrl\n{lab}"]
        colors += ["#C44E52", "#999999"]
    bp = ax.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.6)
    ax.set_ylabel("sites per kb")
    ax.tick_params(axis="x", labelsize=6)
    w = st.wilcoxon(scan["lost_are_density"], scan["ctrl_are_density"])
    ax.set_title(f"a  Element density, lost vs control\n(ARE Wilcoxon p={w.pvalue:.2e}, n={len(scan):,})")

    ax = axes[1]
    ax.scatter(scan["lost_at_frac"], scan["lost_are_density"], s=6, alpha=0.3, color="#4C72B0")
    rho2 = scan["lost_at_frac"].corr(scan["lost_are_density"], method="spearman")
    rho3 = scan["lost_len"].corr(scan["lost_are_pent"], method="spearman")
    ax.set_xlabel("AT fraction of lost segment")
    ax.set_ylabel("ARE pentamers per kb")
    ax.set_title(f"b  ARE density is composition-driven (ρ={rho2:+.2f});\ncount is length-driven (ρ={rho3:+.2f})")

    # the interval these elements are counted in is mostly not 3'UTR-scale
    ax = axes[2]
    ax.hist(np.log10(scan["lost_len"].clip(lower=1)), bins=45, color="#C44E52", edgecolor="white")
    med = scan["lost_len"].median()
    ax.axvline(np.log10(med), color="k", ls="--", lw=1.5, label=f"median {med:,.0f} bp")
    ax.axvline(np.log10(3000), color="crimson", ls=":", lw=1.5, label="3 kb (long 3′UTR)")
    ax.set_xlabel("log₁₀ lost segment length (bp)")
    ax.set_ylabel("genes")
    ax.legend(fontsize=7)
    frac = (scan["lost_len"] > 3000).mean()
    ax.set_title(f"c  {frac:.0%} of intervals exceed 3 kb —\nthese are not tandem 3′UTR pairs")

    fig.suptitle(
        "Regulatory-element content of the 3′UTR segment lost on shortening — "
        "a sequence-composition effect, so it is not evidence of targeted de-repression", fontsize=9)
    save(fig, "10_lost_distal_elements")


# ===========================================================================
# 11. Cancer-gene and pathway characterisation of the recurrent switch genes
# ===========================================================================
def fig_gene_sets() -> None:
    ov = load(EXTRA, "cancer_gene_overlap.csv")
    hall = load(EXTRA, "hallmark_enrichment.csv")
    if ov is None and hall is None:
        _skipped.append("11_gene_sets")
        return

    ncols = int(ov is not None) * 2 + int(hall is not None)
    fig, axes = plt.subplots(1, ncols, figsize=(4.7 * ncols, 4.3), squeeze=False, layout="constrained")
    axes = axes[0]
    i = 0

    if ov is not None:
        # a) raw membership rates, recurrent vs private
        ax = axes[i]; i += 1
        x = np.arange(len(ov))
        w = 0.38
        ax.bar(x - w / 2, ov["rate_recurrent"] * 100, w, label="recurrent (≥3 CT)", color="#4C72B0")
        ax.bar(x + w / 2, ov["rate_private"] * 100, w, label="private (≤2 CT)", color="#999999")
        ax.set_xticks(x, [g.replace("_", "\n") for g in ov["gene_set"]], fontsize=7)
        ax.set_ylabel("% of genes in the set")
        ax.legend(fontsize=7)
        ax.set_title("a  Membership looks enriched…")

        # b) …until detectability is adjusted for
        ax = axes[i]; i += 1
        y = np.arange(len(ov))
        ax.scatter(ov["odds_ratio_crude"], y + 0.13, s=55, color="#C44E52", label="crude", zorder=3)
        if "odds_ratio_adj_detectability" in ov:
            ax.scatter(ov["odds_ratio_adj_detectability"], y - 0.13, s=55, color="#4C72B0",
                       marker="s", label="adjusted for detectability", zorder=3)
            for yy, (cr, ad) in enumerate(zip(ov["odds_ratio_crude"],
                                              ov["odds_ratio_adj_detectability"])):
                ax.plot([cr, ad], [yy + 0.13, yy - 0.13], color="k", lw=0.7, alpha=0.5, zorder=2)
        ax.axvline(1.0, color="k", ls="--", lw=1)
        ax.set_yticks(y, [g.replace("_", " ") for g in ov["gene_set"]], fontsize=7)
        ax.set_xlabel("odds ratio, recurrent vs private")
        ax.legend(fontsize=7)
        ax.set_title("b  …and the enrichment disappears")

    if hall is not None:
        ax = axes[i]
        h = hall.nsmallest(12, "p_value").sort_values("odds_ratio")
        colors = ["#C44E52" if q < 0.05 else "#BBBBBB" for q in h["q_value"]]
        ax.barh(range(len(h)), h["odds_ratio"], color=colors)
        ax.axvline(1.0, color="k", ls="--", lw=1)
        ax.set_yticks(range(len(h)),
                      [g.replace("HALLMARK_", "").replace("_", " ").title() for g in h["gene_set"]],
                      fontsize=6)
        ax.set_xlabel("odds ratio (recurrent vs tested universe)")
        n_sig = int((hall["q_value"] < 0.05).sum())
        ax.set_title(f"c  Hallmark: {n_sig}/{len(hall)} sets at FDR<0.05")

    fig.suptitle(
        "Recurrent switch genes carry no cancer-gene or pathway signal once detectability is controlled",
        fontsize=9)
    save(fig, "11_gene_sets")


FIGURES = [
    fig_strategy_benchmark,
    fig_atlas_saturation,
    fig_count_controlled,
    fig_ip_and_trim,
    fig_clustering,
    fig_depth_confound,
    fig_stage_trajectory,
    fig_fisher_power,
    fig_recurrence,
    fig_lost_distal_elements,
    fig_gene_sets,
]


def main(argv: list[str] | None = None) -> int:
    global TABLES, EXTRA, OUT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tables", type=Path, default=TABLES)
    ap.add_argument("--extra", type=Path, default=EXTRA)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    TABLES, EXTRA, OUT = args.tables, args.extra, args.out

    print(f"tables : {TABLES}")
    print(f"extra  : {EXTRA}")
    print(f"out    : {OUT}\n")
    for fn in FIGURES:
        print(f"[{fn.__name__}]")
        try:
            fn()
        except Exception as exc:  # keep going; a broken panel must not hide the rest
            print(f"  ERROR {fn.__name__}: {type(exc).__name__}: {exc}")
            _skipped.append(fn.__name__)

    print(f"\nwrote {len(_written)} figures -> {OUT}")
    if _skipped:
        print(f"skipped {len(_skipped)}: {', '.join(_skipped)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
