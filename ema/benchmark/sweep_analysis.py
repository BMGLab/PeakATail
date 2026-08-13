"""Cumulative cross-experiment analysis over a finished PeakATail sweep.

This is the SWEEP-level analysis layer on top of the PER-RUN tools already in
``ema/benchmark`` (``metrics.py``, ``switch_diff_compare.py``,
``length_compare.py``, ``visualize.py``). It performs no pipeline
re-computation: it reads finished run-directory artifacts (``pasbed.bed``,
``benchmark_vs_polyasite_v3.json``, ``atlas_stats.json``,
``07_clustering/*/clusters.h5ad``, ``B3_switch/{diff,trend,length}``) across
every run in a sweep and emits sweep-level CSV tables, PNG/SVG figures, and a
cumulative markdown + JSON summary suitable for the reports under
``reports/``.

Covers, per the sweep layout under ``<sweep_root>/runs/``:

1. Peak-calling strategy comparison (``runs/grid/*``): n_PAS, precision /
   recall / F1 @ several distance cutoffs vs PolyASite v3, atlas match rate.
2. Fasta / internal-priming axis (``lg_annotate`` vs ``lg_ip_filter`` vs
   ``lg_ip_off``): n_PAS dropped and % flagged internal-priming.
3. Trim sensitivity (``runs/reannotate/A2_*``): gene-distance / UTR-multiplier
   / include-extended vs n_PAS, n_cells, n_clusters.
4. Clustering comparison (``runs/reannotate/A3_*``): method / resolution /
   n_neighbors vs mean n_clusters, plus cohort GEX concordance
   (``runs/B1_cohort_full/B2_gex_celltyping/concordance.csv``).
5. Switch cumulative (``runs/B1_cohort_full/B3_switch``) — the biological
   headline: shortening vs lengthening across cell types, slope
   distribution, mean-PDUI-by-stage matrix, per-stage-contrast Fisher hit
   counts, and recurrent switch genes (trending in >= N cell types).

CLI::

    python -m ema.benchmark.sweep_analysis <sweep_root> --out <dir>

See ``scripts/harvest_sweep_comparison.py`` for the earlier, narrower
prototype this module folds in and extends (full-cutoff strategy comparison,
the internal-priming axis, and the switch cumulative layer across all cell
types).
"""
from __future__ import annotations

import argparse
import json
import logging
import warnings
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Stage progression used to order celltype x stage matrices / trend json keys.
STAGE_ORDER = ["Normal", "StageI", "IVprimary", "Met"]

#: Distance cutoffs (bp) reported in the strategy comparison table/figures.
DEFAULT_CUTOFFS: tuple[int, ...] = (50, 100, 500, 1000)

#: Fallback axis metadata for ``runs/grid/<run_id>`` when a run doesn't carry
#: its own axis params in an artifact. Mirrors scripts/harvest_sweep_comparison.py.
GRID_META: dict[str, dict] = {
    "lg_annotate": {"peak_strategy": "lambda_gradient", "atlas_mode": "annotate", "ip_mode": "annotate"},
    "lp_annotate": {"peak_strategy": "lambda_poisson", "atlas_mode": "annotate", "ip_mode": "annotate"},
    "si_annotate": {"peak_strategy": "sierra_iterative", "atlas_mode": "annotate", "ip_mode": "annotate"},
    "lg_ip_filter": {"peak_strategy": "lambda_gradient", "atlas_mode": "annotate", "ip_mode": "filter"},
    "lg_ip_off": {"peak_strategy": "lambda_gradient", "atlas_mode": "annotate", "ip_mode": "off"},
}

#: Fallback trim-axis metadata for ``runs/reannotate/A2_*`` branches.
TRIM_META: dict[str, dict] = {
    "A2_trim_d1000_ext": {"max_gene_distance": 1000, "utr_multiplier": 2.0, "include_extended": True},
    "A2_trim_d2000_ext": {"max_gene_distance": 2000, "utr_multiplier": 2.0, "include_extended": True},
    "A2_trim_d3000_ext": {"max_gene_distance": 3000, "utr_multiplier": 2.0, "include_extended": True},
    "A2_trim_d5000_ext": {"max_gene_distance": 5000, "utr_multiplier": 2.0, "include_extended": True},
    "A2_trim_d10000_ext": {"max_gene_distance": 10000, "utr_multiplier": 2.0, "include_extended": True},
    "A2_trim_default": {"max_gene_distance": 5000, "utr_multiplier": 2.0, "include_extended": False},
    "A2_trim_mult1.5": {"max_gene_distance": 5000, "utr_multiplier": 1.5, "include_extended": False},
    "A2_trim_mult3.0": {"max_gene_distance": 5000, "utr_multiplier": 3.0, "include_extended": False},
}

#: Fallback clustering-axis metadata for ``runs/reannotate/A3_*`` branches
#: (``A2_trim_default`` doubles as the clustering-default row).
CLUSTER_META: dict[str, dict] = {
    "A2_trim_default": {"cluster_method": "leiden_tfidf", "resolution": 1.0, "n_neighbors": 30},
    "A3_res0.5": {"cluster_method": "leiden_tfidf", "resolution": 0.5, "n_neighbors": 30},
    "A3_res2.0": {"cluster_method": "leiden_tfidf", "resolution": 2.0, "n_neighbors": 30},
    "A3_nn15": {"cluster_method": "leiden_tfidf", "resolution": 1.0, "n_neighbors": 15},
    "A3_nn50": {"cluster_method": "leiden_tfidf", "resolution": 1.0, "n_neighbors": 50},
    "A3_libsize": {"cluster_method": "leiden_libsize", "resolution": 1.0, "n_neighbors": 30},
}


# ---------------------------------------------------------------------------
# Small IO / formatting helpers
# ---------------------------------------------------------------------------


def _read_json(path: Path) -> dict:
    """Read a JSON file, returning ``{}`` on any missing-file or parse error."""
    path = Path(path)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        log.warning("_read_json: failed to read %s: %s", path, exc)
        return {}


def count_bed_lines(path: Path) -> Optional[int]:
    """Count entries (== n_PAS) in a BED file. Returns ``None`` if missing."""
    path = Path(path)
    if not path.exists():
        return None
    with open(path) as fh:
        return sum(1 for _ in fh)


def _fmt(value, spec: str = "") -> str:
    """Format a possibly-``None``/NaN value for markdown, else ``'n/a'``."""
    if value is None:
        return "n/a"
    if isinstance(value, float) and np.isnan(value):
        return "n/a"
    try:
        return format(value, spec)
    except (ValueError, TypeError):
        return str(value)


def _json_default(o):
    """``json.dump`` default= hook for numpy scalar/array types."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o) if not np.isnan(o) else None
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _leiden_column(obs_columns: Sequence[str]) -> Optional[str]:
    """Find the cluster-label column in an AnnData ``.obs`` (``leiden`` or similar)."""
    if "leiden" in obs_columns:
        return "leiden"
    for c in obs_columns:
        cl = c.lower()
        if "leiden" in cl or cl == "cluster":
            return c
    return None


def _h5ad_stats(h5ad_path: Path) -> dict:
    """n_cells / n_clusters for one ``clusters.h5ad``, read backed (no full load)."""
    import anndata as ad

    try:
        a = ad.read_h5ad(h5ad_path, backed="r")
    except Exception as exc:  # pragma: no cover - defensive, depends on disk state
        log.warning("_h5ad_stats: failed reading %s: %s", h5ad_path, exc)
        return {}
    out: dict = {"n_cells": int(a.n_obs)}
    col = _leiden_column(list(a.obs.columns))
    if col is not None:
        out["n_clusters"] = int(a.obs[col].astype(str).nunique())
    return out


def _branch_cluster_stats(branch_dir: Path) -> dict:
    """Aggregate n_cells/n_clusters across all datasets' clusters.h5ad under a run/branch dir."""
    h5ad_paths = sorted((Path(branch_dir) / "07_clustering").glob("*/clusters.h5ad"))
    n_cells_list: list[int] = []
    n_clusters_list: list[int] = []
    for h in h5ad_paths:
        stats = _h5ad_stats(h)
        if "n_cells" in stats:
            n_cells_list.append(stats["n_cells"])
        if "n_clusters" in stats:
            n_clusters_list.append(stats["n_clusters"])
    return {
        "n_datasets": len(h5ad_paths),
        "n_cells_total": int(sum(n_cells_list)) if n_cells_list else None,
        "n_clusters_mean": round(float(np.mean(n_clusters_list)), 2) if n_clusters_list else None,
        "n_clusters_median": float(np.median(n_clusters_list)) if n_clusters_list else None,
    }


def _list_celltypes(sweep_root: Path) -> list[str]:
    """List cell-type directory names under ``B3_switch/diff``."""
    diff_dir = Path(sweep_root) / "runs" / "B1_cohort_full" / "B3_switch" / "diff"
    if not diff_dir.exists():
        return []
    return sorted(p.name for p in diff_dir.iterdir() if p.is_dir())


# ---------------------------------------------------------------------------
# 1. Peak-calling strategy comparison
# ---------------------------------------------------------------------------


def harvest_strategy_comparison(
    sweep_root: Path, cutoffs: Sequence[int] = DEFAULT_CUTOFFS
) -> pd.DataFrame:
    """Per peak-calling strategy: n_PAS, precision/recall/F1 @ cutoffs, atlas match.

    Reads ``runs/grid/<run>/pasbed.bed``, ``benchmark_vs_polyasite_v3.json``,
    and ``**/atlas_stats.json`` for every run listed in :data:`GRID_META`
    that exists on disk. Ranked descending by F1 at the widest cutoff.
    """
    grid_dir = Path(sweep_root) / "runs" / "grid"
    rows = []
    for run_id, meta in GRID_META.items():
        d = grid_dir / run_id
        if not d.exists():
            continue
        row = {"run": run_id, **meta, "n_pas": count_bed_lines(d / "pasbed.bed")}
        bench = _read_json(d / "benchmark_vs_polyasite_v3.json")
        row["n_reference"] = bench.get("n_reference")
        bench_cutoffs = bench.get("cutoffs") or {}
        for cut in cutoffs:
            c = bench_cutoffs.get(str(cut), {})
            row[f"precision_{cut}"] = c.get("precision")
            row[f"recall_{cut}"] = c.get("recall")
            row[f"f1_{cut}"] = c.get("f1")
        atlas_paths = sorted(d.rglob("atlas_stats.json"))
        if atlas_paths:
            atlas = _read_json(atlas_paths[0])
            row["atlas_match_rate"] = atlas.get("atlas_match_rate")
            row["n_atlas_matched"] = atlas.get("n_atlas_matched")
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    rank_col = f"f1_{max(cutoffs)}"
    if rank_col in df.columns:
        df = df.sort_values(rank_col, ascending=False, na_position="last").reset_index(drop=True)
        df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df


def harvest_ip_filter_axis(strategy_df: pd.DataFrame) -> pd.DataFrame:
    """The lambda_gradient fasta/internal-priming axis: annotate vs filter vs off.

    ``ip_mode="off"`` applies no internal-priming logic (baseline n_PAS);
    ``ip_mode="filter"`` actively drops PAS flagged as internal-priming
    artifacts; ``ip_mode="annotate"`` tags but does not drop. n_PAS dropped
    and % flagged are both computed relative to the ``off`` baseline.
    """
    if strategy_df.empty or "peak_strategy" not in strategy_df.columns:
        return pd.DataFrame()
    sub = strategy_df[strategy_df["peak_strategy"] == "lambda_gradient"].copy()
    if sub.empty or "off" not in set(sub.get("ip_mode", [])):
        return sub.reset_index(drop=True)

    baseline = sub.loc[sub["ip_mode"] == "off", "n_pas"].iloc[0]
    if baseline:
        sub["n_pas_dropped_vs_off"] = baseline - sub["n_pas"]
        sub["pct_flagged_internal_priming"] = (
            (baseline - sub["n_pas"]) / baseline * 100
        ).round(2)
    keep_cols = [c for c in sub.columns if not c.startswith(("precision_", "recall_", "f1_"))]
    return sub[keep_cols].sort_values("ip_mode").reset_index(drop=True)


# ---------------------------------------------------------------------------
# 2. Trim sensitivity
# ---------------------------------------------------------------------------


def harvest_trim_comparison(sweep_root: Path) -> pd.DataFrame:
    """Per trim branch (``runs/reannotate/A2_*``): trim params vs n_PAS/n_cells/n_clusters."""
    reannotate_dir = Path(sweep_root) / "runs" / "reannotate"
    rows = []
    for branch, fallback in TRIM_META.items():
        d = reannotate_dir / branch
        if not d.exists():
            continue
        manifest = _read_json(d / "branch_manifest.json")
        trim_meta = manifest.get("trim") or fallback
        row = {"branch": branch, **trim_meta, "n_pas": count_bed_lines(d / "pasbed.bed")}
        row.update(_branch_cluster_stats(d))
        rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty and "max_gene_distance" in df.columns:
        df = df.sort_values(["max_gene_distance", "utr_multiplier"]).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# 3. Clustering comparison + cohort GEX concordance
# ---------------------------------------------------------------------------


def harvest_clustering_comparison(sweep_root: Path) -> pd.DataFrame:
    """Per clustering branch (``runs/reannotate/A3_*`` + default): params vs mean n_clusters."""
    reannotate_dir = Path(sweep_root) / "runs" / "reannotate"
    rows = []
    for branch, fallback in CLUSTER_META.items():
        d = reannotate_dir / branch
        if not d.exists():
            continue
        manifest = _read_json(d / "branch_manifest.json")
        cl_meta = manifest.get("clustering") or {}
        row = {
            "branch": branch,
            "cluster_method": cl_meta.get("method", fallback.get("cluster_method")),
            "resolution": cl_meta.get("resolution", fallback.get("resolution")),
            "n_neighbors": cl_meta.get("n_neighbors", fallback.get("n_neighbors")),
        }
        row.update(_branch_cluster_stats(d))
        rows.append(row)
    return pd.DataFrame(rows)


def harvest_gex_concordance(sweep_root: Path) -> pd.DataFrame:
    """Per-GSM cohort GEX-vs-PAS-leiden concordance (ARI/AMI)."""
    conc_path = (
        Path(sweep_root) / "runs" / "B1_cohort_full" / "B2_gex_celltyping" / "concordance.csv"
    )
    if not conc_path.exists():
        return pd.DataFrame()
    return pd.read_csv(conc_path)


def summarize_gex_concordance(concordance_df: pd.DataFrame) -> dict:
    """Mean/median/std of each ARI_*/AMI_* column across GSMs."""
    if concordance_df.empty:
        return {}
    cols = [c for c in concordance_df.columns if c.startswith(("ARI_", "AMI_"))]
    summary: dict = {"n_gsm": int(len(concordance_df))}
    for c in cols:
        col = concordance_df[c].dropna()
        if col.empty:
            continue
        summary[c] = {
            "mean": round(float(col.mean()), 4),
            "median": round(float(col.median()), 4),
            "std": round(float(col.std()), 4) if len(col) > 1 else 0.0,
        }
    return summary


# ---------------------------------------------------------------------------
# 4. Switch cumulative — the biological headline
# ---------------------------------------------------------------------------


def harvest_switch_trends(sweep_root: Path) -> pd.DataFrame:
    """Per-celltype length-trend summary: slope/spearman/direction + mean PDUI by stage.

    Reads only ``B3_switch/trend/<celltype>/length_trend.json`` (small, one
    per cell type) — never the multi-GB per-cell length tables.
    """
    trend_dir = Path(sweep_root) / "runs" / "B1_cohort_full" / "B3_switch" / "trend"
    rows = []
    for ct in _list_celltypes(sweep_root):
        d = _read_json(trend_dir / ct / "length_trend.json")
        if not d:
            continue
        row = {
            "celltype": ct,
            "n_stages": d.get("n_stages"),
            "slope": d.get("slope"),
            "spearman": d.get("spearman"),
            "direction": d.get("direction"),
        }
        mean_by_stage = d.get("mean_by_stage") or {}
        for stage in STAGE_ORDER:
            row[f"pdui_{stage}"] = mean_by_stage.get(stage)
        rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("slope", na_position="last").reset_index(drop=True)
    return df


def summarize_shortening_lengthening(trend_df: pd.DataFrame) -> dict:
    """Counts + slope distribution across cell types, by trend direction.

    ``direction`` values follow the trend JSON convention: ``"decreasing"``
    (3'UTR shortening — PDUI falling with stage) and ``"increasing"``
    (lengthening). Anything else (e.g. ``"flat"`` or missing) is bucketed
    into ``n_flat_or_other``.
    """
    if trend_df.empty:
        return {"n_celltypes": 0, "n_shortening": 0, "n_lengthening": 0, "n_flat_or_other": 0}

    counts = trend_df["direction"].value_counts()
    n_shortening = int(counts.get("decreasing", 0))
    n_lengthening = int(counts.get("increasing", 0))
    n_total = int(len(trend_df))
    slopes = trend_df["slope"].dropna()
    return {
        "n_celltypes": n_total,
        "n_shortening": n_shortening,
        "n_lengthening": n_lengthening,
        "n_flat_or_other": n_total - n_shortening - n_lengthening,
        "slope_mean": round(float(slopes.mean()), 6) if len(slopes) else None,
        "slope_median": round(float(slopes.median()), 6) if len(slopes) else None,
        "slope_std": round(float(slopes.std()), 6) if len(slopes) > 1 else None,
    }


def celltype_stage_pdui_matrix(trend_df: pd.DataFrame) -> pd.DataFrame:
    """celltype x stage matrix of mean PDUI, built from the trend jsons' mean_by_stage."""
    if trend_df.empty:
        return pd.DataFrame()
    cols = [f"pdui_{s}" for s in STAGE_ORDER if f"pdui_{s}" in trend_df.columns]
    if not cols:
        return pd.DataFrame()
    matrix = trend_df.set_index("celltype")[cols].copy()
    matrix.columns = [c.replace("pdui_", "") for c in matrix.columns]
    return matrix


def harvest_fisher_hit_counts(sweep_root: Path, fdr: float = 0.05) -> pd.DataFrame:
    """n significant PAS (qvalue < fdr) per celltype x stage-contrast, from Fisher pairwise TSVs.

    Only the ``qvalue`` column is read from each TSV (``usecols``) to keep
    memory bounded; these per-contrast tables are small compared to the raw
    length tables.
    """
    diff_dir = Path(sweep_root) / "runs" / "B1_cohort_full" / "B3_switch" / "diff"
    rows = []
    for ct in _list_celltypes(sweep_root):
        fisher_dir = diff_dir / ct / "fisher" / "differential"
        if not fisher_dir.exists():
            continue
        for tsv in sorted(fisher_dir.glob("fisher_*_vs_*.tsv")):
            contrast = tsv.stem.replace("fisher_", "")
            try:
                df = pd.read_csv(tsv, sep="\t", usecols=["qvalue"])
            except (ValueError, pd.errors.EmptyDataError) as exc:
                log.warning("harvest_fisher_hit_counts: skipping %s: %s", tsv, exc)
                continue
            n_tests = int(len(df))
            n_sig = int((df["qvalue"] < fdr).sum())
            rows.append(
                {
                    "celltype": ct,
                    "contrast": contrast,
                    "n_tests": n_tests,
                    "n_sig": n_sig,
                    "frac_sig": round(n_sig / n_tests, 4) if n_tests else None,
                }
            )
    return pd.DataFrame(rows)


def summarize_fisher_hits(fisher_df: pd.DataFrame) -> pd.DataFrame:
    """Roll up per-(celltype, contrast) Fisher hit counts to per-contrast totals."""
    if fisher_df.empty:
        return pd.DataFrame()
    return (
        fisher_df.groupby("contrast")
        .agg(
            n_celltypes=("celltype", "nunique"),
            total_tests=("n_tests", "sum"),
            total_sig=("n_sig", "sum"),
            mean_frac_sig=("frac_sig", "mean"),
        )
        .reset_index()
        .sort_values("total_sig", ascending=False)
        .reset_index(drop=True)
    )


def harvest_recurrent_switch_genes(
    sweep_root: Path, min_celltypes: int = 3, spearman_thresh: float = 0.8
) -> pd.DataFrame:
    """Genes with a consistent stage trend in >= ``min_celltypes`` cell types.

    A gene "trends" in a cell type when ``|spearman| >= spearman_thresh`` in
    that cell type's ``length_trend_by_gene.tsv``. Reads only these small
    per-gene summary TSVs — never the multi-GB per-cell length tables.

    Returns a DataFrame with one row per recurrent gene, sorted by
    ``n_celltypes_trending`` descending, columns::

        gene_id, n_celltypes_trending, dominant_direction, n_consistent,
        mean_slope, mean_abs_spearman
    """
    trend_dir = Path(sweep_root) / "runs" / "B1_cohort_full" / "B3_switch" / "trend"
    records = []
    for ct in _list_celltypes(sweep_root):
        p = trend_dir / ct / "length_trend_by_gene.tsv"
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, sep="\t")
        except (pd.errors.EmptyDataError, OSError) as exc:
            log.warning("harvest_recurrent_switch_genes: skipping %s: %s", p, exc)
            continue
        if df.empty or "spearman" not in df.columns:
            continue
        trending = df[df["spearman"].abs() >= spearman_thresh]
        for _, r in trending.iterrows():
            records.append(
                {
                    "gene_id": r["gene_id"],
                    "celltype": ct,
                    "direction": r.get("direction"),
                    "slope": r.get("slope"),
                    "spearman": r.get("spearman"),
                }
            )

    if not records:
        return pd.DataFrame(
            columns=[
                "gene_id", "n_celltypes_trending", "dominant_direction",
                "n_consistent", "mean_slope", "mean_abs_spearman",
            ]
        )

    long_df = pd.DataFrame(records)
    grouped = long_df.groupby("gene_id")

    def _dominant(s: pd.Series):
        vc = s.value_counts()
        return vc.idxmax() if len(vc) else None

    agg = pd.DataFrame(
        {
            "n_celltypes_trending": grouped["celltype"].nunique(),
            "mean_slope": grouped["slope"].mean(),
            "mean_abs_spearman": grouped["spearman"].apply(lambda s: s.abs().mean()),
        }
    )
    agg["dominant_direction"] = grouped["direction"].apply(_dominant)
    # n_consistent = how many of the trending calls agree with the dominant direction.
    dom_map = agg["dominant_direction"].to_dict()
    long_df["_matches_dominant"] = long_df.apply(
        lambda r: r["direction"] == dom_map.get(r["gene_id"]), axis=1
    )
    agg["n_consistent"] = long_df.groupby("gene_id")["_matches_dominant"].sum().astype(int)
    agg = agg.reset_index()

    recurrent = (
        agg[agg["n_celltypes_trending"] >= min_celltypes]
        .sort_values(["n_celltypes_trending", "mean_abs_spearman"], ascending=[False, False])
        .reset_index(drop=True)
    )
    return recurrent[
        ["gene_id", "n_celltypes_trending", "dominant_direction", "n_consistent",
         "mean_slope", "mean_abs_spearman"]
    ]


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _savefig(fig, output_path: Path) -> None:
    """Save a figure as both PNG (dpi=300) and SVG, creating parent dirs as needed."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_strategy_precision_recall(
    strategy_df: pd.DataFrame, output_path: Path, cutoffs: Sequence[int] = DEFAULT_CUTOFFS
) -> None:
    """Bar chart of n_PAS per strategy + n_PAS-vs-recall tradeoff scatter."""
    if strategy_df.empty:
        log.warning("plot_strategy_precision_recall: empty df, skipping")
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    x = np.arange(len(strategy_df))
    colors = plt.cm.tab10(np.linspace(0, 1, len(strategy_df)))

    axes[0].bar(x, strategy_df["n_pas"], color=colors, edgecolor="white")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(strategy_df["run"], rotation=30, ha="right")
    axes[0].set_ylabel("n_PAS")
    axes[0].set_title("PAS count by peak-calling strategy")
    for xi, v in zip(x, strategy_df["n_pas"]):
        if pd.notna(v):
            axes[0].text(xi, v, f"{int(v):,}", ha="center", va="bottom", fontsize=9)
    axes[0].grid(True, alpha=0.2, axis="y")

    rc = f"recall_{max(cutoffs)}"
    if rc in strategy_df.columns:
        axes[1].scatter(strategy_df["n_pas"], strategy_df[rc], s=90, c=colors, edgecolor="black")
        for _, row in strategy_df.iterrows():
            if pd.notna(row.get(rc)):
                axes[1].annotate(
                    row["run"], (row["n_pas"], row[rc]), fontsize=8,
                    xytext=(6, 4), textcoords="offset points",
                )
        axes[1].set_xlabel("n_PAS")
        axes[1].set_ylabel(f"recall@{max(cutoffs)}bp")
        axes[1].set_title("n_PAS vs recall tradeoff")
        axes[1].grid(True, alpha=0.2)

    plt.tight_layout()
    _savefig(fig, output_path)


def plot_trim_cluster_sensitivity(
    trim_df: pd.DataFrame, cluster_df: pd.DataFrame, output_path: Path
) -> None:
    """Line plots: n_PAS/mean-n_clusters vs gene-distance, and mean-n_clusters vs resolution."""
    if trim_df.empty and cluster_df.empty:
        log.warning("plot_trim_cluster_sensitivity: no data, skipping")
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    if not trim_df.empty and "max_gene_distance" in trim_df.columns:
        base_mult = trim_df["utr_multiplier"].mode().iloc[0] if "utr_multiplier" in trim_df else None
        sub = trim_df[trim_df["utr_multiplier"] == base_mult] if base_mult is not None else trim_df
        sub = sub.sort_values("max_gene_distance")
        axes[0].plot(sub["max_gene_distance"], sub["n_pas"], "o-", color="#1f77b4", label="n_PAS")
        axes[0].set_xlabel("max_gene_distance (bp)")
        axes[0].set_ylabel("n_PAS", color="#1f77b4")
        if "n_clusters_mean" in sub.columns:
            ax2 = axes[0].twinx()
            ax2.plot(sub["max_gene_distance"], sub["n_clusters_mean"], "s--", color="#ff7f0e",
                      label="mean n_clusters")
            ax2.set_ylabel("mean n_clusters", color="#ff7f0e")
        axes[0].set_title("Trim sensitivity: gene-distance axis")
        axes[0].grid(True, alpha=0.2)

    if not cluster_df.empty and "resolution" in cluster_df.columns:
        sub = cluster_df.dropna(subset=["resolution"]).sort_values("resolution")
        axes[1].plot(sub["resolution"], sub["n_clusters_mean"], "o-", color="#2ca02c")
        axes[1].set_xlabel("Leiden resolution")
        axes[1].set_ylabel("mean n_clusters")
        axes[1].set_title("Clustering sensitivity: resolution axis")
        axes[1].grid(True, alpha=0.2)

    plt.tight_layout()
    _savefig(fig, output_path)


def plot_celltype_stage_heatmap(matrix_df: pd.DataFrame, output_path: Path) -> None:
    """Heatmap of mean PDUI, celltype (rows) x stage (columns)."""
    if matrix_df.empty:
        log.warning("plot_celltype_stage_heatmap: empty matrix, skipping")
        return
    fig, ax = plt.subplots(figsize=(8, max(6, 0.35 * len(matrix_df))))
    im = ax.imshow(matrix_df.values.astype(float), aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(len(matrix_df.columns)))
    ax.set_xticklabels(matrix_df.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(matrix_df.index)))
    ax.set_yticklabels(matrix_df.index, fontsize=7)
    fig.colorbar(im, ax=ax, label="mean PDUI")
    ax.set_title("Mean PDUI by cell type x stage")
    plt.tight_layout()
    _savefig(fig, output_path)


def plot_shortening_lengthening_summary(trend_df: pd.DataFrame, output_path: Path) -> None:
    """Direction-count bar + per-celltype slope bar, colored by shortening/lengthening."""
    if trend_df.empty:
        log.warning("plot_shortening_lengthening_summary: empty, skipping")
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, max(6, 0.3 * len(trend_df))))

    counts = trend_df["direction"].value_counts()
    dir_colors = {"decreasing": "#2ca02c", "increasing": "#d62728"}
    axes[0].bar(counts.index.astype(str), counts.values,
                color=[dir_colors.get(d, "#7f7f7f") for d in counts.index])
    axes[0].set_ylabel("n cell types")
    axes[0].set_title("Shortening (decreasing) vs lengthening (increasing)")
    axes[0].grid(True, alpha=0.2, axis="y")

    sorted_df = trend_df.dropna(subset=["slope"]).sort_values("slope")
    bar_colors = ["#2ca02c" if s < 0 else "#d62728" for s in sorted_df["slope"]]
    axes[1].barh(sorted_df["celltype"], sorted_df["slope"], color=bar_colors)
    axes[1].set_xlabel("length-trend slope (PDUI / stage)")
    axes[1].tick_params(axis="y", labelsize=6)
    axes[1].set_title("Slope by cell type")
    axes[1].axvline(0, color="black", linewidth=0.8)
    axes[1].grid(True, alpha=0.2, axis="x")

    plt.tight_layout()
    _savefig(fig, output_path)


def plot_recurrent_genes(recurrent_df: pd.DataFrame, output_path: Path, top_n: int = 20) -> None:
    """Horizontal bar of the top recurrent switch genes by n_celltypes_trending."""
    if recurrent_df.empty:
        log.warning("plot_recurrent_genes: empty, skipping")
        return
    top = recurrent_df.head(top_n)
    dir_colors = {"decreasing": "#2ca02c", "increasing": "#d62728"}
    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * len(top))))
    colors = [dir_colors.get(d, "#7f7f7f") for d in top["dominant_direction"]]
    ax.barh(top["gene_id"].astype(str), top["n_celltypes_trending"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("n cell types trending")
    ax.set_title(f"Top {len(top)} recurrent switch genes")
    plt.tight_layout()
    _savefig(fig, output_path)


# ---------------------------------------------------------------------------
# Cumulative summary
# ---------------------------------------------------------------------------


def write_markdown_summary(summary: dict, output_path: Path) -> None:
    """Render the cumulative summary dict to a markdown report."""
    lines = ["# PeakATail Cumulative Sweep Analysis", ""]
    lines.append(f"Sweep root: `{summary.get('sweep_root')}`")
    lines.append("")

    lines.append("## 1. Peak-calling strategy comparison")
    strat = summary.get("strategy_comparison") or []
    cutoffs = summary.get("cutoffs") or DEFAULT_CUTOFFS
    if strat:
        header = "| run | strategy | n_PAS | " + " | ".join(f"F1@{c}bp" for c in cutoffs) + " | atlas match |"
        sep = "|---|---|---|" + "---|" * len(cutoffs) + "---|"
        lines.append(header)
        lines.append(sep)
        for r in strat:
            f1s = " | ".join(_fmt(r.get(f"f1_{c}"), ".3f") for c in cutoffs)
            lines.append(
                f"| {r.get('run')} | {r.get('peak_strategy')} | {_fmt(r.get('n_pas'), ',')} "
                f"| {f1s} | {_fmt(r.get('atlas_match_rate'), '.1%')} |"
            )
    else:
        lines.append("_no grid runs found_")
    lines.append("")

    lines.append("## 2. Internal-priming / fasta-filter axis")
    ip = summary.get("ip_filter_axis") or []
    if ip:
        lines.append("| ip_mode | n_PAS | n_PAS dropped vs off | % flagged internal-priming |")
        lines.append("|---|---|---|---|")
        for r in ip:
            lines.append(
                f"| {r.get('ip_mode')} | {_fmt(r.get('n_pas'), ',')} "
                f"| {_fmt(r.get('n_pas_dropped_vs_off'), ',')} "
                f"| {_fmt(r.get('pct_flagged_internal_priming'))}% |"
            )
    else:
        lines.append("_no lambda_gradient ip-axis runs found_")
    lines.append("")

    lines.append("## 3. Trim sensitivity")
    trim = summary.get("trim_comparison") or []
    if trim:
        lines.append(
            "| branch | max_gene_distance | utr_multiplier | include_extended "
            "| n_PAS | n_cells_total | mean n_clusters |"
        )
        lines.append("|---|---|---|---|---|---|---|")
        for r in trim:
            lines.append(
                f"| {r.get('branch')} | {r.get('max_gene_distance')} | {r.get('utr_multiplier')} "
                f"| {r.get('include_extended')} | {_fmt(r.get('n_pas'), ',')} "
                f"| {_fmt(r.get('n_cells_total'), ',')} | {_fmt(r.get('n_clusters_mean'))} |"
            )
    else:
        lines.append("_no reannotate trim branches found_")
    lines.append("")

    lines.append("## 4. Clustering comparison")
    clus = summary.get("clustering_comparison") or []
    if clus:
        lines.append("| branch | method | resolution | n_neighbors | mean n_clusters |")
        lines.append("|---|---|---|---|---|")
        for r in clus:
            lines.append(
                f"| {r.get('branch')} | {r.get('cluster_method')} | {r.get('resolution')} "
                f"| {r.get('n_neighbors')} | {_fmt(r.get('n_clusters_mean'))} |"
            )
    else:
        lines.append("_no reannotate clustering branches found_")
    lines.append("")

    gex = summary.get("gex_concordance_summary") or {}
    if gex:
        lines.append(f"Cohort GEX concordance (n={gex.get('n_gsm')} GSMs):")
        for k, v in gex.items():
            if isinstance(v, dict):
                lines.append(f"- `{k}`: mean={v['mean']}, median={v['median']}, std={v['std']}")
        lines.append("")

    lines.append("## 5. Switch cumulative (biological headline)")
    sw = summary.get("switch") or {}
    sl = sw.get("shortening_vs_lengthening") or {}
    lines.append(
        f"- **{sl.get('n_celltypes', 0)} cell types** analyzed: "
        f"**{sl.get('n_shortening', 0)} shortening** (3'UTR PDUI decreasing with stage), "
        f"**{sl.get('n_lengthening', 0)} lengthening**, "
        f"{sl.get('n_flat_or_other', 0)} flat/other "
        f"(mean slope {_fmt(sl.get('slope_mean'))}, median {_fmt(sl.get('slope_median'))})"
    )
    lines.append("")

    lines.append("### Recurrent switch genes (top 20)")
    rec = sw.get("recurrent_genes_top20") or []
    params = sw.get("recurrent_gene_params") or {}
    lines.append(
        f"_Genes with `|spearman| >= {params.get('spearman_thresh')}` in "
        f">= {params.get('min_celltypes')} cell types._"
    )
    if rec:
        lines.append("| gene_id | n_celltypes_trending | dominant_direction | mean_slope | mean_abs_spearman |")
        lines.append("|---|---|---|---|---|")
        for r in rec:
            lines.append(
                f"| {r.get('gene_id')} | {r.get('n_celltypes_trending')} "
                f"| {r.get('dominant_direction')} | {_fmt(r.get('mean_slope'), '.4f')} "
                f"| {_fmt(r.get('mean_abs_spearman'), '.3f')} |"
            )
    else:
        lines.append("_none met the recurrence threshold_")
    lines.append("")

    lines.append("### Fisher significant-PAS hit counts by stage contrast")
    fc = sw.get("fisher_hits_by_contrast") or []
    if fc:
        lines.append("| contrast | n_celltypes | total_tests | total_sig (FDR<0.05) | mean frac sig |")
        lines.append("|---|---|---|---|---|")
        for r in fc:
            lines.append(
                f"| {r.get('contrast')} | {r.get('n_celltypes')} | {_fmt(r.get('total_tests'), ',')} "
                f"| {_fmt(r.get('total_sig'), ',')} | {_fmt(r.get('mean_frac_sig'), '.4f')} |"
            )
    else:
        lines.append("_no Fisher differential results found_")
    lines.append("")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines))


def analyze_sweep(
    sweep_root: Path,
    out_dir: Path,
    *,
    cutoffs: Sequence[int] = DEFAULT_CUTOFFS,
    fdr: float = 0.05,
    min_recurrent_celltypes: int = 3,
    spearman_thresh: float = 0.8,
) -> dict:
    """Run the full cumulative sweep analysis and write tables/figures/summary.

    Args:
        sweep_root: Root of the finished sweep (contains ``runs/``).
        out_dir: Output directory; ``tables/`` and ``figures/`` subdirs plus
            ``cumulative_summary.{json,md}`` are written here.
        cutoffs: Distance cutoffs (bp) for the strategy comparison table.
        fdr: FDR threshold for Fisher significance counts.
        min_recurrent_celltypes: Minimum cell-type count for a gene to be
            reported as "recurrent" in the switch-gene table.
        spearman_thresh: Minimum ``|spearman|`` for a per-celltype gene trend
            to count toward recurrence.

    Returns:
        The cumulative summary dict (also written to
        ``cumulative_summary.json``).
    """
    sweep_root = Path(sweep_root)
    out_dir = Path(out_dir)
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    log.info("analyze_sweep: 1/6 peak-calling strategy comparison")
    strategy_df = harvest_strategy_comparison(sweep_root, cutoffs)
    strategy_df.to_csv(tables_dir / "strategy_comparison.csv", index=False)

    ip_df = harvest_ip_filter_axis(strategy_df)
    ip_df.to_csv(tables_dir / "ip_filter_axis.csv", index=False)

    log.info("analyze_sweep: 2/6 trim sensitivity")
    trim_df = harvest_trim_comparison(sweep_root)
    trim_df.to_csv(tables_dir / "trim_comparison.csv", index=False)

    log.info("analyze_sweep: 3/6 clustering comparison + GEX concordance")
    cluster_df = harvest_clustering_comparison(sweep_root)
    cluster_df.to_csv(tables_dir / "clustering_comparison.csv", index=False)

    gex_df = harvest_gex_concordance(sweep_root)
    gex_df.to_csv(tables_dir / "gex_concordance.csv", index=False)
    gex_summary = summarize_gex_concordance(gex_df)

    log.info("analyze_sweep: 4/6 switch trend harvesting (%d cell types)", len(_list_celltypes(sweep_root)))
    trend_df = harvest_switch_trends(sweep_root)
    trend_df.to_csv(tables_dir / "switch_trend_by_celltype.csv", index=False)
    shortening_summary = summarize_shortening_lengthening(trend_df)

    pdui_matrix = celltype_stage_pdui_matrix(trend_df)
    pdui_matrix.to_csv(tables_dir / "celltype_stage_pdui_matrix.csv")

    log.info("analyze_sweep: 5/6 fisher hit counts + recurrent genes")
    fisher_df = harvest_fisher_hit_counts(sweep_root, fdr=fdr)
    fisher_df.to_csv(tables_dir / "fisher_hits_by_celltype_contrast.csv", index=False)
    fisher_contrast_summary = summarize_fisher_hits(fisher_df)
    fisher_contrast_summary.to_csv(tables_dir / "fisher_hits_by_contrast_summary.csv", index=False)

    recurrent_df = harvest_recurrent_switch_genes(
        sweep_root, min_celltypes=min_recurrent_celltypes, spearman_thresh=spearman_thresh
    )
    recurrent_df.to_csv(tables_dir / "recurrent_switch_genes.csv", index=False)

    log.info("analyze_sweep: 6/6 figures")
    plot_strategy_precision_recall(strategy_df, figures_dir / "strategy_precision_recall", cutoffs)
    plot_trim_cluster_sensitivity(trim_df, cluster_df, figures_dir / "trim_cluster_sensitivity")
    plot_celltype_stage_heatmap(pdui_matrix, figures_dir / "celltype_stage_pdui_heatmap")
    plot_shortening_lengthening_summary(trend_df, figures_dir / "shortening_lengthening_summary")
    plot_recurrent_genes(recurrent_df, figures_dir / "recurrent_switch_genes")

    summary = {
        "sweep_root": str(sweep_root),
        "cutoffs": list(cutoffs),
        "strategy_comparison": strategy_df.to_dict(orient="records"),
        "ip_filter_axis": ip_df.to_dict(orient="records"),
        "trim_comparison": trim_df.to_dict(orient="records"),
        "clustering_comparison": cluster_df.to_dict(orient="records"),
        "gex_concordance_summary": gex_summary,
        "switch": {
            "shortening_vs_lengthening": shortening_summary,
            "n_celltypes_analyzed": int(len(trend_df)),
            "fisher_hits_by_contrast": fisher_contrast_summary.to_dict(orient="records"),
            "recurrent_genes_top20": recurrent_df.head(20).to_dict(orient="records"),
            "recurrent_gene_params": {
                "min_celltypes": min_recurrent_celltypes,
                "spearman_thresh": spearman_thresh,
            },
        },
    }

    with open(out_dir / "cumulative_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=_json_default)

    write_markdown_summary(summary, out_dir / "cumulative_summary.md")
    log.info("analyze_sweep: done -> %s", out_dir)
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m ema.benchmark.sweep_analysis",
        description="Cumulative cross-experiment analysis over a finished PeakATail sweep.",
    )
    parser.add_argument("sweep_root", type=Path, help="Root of the finished sweep (contains runs/)")
    parser.add_argument(
        "--out", type=Path, default=None,
        help="Output directory (default: <sweep_root>/logs/cumulative_analysis)",
    )
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold for Fisher hit counts")
    parser.add_argument(
        "--min-recurrent-celltypes", type=int, default=3,
        help="Minimum cell-type count for a gene to be reported as recurrent",
    )
    parser.add_argument(
        "--spearman-thresh", type=float, default=0.8,
        help="Minimum |spearman| for a per-celltype gene trend to count toward recurrence",
    )
    parser.add_argument(
        "--cutoffs", type=int, nargs="+", default=list(DEFAULT_CUTOFFS),
        help="Distance cutoffs (bp) for the strategy comparison table",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug-level logging")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    warnings.filterwarnings("ignore")

    out_dir = args.out or (args.sweep_root / "logs" / "cumulative_analysis")
    analyze_sweep(
        args.sweep_root,
        out_dir,
        cutoffs=tuple(args.cutoffs),
        fdr=args.fdr,
        min_recurrent_celltypes=args.min_recurrent_celltypes,
        spearman_thresh=args.spearman_thresh,
    )
    log.info("Cumulative analysis complete: %s", out_dir / "cumulative_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
