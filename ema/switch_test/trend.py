"""A3: ordered-stage APA trend aggregation.

``ema switch length`` computes a per-cluster PDUI (or proportion / entropy)
value, but treats clusters as an unordered set. When the clusters are actually
an ORDERED series — disease stages (Normal → StageI → … → Metastasis),
timepoints, differentiation pseudo-stages — the scientific question is the
*trend* across that order: does the 3'UTR shorten (PDUI = distal fraction
falls) monotonically with progression?

This module is the genuine package gap the Laughney reanalysis surfaced
(previously living only in ``scripts/b3_stage_celltype_switch.py``). It is a
pure aggregation over an existing PDUI long table plus a caller-supplied stage
order — no statistics are re-implemented, no biology is fabricated: the slope
and monotonicity are definitional summaries of the per-stage means.

Direction convention: the reported ``direction`` is generic
(``increasing`` / ``decreasing`` / ``flat``). For a PDUI value column (distal
fraction), ``decreasing`` == 3'UTR **shortening** with stage progression and
``increasing`` == lengthening; callers that know the value is PDUI can relabel.
"""
from __future__ import annotations

from typing import Any

__all__ = [
    "build_stage_rank",
    "trend_from_stage_means",
    "pdui_trend",
    "pdui_trend_by_gene",
]


def build_stage_rank(stage_order: list[str]) -> dict[str, int]:
    """Map each stage label to an integer rank from a caller-ordered list.

    The first stage gets rank 0, the next rank 1, and so on. Duplicate labels
    keep their first-seen rank (so a caller may repeat a label to tie stages,
    e.g. ``["Normal", "StageIA", "StageIB", ...]`` where IA/IB should tie is
    expressed by passing the same label twice — but typically callers pass
    distinct labels). Empty/whitespace labels are ignored.

    Raises:
        ValueError: if ``stage_order`` has no usable labels.
    """
    rank: dict[str, int] = {}
    next_rank = 0
    for label in stage_order:
        key = str(label).strip()
        if not key:
            continue
        if key not in rank:
            rank[key] = next_rank
            next_rank += 1
    if not rank:
        raise ValueError("stage_order contained no usable stage labels")
    return rank


def trend_from_stage_means(mean_by_rank: dict[int, float]) -> dict[str, Any]:
    """Fit a linear trend to (rank -> mean value) points.

    Args:
        mean_by_rank: mapping ordered-stage rank -> mean value at that stage.

    Returns:
        dict with:
            ``n_stages``   number of distinct ranks with a finite mean
            ``slope``      least-squares slope of mean vs rank (NaN if <2 pts)
            ``spearman``   Spearman rho of rank vs mean (NaN if <2 pts or
                           constant); a rank-monotonicity summary
            ``direction``  ``"increasing"`` / ``"decreasing"`` / ``"flat"`` /
                           ``"undetermined"`` (from the sign of slope)
    """
    import numpy as np

    ranks = sorted(r for r, v in mean_by_rank.items() if v is not None and np.isfinite(v))
    vals = [float(mean_by_rank[r]) for r in ranks]
    n = len(ranks)
    out: dict[str, Any] = {
        "n_stages": n, "slope": float("nan"),
        "spearman": float("nan"), "direction": "undetermined",
    }
    if n < 2:
        return out
    x = np.asarray(ranks, dtype=float)
    y = np.asarray(vals, dtype=float)
    slope = float(np.polyfit(x, y, 1)[0])
    out["slope"] = slope
    # Spearman via ranks of x and y (x is already strictly increasing → its
    # rank order is identity). Constant y → undefined correlation → NaN.
    if np.ptp(y) == 0:
        out["spearman"] = float("nan")
    else:
        try:
            from scipy.stats import spearmanr
            rho, _ = spearmanr(x, y)
            out["spearman"] = float(rho)
        except Exception:
            out["spearman"] = float("nan")
    eps = 1e-12
    out["direction"] = ("increasing" if slope > eps
                        else "decreasing" if slope < -eps else "flat")
    return out


def _stage_means(df, stage_order, stage_col, value_col):
    """Return (mean_by_label, mean_by_rank) restricted to known ordered stages."""
    import numpy as np
    import pandas as pd

    rank = build_stage_rank(stage_order)
    if stage_col not in df.columns or value_col not in df.columns:
        raise ValueError(
            f"table missing required columns: need {stage_col!r} and {value_col!r}, "
            f"have {list(df.columns)}"
        )
    sub = df[[stage_col, value_col]].copy()
    sub[value_col] = pd.to_numeric(sub[value_col], errors="coerce")
    sub = sub.dropna(subset=[value_col])
    sub = sub[sub[stage_col].astype(str).str.strip().isin(rank)]
    mean_by_label = sub.groupby(sub[stage_col].astype(str).str.strip())[value_col].mean()
    mean_by_rank: dict[int, float] = {}
    for label, mean in mean_by_label.items():
        r = rank[label]
        # If two labels tie on a rank, average their means.
        if r in mean_by_rank:
            mean_by_rank[r] = (mean_by_rank[r] + float(mean)) / 2.0
        else:
            mean_by_rank[r] = float(mean)
    return mean_by_label.round(6).to_dict(), mean_by_rank


def pdui_trend(
    df,
    stage_order: list[str],
    *,
    stage_col: str = "cluster",
    value_col: str = "pdui",
) -> dict[str, Any]:
    """Overall ordered-stage trend of ``value_col`` across ``stage_order``.

    Args:
        df: long PDUI table (e.g. one row per cell/gene) with a stage column
            and a numeric value column.
        stage_order: ordered stage labels defining the progression axis. Rows
            whose stage is not in this list are ignored (with no error).
        stage_col: name of the stage/cluster column (default ``"cluster"``).
        value_col: name of the value column (default ``"pdui"``).

    Returns:
        dict from :func:`trend_from_stage_means` plus ``mean_by_stage``
        (label -> mean) and echoed ``value_col``.
    """
    mean_by_label, mean_by_rank = _stage_means(df, stage_order, stage_col, value_col)
    result = trend_from_stage_means(mean_by_rank)
    result["mean_by_stage"] = mean_by_label
    result["value_col"] = value_col
    return result


def pdui_trend_by_gene(
    df,
    stage_order: list[str],
    *,
    stage_col: str = "cluster",
    value_col: str = "pdui",
    gene_col: str = "gene_id",
):
    """Per-gene ordered-stage trend table.

    Returns a DataFrame indexed by gene with columns
    ``[n_stages, slope, spearman, direction]`` — one linear trend per gene.
    Genes with fewer than 2 covered stages get NaN slope / ``undetermined``.

    Raises:
        ValueError: if ``gene_col`` is absent (use :func:`pdui_trend` instead).
    """
    import pandas as pd

    if gene_col not in df.columns:
        raise ValueError(
            f"pdui_trend_by_gene needs a {gene_col!r} column; have {list(df.columns)}"
        )
    rows = []
    for gene, sub in df.groupby(df[gene_col].astype(str)):
        if not gene:
            continue
        res = pdui_trend(sub, stage_order, stage_col=stage_col, value_col=value_col)
        rows.append({
            "gene_id": gene,
            "n_stages": res["n_stages"],
            "slope": res["slope"],
            "spearman": res["spearman"],
            "direction": res["direction"],
        })
    if not rows:
        return pd.DataFrame(columns=["gene_id", "n_stages", "slope", "spearman", "direction"])
    return pd.DataFrame(rows).set_index("gene_id").sort_values("slope")
