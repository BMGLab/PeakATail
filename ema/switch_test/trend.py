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
    "distal_proportion_trend",
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


def distal_proportion_trend(
    df,
    stage_order: list[str],
    *,
    stage_col: str = "cluster",
    value_col: str = "proportion",
    gene_col: str = "gene_id",
    rank_col: str = "rank",
    cell_col: str = "cell",
    min_pas: int = 2,
) -> dict[str, Any]:
    """Ordered-stage trend of DISTAL-PAS usage, restricted to multi-PAS genes.

    This is the scientifically-meaningful trend for the per-PAS
    ``proportion`` strategy. A single-PAS gene has ``proportion`` == 1.0 in
    every covered cell by construction — that's not APA usage information,
    it's just "the only PAS got all the reads because there was only one
    PAS". Feeding those rows into a naive whole-table mean floods the trend
    with a trivial constant and hides any real multi-PAS signal (this is
    the mechanism behind the previously-reported "flat proportion trend"
    finding, compounded by the writer padding bug — see
    :mod:`ema.quantification.strategies.proportion`).

    This function:
      1. Restricts to genes with >= *min_pas* distinct PAS (``n_PAS`` is
         read off ``rank_col``'s per-gene max, i.e. the same rank column the
         corrected proportion writer now populates).
      2. Within each remaining (gene, cell), keeps only the row for the
         most-DISTAL PAS (max ``rank_col`` within the gene) as the per-cell
         APA readout — a PDUI-like "how much of this multi-PAS gene's usage
         sits at the distal end" index.
      3. Delegates the actual trend fit to :func:`pdui_trend` (linear slope
         + Spearman rho across ``stage_order``). ``spearman`` is then forced
         to NaN whenever ``n_stages < 3`` — a 2-point Spearman rho is always
         exactly +/-1 by construction (there are only two possible
         orderings of two points), so it carries no information beyond the
         sign of the slope; this is the report-analyst convention for this
         trend. ``n_stages`` is always reported regardless, so a caller can
         see exactly how many stages backed the (possibly NaN) rho.

    Args:
        df: Long-format ``proportion`` table (one row per gene/pas/cell) —
            e.g. the direct output of ``ema switch length --strategy
            proportion``, after the padding fix (zero-coverage rows already
            omitted; a stale/pre-fix table with padded rows would silently
            include the padding in the "distal" mean, since a padded row's
            ``rank`` was always the same broken sentinel — this function
            does not itself detect padding, it composes with the corrected
            writer).
        stage_order: Ordered stage labels (see :func:`build_stage_rank`).
        stage_col: Stage/cluster column name.
        value_col: The proportion (or any per-PAS numeric) column name.
        gene_col: Gene identifier column name.
        rank_col: Proximal(low)->distal(high) rank column name.
        cell_col: Cell identifier column name (used only to report
            ``n_cells``).
        min_pas: Minimum distinct PAS per gene to be considered informative
            (default 2 — excludes trivial single-PAS genes).

    Returns:
        Same dict shape as :func:`pdui_trend` (``n_stages``, ``slope``,
        ``spearman``, ``direction``, ``mean_by_stage``, ``value_col``),
        plus ``n_genes`` (distinct multi-PAS genes contributing to the
        distal-usage series) and ``n_cells`` (distinct covered cells
        contributing).

    Raises:
        ValueError: if any required column is missing.
    """
    import pandas as pd

    required = {gene_col, rank_col, cell_col, stage_col, value_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"distal_proportion_trend missing columns: {sorted(missing)}; "
            f"have {list(df.columns)}"
        )

    work = df[[gene_col, rank_col, cell_col, stage_col, value_col]].copy()
    work[rank_col] = pd.to_numeric(work[rank_col], errors="coerce")
    work = work.dropna(subset=[rank_col])

    def _empty_result() -> dict[str, Any]:
        result = trend_from_stage_means({})
        result["mean_by_stage"] = {}
        result["value_col"] = value_col
        result["n_genes"] = 0
        result["n_cells"] = 0
        return result

    if work.empty:
        return _empty_result()

    max_rank_by_gene = work.groupby(gene_col)[rank_col].transform("max")
    work = work[max_rank_by_gene >= min_pas]
    if work.empty:
        return _empty_result()

    is_distal = work[rank_col] == work.groupby(gene_col)[rank_col].transform("max")
    distal = work[is_distal]
    if distal.empty:
        return _empty_result()

    result = pdui_trend(distal, stage_order, stage_col=stage_col, value_col=value_col)
    result["n_genes"] = int(distal[gene_col].nunique())
    result["n_cells"] = int(distal[cell_col].nunique())
    # A 2-point Spearman rho is +/-1 by construction (only two possible
    # orderings) and carries no monotonicity information — null it out,
    # per the report-analyst convention, but keep n_stages so callers can
    # see why.
    if result["n_stages"] < 3:
        result["spearman"] = float("nan")
    return result
