"""A3: ordered-stage APA trend aggregation.

``peakatail switch length`` computes a per-cluster PDUI (or proportion / entropy)
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
    "derive_distal_fraction",
    "consensus_trend_by_gene",
    "consensus_trend_summary",
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


def derive_distal_fraction(
    df,
    *,
    gene_col: str = "gene_id",
    cell_col: str = "cell",
    stage_col: str = "cluster",
    rank_col: str = "rank",
    proportion_col: str = "proportion",
):
    """Derive a per-(gene, cell) DISTAL-USAGE scalar from a per-PAS proportion table.

    ``peakatail switch length --strategy proportion`` reports one row per (gene,
    PAS, cell) with the fraction of that gene's reads landing on that PAS
    (``proportion``) and the PAS's rank along the transcript (``1`` = most
    proximal ... max rank = most distal). This collapses that per-PAS table
    to one row per (gene, cell): the proportion carried by the DISTAL-most
    PAS (max ``rank`` within the group), named ``value`` in the returned
    frame so it mirrors PDUI's convention — a *falling* ``value`` across
    stages means less distal usage, i.e. 3'UTR **shortening**, exactly like
    a falling PDUI.

    Ties on the max rank within a (gene, cell) group are averaged.

    Args:
        df: long per-PAS proportion table with ``gene_col``, ``cell_col``,
            ``stage_col``, ``rank_col`` and ``proportion_col`` columns.
        gene_col: gene identifier column.
        cell_col: cell identifier column.
        stage_col: stage/cluster column, carried through unchanged (all rows
            for one cell are expected to share one stage).
        rank_col: PAS rank-along-transcript column (larger = more distal).
        proportion_col: per-PAS fraction-of-gene-reads column.

    Returns:
        DataFrame with columns ``[gene_id, cell, cluster, value]``, one row
        per (gene, cell) present in ``df`` after coercion/NaN-dropping.

    Raises:
        ValueError: if any required input column is missing.
    """
    import pandas as pd

    required = [gene_col, cell_col, stage_col, rank_col, proportion_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"derive_distal_fraction needs columns {required}; missing {missing}, "
            f"have {list(df.columns)}"
        )
    sub = df[required].copy()
    sub[rank_col] = pd.to_numeric(sub[rank_col], errors="coerce")
    sub[proportion_col] = pd.to_numeric(sub[proportion_col], errors="coerce")
    sub = sub.dropna(subset=[rank_col, proportion_col])
    if sub.empty:
        return pd.DataFrame(columns=["gene_id", "cell", "cluster", "value"])

    rows = []
    for (gene, cell), g in sub.groupby([gene_col, cell_col], sort=False):
        max_rank = g[rank_col].max()
        at_max = g[g[rank_col] == max_rank]
        rows.append({
            "gene_id": gene,
            "cell": cell,
            "cluster": at_max[stage_col].iloc[0],
            "value": float(at_max[proportion_col].mean()),
        })
    return pd.DataFrame(rows)


def consensus_trend_by_gene(
    per_gene_by_metric: dict[str, Any],
    *,
    vote_metrics: list[str],
    min_agree: int,
) -> Any:
    """Consensus per-gene trend across multiple length-strategy metrics.

    Args:
        per_gene_by_metric: metric name -> per-gene trend DataFrame as
            returned by :func:`pdui_trend_by_gene` (indexed by ``gene_id``,
            with ``slope``/``direction`` columns), one entry per metric
            (e.g. ``{"pdui": ..., "proportion": ..., "entropy": ...}``).
        vote_metrics: subset of ``per_gene_by_metric`` keys whose direction
            counts toward the shortening/lengthening vote. By convention
            (decided by the caller, not this function) ``"entropy"`` is
            usually excluded — rising entropy means more even PAS usage, not
            cleanly "lengthening" like a rising PDUI/distal-fraction.
        min_agree: minimum number of *voting* metrics that must agree on
            "increasing" or "decreasing" (and outnumber the opposite
            direction) for a gene's consensus to be called; otherwise
            ``"ambiguous"``.

    Returns:
        DataFrame indexed by ``gene_id`` (union of genes across all metrics)
        with ``<metric>_direction`` / ``<metric>_slope`` per metric, plus
        ``n_voting`` (voting metrics with a determined direction for this
        gene), ``n_agree`` (size of the majority bucket among those),
        ``mean_slope`` (mean of all available per-metric slopes), and
        ``consensus_direction`` (``"shortening"`` / ``"lengthening"`` /
        ``"ambiguous"``), ``agree`` (bool, ``consensus_direction`` != ambiguous).

    Raises:
        ValueError: if ``per_gene_by_metric`` is empty.
    """
    import numpy as np
    import pandas as pd

    if not per_gene_by_metric:
        raise ValueError("consensus_trend_by_gene needs at least one metric")

    all_genes: set = set()
    for gdf in per_gene_by_metric.values():
        all_genes.update(gdf.index)

    rows = []
    for gene in sorted(all_genes):
        row: dict[str, Any] = {"gene_id": gene}
        slopes: list[float] = []
        inc = dec = 0
        for metric, gdf in per_gene_by_metric.items():
            if gene in gdf.index:
                direction = gdf.loc[gene, "direction"]
                slope = gdf.loc[gene, "slope"]
            else:
                direction = "undetermined"
                slope = float("nan")
            row[f"{metric}_direction"] = direction
            row[f"{metric}_slope"] = float(slope) if slope is not None else float("nan")
            if slope is not None and np.isfinite(slope):
                slopes.append(float(slope))
            if metric in vote_metrics:
                if direction == "increasing":
                    inc += 1
                elif direction == "decreasing":
                    dec += 1
        if dec >= min_agree and dec > inc:
            consensus = "shortening"
        elif inc >= min_agree and inc > dec:
            consensus = "lengthening"
        else:
            consensus = "ambiguous"
        row["n_voting"] = inc + dec
        row["n_agree"] = max(inc, dec)
        row["mean_slope"] = float(np.mean(slopes)) if slopes else float("nan")
        row["consensus_direction"] = consensus
        row["agree"] = consensus != "ambiguous"
        rows.append(row)

    cols = ["gene_id"]
    for metric in per_gene_by_metric:
        cols += [f"{metric}_direction", f"{metric}_slope"]
    cols += ["n_voting", "n_agree", "mean_slope", "consensus_direction", "agree"]
    if not rows:
        return pd.DataFrame(columns=cols).set_index("gene_id")
    return pd.DataFrame(rows)[cols].set_index("gene_id")


def consensus_trend_summary(
    combined,
    *,
    metrics: list[str],
    vote_metrics: list[str],
    min_agree: int,
) -> dict[str, Any]:
    """Cohort-level summary of a :func:`consensus_trend_by_gene` table.

    Args:
        combined: the per-gene consensus DataFrame from
            :func:`consensus_trend_by_gene`.
        metrics: all metric names present as ``<metric>_direction`` columns
            (used to compute every pairwise agreement rate, including
            non-voting metrics like entropy, for visibility).
        vote_metrics: the metrics that counted toward the consensus vote.
        min_agree: the agreement threshold used to build ``combined``.

    Returns:
        dict with ``n_genes``, ``n_shortening``, ``n_lengthening``,
        ``n_ambiguous``, echoed ``metrics``/``vote_metrics``/
        ``combine_min_agree``, and ``pairwise_agreement``: for every metric
        pair, ``{"n_genes_compared": int, "agreement_rate": float}`` — the
        rate, among genes where BOTH metrics have a determined
        (increasing/decreasing) direction, that the two directions match
        (NaN if no such genes).
    """
    import itertools

    counts = combined["consensus_direction"].value_counts().to_dict() if len(combined) else {}
    out: dict[str, Any] = {
        "n_genes": int(len(combined)),
        "n_shortening": int(counts.get("shortening", 0)),
        "n_lengthening": int(counts.get("lengthening", 0)),
        "n_ambiguous": int(counts.get("ambiguous", 0)),
        "metrics": list(metrics),
        "vote_metrics": list(vote_metrics),
        "combine_min_agree": int(min_agree),
        "pairwise_agreement": {},
    }
    determined = {"increasing", "decreasing"}
    for m1, m2 in itertools.combinations(metrics, 2):
        c1 = combined.get(f"{m1}_direction")
        c2 = combined.get(f"{m2}_direction")
        if c1 is None or c2 is None:
            continue
        both = c1.isin(determined) & c2.isin(determined)
        n_both = int(both.sum())
        rate = float((c1[both] == c2[both]).mean()) if n_both else float("nan")
        out["pairwise_agreement"][f"{m1}_vs_{m2}"] = {
            "n_genes_compared": n_both,
            "agreement_rate": rate,
        }
    return out
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
            e.g. the direct output of ``peakatail switch length --strategy
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
