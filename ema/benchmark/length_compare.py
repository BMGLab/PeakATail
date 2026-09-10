"""Compare multiple 3'UTR length quantification strategies on a single h5ad.

Production module for benchmarking :mod:`ema.quantification` strategies
(``classic``, ``proportion``, ``shannon``) with per-gene and per-isoform
aggregation against a clustered AnnData, emitting per-strategy summary
metrics.

Typical usage::

    from pathlib import Path
    from ema.benchmark.length_compare import compare_length_strategies

    summary = compare_length_strategies(
        h5ad_path=Path("clusters.h5ad"),
        pasbed_path=Path("pasbed.bed"),
        gtf_path=Path("Homo_sapiens.GRCh38.99.gtf"),
        output_dir=Path("reports/length_compare"),
        strategies_and_aggs=[
            ("classic",    "per_gene"),
            ("proportion", "per_gene"),
            ("proportion", "per_isoform"),
            ("shannon",    "per_gene"),
        ],
        threads=4,
    )
    print(summary)

Notes
-----
- The h5ad and the ``pas_isoform_map`` are built once per ``isoform_agg``
  level (not once per strategy) to avoid re-parsing the GTF for every
  strategy with the same aggregation level.
- ``per_isoform`` aggregation requires a GTF; if ``gtf_path`` is ``None``
  or does not exist, that aggregation level is skipped with a warning.
- The ``n_units_with_inter_cluster_shift`` metric counts genes/transcripts
  whose max score across clusters minus min score across clusters exceeds 0.1
  in the strategy's native score column (``pdui``, ``proportion``, or
  ``entropy``).
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _to_count_matrix(adata: ad.AnnData) -> pd.DataFrame:
    """Return a PAS × cells DataFrame of integer counts.

    Args:
        adata: Clustered AnnData with sparse or dense ``X``.

    Returns:
        DataFrame with shape ``(n_pas, n_cells)``; index is integer PAS IDs
        where possible, falling back to raw var-names.
    """
    X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    try:
        pas_index = [int(v) for v in adata.var_names]
    except (TypeError, ValueError):
        pas_index = list(adata.var_names)
    return pd.DataFrame(X.T, index=pas_index, columns=list(adata.obs_names))


def _build_per_gene_map(
    adata: ad.AnnData,
) -> dict[int, list[tuple[str, str, int, int, int]]]:
    """Build a minimal ``pas_isoform_map`` from ``adata.var['gene_id']``.

    Each PAS maps to a single ``(gene_id, "_gene_", rank=0, 1, 1)`` entry —
    sufficient for per_gene aggregation which does not use transcript
    information.

    Args:
        adata: AnnData with a ``gene_id`` column in ``var``.

    Returns:
        ``{pas_id_int: [(gene_id, "_gene_", 0, rank, 1)]}`` mapping, where
        rank is assigned by occurrence order within the gene (0-indexed).
        Returns an empty dict if ``gene_id`` is absent from ``adata.var``.
    """
    if "gene_id" not in adata.var.columns:
        log.warning("adata.var has no 'gene_id'; per_gene map will be empty")
        return {}

    gene_id_col = adata.var["gene_id"]
    gene_counter: dict[str, int] = {}
    result: dict[int, list[tuple[str, str, int, int, int]]] = {}
    for pas_id_str, gene_id in gene_id_col.items():
        if not gene_id or (isinstance(gene_id, float) and pd.isna(gene_id)):
            continue
        gene_id_str = str(gene_id)
        try:
            pas_id_int = int(pas_id_str)
        except (TypeError, ValueError):
            continue
        rank = gene_counter.get(gene_id_str, 0)
        gene_counter[gene_id_str] = rank + 1
        result[pas_id_int] = [(gene_id_str, "_gene_", 0, rank, 1)]
    log.info("_build_per_gene_map: %d PAS mapped to %d genes", len(result), len(gene_counter))
    return result


def _inter_cluster_shift(
    df: pd.DataFrame,
    score_col: str,
    unit_cols: list[str],
    cluster_labels: pd.Series,
) -> int:
    """Count quantification units with large between-cluster score variation.

    For each unique ``(unit_cols)`` combination, the per-cluster mean score
    is computed.  A unit is counted if ``max_cluster_mean - min_cluster_mean
    > 0.1``.

    Args:
        df: Long-format strategy output containing a ``cell`` column and
            the ``score_col``.
        score_col: Name of the numeric score column (``pdui``,
            ``proportion``, ``entropy``).
        unit_cols: Column(s) that identify a unique quantification unit
            (e.g. ``["gene_id"]`` or ``["gene_id", "transcript_id"]``).
        cluster_labels: Series mapping cell barcodes to cluster labels
            (aligned to ``df["cell"]``).

    Returns:
        Number of units with inter-cluster range > 0.1.
    """
    if score_col not in df.columns or "cell" not in df.columns:
        return 0
    df2 = df[unit_cols + ["cell", score_col]].dropna(subset=[score_col]).copy()
    if df2.empty:
        return 0
    df2["cluster"] = df2["cell"].map(cluster_labels)
    df2 = df2.dropna(subset=["cluster"])
    grp = (
        df2.groupby(unit_cols + ["cluster"])[score_col]
        .mean()
        .reset_index()
        .groupby(unit_cols)[score_col]
        .agg(lambda x: x.max() - x.min())
    )
    return int((grp > 0.1).sum())


# ---------------------------------------------------------------------------
# Core public function
# ---------------------------------------------------------------------------


def compare_length_strategies(
    h5ad_path: Path,
    pasbed_path: Path | None,
    gtf_path: Path | None,
    output_dir: Path,
    *,
    strategies_and_aggs: list[tuple[str, str]] | None = None,
    cluster_pairs: list[tuple[str, str]] | None = None,  # reserved
    threads: int = 4,
) -> pd.DataFrame:
    """Run length strategies + isoform aggregations and emit metrics.

    Loads ``h5ad_path`` once; builds the ``pas_isoform_map`` separately for
    each unique ``isoform_agg`` level encountered in ``strategies_and_aggs``
    so the GTF is only parsed once per aggregation type.

    Output files written to ``output_dir``:
        - ``<strategy>__<agg>/*.tsv`` — per-cell or per-unit raw scores.
        - ``summary.tsv`` — the returned DataFrame serialised to disk.

    Args:
        h5ad_path: Path to a clustered AnnData h5ad.
        pasbed_path: PAS BED file (for future isoform-level PAS coordinate
            mapping; currently read by the isoform branch of ``run_length``).
            Pass ``None`` to skip coordinate annotation.
        gtf_path: GTF file required for ``per_isoform`` aggregation.  If
            ``None`` or non-existent, any ``per_isoform`` entry in
            ``strategies_and_aggs`` is silently skipped.
        output_dir: Root directory for all outputs.
        strategies_and_aggs: List of ``(strategy_name, isoform_agg)`` tuples
            to benchmark.  Defaults to::

                [
                    ("classic",    "per_gene"),
                    ("proportion", "per_gene"),
                    ("proportion", "per_isoform"),
                    ("shannon",    "per_gene"),
                ]

        cluster_pairs: Reserved for future per-pair length scoring; currently
            unused.  Accepted to maintain interface symmetry with
            :func:`compare_diff_strategies`.
        threads: Maximum parallel workers forwarded to strategy workers
            via ``ResourceManager``.

    Returns:
        DataFrame indexed by ``"<strategy>__<agg>"`` with columns::

            n_quantified_units, mean_score, score_std, runtime_seconds,
            n_units_with_inter_cluster_shift

        ``n_quantified_units`` = number of unique (gene_id, transcript_id)
        pairs in the output.
    """
    if strategies_and_aggs is None:
        strategies_and_aggs = [
            ("classic",    "per_gene"),
            ("proportion", "per_gene"),
            ("proportion", "per_isoform"),
            ("shannon",    "per_gene"),
        ]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load h5ad once ---
    log.info("compare_length_strategies: loading %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    count_matrix = _to_count_matrix(adata)
    cluster_labels = pd.Series(
        adata.obs["leiden"].values, index=adata.obs_names, name="leiden"
    )
    log.info(
        "compare_length_strategies: %d cells × %d PAS, %d clusters",
        *adata.shape,
        cluster_labels.nunique(),
    )

    # Pre-build per_gene map once (from adata.var); per_isoform needs GTF.
    per_gene_map = _build_per_gene_map(adata)

    # Cache per_isoform map to avoid re-parsing GTF for multiple strategies.
    per_isoform_map: dict[int, list[tuple[str, str, int, int, int]]] | None = None

    from ema.quantification.strategies import get_pdui_strategy

    summary_rows: list[dict] = []

    for strategy_name, isoform_agg in strategies_and_aggs:
        run_key = f"{strategy_name}__{isoform_agg}"
        log.info("compare_length_strategies: running %s", run_key)

        # Resolve the isoform map for this aggregation level.
        if isoform_agg == "per_isoform":
            if per_isoform_map is None:
                if gtf_path is None or not Path(gtf_path).exists():
                    log.warning(
                        "compare_length_strategies: per_isoform requires --gtf "
                        "but gtf_path=%r is missing; skipping %s",
                        str(gtf_path), run_key,
                    )
                    continue
                # Parse isoform UTRs.  This is the slow step (~30 s for human GTF).
                try:
                    from ema.annotate.gtf2isoform_utr import parse_isoform_utrs
                    from ema.quantification.pas_to_isoform import map_pas_to_isoforms
                    _cache_dir = Path(h5ad_path).resolve().parent.parent / "gtf_cache"
                    _cache_dir.mkdir(parents=True, exist_ok=True)
                    isoform_utrs = parse_isoform_utrs(
                        Path(gtf_path),
                        cache_dir=_cache_dir,
                        n_workers=min(threads, 4),
                    )
                    if pasbed_path is not None and Path(pasbed_path).exists():
                        per_isoform_map = map_pas_to_isoforms(pasbed_path, isoform_utrs)
                    else:
                        log.warning(
                            "compare_length_strategies: per_isoform needs pasbed; "
                            "pasbed_path=%r missing; skipping %s",
                            str(pasbed_path), run_key,
                        )
                        continue
                except Exception as exc:
                    log.warning(
                        "compare_length_strategies: failed to build per_isoform map: %s; "
                        "skipping %s",
                        exc, run_key,
                    )
                    continue
            pas_map = per_isoform_map
        else:
            pas_map = per_gene_map

        if not pas_map:
            log.warning("compare_length_strategies: empty pas_map for %s; skipping", run_key)
            continue

        strat = get_pdui_strategy(strategy_name)
        t0 = time.perf_counter()
        try:
            df = strat.compute(
                count_matrix=count_matrix,
                pas_isoform_map=pas_map,
                aggregation=isoform_agg,
                isoform_collapse="none",
                pseudocount=0.0,
            )
        except Exception as exc:
            log.warning("compare_length_strategies: strategy %s failed: %s", run_key, exc)
            summary_rows.append({
                "run_key": run_key,
                "n_quantified_units": 0,
                "mean_score": float("nan"),
                "score_std": float("nan"),
                "runtime_seconds": round(time.perf_counter() - t0, 2),
                "n_units_with_inter_cluster_shift": 0,
            })
            continue

        elapsed = time.perf_counter() - t0

        # Write output TSV
        run_dir = output_dir / run_key
        run_dir.mkdir(parents=True, exist_ok=True)
        out_filename = getattr(strat, "output_filename", f"{run_key}.tsv")
        out_path = run_dir / out_filename
        df.to_csv(out_path, sep="\t", index=False)
        log.info(
            "compare_length_strategies: %s -> %d rows written to %s",
            run_key, len(df), out_path,
        )

        # Summary metrics
        score_col = getattr(strat, "score_column", "pdui")
        unit_cols = (
            ["gene_id", "transcript_id"]
            if "transcript_id" in df.columns
            else ["gene_id"]
        )
        n_units = (
            df[unit_cols].drop_duplicates().shape[0]
            if unit_cols and all(c in df.columns for c in unit_cols)
            else 0
        )

        scores = df[score_col].dropna() if score_col in df.columns else pd.Series(dtype=float)
        mean_score = float(scores.mean()) if len(scores) > 0 else float("nan")
        score_std = float(scores.std()) if len(scores) > 1 else float("nan")

        n_shift = _inter_cluster_shift(df, score_col, ["gene_id"], cluster_labels)

        summary_rows.append({
            "run_key": run_key,
            "n_quantified_units": n_units,
            "mean_score": round(mean_score, 4),
            "score_std": round(score_std, 4),
            "runtime_seconds": round(elapsed, 2),
            "n_units_with_inter_cluster_shift": n_shift,
        })

    summary = pd.DataFrame(summary_rows).set_index("run_key")
    summary_path = output_dir / "summary.tsv"
    summary.to_csv(summary_path, sep="\t")
    log.info("compare_length_strategies: summary written to %s\n%s", summary_path, summary.to_string())
    return summary
