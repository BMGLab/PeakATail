"""Compare multiple differential APA strategies on a single h5ad.

Production module for benchmarking :mod:`ema.switch_test` strategies
(``fisher``, ``nb_pairwise``, ``nb_multi``) against a clustered AnnData
and emitting per-strategy summary metrics plus concordance statistics.

Typical usage::

    from pathlib import Path
    from ema.benchmark.switch_diff_compare import compare_diff_strategies

    summary = compare_diff_strategies(
        h5ad_path=Path("clusters.h5ad"),
        pasbed_path=Path("pasbed.bed"),
        gtf_path=None,
        output_dir=Path("reports/diff_compare"),
        strategies=["fisher", "nb_pairwise", "nb_multi"],
        cluster_pairs=[("0", "1")],
        fdr=0.05,
        threads=4,
    )
    print(summary)

Notes
-----
- The h5ad is loaded **once** and reused across all strategies to minimise
  I/O and memory overhead on large datasets.
- ``nb_multi`` is an omnibus test (one stat per PAS across all clusters)
  and therefore lacks ``log2fc`` / ``delta_proportion``; its sig-count
  definition falls back to ``qvalue < fdr`` without an effect-size filter.
- Concordance Jaccard is computed only across per-pair strategies (``fisher``
  and ``nb_pairwise``) that share the same ``(c1, c2)`` namespace.  The
  ``nb_multi`` omnibus result is excluded from pairwise Jaccard but still
  appears in the summary row.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _to_dense(adata: ad.AnnData) -> pd.DataFrame:
    """Return a cells × PAS DataFrame of integer counts.

    Args:
        adata: Clustered AnnData with sparse or dense ``X``.

    Returns:
        DataFrame with shape ``(n_cells, n_pas)``; columns are PAS IDs
        converted to ``int`` where possible, falling back to the raw string
        var-name otherwise.
    """
    X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    try:
        pas_index = [int(v) for v in adata.var_names]
    except (TypeError, ValueError):
        pas_index = list(adata.var_names)
    return pd.DataFrame(X, index=list(adata.obs_names), columns=pas_index)


def _build_pas_gene_map(adata: ad.AnnData) -> dict[str, str] | None:
    """Build a ``{pas_id_str: gene_id_str}`` mapping from ``adata.var``.

    Args:
        adata: AnnData whose ``var`` DataFrame should contain a ``gene_id``
            column.

    Returns:
        Mapping dict, or ``None`` if ``gene_id`` is absent from ``adata.var``.
    """
    if "gene_id" not in adata.var.columns:
        log.warning("adata.var has no 'gene_id' column; Fisher will use global framing")
        return None
    return {
        str(k): str(v)
        for k, v in adata.var["gene_id"].dropna().items()
        if str(v).strip()
    }


def _sig_set(
    df: pd.DataFrame,
    fdr: float,
    delta_col: str | None = "delta_proportion",
    delta_thresh: float = 0.2,
) -> set:
    """Return the set of significant PAS IDs from a per-pair result DataFrame.

    A PAS is significant when ``qvalue < fdr``.  An *strongly* significant set
    additionally requires ``|delta_col| > delta_thresh``; however, for the
    Jaccard concordance the simple ``qvalue < fdr`` set is used so that all
    strategies are comparable (``nb_pairwise`` has ``log2fc`` but no
    ``delta_proportion``).

    Args:
        df: Differential result DataFrame indexed by pas_id with a ``qvalue``
            column (and optionally a ``delta_proportion`` / ``log2fc`` column).
        fdr: FDR threshold.
        delta_col: Optional effect-size column for the "strong" filter.
        delta_thresh: Absolute threshold on ``delta_col``.

    Returns:
        Frozenset of PAS IDs (index values) that pass the basic ``qvalue < fdr``
        filter.
    """
    if "qvalue" not in df.columns:
        return set()
    return set(df.index[df["qvalue"] < fdr].tolist())


def _jaccard(a: set, b: set) -> float:
    """Compute Jaccard similarity between two sets.

    Args:
        a: First set.
        b: Second set.

    Returns:
        ``|a ∩ b| / |a ∪ b|``, or ``0.0`` when both sets are empty.
    """
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


# ---------------------------------------------------------------------------
# Core public function
# ---------------------------------------------------------------------------


def compare_diff_strategies(
    h5ad_path: Path,
    pasbed_path: Path | None,
    gtf_path: Path | None,
    output_dir: Path,
    *,
    strategies: list[str] | None = None,
    cluster_pairs: list[tuple[str, str]] | None = None,
    fdr: float = 0.05,
    threads: int = 4,
) -> pd.DataFrame:
    """Run multiple differential APA strategies and emit per-strategy + concordance metrics.

    Loads ``h5ad_path`` once and passes the same in-memory AnnData to each
    strategy runner, avoiding repeated I/O overhead on large datasets.

    Output files written to ``output_dir``:
        - ``<strategy>/differential/*.tsv`` — per-pair raw stats.
        - ``summary.tsv`` — the returned DataFrame serialised to disk.
        - ``concordance_jaccard.tsv`` — pairwise Jaccard on significant PAS
          sets (rows × cols = pairwise strategies only; nb_multi is excluded).
        - ``concordance_overlap.tsv`` — raw count of overlapping significant PAS.

    Args:
        h5ad_path: Path to a clustered AnnData h5ad produced by
            ``ema run`` / ``ema cluster``.
        pasbed_path: Optional PAS BED file for coordinate annotation.  When
            supplied, genomic coordinates are joined onto output TSVs.
        gtf_path: Optional GTF path (reserved for future isoform-aware
            annotation; currently unused by diff strategies).
        output_dir: Directory under which all output subdirectories and
            summary files are written.  Created if absent.
        strategies: Strategy names to benchmark.  Defaults to
            ``["fisher", "nb_pairwise", "nb_multi"]``.
        cluster_pairs: Cluster label pairs to test (pairwise strategies).
            Defaults to ``[("0", "1")]``.
        fdr: Benjamini-Hochberg FDR threshold for significance calls.
        threads: Maximum parallel worker threads forwarded to each strategy's
            internal ``ResourceManager``.

    Returns:
        DataFrame indexed by strategy name with columns::

            n_tests, n_sig, n_sig_strong, median_qvalue,
            mean_log2fc_abs, mean_delta_proportion_abs, runtime_seconds

        ``n_sig_strong`` applies ``|delta_proportion| > 0.2`` (or
        ``|log2fc| > 1`` for ``nb_pairwise``) on top of the FDR filter.
        Columns that are not meaningful for a strategy (e.g. ``log2fc`` for
        ``nb_multi``) are filled with ``NaN``.
    """
    strategies = strategies or ["fisher", "nb_pairwise", "nb_multi"]
    cluster_pairs = cluster_pairs or [("0", "1")]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load h5ad once ---
    log.info("compare_diff_strategies: loading %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    diff_df = _to_dense(adata)
    cluster_labels = pd.Series(
        adata.obs["leiden"].values, index=adata.obs_names, name="leiden"
    )
    pas_gene_map = _build_pas_gene_map(adata)

    log.info(
        "compare_diff_strategies: %d cells × %d PAS, %d clusters",
        *adata.shape,
        cluster_labels.nunique(),
    )

    # Import strategy registry and runner
    from ema.switch_test.strategies import get_diff_strategy
    from ema.switch_test.pair_runner import run_one_pair

    summary_rows: list[dict] = []
    # Collect per-pair significant sets per pairwise strategy for concordance.
    # Structure: strategy_name -> pair_key -> set[pas_id]
    pairwise_sig_sets: dict[str, dict[tuple[str, str], set]] = {}

    for strategy_name in strategies:
        log.info("compare_diff_strategies: running strategy=%s", strategy_name)
        strat = get_diff_strategy(strategy_name)
        strat_dir = output_dir / strategy_name
        strat_dir.mkdir(parents=True, exist_ok=True)
        diff_dir = strat_dir / "differential"
        diff_dir.mkdir(exist_ok=True)

        t0 = time.perf_counter()

        if strat.supports_multi_condition:
            # Omnibus test: one call across all clusters.
            result_df = strat.test(
                count_matrix=diff_df,
                cluster_labels=cluster_labels,
                min_cells_per_group=10,
                n_jobs=threads,
            )
            elapsed = time.perf_counter() - t0
            out_path = diff_dir / f"{strategy_name}_omnibus.tsv"
            result_df.to_csv(out_path, sep="\t")
            log.info(
                "compare_diff_strategies: %s omnibus -> %d rows, written to %s",
                strategy_name, len(result_df), out_path,
            )

            # Summary metrics
            n_tests = int(len(result_df))
            sig_mask = result_df["qvalue"] < fdr if "qvalue" in result_df.columns else pd.Series([], dtype=bool)
            n_sig = int(sig_mask.sum()) if len(sig_mask) else 0
            # nb_multi has no log2fc or delta_proportion
            n_sig_strong = n_sig  # no additional effect-size filter for omnibus
            med_q = float(result_df["qvalue"].median()) if "qvalue" in result_df.columns and n_tests > 0 else float("nan")
            summary_rows.append({
                "strategy": strategy_name,
                "n_tests": n_tests,
                "n_sig": n_sig,
                "n_sig_strong": n_sig_strong,
                "median_qvalue": med_q,
                "mean_log2fc_abs": float("nan"),
                "mean_delta_proportion_abs": float("nan"),
                "runtime_seconds": round(elapsed, 2),
            })

        else:
            # Pairwise test across specified cluster pairs.
            pair_dfs: list[pd.DataFrame] = []
            pairwise_sig_sets[strategy_name] = {}

            for c1, c2 in cluster_pairs:
                _, _, df = run_one_pair(
                    strategy_name,
                    diff_df,
                    cluster_labels,
                    c1,
                    c2,
                    n_jobs_inner=threads,
                    min_cells_per_group=10,
                    pas_gene_map=pas_gene_map,
                )
                out_path = diff_dir / f"{strategy_name}_{c1}_vs_{c2}.tsv"
                df.to_csv(out_path, sep="\t")
                log.info(
                    "compare_diff_strategies: %s %s vs %s -> %d rows",
                    strategy_name, c1, c2, len(df),
                )
                pair_dfs.append(df)
                pairwise_sig_sets[strategy_name][(c1, c2)] = _sig_set(df, fdr)

            elapsed = time.perf_counter() - t0

            if pair_dfs:
                combined = pd.concat(pair_dfs, axis=0)
            else:
                combined = pd.DataFrame()

            n_tests = int(len(combined))
            sig_mask = (combined["qvalue"] < fdr) if "qvalue" in combined.columns else pd.Series(dtype=bool)
            n_sig = int(sig_mask.sum())

            # Strong: qvalue < fdr AND effect size threshold
            if "delta_proportion" in combined.columns:
                strong_mask = sig_mask & (combined["delta_proportion"].abs() > 0.2)
                mean_delta = float(combined.loc[sig_mask, "delta_proportion"].abs().mean()) if n_sig > 0 else float("nan")
            else:
                strong_mask = sig_mask
                mean_delta = float("nan")

            if "log2fc" in combined.columns and n_sig > 0:
                mean_l2fc = float(combined.loc[sig_mask, "log2fc"].abs().mean())
                strong_mask = sig_mask & (combined["log2fc"].abs() > 1.0)
            else:
                mean_l2fc = float("nan")

            n_sig_strong = int(strong_mask.sum())
            med_q = float(combined["qvalue"].median()) if "qvalue" in combined.columns and n_tests > 0 else float("nan")

            summary_rows.append({
                "strategy": strategy_name,
                "n_tests": n_tests,
                "n_sig": n_sig,
                "n_sig_strong": n_sig_strong,
                "median_qvalue": med_q,
                "mean_log2fc_abs": mean_l2fc,
                "mean_delta_proportion_abs": mean_delta,
                "runtime_seconds": round(elapsed, 2),
            })

    # --- Build summary DataFrame ---
    summary = pd.DataFrame(summary_rows).set_index("strategy")
    summary_path = output_dir / "summary.tsv"
    summary.to_csv(summary_path, sep="\t")
    log.info("compare_diff_strategies: summary written to %s", summary_path)

    # --- Concordance: Jaccard + overlap for pairwise strategies ---
    pairwise_strats = [s for s in strategies if s in pairwise_sig_sets]
    if len(pairwise_strats) >= 2:
        # Union significant sets across cluster pairs per strategy
        union_sig: dict[str, set] = {}
        for s in pairwise_strats:
            all_sig: set = set()
            for pair_set in pairwise_sig_sets[s].values():
                all_sig |= pair_set
            union_sig[s] = all_sig

        jaccard_mat = pd.DataFrame(index=pairwise_strats, columns=pairwise_strats, dtype=float)
        overlap_mat = pd.DataFrame(index=pairwise_strats, columns=pairwise_strats, dtype=int)
        for s1 in pairwise_strats:
            for s2 in pairwise_strats:
                jaccard_mat.loc[s1, s2] = _jaccard(union_sig[s1], union_sig[s2])
                overlap_mat.loc[s1, s2] = len(union_sig[s1] & union_sig[s2])

        jaccard_mat.to_csv(output_dir / "concordance_jaccard.tsv", sep="\t")
        overlap_mat.to_csv(output_dir / "concordance_overlap.tsv", sep="\t")
        log.info(
            "compare_diff_strategies: concordance Jaccard written to %s",
            output_dir / "concordance_jaccard.tsv",
        )
    else:
        log.info("compare_diff_strategies: only one pairwise strategy; skipping concordance")

    log.info("compare_diff_strategies: done.\n%s", summary.to_string())
    return summary
