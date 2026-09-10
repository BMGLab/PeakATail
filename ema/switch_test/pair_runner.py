"""Per-pair worker function for parallel differential APA testing.

This module provides :func:`run_one_pair`, a top-level function that is
pickle-safe and suitable for dispatch via ``multiprocessing.Pool``.

Design notes
------------
- Top-level (module-level) function: required for ``spawn``-based
  ``multiprocessing.Pool`` to pickle correctly.
- Self-contained: imports the strategy registry internally so each worker
  subprocess performs its own registration on import.
- Stateless: receives all inputs by value (DataFrames passed directly, not
  by path) so no shared memory is required.
- ``n_jobs_inner`` lets the caller (cli.py) pass the inner-parallelism
  budget that ResourceManager.split_jobs() computed, avoiding CPU
  over-subscription when many pairs run in parallel.
"""

from __future__ import annotations

import pandas as pd


def run_one_pair(
    strategy_name: str,
    diff_df: pd.DataFrame,
    cluster_labels: pd.Series,
    c1: str,
    c2: str,
    n_jobs_inner: int = 1,
    min_cells_per_group: int = 10,
    pas_gene_map: dict[str, str] | None = None,
    count_mode: str = "cells",
    full_count_matrix: pd.DataFrame | None = None,
) -> tuple[str, str, pd.DataFrame]:
    """Run a single cluster-pair differential APA test.

    This function is intentionally a top-level module function (not a method
    or closure) so that ``multiprocessing.Pool`` with ``spawn`` context can
    pickle it without error.

    Args:
        strategy_name: Registry key for the differential strategy, e.g.
            ``"fisher"`` or ``"nb_pairwise"``.
        diff_df: Count matrix, shape ``(n_cells, n_pas)``.
        cluster_labels: Series indexed by cell barcode with cluster labels.
        c1: Label of the first cluster.
        c2: Label of the second cluster.
        n_jobs_inner: Number of inner parallel workers for NB-style per-PAS
            fits.  Controlled by ``ResourceManager.split_jobs()`` to avoid
            over-subscription.
        min_cells_per_group: Minimum cells per group to include a PAS.
        full_count_matrix: Unrestricted count matrix used ONLY for the
            within-gene denominator when ``diff_df`` has been column-restricted
            by ``--marker-top-n`` pre-selection (issue #94).  ``None`` when no
            restriction is active.

    Returns:
        Tuple of ``(c1, c2, result_df)`` so callers can reconstruct the
        keyed result dict after ``pool.imap_unordered`` reorders outputs.
    """
    # Import inside the worker so each subprocess registers strategies fresh.
    from ema.switch_test.strategies import get_diff_strategy

    strategy = get_diff_strategy(strategy_name)
    # pas_gene_map is keyword-only on the strategy interface; pass it through
    # only when the caller supplied one so strategies that ignore it (NB)
    # don't see an unexpected None in their **kwargs path.
    extra: dict = {"pas_gene_map": pas_gene_map} if pas_gene_map is not None else {}
    # Issue #94: only fisher consumes full_count_matrix (the NB strategies
    # accept-and-drop via **_ignored), and it is only ever set when the caller
    # restricted the tested PAS -- so forward it only when it exists.
    if full_count_matrix is not None:
        extra["full_count_matrix"] = full_count_matrix
    # count_mode is a fisher-specific knob (reads vs de-pseudoreplicated
    # cells); the NB strategies accept-and-drop it via **_ignored.
    result_df = strategy.test(
        count_matrix=diff_df,
        cluster_labels=cluster_labels,
        cluster1=c1,
        cluster2=c2,
        min_cells_per_group=min_cells_per_group,
        n_jobs=n_jobs_inner,
        count_mode=count_mode,
        **extra,
    )
    return (c1, c2, result_df)
