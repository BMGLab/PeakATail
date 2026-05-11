"""Negative-Binomial GLM pairwise differential APA strategy.

For each PAS, fits a NB GLM::

    count ~ cluster_indicator + offset(log(library_size))

using :mod:`statsmodels`.  Dispersion (alpha) is estimated per-PAS via MLE;
if the solver fails to converge, a method-of-moments estimate is used as
fallback.  p-values come from the Wald z-test on the cluster coefficient.
BH FDR correction is applied across all tested PAS.

Parallelisation
---------------
Per-PAS GLM fits are embarrassingly parallel.  The implementation uses
``joblib.Parallel(backend="loky")`` over PAS batches.  Dense numpy slices
(not sparse matrices) are pre-extracted before dispatch so workers avoid
pickling large sparse objects.
"""

from __future__ import annotations

import math
import warnings
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import false_discovery_control
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from ema.switch_test.strategies import register_diff_strategy
from ema.switch_test.strategies.base import DiffAPAStrategy

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _mom_dispersion(counts: np.ndarray) -> float:
    """Method-of-moments dispersion (alpha) estimate for a 1-D count vector.

    For NB parameterisation ``Var = mu + alpha * mu^2``:

    .. math::
        \\hat{\\alpha} = \\frac{s^2 - \\bar{x}}{\\bar{x}^2}

    Clipped to ``[1e-4, 10]`` to keep the GLM numerically stable.
    """
    mu = counts.mean()
    if mu <= 0:
        return 1.0
    var = counts.var(ddof=1)
    alpha = (var - mu) / (mu ** 2)
    return float(np.clip(alpha, 1e-4, 10.0))


def _fit_single_pas(
    pas_id: Any,
    y: np.ndarray,           # shape (n_cells,) — raw counts for this PAS
    group_indicator: np.ndarray,  # shape (n_cells,) — 0/1 cluster indicator
    log_lib_size: np.ndarray,     # shape (n_cells,) — offset
    n_c1: int,
    n_c2: int,
    min_cells_per_group: int,
) -> dict | None:
    """Fit NB GLM for a single PAS and return result dict or None if skipped."""
    # Cell-level filter: require enough cells with >0 counts in each group
    mask1 = group_indicator == 0
    mask2 = group_indicator == 1
    nonzero1 = int((y[mask1] > 0).sum())
    nonzero2 = int((y[mask2] > 0).sum())
    if nonzero1 < min_cells_per_group or nonzero2 < min_cells_per_group:
        return None

    # Build design matrix [intercept, cluster_indicator]
    X = np.column_stack([np.ones(len(y), dtype=np.float64),
                         group_indicator.astype(np.float64)])

    # --- Step 1: estimate dispersion via NegativeBinomial MLE ---
    import statsmodels.api as sm
    from statsmodels.discrete.discrete_model import NegativeBinomial

    alpha = None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            nb_model = NegativeBinomial(y.astype(np.float64), X,
                                        loglike_method="nb2",
                                        offset=log_lib_size)
            nb_res = nb_model.fit(
                method="bfgs",
                maxiter=200,
                disp=False,
                skip_hessian=True,
            )
            alpha = float(nb_res.params[-1])  # last param is ln(alpha) or alpha
            # statsmodels NegativeBinomial reports alpha directly
            alpha = max(1e-4, min(alpha, 10.0))
    except Exception:
        alpha = None

    if alpha is None or not np.isfinite(alpha):
        alpha = _mom_dispersion(y)

    # --- Step 2: GLM with fixed dispersion ---
    from statsmodels.genmod.families import NegativeBinomial as NBFamily

    fallback_used = False
    coef = np.nan
    pvalue = np.nan
    test_stat = np.nan

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            glm = sm.GLM(
                y.astype(np.float64),
                X,
                family=NBFamily(alpha=alpha),
                offset=log_lib_size,
            )
            glm_res = glm.fit(method="irls", maxiter=200, disp=False)

        if not glm_res.converged:
            fallback_used = True

        coef = float(glm_res.params[1])          # cluster indicator coefficient
        pvalue = float(glm_res.pvalues[1])
        test_stat = float(glm_res.tvalues[1])
    except Exception:
        fallback_used = True

    if fallback_used or not np.isfinite(pvalue):
        # Fallback: Wald test with method-of-moments dispersion
        alpha = _mom_dispersion(y)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", ConvergenceWarning)
                glm = sm.GLM(
                    y.astype(np.float64),
                    X,
                    family=NBFamily(alpha=alpha),
                    offset=log_lib_size,
                )
                glm_res = glm.fit(method="irls", maxiter=100, disp=False)
            coef = float(glm_res.params[1])
            pvalue = float(glm_res.pvalues[1])
            test_stat = float(glm_res.tvalues[1])
            fallback_used = True
        except Exception:
            # Ultimate fallback: return p=1.0 so it is never significant
            return {
                "pas_id": pas_id,
                "pvalue": 1.0,
                "log2fc": 0.0,
                "dispersion": alpha,
                "n_cells": n_c1 + n_c2,
                "test_stat": 0.0,
                "fallback_used": True,
            }

    log2fc = coef / math.log(2) if np.isfinite(coef) else 0.0
    pvalue = float(np.clip(pvalue, 0.0, 1.0))

    return {
        "pas_id": pas_id,
        "pvalue": pvalue,
        "log2fc": log2fc,
        "dispersion": alpha,
        "n_cells": n_c1 + n_c2,
        "test_stat": test_stat,
        "fallback_used": fallback_used,
    }


def _fit_pas_batch(
    batch: list[tuple[Any, np.ndarray]],
    group_indicator: np.ndarray,
    log_lib_size: np.ndarray,
    n_c1: int,
    n_c2: int,
    min_cells_per_group: int,
) -> list[dict]:
    """Fit a batch of (pas_id, count_vector) pairs; returns non-None results."""
    out = []
    for pas_id, y in batch:
        result = _fit_single_pas(
            pas_id, y, group_indicator, log_lib_size, n_c1, n_c2,
            min_cells_per_group,
        )
        if result is not None:
            out.append(result)
    return out


# ---------------------------------------------------------------------------
# Strategy class
# ---------------------------------------------------------------------------

@register_diff_strategy
class NbPairwiseStrategy(DiffAPAStrategy):
    """Per-PAS NB GLM pairwise differential APA test.

    Fits ``count ~ cluster_indicator + offset(log(library_size))`` with
    family=NegativeBinomial per PAS.  Dispersion is estimated per-PAS via
    MLE (NegativeBinomial.fit), falling back to method-of-moments if the
    solver fails.  p-values are Wald z-statistics; BH FDR is applied across
    all tested PAS.

    Parallelised with ``joblib.Parallel(backend="loky")`` over PAS batches.
    Dense numpy arrays are sliced before dispatch so workers never need to
    pickle sparse matrices.
    """

    name: str = "nb_pairwise"
    supports_multi_condition: bool = False

    def test(
        self,
        count_matrix: pd.DataFrame,
        cluster_labels: pd.Series,
        cluster1: str | None = None,
        cluster2: str | None = None,
        min_cells_per_group: int = 10,
        n_jobs: int = -1,
    ) -> pd.DataFrame:
        """Run NB pairwise test.

        Parameters
        ----------
        count_matrix:
            Shape ``(n_cells, n_pas)``; integer UMI counts.
        cluster_labels:
            Per-cell cluster labels aligned to ``count_matrix`` rows.
        cluster1, cluster2:
            The two cluster labels to compare.  Both required.
        min_cells_per_group:
            Minimum cells with non-zero counts in each group per PAS.
        n_jobs:
            Parallel workers (joblib convention; -1 = all CPUs).

        Returns
        -------
        pd.DataFrame
            Indexed by ``pas_id`` with columns
            ``[pvalue, qvalue, log2fc, dispersion, n_cells, test_stat]``.
        """
        if cluster1 is None or cluster2 is None:
            raise ValueError(
                "NbPairwiseStrategy is pairwise: cluster1 and cluster2 must "
                "both be provided."
            )

        # --- restrict to the two clusters ---
        mask = cluster_labels.isin([cluster1, cluster2])
        labels_sub = cluster_labels[mask]
        mat_sub = count_matrix.loc[mask]  # (n_cells_sub, n_pas)

        mask1 = labels_sub == cluster1
        mask2 = labels_sub == cluster2
        n_c1 = int(mask1.sum())
        n_c2 = int(mask2.sum())

        if n_c1 < min_cells_per_group or n_c2 < min_cells_per_group:
            return pd.DataFrame(
                columns=["pvalue", "qvalue", "log2fc", "dispersion",
                         "n_cells", "test_stat"]
            )

        # --- pre-compute shared arrays (avoid re-pickling inside workers) ---
        # group_indicator: 0 = cluster1, 1 = cluster2
        group_indicator = np.where(labels_sub == cluster2, 1, 0)  # shape (n_cells_sub,)

        # library size = total counts per cell across all PAS
        lib_sizes = mat_sub.values.sum(axis=1).astype(np.float64)
        lib_sizes = np.where(lib_sizes < 1, 1.0, lib_sizes)  # avoid log(0)
        log_lib_size = np.log(lib_sizes)

        # Convert to dense numpy array once (avoids repeated sparse -> dense in workers)
        dense = mat_sub.values  # shape (n_cells_sub, n_pas)
        pas_ids = mat_sub.columns.tolist()
        n_pas = len(pas_ids)

        # --- build batches ---
        # When n_jobs=-1, delegate to ResourceManager so the user's --threads
        # ceiling is respected.  Explicit n_jobs=N is honoured as-is.
        actual_jobs: int
        if n_jobs == -1:
            from ema.utils import get_resource_manager
            actual_jobs = get_resource_manager().get_n_jobs(
                per_worker_mb=300, stage="nb_pairwise"
            )
        else:
            actual_jobs = max(1, n_jobs)

        batch_size = max(1, math.ceil(n_pas / actual_jobs))
        batches: list[list[tuple]] = []
        for start in range(0, n_pas, batch_size):
            end = min(start + batch_size, n_pas)
            batches.append(
                [(pas_ids[i], dense[:, i]) for i in range(start, end)]
            )

        # --- parallel execution ---
        raw_results = Parallel(n_jobs=actual_jobs, backend="loky")(
            delayed(_fit_pas_batch)(
                batch, group_indicator, log_lib_size,
                n_c1, n_c2, min_cells_per_group,
            )
            for batch in batches
        )

        # Flatten list-of-lists
        flat: list[dict] = [item for sublist in raw_results for item in sublist]

        if not flat:
            return pd.DataFrame(
                columns=["pvalue", "qvalue", "log2fc", "dispersion",
                         "n_cells", "test_stat"]
            )

        df = pd.DataFrame(flat).set_index("pas_id")

        # BH FDR correction across all tested PAS
        qvalues = false_discovery_control(df["pvalue"].values, method="bh")
        df["qvalue"] = qvalues

        # Drop internal bookkeeping column
        df = df.drop(columns=["fallback_used"], errors="ignore")

        return df[["pvalue", "qvalue", "log2fc", "dispersion",
                   "n_cells", "test_stat"]].sort_values("qvalue")
