"""Negative-Binomial GLM pairwise differential APA strategy.

For each PAS, fits a NB GLM::

    count ~ cluster_indicator + offset(log(library_size))

using :mod:`statsmodels`.  Dispersion (alpha) is estimated per-PAS via MLE;
if the solver fails to converge, a method-of-moments estimate is used as
fallback.  p-values come from the Wald z-test on the cluster coefficient.
BH FDR correction is applied across all tested PAS, except that PAS whose
dispersion collapsed to the numerical floor get no q-value (see
``ALPHA_FLOOR`` and issue #94).

Parallelisation
---------------
Per-PAS GLM fits are embarrassingly parallel.  The implementation uses
``joblib.Parallel(backend="loky")`` over PAS batches.  Dense numpy slices
(not sparse matrices) are pre-extracted before dispatch so workers avoid
pickling large sparse objects.
"""

from __future__ import annotations

import logging
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

log = logging.getLogger(__name__)

# Numerical clip on the per-PAS NB dispersion alpha.  Hitting the LOWER bound
# collapses the GLM to Poisson and the Wald standard error with it, which makes
# such a test anti-conservative (issue #94: 8.5% of null tests hit the floor but
# produced 67% of the false q<0.05 calls).  Rows that hit it are flagged and
# their q-value is withheld -- see ``NbPairwiseStrategy.test``.
ALPHA_FLOOR = 1e-4
ALPHA_CEIL = 10.0

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
    return float(np.clip(alpha, ALPHA_FLOOR, ALPHA_CEIL))


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
            alpha = max(ALPHA_FLOOR, min(alpha, ALPHA_CEIL))
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
                "dispersion_floored": bool(alpha <= ALPHA_FLOOR),
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
        "dispersion_floored": bool(alpha <= ALPHA_FLOOR),
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

    Tunable hyperparameters:
        min_cells_per_group (default 10): Minimum cells with *non-zero* counts
            in each cluster for a PAS to be tested.  Raise to reduce noisy
            low-count tests; lower only if your clusters are very small.
            CLI: ``--min-cells-per-group`` / YAML: ``min_cells_per_group``.
        fdr (default 0.05): Benjamini-Hochberg FDR threshold applied after
            testing.  CLI: ``--fdr`` / YAML: ``fdr``.

    Rows whose dispersion collapsed to the lower clip (``ALPHA_FLOOR``) are
    marked ``dispersion_floored=True`` and get ``qvalue=NaN``: a floored
    dispersion makes the fit Poisson and the Wald p-value anti-conservative,
    so such a test is never reported as significant.  Even for the remaining
    rows the q-values are not permutation-calibrated (issue #94).

    Deliberately left hardcoded (internal numerics, not researcher-facing):
        alpha dispersion clip: [1e-4, 10.0] — keeps GLM numerically stable.
        maxiter (MLE): 200 iterations; fallback to 100 for the MoM path.
        per_worker_mb: 300 MB — used by ResourceManager to cap parallelism.
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
        **_ignored,  # pas_gene_map is fisher-only — accept & drop
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
            ``[pvalue, qvalue, log2fc, dispersion, dispersion_floored,
            n_cells, test_stat]``.  ``qvalue`` is NaN wherever
            ``dispersion_floored`` is True (issue #94).
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
                         "dispersion_floored", "n_cells", "test_stat"]
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
                         "dispersion_floored", "n_cells", "test_stat"]
            )

        df = pd.DataFrame(flat).set_index("pas_id")

        # BH FDR correction across all tested PAS.  Dispersion-floored rows
        # stay in the BH input so the family size m is unchanged (dropping them
        # would shrink m and make every OTHER q-value less conservative); only
        # their reported q is withheld below.
        qvalues = false_discovery_control(df["pvalue"].values, method="bh")
        df["qvalue"] = qvalues

        # --- withhold q for dispersion-floored tests (issue #94) ---
        # When the per-PAS dispersion collapses to ALPHA_FLOOR the GLM is a
        # Poisson fit and its Wald SE is a lower bound, so the p-value is
        # anti-conservative.  The dispersion is also estimated under the FULL
        # model, so noise that mimics a group difference is absorbed into the
        # fit and pushes alpha down -- exactly the tests that then look most
        # significant.  Under a label-permutation null those rows produced 67%
        # of the false q<0.05 calls, so their q is reported as NaN (never
        # significant) while the raw p-value is kept for inspection.
        floored = df["dispersion_floored"].astype(bool)
        n_floored = int(floored.sum())
        if n_floored:
            df.loc[floored, "qvalue"] = np.nan
            log.warning(
                "NbPairwiseStrategy: %d/%d tested PAS (%.1f%%) hit the "
                "dispersion floor (alpha=%g); their q-values are withheld "
                "(NaN) because a floored dispersion collapses the GLM to "
                "Poisson and makes the Wald p-value anti-conservative.  The "
                "raw p-values are still reported in the 'pvalue' column and "
                "the rows are marked in 'dispersion_floored'.  nb_pairwise "
                "q-values are NOT permutation-calibrated -- see issue #94.",
                n_floored, len(df), 100.0 * n_floored / len(df), ALPHA_FLOOR,
            )

        # Drop internal bookkeeping column
        df = df.drop(columns=["fallback_used"], errors="ignore")

        return df[["pvalue", "qvalue", "log2fc", "dispersion",
                   "dispersion_floored", "n_cells", "test_stat"]].sort_values(
            "qvalue")
