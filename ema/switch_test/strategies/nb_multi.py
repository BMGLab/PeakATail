"""Negative-Binomial GLM multi-condition differential APA strategy.

For each PAS, fits two NB GLMs:

* **Full model**: ``count ~ C(cluster) + offset(log(library_size))``
* **Null model**: ``count ~ 1 + offset(log(library_size))``

The likelihood ratio test (LRT) statistic is::

    LRT = 2 * (LL_full - LL_null)

which is chi-squared distributed with ``df = K - 1`` degrees of freedom,
where ``K`` is the number of distinct clusters.  This gives an omnibus
p-value: significant PAS are those where *any* cluster shows differential
usage, regardless of direction.

No ``log2fc`` column is produced because the test is multi-level; the caller
must inspect cluster-level coefficients from the full model if contrasts
are needed.

Parallelisation
---------------
Same joblib / loky pattern as :mod:`nb_pairwise`:  dense numpy slices are
pre-extracted and chunked into batches before dispatch.
"""

from __future__ import annotations

import math
import warnings
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import chi2, false_discovery_control
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from ema.switch_test.strategies import register_diff_strategy
from ema.switch_test.strategies.base import DiffAPAStrategy

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _mom_dispersion(counts: np.ndarray) -> float:
    """Method-of-moments dispersion estimate — identical to nb_pairwise."""
    mu = counts.mean()
    if mu <= 0:
        return 1.0
    var = counts.var(ddof=1)
    alpha = (var - mu) / (mu ** 2)
    return float(np.clip(alpha, 1e-4, 10.0))


def _lrt_single_pas(
    pas_id: Any,
    y: np.ndarray,           # (n_cells,) raw counts
    X_full: np.ndarray,      # (n_cells, K) one-hot / dummies for clusters
    X_null: np.ndarray,      # (n_cells, 1) intercept only
    log_lib_size: np.ndarray,
    df_lrt: int,
    n_cells: int,
    min_cells_per_group: int,
    cluster_indicator: np.ndarray,  # (n_cells,) integer cluster index 0..K-1
    n_clusters: int,
) -> dict | None:
    """Fit full + null NB GLM, compute LRT p-value for one PAS."""
    # PAS-level filter (scRNA convention): require at least min_cells_per_group
    # cells with nonzero counts TOTAL, AND signal in at least 2 clusters.
    total_nonzero = int((y > 0).sum())
    if total_nonzero < min_cells_per_group:
        return None
    clusters_with_signal = sum(
        1 for k in range(n_clusters)
        if (y[cluster_indicator == k] > 0).any()
    )
    if clusters_with_signal < 2:
        return None

    import statsmodels.api as sm
    from statsmodels.discrete.discrete_model import NegativeBinomial
    from statsmodels.genmod.families import NegativeBinomial as NBFamily

    # --- dispersion estimate (on full model) ---
    alpha: float
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            nb_disp = NegativeBinomial(y.astype(np.float64), X_full,
                                       loglike_method="nb2",
                                       offset=log_lib_size)
            nb_res = nb_disp.fit(
                method="bfgs", maxiter=200, disp=False, skip_hessian=True
            )
            alpha_raw = float(nb_res.params[-1])
            alpha = float(np.clip(alpha_raw, 1e-4, 10.0))
    except Exception:
        alpha = _mom_dispersion(y)

    if not np.isfinite(alpha):
        alpha = _mom_dispersion(y)

    family = NBFamily(alpha=alpha)

    # --- fit full model ---
    ll_full: float
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            glm_full = sm.GLM(
                y.astype(np.float64), X_full,
                family=family, offset=log_lib_size,
            )
            res_full = glm_full.fit(method="irls", maxiter=200, disp=False)
        ll_full = float(res_full.llf)
    except Exception:
        return {
            "pas_id": pas_id,
            "pvalue": 1.0,
            "test_stat": 0.0,
            "df": df_lrt,
            "dispersion": alpha,
            "n_cells": n_cells,
        }

    # --- fit null model ---
    ll_null: float
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            glm_null = sm.GLM(
                y.astype(np.float64), X_null,
                family=family, offset=log_lib_size,
            )
            res_null = glm_null.fit(method="irls", maxiter=200, disp=False)
        ll_null = float(res_null.llf)
    except Exception:
        return {
            "pas_id": pas_id,
            "pvalue": 1.0,
            "test_stat": 0.0,
            "df": df_lrt,
            "dispersion": alpha,
            "n_cells": n_cells,
        }

    # --- LRT ---
    lrt_stat = max(0.0, 2.0 * (ll_full - ll_null))
    pvalue = float(chi2.sf(lrt_stat, df=df_lrt))
    pvalue = float(np.clip(pvalue, 0.0, 1.0))

    return {
        "pas_id": pas_id,
        "pvalue": pvalue,
        "test_stat": lrt_stat,
        "df": df_lrt,
        "dispersion": alpha,
        "n_cells": n_cells,
    }


def _lrt_pas_batch(
    batch: list[tuple[Any, np.ndarray]],
    X_full: np.ndarray,
    X_null: np.ndarray,
    log_lib_size: np.ndarray,
    df_lrt: int,
    n_cells: int,
    min_cells_per_group: int,
    cluster_indicator: np.ndarray,
    n_clusters: int,
) -> list[dict]:
    """Fit LRT for a batch of (pas_id, count_vector) pairs."""
    out = []
    for pas_id, y in batch:
        result = _lrt_single_pas(
            pas_id, y, X_full, X_null, log_lib_size, df_lrt,
            n_cells, min_cells_per_group, cluster_indicator, n_clusters,
        )
        if result is not None:
            out.append(result)
    return out


# ---------------------------------------------------------------------------
# Strategy class
# ---------------------------------------------------------------------------


@register_diff_strategy
class NbMultiStrategy(DiffAPAStrategy):
    """Per-PAS NB GLM omnibus test across all clusters.

    Fits ``count ~ C(cluster) + offset(log(library_size))`` (full model)
    and ``count ~ 1 + offset(log(library_size))`` (null model) per PAS,
    then performs a likelihood ratio test.  The LRT statistic follows
    chi-squared with ``K-1`` degrees of freedom where ``K`` is the number
    of distinct clusters.

    Returns one row per PAS with the omnibus p-value; there is no
    ``log2fc`` column because the test is multi-level.

    Parallelised with ``joblib.Parallel(backend="loky")`` over PAS batches.

    Tunable hyperparameters:
        min_cells_per_group (default 10): Minimum cells with *non-zero* counts
            across all clusters for a PAS to pass the omnibus filter.  Also
            requires signal (non-zero counts) in at least 2 clusters.
            CLI: ``--min-cells-per-group`` / YAML: ``min_cells_per_group``.
        fdr (default 0.05): Benjamini-Hochberg FDR threshold applied after
            testing.  CLI: ``--fdr`` / YAML: ``fdr``.

    Deliberately left hardcoded (internal numerics, not researcher-facing):
        alpha dispersion clip: [1e-4, 10.0] — keeps GLM numerically stable.
        maxiter: 200 iterations for both full and null model fitting.
        per_worker_mb: 300 MB — used by ResourceManager to cap parallelism.
    """

    name: str = "nb_multi"
    supports_multi_condition: bool = True

    def test(
        self,
        count_matrix: pd.DataFrame,
        cluster_labels: pd.Series,
        cluster1: str | None = None,   # ignored for multi-condition
        cluster2: str | None = None,   # ignored for multi-condition
        min_cells_per_group: int = 10,
        n_jobs: int = -1,
        sample_split: bool = False,
        split_seed: int = 0,
        full_count_matrix: pd.DataFrame | None = None,
        **_ignored,  # pas_gene_map is fisher-only
    ) -> pd.DataFrame:
        """Run NB omnibus test across all cluster levels.

        Parameters
        ----------
        count_matrix:
            Shape ``(n_cells, n_pas)``; integer UMI counts.
        cluster_labels:
            Per-cell cluster labels.  All unique values define the groups.
        cluster1, cluster2:
            Ignored; accepted for interface compatibility.
        min_cells_per_group:
            Minimum cells with non-zero counts per cluster per PAS.
        n_jobs:
            Parallel workers.

        Returns
        -------
        pd.DataFrame
            Indexed by ``pas_id`` with columns
            ``[pvalue, qvalue, test_stat, df, dispersion, n_cells]``.
        """
        # Align
        common_idx = count_matrix.index.intersection(cluster_labels.index)
        mat = count_matrix.loc[common_idx]
        labs = cluster_labels.loc[common_idx]

        # D5: the omnibus filter (which PAS are "testable") is decided on the
        # SAME cells the LRT then tests — a selection/inference double-dip that
        # inflates significance toward ~100%. ``sample_split=True`` removes it:
        # partition cells 50/50 (deterministic, seeded), SELECT testable PAS on
        # half A, and run the LRT (INFERENCE) on the disjoint half B. Default
        # False (unchanged behaviour); enabling it changes reported numbers and
        # halves power — characterise on the no-atlas re-run (flagged). Note this
        # does NOT remove the upstream circularity of clusters being defined on
        # the same counts; that must be split at clustering time.
        if sample_split:
            n = len(common_idx)
            if n < 4:
                raise ValueError("sample_split needs >=4 cells")
            rng = np.random.RandomState(split_seed)
            perm = rng.permutation(n)
            half = n // 2
            idx_A = common_idx[perm[:half]]
            idx_B = common_idx[perm[half:]]
            mat_A, labs_A = mat.loc[idx_A], labs.loc[idx_A]
            # Selection on A: PAS with >= min_cells_per_group nonzero cells AND
            # nonzero signal in >= 2 clusters (mirrors the per-PAS LRT filter).
            selected: list = []
            for p in mat.columns:
                nz = mat_A[p].values > 0
                if int(nz.sum()) < min_cells_per_group:
                    continue
                if labs_A[nz].nunique() < 2:
                    continue
                selected.append(p)
            if not selected:
                return pd.DataFrame(
                    columns=["pvalue", "qvalue", "test_stat", "df",
                             "dispersion", "n_cells"]
                )
            # Inference on the disjoint half B, restricted to selected PAS.
            return self.test(
                mat.loc[idx_B, selected], labs.loc[idx_B],
                cluster1=cluster1, cluster2=cluster2,
                min_cells_per_group=min_cells_per_group, n_jobs=n_jobs,
                sample_split=False,
            )

        cluster_cats = sorted(labs.unique().tolist())
        K = len(cluster_cats)
        if K < 2:
            raise ValueError(
                "nb_multi requires at least 2 cluster levels; "
                f"found {K}: {cluster_cats}"
            )

        # Integer cluster indicator (0-indexed)
        cat_to_idx = {c: i for i, c in enumerate(cluster_cats)}
        cluster_indicator = np.array([cat_to_idx[c] for c in labs], dtype=np.int32)

        # Design matrices
        # Full: one-hot dummies without intercept (statsmodels needs add_constant separately)
        import statsmodels.api as sm

        # Use pd.get_dummies to build the full design matrix
        dummies = pd.get_dummies(labs, prefix="cl", drop_first=False).astype(np.float64)
        # Drop first column to avoid collinearity; add intercept manually
        dummies_reduced = dummies.iloc[:, 1:]  # K-1 columns
        X_full = sm.add_constant(dummies_reduced.values, prepend=True)  # (n_cells, K)
        X_null = sm.add_constant(
            np.zeros((len(mat), 1), dtype=np.float64), prepend=False
        )  # (n_cells, 1) — just the intercept
        X_null = np.ones((len(mat), 1), dtype=np.float64)

        df_lrt = K - 1
        n_cells = len(mat)

        # Library size must be the cell's SEQUENCING DEPTH, not the depth of
        # whichever PAS survived a pre-selection: an offset that depends on
        # which hypotheses you chose to test is not a depth proxy, and it made
        # p-values shift by up to two orders of magnitude between an
        # unrestricted run and a --marker-top-n / --prefilter-min-cells run
        # (issue #94). ``full_count_matrix`` carries the unrestricted matrix
        # when the caller narrowed the columns; the grouped UTR paths do NOT
        # pass it, because there the per-group scoping IS the analysis unit.
        _depth_src = mat if full_count_matrix is None else \
            full_count_matrix.reindex(index=mat.index).fillna(0.0)
        lib_sizes = _depth_src.values.sum(axis=1).astype(np.float64)
        lib_sizes = np.where(lib_sizes < 1, 1.0, lib_sizes)
        log_lib_size = np.log(lib_sizes)

        # Dense array and pas list
        dense = mat.values  # (n_cells, n_pas)
        pas_ids = mat.columns.tolist()
        n_pas = len(pas_ids)

        # Batching
        # When n_jobs=-1, delegate to ResourceManager so the user's --threads
        # ceiling is respected.  Explicit n_jobs=N is honoured as-is.
        if n_jobs == -1:
            from ema.utils import get_resource_manager
            actual_jobs = get_resource_manager().get_n_jobs(
                per_worker_mb=300, stage="nb_multi"
            )
        else:
            actual_jobs = max(1, n_jobs)
        batch_size = max(1, math.ceil(n_pas / actual_jobs))
        batches: list[list[tuple]] = []
        for start in range(0, n_pas, batch_size):
            end = min(start + batch_size, n_pas)
            batches.append([(pas_ids[i], dense[:, i]) for i in range(start, end)])

        # Parallel execution
        raw_results = Parallel(n_jobs=actual_jobs, backend="loky")(
            delayed(_lrt_pas_batch)(
                batch, X_full, X_null, log_lib_size, df_lrt,
                n_cells, min_cells_per_group, cluster_indicator, K,
            )
            for batch in batches
        )

        flat: list[dict] = [item for sublist in raw_results for item in sublist]

        if not flat:
            return pd.DataFrame(
                columns=["pvalue", "qvalue", "test_stat", "df",
                         "dispersion", "n_cells"]
            )

        df = pd.DataFrame(flat).set_index("pas_id")

        qvalues = false_discovery_control(df["pvalue"].values, method="bh")
        df["qvalue"] = qvalues

        return df[["pvalue", "qvalue", "test_stat", "df",
                   "dispersion", "n_cells"]].sort_values("qvalue")
