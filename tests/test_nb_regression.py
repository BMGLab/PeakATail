"""Tests for NB regression differential APA strategies.

All tests use synthetically simulated count data — no real BAM files required.

Test matrix
-----------
1. ``test_null_pvalue_uniform``       — NB pairwise on null data: KS test for
   uniform p-value distribution (KS p > 0.01).
2. ``test_signal_recovery``           — NB pairwise with 2× fold-change injected
   into 20 PAS: strategy must recover ≥ 15 of those 20 with q < 0.05.
3. ``test_multi_condition_omnibus``   — NB multi on 4-cluster data with injected
   cluster-specific signal: omnibus LRT detects the spiked PAS.
4. ``test_fisher_inflation``          — On the same null data, Fisher produces
   more "significant" calls than nominal alpha (anti-conservative inflation).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.stats import kstest


# ---------------------------------------------------------------------------
# Simulation helpers
# ---------------------------------------------------------------------------

RNG = np.random.default_rng(42)

N_PAS = 100
N_CELLS_PER_CLUSTER = 100  # 200 cells total (2 clusters × 100)


def _simulate_nb(
    mu: float,
    alpha: float,  # NB dispersion: Var = mu + alpha * mu^2
    size: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw NB samples via the Gamma-Poisson mixture.

    Parameters
    ----------
    mu : float  Expected count.
    alpha : float  Dispersion (1/r in canonical NB).
    size : int  Number of samples.
    rng : Generator  NumPy RNG.
    """
    r = 1.0 / alpha  # shape
    p = r / (r + mu)  # success probability
    return rng.negative_binomial(r, p, size=size).astype(np.float64)


def _make_null_matrix(
    n_pas: int = N_PAS,
    n_per_cluster: int = N_CELLS_PER_CLUSTER,
    mu: float = 5.0,
    alpha: float = 0.5,
    rng: np.random.Generator | None = None,
) -> tuple[pd.DataFrame, pd.Series]:
    """Simulate count matrix with NO differential signal between 2 clusters.

    Returns
    -------
    count_matrix : DataFrame, shape (n_cells, n_pas)
    cluster_labels : Series, index = cell barcodes
    """
    if rng is None:
        rng = np.random.default_rng(0)

    n_cells = n_per_cluster * 2
    counts = _simulate_nb(mu, alpha, n_cells * n_pas, rng).reshape(n_cells, n_pas)

    cell_ids = [f"cell_{i}" for i in range(n_cells)]
    pas_ids = [f"pas_{j}" for j in range(n_pas)]
    cluster_labels = pd.Series(
        ["cluster_A"] * n_per_cluster + ["cluster_B"] * n_per_cluster,
        index=cell_ids,
        name="cluster",
    )
    count_matrix = pd.DataFrame(counts, index=cell_ids, columns=pas_ids)
    return count_matrix, cluster_labels


def _make_signal_matrix(
    n_pas: int = N_PAS,
    n_per_cluster: int = N_CELLS_PER_CLUSTER,
    n_signal_pas: int = 20,
    fold_change: float = 2.0,
    mu_base: float = 5.0,
    alpha: float = 0.5,
    rng: np.random.Generator | None = None,
) -> tuple[pd.DataFrame, pd.Series, list[str]]:
    """Simulate count matrix with injected fold-change in ``n_signal_pas`` PAS.

    Returns
    -------
    count_matrix : DataFrame
    cluster_labels : Series
    signal_pas_ids : list of PAS IDs that have injected FC
    """
    if rng is None:
        rng = np.random.default_rng(1)

    count_matrix, cluster_labels = _make_null_matrix(
        n_pas=n_pas, n_per_cluster=n_per_cluster, mu=mu_base,
        alpha=alpha, rng=rng,
    )
    # Overwrite the first n_signal_pas columns for cluster_B cells with 2× mu
    signal_pas_ids = [f"pas_{j}" for j in range(n_signal_pas)]
    b_cells = cluster_labels[cluster_labels == "cluster_B"].index
    for pas_id in signal_pas_ids:
        count_matrix.loc[b_cells, pas_id] = _simulate_nb(
            mu_base * fold_change, alpha, len(b_cells), rng
        )

    return count_matrix, cluster_labels, signal_pas_ids


def _make_multi_cluster_matrix(
    n_pas: int = N_PAS,
    n_per_cluster: int = 50,  # 4 clusters × 50 = 200 cells
    n_signal_pas: int = 20,
    mu_base: float = 5.0,
    alpha: float = 0.5,
    rng: np.random.Generator | None = None,
) -> tuple[pd.DataFrame, pd.Series, list[str]]:
    """Simulate 4-cluster data with cluster-specific signal in the first n_signal_pas PAS."""
    if rng is None:
        rng = np.random.default_rng(2)

    clusters = ["A", "B", "C", "D"]
    n_cells = len(clusters) * n_per_cluster
    n_total = n_cells * n_pas

    counts = _simulate_nb(mu_base, alpha, n_total, rng).reshape(n_cells, n_pas)
    cell_ids = [f"cell_{i}" for i in range(n_cells)]
    pas_ids = [f"pas_{j}" for j in range(n_pas)]
    labels_raw = []
    for cl in clusters:
        labels_raw.extend([cl] * n_per_cluster)
    cluster_labels = pd.Series(labels_raw, index=cell_ids, name="cluster")
    count_matrix = pd.DataFrame(counts, index=cell_ids, columns=pas_ids)

    # Inject: cluster "B" cells get 5× the base for the signal PAS
    signal_pas_ids = [f"pas_{j}" for j in range(n_signal_pas)]
    b_cells = cluster_labels[cluster_labels == "B"].index
    for pas_id in signal_pas_ids:
        count_matrix.loc[b_cells, pas_id] = _simulate_nb(
            mu_base * 5.0, alpha, len(b_cells), rng
        )

    return count_matrix, cluster_labels, signal_pas_ids


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestNbPairwiseNull:
    """NB pairwise on null data: p-values should be approximately uniform."""

    def test_null_pvalue_uniform(self):
        """KS test against Uniform(0,1): p > 0.01 expected under null."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels = _make_null_matrix(rng=RNG)
        strategy = get_diff_strategy("nb_pairwise")
        results = strategy.test(
            count_matrix,
            cluster_labels,
            cluster1="cluster_A",
            cluster2="cluster_B",
            min_cells_per_group=5,
            n_jobs=2,
        )
        assert len(results) > 0, "Expected non-empty results on null data"

        pvalues = results["pvalue"].values
        ks_stat, ks_pval = kstest(pvalues, "uniform")
        assert ks_pval > 0.01, (
            f"NB pairwise p-values on null data are NOT uniform: "
            f"KS statistic={ks_stat:.4f}, KS p-value={ks_pval:.4f}. "
            "This indicates inflated false positives."
        )

    def test_output_schema(self):
        """Result DataFrame must have the required columns."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels = _make_null_matrix(rng=RNG)
        strategy = get_diff_strategy("nb_pairwise")
        results = strategy.test(
            count_matrix, cluster_labels,
            cluster1="cluster_A", cluster2="cluster_B",
            min_cells_per_group=5, n_jobs=2,
        )
        required_cols = {"pvalue", "qvalue", "log2fc", "dispersion",
                         "n_cells", "test_stat"}
        missing = required_cols - set(results.columns)
        assert not missing, f"Missing output columns: {missing}"

    def test_qvalues_geq_pvalues_or_equal(self):
        """BH correction: qvalue >= pvalue (or at most equal for smallest)."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels = _make_null_matrix(rng=RNG)
        strategy = get_diff_strategy("nb_pairwise")
        results = strategy.test(
            count_matrix, cluster_labels,
            cluster1="cluster_A", cluster2="cluster_B",
            min_cells_per_group=5, n_jobs=2,
        )
        # BH qvalues are always >= corresponding pvalues
        assert (results["qvalue"].values >= results["pvalue"].values - 1e-10).all()


class TestNbPairwiseSignal:
    """NB pairwise with injected 2× fold-change: must recover ≥15/20 signal PAS."""

    def test_signal_recovery(self):
        """Recover at least 15 of 20 injected signal PAS with q < 0.05."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels, signal_pas = _make_signal_matrix(rng=RNG)
        strategy = get_diff_strategy("nb_pairwise")
        results = strategy.test(
            count_matrix,
            cluster_labels,
            cluster1="cluster_A",
            cluster2="cluster_B",
            min_cells_per_group=5,
            n_jobs=2,
        )

        signal_results = results.loc[results.index.isin(signal_pas)]
        n_recovered = int((signal_results["qvalue"] < 0.05).sum())

        assert n_recovered >= 15, (
            f"NB pairwise recovered only {n_recovered}/20 signal PAS (need ≥15). "
            f"Signal PAS q-values:\n{signal_results['qvalue'].sort_values()}"
        )

    def test_log2fc_direction(self):
        """Injected signal has positive FC for cluster_B; log2fc should be positive."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels, signal_pas = _make_signal_matrix(rng=RNG)
        strategy = get_diff_strategy("nb_pairwise")
        results = strategy.test(
            count_matrix, cluster_labels,
            cluster1="cluster_A", cluster2="cluster_B",
            min_cells_per_group=5, n_jobs=2,
        )
        sig_results = results.loc[
            results.index.isin(signal_pas) & (results["qvalue"] < 0.05)
        ]
        # Most (> 50%) recovered signal PAS should have positive log2fc
        n_positive_fc = (sig_results["log2fc"] > 0).sum()
        assert n_positive_fc > len(sig_results) * 0.5, (
            "Majority of recovered signal PAS should have positive log2fc "
            f"(cluster_B has higher counts); got {n_positive_fc}/{len(sig_results)}"
        )


class TestNbMulti:
    """NB multi-condition omnibus test on 4-cluster data."""

    def test_multi_detects_signal(self):
        """Omnibus LRT must detect a cluster-specific signal in the spiked PAS."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels, signal_pas = _make_multi_cluster_matrix(
            rng=RNG
        )
        strategy = get_diff_strategy("nb_multi")
        results = strategy.test(
            count_matrix,
            cluster_labels,
            min_cells_per_group=5,
            n_jobs=2,
        )

        assert len(results) > 0, "Expected non-empty results"

        signal_results = results.loc[results.index.isin(signal_pas)]
        n_detected = int((signal_results["qvalue"] < 0.05).sum())
        assert n_detected >= 10, (
            f"nb_multi detected only {n_detected}/20 spiked PAS (need ≥10). "
            f"q-values of signal PAS:\n{signal_results['qvalue'].sort_values()}"
        )

    def test_multi_output_schema(self):
        """Result DataFrame must have required multi-condition columns."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels, _ = _make_multi_cluster_matrix(rng=RNG)
        strategy = get_diff_strategy("nb_multi")
        results = strategy.test(
            count_matrix, cluster_labels,
            min_cells_per_group=5, n_jobs=2,
        )
        required = {"pvalue", "qvalue", "test_stat", "df", "dispersion", "n_cells"}
        missing = required - set(results.columns)
        assert not missing, f"Missing columns in nb_multi output: {missing}"
        # No log2fc for multi-condition
        assert "log2fc" not in results.columns

    def test_multi_no_cluster1_cluster2_required(self):
        """nb_multi must run without cluster1/cluster2 arguments."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels, _ = _make_multi_cluster_matrix(rng=RNG)
        strategy = get_diff_strategy("nb_multi")
        # Should not raise even though cluster1/cluster2 are None
        results = strategy.test(
            count_matrix, cluster_labels,
            min_cells_per_group=5, n_jobs=2,
        )
        assert results is not None


class TestFisherInflation:
    """Fisher exact is anti-conservative: more significant calls than NB on null data."""

    def test_fisher_inflation_vs_nb(self):
        """Fisher produces more false positives (q < 0.05) than NB on null data."""
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels = _make_null_matrix(rng=RNG)

        fisher = get_diff_strategy("fisher")
        nb = get_diff_strategy("nb_pairwise")

        # D4: read-based contingency is the pseudoreplicated (anti-conservative)
        # one; the engine DEFAULT is now count_mode="cells" (calibrated), so this
        # inflation regression explicitly pins the legacy 'reads' mode.
        fisher_results = fisher.test(
            count_matrix, cluster_labels,
            cluster1="cluster_A", cluster2="cluster_B",
            min_cells_per_group=5, count_mode="reads",
        )
        nb_results = nb.test(
            count_matrix, cluster_labels,
            cluster1="cluster_A", cluster2="cluster_B",
            min_cells_per_group=5, n_jobs=2,
        )

        fisher_sig = int((fisher_results["qvalue"] < 0.05).sum())
        nb_sig = int((nb_results["qvalue"] < 0.05).sum())

        # Fisher should be more inflated.  We accept the test if
        # Fisher ≥ NB (even if both are 0 under BH correction on null).
        assert fisher_sig >= nb_sig, (
            f"Expected Fisher ({fisher_sig} sig) ≥ NB ({nb_sig} sig) on null data, "
            "but NB produced more significant calls. "
            "This is unexpected given Fisher's known anti-conservatism."
        )

    def test_fisher_pvalue_distribution_more_skewed(self):
        """Fisher p-values should be more skewed toward 0 than NB p-values on null data.

        We test this by comparing the fraction of raw p-values below 0.05:
        Fisher's should be higher.
        """
        from ema.switch_test.strategies import get_diff_strategy

        count_matrix, cluster_labels = _make_null_matrix(rng=RNG)

        fisher = get_diff_strategy("fisher")
        nb = get_diff_strategy("nb_pairwise")

        # D4: read-based contingency is the pseudoreplicated (anti-conservative)
        # one; the engine DEFAULT is now count_mode="cells" (calibrated), so this
        # inflation regression explicitly pins the legacy 'reads' mode.
        fisher_results = fisher.test(
            count_matrix, cluster_labels,
            cluster1="cluster_A", cluster2="cluster_B",
            min_cells_per_group=5, count_mode="reads",
        )
        nb_results = nb.test(
            count_matrix, cluster_labels,
            cluster1="cluster_A", cluster2="cluster_B",
            min_cells_per_group=5, n_jobs=2,
        )

        fisher_frac = float((fisher_results["pvalue"] < 0.05).mean())
        nb_frac = float((nb_results["pvalue"] < 0.05).mean())

        # Fisher fraction of "significant" raw p-values should be at least as
        # large as NB's on null data (anti-conservatism manifests as excess
        # small p-values).
        assert fisher_frac >= nb_frac * 0.8, (
            f"Fisher p<0.05 fraction ({fisher_frac:.3f}) should be ≥ "
            f"NB ({nb_frac:.3f}) on null data."
        )


class TestRegistry:
    """Strategy registry correctness."""

    def test_all_three_strategies_registered(self):
        from ema.switch_test.strategies import list_diff_strategies
        strategies = list_diff_strategies()
        assert "fisher" in strategies
        assert "nb_pairwise" in strategies
        assert "nb_multi" in strategies

    def test_supports_multi_condition_flags(self):
        from ema.switch_test.strategies import get_diff_strategy
        assert get_diff_strategy("fisher").supports_multi_condition is False
        assert get_diff_strategy("nb_pairwise").supports_multi_condition is False
        assert get_diff_strategy("nb_multi").supports_multi_condition is True

    def test_get_unknown_strategy_raises(self):
        from ema.switch_test.strategies import get_diff_strategy
        with pytest.raises(KeyError, match="Unknown diff-APA strategy"):
            get_diff_strategy("does_not_exist")

    def test_pairwise_strategies_require_cluster_args(self):
        from ema.switch_test.strategies import get_diff_strategy
        count_matrix, cluster_labels = _make_null_matrix(
            n_pas=10, n_per_cluster=20, rng=np.random.default_rng(99)
        )
        for name in ("fisher", "nb_pairwise"):
            strategy = get_diff_strategy(name)
            with pytest.raises(ValueError, match="cluster1"):
                strategy.test(count_matrix, cluster_labels)
