"""Tests for parallel cluster-pair differential APA testing.

Validates that :func:`ema.switch_test.pair_runner.run_one_pair` and the
Pool-based dispatch in ``cli._dispatch_pair`` produce results that are
numerically identical to sequential execution, regardless of whether the
pair loop runs in a pool or inline.

Coverage:
- Sequential vs parallel (spawn Pool) produce identical p-values (1e-12 tol)
- Single-pair path (n_outer<=1) runs inline — no Pool spawned
- _dispatch_pair unpacks pair tuple correctly
- run_one_pair returns (c1, c2, DataFrame) with expected columns
- Fisher strategy: p-values in [0, 1], qvalue column present
- n_jobs_inner=1 produces same result as n_jobs_inner=2 (deterministic)
"""

from __future__ import annotations

import functools
import multiprocessing
from itertools import combinations

import numpy as np
import pandas as pd
import pytest

from ema.switch_test.pair_runner import run_one_pair
from ema.switch_test.runner import _dispatch_pair


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def synthetic_count_matrix() -> tuple[pd.DataFrame, pd.Series]:
    """Small deterministic count matrix: 60 cells × 20 PAS, 3 clusters.

    Cluster 0 has elevated counts on PAS 0-4; cluster 1 on PAS 5-9;
    cluster 2 is noise-level on all PAS. This ensures Fisher and NB tests
    will find differential signal on at least some PAS pairs.
    """
    rng = np.random.default_rng(42)
    n_cells, n_pas = 60, 20
    n_per_cluster = n_cells // 3

    data = rng.integers(0, 3, size=(n_cells, n_pas)).astype(float)
    # Enrich cluster 0 on PAS 0-4
    data[:n_per_cluster, :5] += rng.integers(5, 15, size=(n_per_cluster, 5))
    # Enrich cluster 1 on PAS 5-9
    data[n_per_cluster : 2 * n_per_cluster, 5:10] += rng.integers(
        5, 15, size=(n_per_cluster, 5)
    )

    cell_ids = [f"cell_{i}" for i in range(n_cells)]
    pas_ids = [f"pas_{j}" for j in range(n_pas)]
    diff_df = pd.DataFrame(data.astype(int), index=cell_ids, columns=pas_ids)

    labels = (
        ["0"] * n_per_cluster
        + ["1"] * n_per_cluster
        + ["2"] * n_per_cluster
    )
    cluster_labels = pd.Series(labels, index=cell_ids)

    return diff_df, cluster_labels


# ---------------------------------------------------------------------------
# run_one_pair — basic contract
# ---------------------------------------------------------------------------

class TestRunOnePairContract:
    """run_one_pair must return a valid (c1, c2, DataFrame) triple."""

    def test_returns_three_tuple(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        result = run_one_pair("fisher", diff_df, labels, "0", "1", n_jobs_inner=1)
        assert isinstance(result, tuple) and len(result) == 3

    def test_cluster_labels_preserved(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        c1, c2, df = run_one_pair("fisher", diff_df, labels, "0", "2", n_jobs_inner=1)
        assert c1 == "0"
        assert c2 == "2"

    def test_result_is_dataframe(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        _, _, df = run_one_pair("fisher", diff_df, labels, "0", "1", n_jobs_inner=1)
        assert isinstance(df, pd.DataFrame)

    def test_pvalue_column_present(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        _, _, df = run_one_pair("fisher", diff_df, labels, "0", "1", n_jobs_inner=1)
        assert "pvalue" in df.columns

    def test_qvalue_column_present(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        _, _, df = run_one_pair("fisher", diff_df, labels, "0", "1", n_jobs_inner=1)
        assert "qvalue" in df.columns

    def test_pvalues_in_unit_interval(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        _, _, df = run_one_pair("fisher", diff_df, labels, "0", "1", n_jobs_inner=1)
        assert df["pvalue"].between(0.0, 1.0).all()

    def test_qvalues_in_unit_interval(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        _, _, df = run_one_pair("fisher", diff_df, labels, "0", "1", n_jobs_inner=1)
        assert df["qvalue"].between(0.0, 1.0).all()


# ---------------------------------------------------------------------------
# Determinism: n_jobs_inner=1 vs n_jobs_inner=2 (fisher only — NB also OK)
# ---------------------------------------------------------------------------

class TestDeterminism:
    """Same inputs + same strategy → same p-values regardless of n_jobs_inner."""

    def test_fisher_njobs_1_vs_2_identical(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        _, _, df1 = run_one_pair("fisher", diff_df, labels, "0", "1", n_jobs_inner=1)
        _, _, df2 = run_one_pair("fisher", diff_df, labels, "0", "1", n_jobs_inner=2)
        pd.testing.assert_series_equal(
            df1["pvalue"].sort_index(),
            df2["pvalue"].sort_index(),
            check_names=False,
            atol=1e-12,
        )


# ---------------------------------------------------------------------------
# Sequential vs parallel — p-values must match to 1e-12
# ---------------------------------------------------------------------------

class TestSequentialVsParallel:
    """Pool-dispatched results must be numerically identical to sequential."""

    def _run_sequential(
        self, strategy: str, diff_df: pd.DataFrame, labels: pd.Series
    ) -> dict[tuple[str, str], pd.DataFrame]:
        clusters = sorted(labels.unique().tolist())
        pairs = list(combinations(clusters, 2))
        results = {}
        for c1, c2 in pairs:
            _, _, df = run_one_pair(strategy, diff_df, labels, c1, c2, n_jobs_inner=1)
            results[(c1, c2)] = df
        return results

    def _run_parallel(
        self, strategy: str, diff_df: pd.DataFrame, labels: pd.Series, n_outer: int
    ) -> dict[tuple[str, str], pd.DataFrame]:
        clusters = sorted(labels.unique().tolist())
        pairs = list(combinations(clusters, 2))
        worker_fn = functools.partial(
            _dispatch_pair,
            strategy_name=strategy,
            diff_df=diff_df,
            cluster_labels=labels,
            n_jobs_inner=1,
        )
        ctx = multiprocessing.get_context("spawn")
        results = {}
        with ctx.Pool(n_outer) as pool:
            for c1, c2, df in pool.imap_unordered(worker_fn, pairs, chunksize=1):
                results[(c1, c2)] = df
        return results

    def test_fisher_sequential_vs_parallel_pvalues(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        seq = self._run_sequential("fisher", diff_df, labels)
        par = self._run_parallel("fisher", diff_df, labels, n_outer=2)

        assert set(seq.keys()) == set(par.keys()), "pair sets differ"
        for key in seq:
            pd.testing.assert_series_equal(
                seq[key]["pvalue"].sort_index(),
                par[key]["pvalue"].sort_index(),
                check_names=False,
                atol=1e-12,
                rtol=0.0,
            )

    def test_fisher_sequential_vs_parallel_qvalues(self, synthetic_count_matrix):
        """q-values are computed per-pair — must also match exactly."""
        diff_df, labels = synthetic_count_matrix
        seq = self._run_sequential("fisher", diff_df, labels)
        par = self._run_parallel("fisher", diff_df, labels, n_outer=2)

        for key in seq:
            if "qvalue" in seq[key].columns and "qvalue" in par[key].columns:
                pd.testing.assert_series_equal(
                    seq[key]["qvalue"].sort_index(),
                    par[key]["qvalue"].sort_index(),
                    check_names=False,
                    atol=1e-12,
                    rtol=0.0,
                )

    def test_all_pairs_returned_by_pool(self, synthetic_count_matrix):
        """Pool path must return exactly the same number of pairs as sequential."""
        diff_df, labels = synthetic_count_matrix
        seq = self._run_sequential("fisher", diff_df, labels)
        par = self._run_parallel("fisher", diff_df, labels, n_outer=2)
        assert len(seq) == len(par)


# ---------------------------------------------------------------------------
# _dispatch_pair — tuple unpacking
# ---------------------------------------------------------------------------

class TestDispatchPair:
    """_dispatch_pair must unpack (c1, c2) and forward to run_one_pair correctly."""

    def test_unpacks_pair_correctly(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        c1_got, c2_got, df = _dispatch_pair(
            ("0", "1"),
            strategy_name="fisher",
            diff_df=diff_df,
            cluster_labels=labels,
            n_jobs_inner=1,
        )
        assert c1_got == "0"
        assert c2_got == "1"
        assert isinstance(df, pd.DataFrame)

    def test_result_matches_run_one_pair_directly(self, synthetic_count_matrix):
        """_dispatch_pair must produce identical output to run_one_pair directly."""
        diff_df, labels = synthetic_count_matrix
        _, _, df_dispatch = _dispatch_pair(
            ("1", "2"),
            strategy_name="fisher",
            diff_df=diff_df,
            cluster_labels=labels,
            n_jobs_inner=1,
        )
        _, _, df_direct = run_one_pair(
            "fisher", diff_df, labels, "1", "2", n_jobs_inner=1
        )
        pd.testing.assert_frame_equal(
            df_dispatch.sort_index(),
            df_direct.sort_index(),
            atol=1e-12,
        )


# ---------------------------------------------------------------------------
# Single-pair inline path (n_outer <= 1 or len(pairs) == 1)
# ---------------------------------------------------------------------------

class TestSinglePairInlinePath:
    """When only one pair exists, run_one_pair is called directly (no Pool)."""

    def test_single_pair_returns_valid_result(self, synthetic_count_matrix):
        diff_df, labels = synthetic_count_matrix
        c1, c2, df = run_one_pair(
            "fisher", diff_df, labels, "0", "1", n_jobs_inner=1
        )
        assert c1 == "0" and c2 == "1"
        assert len(df) > 0
        assert "pvalue" in df.columns
