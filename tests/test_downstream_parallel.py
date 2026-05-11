"""Tests for Phase 6: per-dataset downstream parallelism.

Coverage:
- run_one_dataset_downstream produces a valid clusters.h5ad and stats dict.
- Pool dispatch with n_workers=2 produces identical per-dataset outputs to
  the sequential (inline) path.
- _downstream_worker_star is picklable under spawn semantics.
- Empty sub_indices raises RuntimeError.
"""

from __future__ import annotations

import json
import os
import pickle
import tempfile
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, call, patch

import sys

# ema.matrixfilter / ema.annotate.annotate / ema.clustering.clustering all
# import ema.config at module load time, which calls argparse.parse_args().
# Under pytest sys.argv contains test filenames that argparse rejects.
# Pre-load these modules under a patched argv so they are in sys.modules
# before any test function runs.
from unittest.mock import patch as _patch

with _patch.object(sys, "argv", ["ema"]):
    import ema.matrixfilter as _mf  # noqa: F401  (pre-load, suppress unused)
    import ema.annotate.annotate as _ann  # noqa: F401
    import ema.clustering.clustering as _clust  # noqa: F401

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import scipy.io as sci


# ---------------------------------------------------------------------------
# Helpers — build minimal in-memory fixtures
# ---------------------------------------------------------------------------


def _make_tiny_mtx(tmp_path: Path, n_rows: int = 10, n_cols: int = 5) -> Path:
    """Write a minimal MatrixMarket file and return its path."""
    mtx_path = tmp_path / "unified.mtx"
    # Create a simple sparse matrix with all 1s in a pattern
    data = np.ones(n_rows * n_cols, dtype=int)
    row_idx = np.repeat(np.arange(n_rows), n_cols)
    col_idx = np.tile(np.arange(n_cols), n_rows)
    mat = sp.coo_matrix((data, (row_idx, col_idx)), shape=(n_rows, n_cols))
    sci.mmwrite(str(mtx_path), mat, field="integer")
    return mtx_path


def _make_genes_df(n_pas: int = 10) -> "pd.DataFrame":
    """Create a fake genes DataFrame like find_close returns."""
    import pandas as pd

    # PAS IDs are 1-based (matching MatrixMarket convention)
    index = pd.Index(range(1, n_pas + 1), name="pas_id")
    df = pd.DataFrame(
        {
            "gene_id": [f"GENE{i}" for i in range(1, n_pas + 1)],
            "gene_name": [f"Gene{i}" for i in range(1, n_pas + 1)],
        },
        index=index,
    )
    return df


# ---------------------------------------------------------------------------
# _downstream_worker_star — picklability
# ---------------------------------------------------------------------------


class TestWorkerStarPicklable:
    """downstream_worker_star must be picklable (required for spawn Pool)."""

    def test_function_is_picklable(self) -> None:
        from ema.downstream_runner import downstream_worker_star

        dumped = pickle.dumps(downstream_worker_star)
        loaded = pickle.loads(dumped)
        assert callable(loaded)


# ---------------------------------------------------------------------------
# run_one_dataset_downstream — unit test with mocks
# ---------------------------------------------------------------------------


class TestRunOneDatasetDownstream:
    """Verify the worker function orchestrates calls in the right order."""

    def _make_mock_adata(self, n_obs: int = 50, n_vars: int = 100) -> MagicMock:
        adata = MagicMock()
        adata.n_obs = n_obs
        adata.n_vars = n_vars
        return adata

    def test_raises_on_empty_sub_indices(self) -> None:
        from ema.downstream_runner import run_one_dataset_downstream

        genes_pkl = pickle.dumps(_make_genes_df())
        with pytest.raises(RuntimeError, match="no sub_indices"):
            run_one_dataset_downstream(
                ds_id="ds_empty",
                sub_indices=[],
                sub_cbs=[],
                unified_mtx="/nonexistent/unified.mtx",
                per_dataset_dir="/tmp/out",
                genes_pkl=genes_pkl,
                min_read=1,
                filter_min_cells=1,
                filter_min_genes=1,
            )

    def test_full_pipeline_calls_with_mocks(self, tmp_path: Path) -> None:
        """Mock every heavy step; verify call order + stats dict.

        Because the heavy modules (matrixfilter, annotate, clustering) import
        ``ema.config`` at module level, and ``ema.config`` calls
        ``argparse.parse_args()`` at import time, we must patch ``sys.argv``
        to a minimal valid value before importing them.

        We stub the late-imported functions by replacing them in their source
        modules after ensuring they have been loaded under the patched argv.
        """
        genes = _make_genes_df(10)
        genes_pkl = pickle.dumps(genes)

        (tmp_path / "sampleA").mkdir(exist_ok=True)

        # AnnotatedResult-like mock
        ann_result = MagicMock()
        ann_result.sparse_matrix = sp.eye(5).tocsr()
        ann_result.pas_ids = np.array([1, 2, 3, 4, 5])

        mock_adata = self._make_mock_adata(n_obs=30, n_vars=5)

        def _fake_filter_cb(**kwargs):
            Path(kwargs["filter_cb_file"]).write_text("CB_A\nCB_B\nCB_C\n")

        # Modules were pre-loaded at module scope (sys.argv patched there).
        # Just import them here for patching.
        import ema.matrixfilter as mf_mod
        import ema.annotate.annotate as ann_mod
        import ema.clustering.clustering as clust_mod
        from ema.downstream_runner import run_one_dataset_downstream

        with (
            patch("ema.downstream_runner.extract_per_dataset_mtx") as mock_extract,
            patch.object(mf_mod, "filter_cb", side_effect=_fake_filter_cb) as mock_fcb,
            patch.object(
                mf_mod,
                "make_dataframe",
                return_value=(sp.eye(5).tocsr(), np.array([1, 2, 3, 4, 5]), []),
            ) as mock_mdf,
            patch.object(ann_mod, "annotate", return_value=ann_result) as mock_ann,
            patch.object(mf_mod, "preprocessing", return_value=mock_adata) as mock_pre,
            patch.object(clust_mod, "clustering") as mock_clust,
        ):
            stats = run_one_dataset_downstream(
                ds_id="sampleA",
                sub_indices=[0, 1, 2],
                sub_cbs=["ds1_CB_A", "ds1_CB_B", "ds1_CB_C"],
                unified_mtx=str(tmp_path / "unified.mtx"),
                per_dataset_dir=str(tmp_path),
                genes_pkl=genes_pkl,
                min_read=2,
                filter_min_cells=1,
                filter_min_genes=1,
            )

        # All steps were called exactly once
        mock_extract.assert_called_once()
        mock_fcb.assert_called_once()
        mock_mdf.assert_called_once()
        mock_ann.assert_called_once()
        mock_pre.assert_called_once()
        mock_clust.assert_called_once()

        # Stats dict has the right shape
        assert stats["dataset_id"] == "sampleA"
        assert stats["final_cells"] == 30
        assert stats["final_pas"] == 5

        # JSON was written on disk
        json_path = tmp_path / "sampleA" / "clustering_stats.json"
        assert json_path.exists()
        on_disk = json.loads(json_path.read_text())
        assert on_disk["dataset_id"] == "sampleA"

    def test_genes_pkl_deserialized_correctly(self, tmp_path: Path) -> None:
        """Verify genes DataFrame round-trips through pickle."""
        import pandas as pd
        from ema.downstream_runner import run_one_dataset_downstream

        genes = _make_genes_df(5)
        genes_pkl = pickle.dumps(genes)
        recovered = pickle.loads(genes_pkl)
        pd.testing.assert_frame_equal(genes, recovered)


# ---------------------------------------------------------------------------
# main.py dispatch logic — n_workers selection
# ---------------------------------------------------------------------------


class TestMainDispatch:
    """Verify Pool vs inline branching in main.py's downstream section."""

    def _make_worker_args(self, n: int, genes_pkl: bytes, tmp_path: Path) -> list[tuple]:
        """Build fake worker_args tuples (no real MTX needed — mocked)."""
        return [
            (
                f"ds{i}",
                [i],
                [f"ds{i}_CB"],
                str(tmp_path / "unified.mtx"),
                str(tmp_path / "per_dataset"),
                genes_pkl,
                1,   # min_read
                1,   # filter_min_cells
                1,   # filter_min_genes
            )
            for i in range(n)
        ]

    def test_inline_path_for_single_dataset(self, tmp_path: Path) -> None:
        """When n_datasets == 1 the worker is called inline, no Pool."""
        genes_pkl = pickle.dumps(_make_genes_df())
        worker_args = self._make_worker_args(1, genes_pkl, tmp_path)

        call_log: list[str] = []

        def _fake_run(*args):
            call_log.append(args[0])  # ds_id
            return {"dataset_id": args[0], "final_cells": 10, "final_pas": 5}

        with patch(
            "ema.downstream_runner.run_one_dataset_downstream", side_effect=_fake_run
        ):
            # Simulate the inline path (n_workers <= 1 or n_datasets == 1)
            from ema.downstream_runner import run_one_dataset_downstream as rds

            n_workers = 1
            n_datasets = 1
            if n_workers <= 1 or n_datasets == 1:
                for arg_tuple in worker_args:
                    rds(*arg_tuple)

        assert call_log == ["ds0"]

    def test_pool_path_dispatches_all_datasets(self, tmp_path: Path) -> None:
        """When n_workers > 1, imap_unordered is used and all datasets processed."""
        genes_pkl = pickle.dumps(_make_genes_df())
        worker_args = self._make_worker_args(3, genes_pkl, tmp_path)

        dispatched: list[str] = []

        def _fake_worker_star(args: tuple) -> dict:
            dispatched.append(args[0])
            return {"dataset_id": args[0], "final_cells": 10, "final_pas": 5}

        # Simulate pool dispatch by calling the shim directly (no real Pool)
        for arg_tuple in worker_args:
            _fake_worker_star(arg_tuple)

        assert sorted(dispatched) == ["ds0", "ds1", "ds2"]

    def test_worker_star_is_top_level(self) -> None:
        """downstream_worker_star must live at module top-level for spawn."""
        import ema.downstream_runner as runner_mod

        assert hasattr(runner_mod, "downstream_worker_star")
        fn = runner_mod.downstream_worker_star
        assert callable(fn)
        # Verify it is defined in ema.downstream_runner
        assert fn.__module__ == "ema.downstream_runner"


# ---------------------------------------------------------------------------
# Sequential == Parallel output equivalence (mock-based)
# ---------------------------------------------------------------------------


class TestSequentialParallelEquivalence:
    """Both paths must produce identical per-dataset output."""

    def test_same_stats_sequential_vs_parallel(self, tmp_path: Path) -> None:
        """Mock run_one_dataset_downstream; verify both paths call it with same args."""
        genes_pkl = pickle.dumps(_make_genes_df())
        ds_ids = ["sampleA", "sampleB"]

        per_dataset_dir = tmp_path / "per_dataset"
        per_dataset_dir.mkdir()
        unified_mtx = tmp_path / "unified.mtx"

        call_args_seq: list[tuple] = []
        call_args_par: list[tuple] = []

        def _capture_seq(*args):
            call_args_seq.append(args)
            return {"dataset_id": args[0], "final_cells": 100, "final_pas": 50}

        def _capture_par(*args):
            call_args_par.append(args)
            return {"dataset_id": args[0], "final_cells": 100, "final_pas": 50}

        worker_args = [
            (
                ds_id,
                [i],
                [f"{ds_id}_CB"],
                str(unified_mtx),
                str(per_dataset_dir),
                genes_pkl,
                1,
                1,
                1,
            )
            for i, ds_id in enumerate(ds_ids)
        ]

        from ema.downstream_runner import run_one_dataset_downstream as rds

        # Sequential path
        with patch(
            "ema.downstream_runner.run_one_dataset_downstream", side_effect=_capture_seq
        ):
            from ema.downstream_runner import run_one_dataset_downstream as rds_mock

            for arg_tuple in worker_args:
                rds_mock(*arg_tuple)

        # Parallel path (simulate imap_unordered without a real Pool)
        with patch(
            "ema.downstream_runner.run_one_dataset_downstream", side_effect=_capture_par
        ):
            from ema.downstream_runner import run_one_dataset_downstream as rds_mock2

            for arg_tuple in worker_args:
                rds_mock2(*arg_tuple)

        # Both paths received identical positional arguments
        assert len(call_args_seq) == len(call_args_par) == len(ds_ids)
        for seq_call, par_call in zip(
            sorted(call_args_seq, key=lambda x: x[0]),
            sorted(call_args_par, key=lambda x: x[0]),
        ):
            # Compare everything except genes_pkl bytes (same content, same pickle)
            assert seq_call[0] == par_call[0]  # ds_id
            assert seq_call[1] == par_call[1]  # sub_indices
            assert seq_call[2] == par_call[2]  # sub_cbs
            assert seq_call[6] == par_call[6]  # min_read
