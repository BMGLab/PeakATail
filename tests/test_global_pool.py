"""Tests for the global flat-pool dispatch (Phase 3).

Coverage:
- JobSpec is a frozen dataclass (hashable, picklable)
- build_job_specs produces correct count across datasets × chroms × tiles × directions
- build_job_specs respects per_bam_tile_sizes
- build_job_specs assigns monotonically increasing job_ids
- run_all_jobs groups results by (dataset_id, direction)
- run_all_jobs handles empty job list gracefully
- run_all_jobs result total equals input job count
"""

from __future__ import annotations

import pickle
from dataclasses import FrozenInstanceError
from unittest.mock import MagicMock, patch

import pytest

from ema.countmatrix.tile_runner import (
    JobSpec,
    build_job_specs,
    run_all_jobs,
)


# ---------------------------------------------------------------------------
# JobSpec — dataclass contract
# ---------------------------------------------------------------------------


class TestJobSpec:
    """JobSpec must be frozen (immutable), hashable, and pickle-safe."""

    def _make(self, **kwargs) -> JobSpec:
        defaults = dict(
            job_id=0,
            dataset_id="ds1",
            bam_path="/fake/sample.bam",
            chrom="chr1",
            tile_start=0,
            tile_end=5_000_000,
            fetch_start=0,
            fetch_end=5_010_000,
            direction=False,
        )
        defaults.update(kwargs)
        return JobSpec(**defaults)

    def test_is_hashable(self):
        spec = self._make()
        assert hash(spec) is not None
        s = {spec}
        assert spec in s

    def test_is_picklable(self):
        spec = self._make(job_id=7, dataset_id="sample_A", direction=True)
        round_tripped = pickle.loads(pickle.dumps(spec))
        assert round_tripped == spec

    def test_is_frozen(self):
        spec = self._make()
        with pytest.raises((FrozenInstanceError, AttributeError)):
            spec.job_id = 99  # type: ignore[misc]

    def test_fields_accessible(self):
        spec = self._make(
            job_id=3,
            dataset_id="ds2",
            bam_path="/tmp/b.bam",
            chrom="chr22",
            tile_start=1_000_000,
            tile_end=6_000_000,
            fetch_start=990_000,
            fetch_end=6_010_000,
            direction=True,
            default_threshold=7,
            strategy_name="lambda_gradient",
        )
        assert spec.job_id == 3
        assert spec.dataset_id == "ds2"
        assert spec.chrom == "chr22"
        assert spec.direction is True
        assert spec.default_threshold == 7
        assert spec.strategy_name == "lambda_gradient"


# ---------------------------------------------------------------------------
# build_job_specs
# ---------------------------------------------------------------------------


def _fake_get_chromosomes(bam_path: str) -> list[tuple[str, int]]:
    """Return 2 fake chromosomes of 10 Mb each."""
    return [("chr1", 10_000_000), ("chr2", 10_000_000)]


class TestBuildJobSpecs:
    """build_job_specs must enumerate all (dataset × chrom × tile × direction) combos."""

    def test_single_dataset_single_tile_per_chrom(self):
        """tile_size > chrom_length → 1 tile per chrom × 2 dirs × 1 dataset × 2 chroms = 4 jobs."""
        bam_list = [("ds1", "/fake/ds1.bam")]
        with patch(
            "ema.countmatrix.tile_runner.get_chromosomes",
            side_effect=_fake_get_chromosomes,
        ):
            specs = build_job_specs(
                bam_list=bam_list,
                tile_size=50_000_000,  # larger than 10 Mb chroms → 1 tile each
                directions=[False, True],
            )
        assert len(specs) == 4  # 1 ds × 2 chroms × 1 tile × 2 dirs

    def test_single_dataset_two_tiles_per_chrom(self):
        """tile_size = 5 Mb, chrom = 10 Mb → 2 tiles × 2 chroms × 2 dirs = 8 jobs."""
        bam_list = [("ds1", "/fake/ds1.bam")]
        with patch(
            "ema.countmatrix.tile_runner.get_chromosomes",
            side_effect=_fake_get_chromosomes,
        ):
            specs = build_job_specs(
                bam_list=bam_list,
                tile_size=5_000_000,
                directions=[False, True],
            )
        assert len(specs) == 8  # 1 ds × 2 chroms × 2 tiles × 2 dirs

    def test_two_datasets_job_count(self):
        """2 datasets × 2 chroms × 1 tile × 2 dirs = 8 jobs."""
        bam_list = [("ds1", "/fake/ds1.bam"), ("ds2", "/fake/ds2.bam")]
        with patch(
            "ema.countmatrix.tile_runner.get_chromosomes",
            side_effect=_fake_get_chromosomes,
        ):
            specs = build_job_specs(
                bam_list=bam_list,
                tile_size=50_000_000,
                directions=[False, True],
            )
        assert len(specs) == 8  # 2 ds × 2 chroms × 1 tile × 2 dirs

    def test_job_ids_are_monotonically_increasing(self):
        bam_list = [("ds1", "/fake/ds1.bam"), ("ds2", "/fake/ds2.bam")]
        with patch(
            "ema.countmatrix.tile_runner.get_chromosomes",
            side_effect=_fake_get_chromosomes,
        ):
            specs = build_job_specs(
                bam_list=bam_list,
                tile_size=50_000_000,
                directions=[False, True],
            )
        ids = [s.job_id for s in specs]
        assert ids == list(range(len(specs)))

    def test_dataset_ids_correct(self):
        bam_list = [("alpha", "/fake/a.bam"), ("beta", "/fake/b.bam")]
        with patch(
            "ema.countmatrix.tile_runner.get_chromosomes",
            side_effect=_fake_get_chromosomes,
        ):
            specs = build_job_specs(
                bam_list=bam_list,
                tile_size=50_000_000,
                directions=[False],
            )
        ds_ids = {s.dataset_id for s in specs}
        assert ds_ids == {"alpha", "beta"}

    def test_per_bam_tile_sizes_override(self):
        """When per_bam_tile_sizes is given, the override is used instead of tile_size."""
        bam_list = [("ds1", "/fake/a.bam")]
        per_bam = {"/fake/a.bam": 5_000_000}  # 2 tiles per 10 Mb chrom

        with patch(
            "ema.countmatrix.tile_runner.get_chromosomes",
            side_effect=_fake_get_chromosomes,
        ):
            specs = build_job_specs(
                bam_list=bam_list,
                tile_size=50_000_000,  # would give 1 tile without override
                directions=[False, True],
                per_bam_tile_sizes=per_bam,
            )
        # 2 chroms × 2 tiles × 2 dirs = 8
        assert len(specs) == 8

    def test_single_direction_only(self):
        """Passing directions=[False] halves the job count."""
        bam_list = [("ds1", "/fake/a.bam")]
        with patch(
            "ema.countmatrix.tile_runner.get_chromosomes",
            side_effect=_fake_get_chromosomes,
        ):
            both = build_job_specs(bam_list, tile_size=50_000_000, directions=[False, True])
            one = build_job_specs(bam_list, tile_size=50_000_000, directions=[False])
        assert len(both) == 2 * len(one)

    def test_empty_bam_list_returns_empty(self):
        specs = build_job_specs(bam_list=[], tile_size=25_000_000)
        assert specs == []


# ---------------------------------------------------------------------------
# run_all_jobs
# ---------------------------------------------------------------------------


def _make_fake_tile_worker_result(spec: JobSpec) -> dict:
    """Produce a synthetic tile_worker result for a JobSpec."""
    return {
        "tile_id": spec.job_id,
        "dataset_id": spec.dataset_id,
        "direction": spec.direction,
        "chrom": spec.chrom,
        "tile_start": spec.tile_start,
        "tile_end": spec.tile_end,
        "bed_path": "/tmp/fake.bed",
        "mtx_path": "/tmp/fake.mtx",
        "cb_path": "/tmp/fake.cb.tsv",
        "workdir": "/tmp/fake_workdir",
    }


class TestRunAllJobs:
    """run_all_jobs must group by (dataset_id, direction) and cover all jobs."""

    def _build_specs(
        self,
        n_datasets: int = 2,
        tile_size: int = 50_000_000,
        directions: list[bool] | None = None,
    ) -> list[JobSpec]:
        bam_list = [(f"ds{i}", f"/fake/ds{i}.bam") for i in range(n_datasets)]
        with patch(
            "ema.countmatrix.tile_runner.get_chromosomes",
            side_effect=_fake_get_chromosomes,
        ):
            return build_job_specs(
                bam_list=bam_list,
                tile_size=tile_size,
                directions=directions or [False, True],
            )

    def test_empty_jobs_returns_empty_dict(self):
        result = run_all_jobs([], n_workers=2)
        assert result == {}

    def test_grouping_keys_match_dataset_direction_combos(self):
        specs = self._build_specs(n_datasets=2, directions=[False, True])
        # Patch tile_worker so we don't actually spawn processes
        with patch(
            "ema.countmatrix.tile_runner.tile_worker",
            side_effect=_make_fake_tile_worker_result,
        ):
            with patch("multiprocessing.get_context") as mock_ctx:
                # Build a mock pool whose imap_unordered applies the real side_effect
                mock_pool = MagicMock()
                mock_pool.__enter__ = lambda s: s
                mock_pool.__exit__ = MagicMock(return_value=False)
                mock_pool.imap_unordered = lambda fn, jobs, chunksize=1: (
                    _make_fake_tile_worker_result(j) for j in jobs
                )
                mock_ctx.return_value.Pool.return_value = mock_pool

                grouped = run_all_jobs(specs, n_workers=4)

        expected_keys = {
            ("ds0", False), ("ds0", True),
            ("ds1", False), ("ds1", True),
        }
        assert set(grouped.keys()) == expected_keys

    def test_total_results_equals_job_count(self):
        specs = self._build_specs(n_datasets=2, tile_size=5_000_000, directions=[False, True])
        with patch("multiprocessing.get_context") as mock_ctx:
            mock_pool = MagicMock()
            mock_pool.__enter__ = lambda s: s
            mock_pool.__exit__ = MagicMock(return_value=False)
            mock_pool.imap_unordered = lambda fn, jobs, chunksize=1: (
                _make_fake_tile_worker_result(j) for j in jobs
            )
            mock_ctx.return_value.Pool.return_value = mock_pool

            grouped = run_all_jobs(specs, n_workers=4)

        total_results = sum(len(v) for v in grouped.values())
        assert total_results == len(specs)

    def test_n_workers_capped_at_job_count(self):
        """Pool size must not exceed number of jobs."""
        specs = self._build_specs(n_datasets=1, tile_size=50_000_000, directions=[False])
        # 1 ds × 2 chroms × 1 tile × 1 dir = 2 jobs

        with patch("multiprocessing.get_context") as mock_ctx:
            mock_pool = MagicMock()
            mock_pool.__enter__ = lambda s: s
            mock_pool.__exit__ = MagicMock(return_value=False)
            mock_pool.imap_unordered = lambda fn, jobs, chunksize=1: (
                _make_fake_tile_worker_result(j) for j in jobs
            )
            mock_ctx.return_value.Pool.return_value = mock_pool

            run_all_jobs(specs, n_workers=100)

        # Pool was called with min(100, len(specs)) = 2
        call_kwargs = mock_ctx.return_value.Pool.call_args
        assert call_kwargs.kwargs["processes"] == len(specs)
