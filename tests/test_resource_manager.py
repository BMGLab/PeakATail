"""Tests for ResourceManager and the get_resource_manager() singleton accessor.

Coverage:
- get_n_jobs() respects user_max_threads ceiling
- get_n_jobs() with no ceiling returns sensible value (>=1, <=physical_cpu)
- split_jobs() obeys constraint outer × inner <= total_budget
- split_jobs() with small n_outer works correctly
- Singleton identity and reset behaviour
- Singleton picks up user_max_threads after reset + mock
"""

from __future__ import annotations

import sys
import types
from unittest.mock import MagicMock, patch

import pytest

from ema.utils.resource_manager import ResourceManager
from ema.utils import get_resource_manager, reset_resource_manager


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rm(user_max_threads: int | None = None, **kwargs) -> ResourceManager:
    """Construct a ResourceManager with a stable RAM/CPU environment."""
    rm = ResourceManager(user_max_threads=user_max_threads, **kwargs)
    return rm


# ---------------------------------------------------------------------------
# get_n_jobs — user_max_threads ceiling
# ---------------------------------------------------------------------------

class TestGetNJobsWithCeiling:
    """ResourceManager.get_n_jobs() must never exceed user_max_threads."""

    def test_ceiling_of_1_returns_1(self):
        rm = _rm(user_max_threads=1)
        assert rm.get_n_jobs(per_worker_mb=200) == 1

    def test_ceiling_of_4_returns_at_most_4(self):
        rm = _rm(user_max_threads=4)
        assert rm.get_n_jobs(per_worker_mb=200) <= 4

    def test_ceiling_of_4_returns_at_least_1(self):
        rm = _rm(user_max_threads=4)
        assert rm.get_n_jobs(per_worker_mb=200) >= 1

    def test_ceiling_applied_when_system_would_allow_more(self):
        """Patch free_ram_mb and physical_cpu_count to simulate a large machine."""
        rm = _rm(user_max_threads=3)
        with (
            patch.object(rm, "free_ram_mb", return_value=128_000),
            patch.object(rm, "physical_cpu_count", return_value=64),
        ):
            result = rm.get_n_jobs(per_worker_mb=200)
        assert result == 3

    def test_ceiling_not_exceeded_under_various_per_worker_mb(self):
        for mb in [50, 200, 500, 1000]:
            rm = _rm(user_max_threads=2)
            assert rm.get_n_jobs(per_worker_mb=mb) <= 2


# ---------------------------------------------------------------------------
# get_n_jobs — no user ceiling
# ---------------------------------------------------------------------------

class TestGetNJobsNoCeiling:
    """With no user ceiling, get_n_jobs() should return a sensible value."""

    def test_no_ceiling_returns_at_least_1(self):
        rm = _rm(user_max_threads=None)
        assert rm.get_n_jobs(per_worker_mb=200) >= 1

    def test_no_ceiling_returns_at_most_physical_cpu(self):
        rm = _rm(user_max_threads=None)
        n = rm.get_n_jobs(per_worker_mb=200)
        # Can't exceed the hard cap or the CPU count (whichever is smaller)
        assert n <= max(1, rm.physical_cpu_count())

    def test_no_ceiling_respects_hard_cap(self):
        rm = ResourceManager(user_max_threads=None, hard_n_jobs_cap=4)
        with (
            patch.object(rm, "free_ram_mb", return_value=512_000),
            patch.object(rm, "physical_cpu_count", return_value=128),
        ):
            result = rm.get_n_jobs(per_worker_mb=200)
        assert result <= 4

    def test_stage_argument_accepted(self):
        """stage= is a label only; should not change the numeric result."""
        rm = _rm(user_max_threads=None)
        n1 = rm.get_n_jobs(per_worker_mb=200, stage="nb_pairwise")
        n2 = rm.get_n_jobs(per_worker_mb=200, stage="pdui")
        assert n1 == n2  # same budget, different label

    def test_low_ram_gives_fewer_workers(self):
        rm = _rm(user_max_threads=None, cpu_reserve=0, hard_n_jobs_cap=32)
        with (
            patch.object(rm, "free_ram_mb", return_value=400),
            patch.object(rm, "physical_cpu_count", return_value=32),
        ):
            # 400 MB * 0.7 = 280 MB budget; 280 / 200 per worker = 1 worker
            result = rm.get_n_jobs(per_worker_mb=200)
        assert result == 1


# ---------------------------------------------------------------------------
# split_jobs
# ---------------------------------------------------------------------------

class TestSplitJobs:
    """ResourceManager.split_jobs() must satisfy outer × inner <= total."""

    def _total(self, rm: ResourceManager, per_inner_mb: int) -> int:
        return rm.get_n_jobs(per_inner_mb)

    def test_constraint_holds_large_n_outer(self):
        rm = _rm(user_max_threads=8)
        with (
            patch.object(rm, "free_ram_mb", return_value=128_000),
            patch.object(rm, "physical_cpu_count", return_value=64),
        ):
            outer, inner = rm.split_jobs(n_outer=10, per_inner_mb=300)
            total = rm.get_n_jobs(per_worker_mb=300)
        assert outer * inner <= total
        assert outer >= 1
        assert inner >= 1

    def test_constraint_holds_small_n_outer(self):
        rm = _rm(user_max_threads=8)
        with (
            patch.object(rm, "free_ram_mb", return_value=128_000),
            patch.object(rm, "physical_cpu_count", return_value=64),
        ):
            outer, inner = rm.split_jobs(n_outer=2, per_inner_mb=300)
            total = rm.get_n_jobs(per_worker_mb=300)
        assert outer <= 2
        assert inner >= 1
        assert outer * inner <= total

    def test_n_outer_equals_1(self):
        rm = _rm(user_max_threads=4)
        outer, inner = rm.split_jobs(n_outer=1, per_inner_mb=200)
        assert outer == 1
        assert inner >= 1

    def test_n_outer_exceeds_budget_clamps_outer(self):
        """When n_outer >= total, outer = total and inner = 1."""
        rm = ResourceManager(user_max_threads=2, hard_n_jobs_cap=2)
        with (
            patch.object(rm, "free_ram_mb", return_value=128_000),
            patch.object(rm, "physical_cpu_count", return_value=64),
        ):
            outer, inner = rm.split_jobs(n_outer=100, per_inner_mb=200)
        assert outer == 2
        assert inner == 1

    def test_returns_tuple_of_two_positive_ints(self):
        rm = _rm(user_max_threads=4)
        result = rm.split_jobs(n_outer=3, per_inner_mb=200)
        assert len(result) == 2
        assert all(isinstance(v, int) for v in result)
        assert all(v >= 1 for v in result)


# ---------------------------------------------------------------------------
# Singleton: identity and reset
# ---------------------------------------------------------------------------

class TestSingleton:
    """get_resource_manager() must return the same instance across calls;
    reset_resource_manager() must cause the next call to produce a fresh one."""

    def setup_method(self):
        """Start each test with a clean singleton state."""
        reset_resource_manager()

    def teardown_method(self):
        """Leave state clean for subsequent tests."""
        reset_resource_manager()

    def test_same_instance_on_repeated_calls(self):
        rm1 = get_resource_manager()
        rm2 = get_resource_manager()
        assert rm1 is rm2

    def test_reset_causes_new_instance(self):
        rm1 = get_resource_manager()
        reset_resource_manager()
        rm2 = get_resource_manager()
        assert rm1 is not rm2

    def test_reset_then_get_returns_valid_rm(self):
        reset_resource_manager()
        rm = get_resource_manager()
        assert isinstance(rm, ResourceManager)
        assert rm.get_n_jobs() >= 1


# ---------------------------------------------------------------------------
# Singleton picks up user_max_threads from CLI args mock
# ---------------------------------------------------------------------------

class TestSingletonPicksUpCliArgs:
    """After reset, the new singleton must read user_max_threads from ema.config."""

    def setup_method(self):
        reset_resource_manager()

    def teardown_method(self):
        reset_resource_manager()
        # Remove the mock config module if we injected it
        sys.modules.pop("ema.config", None)

    def test_picks_up_threads_from_ema_config(self):
        """Inject a mock ema.config.args with threads=3 and verify the RM respects it."""
        mock_args = MagicMock()
        mock_args.threads = 3

        mock_config = types.ModuleType("ema.config")
        mock_config.args = mock_args
        sys.modules["ema.config"] = mock_config

        reset_resource_manager()
        rm = get_resource_manager()

        assert rm.user_max_threads == 3
        assert rm.get_n_jobs(per_worker_mb=200) <= 3

    def test_no_ema_config_module_defaults_to_none(self):
        """If ema.config is absent, user_max_threads should be None (no ceiling)."""
        sys.modules.pop("ema.config", None)
        reset_resource_manager()
        rm = get_resource_manager()
        assert rm.user_max_threads is None

    def test_ema_config_without_threads_attr_defaults_to_none(self):
        """If ema.config.args has no 'threads' attribute, user_max_threads is None."""
        mock_args = MagicMock(spec=[])  # no attributes
        mock_config = types.ModuleType("ema.config")
        mock_config.args = mock_args
        sys.modules["ema.config"] = mock_config

        reset_resource_manager()
        rm = get_resource_manager()
        assert rm.user_max_threads is None


# ---------------------------------------------------------------------------
# ResourceManager.report() exposes user_max_threads
# ---------------------------------------------------------------------------

class TestReport:
    def test_report_includes_user_max_threads(self):
        rm = ResourceManager(user_max_threads=6)
        report = rm.report()
        assert "user_max_threads" in report
        assert report["user_max_threads"] == 6

    def test_report_user_max_threads_none_when_not_set(self):
        rm = ResourceManager(user_max_threads=None)
        report = rm.report()
        assert report["user_max_threads"] is None
