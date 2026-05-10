"""Tests for ResourceManager.get_tile_size() (Phase 4).

Coverage:
- get_tile_size returns value in [min_tile_bp, max_tile_bp]
- get_tile_size falls back to min_tile_bp on BAM open failure
- get_tile_size falls back to min_tile_bp for a BAM with no references
- density=0 (empty region) returns max_tile_bp
- high density produces smaller tiles than low density
- more workers means smaller tile (budget is divided by n_workers)
- per-worker budget cap from free_ram is respected
- real BAM (data/chr22.bam) returns a sane value when available
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator
from unittest.mock import MagicMock, patch

import pytest

from ema.utils.resource_manager import ResourceManager


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _rm(**kwargs) -> ResourceManager:
    """Create a ResourceManager with predictable RAM/CPU for testing."""
    return ResourceManager(
        user_max_threads=kwargs.pop("user_max_threads", None),
        ram_safety_fraction=kwargs.pop("ram_safety_fraction", 0.7),
        **kwargs,
    )


@contextmanager
def _mock_bam(n_references: int = 2, reads_in_1mb: int = 1000) -> Iterator[None]:
    """Context manager that patches pysam.AlignmentFile for get_tile_size tests."""
    mock_read = MagicMock()
    mock_bam = MagicMock()
    mock_bam.__enter__ = lambda s: s
    mock_bam.__exit__ = MagicMock(return_value=False)
    mock_bam.nreferences = n_references
    mock_bam.get_reference_name = MagicMock(return_value="chr1")
    mock_bam.get_reference_length = MagicMock(return_value=50_000_000)
    mock_bam.fetch = MagicMock(return_value=iter([mock_read] * reads_in_1mb))

    with patch("pysam.AlignmentFile", return_value=mock_bam):
        yield


# ---------------------------------------------------------------------------
# Basic contract
# ---------------------------------------------------------------------------


class TestGetTileSizeBasicContract:
    """Return value must always be in [min_tile_bp, max_tile_bp]."""

    def test_return_in_bounds_typical_density(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=8192):
            with _mock_bam(reads_in_1mb=1000):
                result = rm.get_tile_size("/fake/a.bam", n_workers=4)
        assert 5_000_000 <= result <= 100_000_000

    def test_return_in_bounds_high_density(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=8192):
            with _mock_bam(reads_in_1mb=1_000_000):
                result = rm.get_tile_size("/fake/a.bam", n_workers=4)
        assert 5_000_000 <= result <= 100_000_000

    def test_return_in_bounds_low_density(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=8192):
            with _mock_bam(reads_in_1mb=1):
                result = rm.get_tile_size("/fake/a.bam", n_workers=4)
        assert 5_000_000 <= result <= 100_000_000

    def test_returns_int(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=8192):
            with _mock_bam(reads_in_1mb=500):
                result = rm.get_tile_size("/fake/a.bam", n_workers=2)
        assert isinstance(result, int)


# ---------------------------------------------------------------------------
# Fallback behaviour on errors
# ---------------------------------------------------------------------------


class TestGetTileSizeFallbacks:
    """get_tile_size must degrade gracefully and return min_tile_bp on failure."""

    def test_bam_open_failure_returns_min(self):
        rm = _rm()
        with patch("pysam.AlignmentFile", side_effect=OSError("file not found")):
            result = rm.get_tile_size("/nonexistent.bam", n_workers=2)
        assert result == 5_000_000

    def test_no_references_returns_min(self):
        rm = _rm()
        with _mock_bam(n_references=0):
            result = rm.get_tile_size("/fake/empty.bam", n_workers=2)
        assert result == 5_000_000

    def test_zero_density_returns_max(self):
        """fetch returns an empty iterator → density=0 → tile_size = max_tile_bp."""
        rm = _rm()
        with _mock_bam(reads_in_1mb=0):
            with patch.object(rm, "free_ram_mb", return_value=8192):
                result = rm.get_tile_size("/fake/a.bam", n_workers=1)
        assert result == 100_000_000  # clamped to max


# ---------------------------------------------------------------------------
# Density sensitivity
# ---------------------------------------------------------------------------


class TestGetTileSizeDensitySensitivity:
    """Higher read density must yield smaller tiles (more reads per bp → less bp per worker)."""

    def test_high_density_smaller_than_low_density(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=16384):
            with _mock_bam(reads_in_1mb=100_000):
                high = rm.get_tile_size("/fake/a.bam", n_workers=4)
            with _mock_bam(reads_in_1mb=1_000):
                low = rm.get_tile_size("/fake/a.bam", n_workers=4)
        assert high <= low


# ---------------------------------------------------------------------------
# Worker count effect
# ---------------------------------------------------------------------------


class TestGetTileSizeWorkerEffect:
    """More workers → smaller per-worker budget → smaller tile (or clamped to min)."""

    def test_more_workers_smaller_or_equal_tile(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=8192):
            with _mock_bam(reads_in_1mb=500):
                few = rm.get_tile_size("/fake/a.bam", n_workers=1)
                many = rm.get_tile_size("/fake/a.bam", n_workers=16)
        # With more workers the per-worker budget shrinks → tile ≤ single-worker tile
        assert many <= few


# ---------------------------------------------------------------------------
# Budget cap from free RAM
# ---------------------------------------------------------------------------


class TestGetTileSizeRamBudget:
    """Very low free RAM should produce tiles at or near min_tile_bp."""

    def test_very_low_ram_very_high_density_produces_min_tile(self):
        # 100 MB free × 0.7 / 4 = 17.5 MB budget; overhead 92 MB → usable = 1 MB
        # density = 1_000_000 reads / 1 MB = 1 read/bp
        # tile_size = 1*1024*1024 / 50 / 1.0 = 20,971 bp → clamps to min_tile_bp
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=100):
            with _mock_bam(reads_in_1mb=1_000_000):
                result = rm.get_tile_size("/fake/a.bam", n_workers=4)
        assert result == 5_000_000  # clamped to min

    def test_high_ram_allows_larger_tile(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=128_000):
            with _mock_bam(reads_in_1mb=10):
                result = rm.get_tile_size("/fake/a.bam", n_workers=1)
        # Large RAM + very low density → result should reach max
        assert result == 100_000_000


# ---------------------------------------------------------------------------
# Custom clamp bounds
# ---------------------------------------------------------------------------


class TestGetTileSizeCustomBounds:
    def test_custom_min_respected(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=200):
            with _mock_bam(reads_in_1mb=500):
                result = rm.get_tile_size(
                    "/fake/a.bam",
                    n_workers=4,
                    min_tile_bp=10_000_000,
                    max_tile_bp=100_000_000,
                )
        assert result >= 10_000_000

    def test_custom_max_respected(self):
        rm = _rm()
        with patch.object(rm, "free_ram_mb", return_value=128_000):
            with _mock_bam(reads_in_1mb=1):
                result = rm.get_tile_size(
                    "/fake/a.bam",
                    n_workers=1,
                    min_tile_bp=5_000_000,
                    max_tile_bp=20_000_000,
                )
        assert result <= 20_000_000


# ---------------------------------------------------------------------------
# Integration: real BAM (skip if not present)
# ---------------------------------------------------------------------------


class TestGetTileSizeRealBam:
    """Smoke test against an actual indexed BAM in the repo data/ directory."""

    @pytest.fixture
    def chr22_bam(self) -> str:
        """Locate any indexed BAM in the repo for tile-sizing smoke tests.

        Tries (in order): the small chr22 test BAM under ``test_run/``,
        the legacy hardcoded ``data/chr22.bam`` path, and finally any
        ``data/*.bam`` with a ``.bai`` sibling. Skips if nothing found.
        """
        import os, glob
        candidates = [
            "/home/user/PeakATail/test_run/chr22.bam",
            "/home/user/PeakATail/data/chr22.bam",
        ]
        candidates += sorted(glob.glob("/home/user/PeakATail/data/*.bam"))
        for bam in candidates:
            if os.path.exists(bam) and os.path.exists(bam + ".bai"):
                return bam
        pytest.skip("no indexed BAM found in test_run/ or data/")

    def test_real_bam_returns_valid_tile_size(self, chr22_bam: str):
        rm = _rm()
        result = rm.get_tile_size(chr22_bam, n_workers=4)
        assert 5_000_000 <= result <= 100_000_000
        assert isinstance(result, int)

    def test_real_bam_single_worker_larger_than_many_workers(self, chr22_bam: str):
        rm = _rm()
        single = rm.get_tile_size(chr22_bam, n_workers=1)
        multi = rm.get_tile_size(chr22_bam, n_workers=8)
        # single-worker has bigger budget → tile ≥ multi-worker tile (or both clamped)
        assert single >= multi
