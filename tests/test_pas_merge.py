"""Unit tests for the strategy-agnostic post-detection PAS merger.

Covers both tiers and the strategy-driven valley threshold:
    Tier 1: ``min_pas_spacing`` distance floor
    Tier 2: valley vs ``strategy.valley_threshold(peak, min_pas_prominence)``
"""
from __future__ import annotations

import pytest

from ema.strategies.utils import merge_close_or_low_prominence


class _StubPeak:
    """Minimal Peak stand-in exposing only what the merger needs."""

    def __init__(self, peak_list):
        self.peak_list = peak_list


class _StaticStrategy:
    """Non-lambda strategy: valley threshold == user setting."""

    def valley_threshold(self, peak, user_min_pas_prominence):
        return float(user_min_pas_prominence)


class _LambdaLikeStrategy:
    """Lambda strategy: valley threshold == fixed lambda value (here: 10)."""

    def __init__(self, lambda_value: float = 10.0):
        self.lambda_value = lambda_value

    def valley_threshold(self, peak, user_min_pas_prominence):
        return float(self.lambda_value)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------


def test_empty_input_passes_through():
    peak = _StubPeak([])
    assert merge_close_or_low_prominence([], peak, _StaticStrategy(), 100, 5.0) == []


def test_single_region_passes_through():
    peak = _StubPeak([(100, 50)])
    assert merge_close_or_low_prominence(
        [(50, 150)], peak, _StaticStrategy(), 100, 5.0
    ) == [(50, 150)]


def test_both_tiers_disabled_returns_input_unchanged():
    """min_pas_spacing=0 and min_pas_prominence<0 → full no-op (baseline)."""
    peak = _StubPeak([(p, 30) for p in range(0, 200)])
    regions = [(0, 10), (11, 20), (50, 60)]
    out = merge_close_or_low_prominence(regions, peak, _StaticStrategy(), 0, -1.0)
    assert out == regions


# ---------------------------------------------------------------------------
# Tier 1: distance floor
# ---------------------------------------------------------------------------


def test_tier1_merges_close_pair_unconditionally():
    """1-bp gap < min_pas_spacing=50 → merge even if valley would not trigger Tier 2."""
    # Heights are all 100, so any valley would be 100 — well above threshold.
    peak = _StubPeak([(p, 100) for p in range(0, 200)])
    regions = [(50, 100), (101, 150)]
    out = merge_close_or_low_prominence(
        regions, peak, _StaticStrategy(), min_pas_spacing=50, min_pas_prominence=5.0,
    )
    assert out == [(50, 150)]


def test_tier1_keeps_far_pair_when_tier2_passes():
    """Gap of 100 bp ≥ spacing=50; valley is high → Tier 2 passes → keep both."""
    peak = _StubPeak([(p, 100) for p in range(0, 300)])
    regions = [(50, 100), (200, 250)]
    out = merge_close_or_low_prominence(
        regions, peak, _StaticStrategy(), min_pas_spacing=50, min_pas_prominence=5.0,
    )
    assert out == regions


def test_tier1_disabled_when_spacing_zero():
    """spacing=0 → Tier 1 inactive, Tier 2 still applies."""
    peak = _StubPeak([(p, 100) for p in range(0, 200)])
    regions = [(50, 100), (101, 150)]
    out = merge_close_or_low_prominence(
        regions, peak, _StaticStrategy(), min_pas_spacing=0, min_pas_prominence=5.0,
    )
    # Tier 1 off; valley is empty (no positions strictly between 100 and 101)
    # → _min_h_between returns 0 → 0 < 5 → merge by Tier 2.
    assert out == [(50, 150)]


# ---------------------------------------------------------------------------
# Tier 2: static threshold (non-lambda strategy)
# ---------------------------------------------------------------------------


def test_tier2_static_merges_shallow_valley():
    """Valley = 3, threshold = 5 → merge."""
    heights = [(p, 100) for p in range(0, 60)] \
              + [(p, 3) for p in range(60, 70)] \
              + [(p, 100) for p in range(70, 200)]
    peak = _StubPeak(heights)
    regions = [(0, 60), (70, 150)]
    out = merge_close_or_low_prominence(
        regions, peak, _StaticStrategy(), min_pas_spacing=0, min_pas_prominence=5.0,
    )
    assert out == [(0, 150)]


def test_tier2_static_keeps_deep_valley():
    """Valley = 20, threshold = 5 → keep separate."""
    heights = [(p, 100) for p in range(0, 60)] \
              + [(p, 20) for p in range(60, 70)] \
              + [(p, 100) for p in range(70, 200)]
    peak = _StubPeak(heights)
    regions = [(0, 60), (70, 150)]
    out = merge_close_or_low_prominence(
        regions, peak, _StaticStrategy(), min_pas_spacing=0, min_pas_prominence=5.0,
    )
    assert out == regions


def test_tier2_disabled_when_prominence_negative():
    """min_pas_prominence=-1 → Tier 2 off."""
    heights = [(p, 100) for p in range(0, 200)]
    peak = _StubPeak(heights)
    regions = [(0, 60), (70, 150)]
    out = merge_close_or_low_prominence(
        regions, peak, _StaticStrategy(), min_pas_spacing=200, min_pas_prominence=-1.0,
    )
    # Tier 1 fires because gap 10 < 200 → merge.
    assert out == [(0, 150)]


# ---------------------------------------------------------------------------
# Tier 2: dynamic threshold (lambda-like strategy)
# ---------------------------------------------------------------------------


def test_tier2_lambda_uses_strategy_threshold_not_user_setting():
    """Lambda strategy returns 10; valley = 8 → merge (8 < 10) regardless of user value."""
    heights = [(p, 100) for p in range(0, 60)] \
              + [(p, 8) for p in range(60, 70)] \
              + [(p, 100) for p in range(70, 200)]
    peak = _StubPeak(heights)
    regions = [(0, 60), (70, 150)]
    # User sets min_pas_prominence=5 (would not merge with static), but the
    # lambda strategy overrides to 10 → valley 8 < 10 → merge.
    out = merge_close_or_low_prominence(
        regions, peak, _LambdaLikeStrategy(10.0),
        min_pas_spacing=0, min_pas_prominence=5.0,
    )
    assert out == [(0, 150)]


def test_tier2_lambda_keeps_valley_above_lambda():
    """Lambda = 10; valley = 50 → keep."""
    heights = [(p, 100) for p in range(0, 60)] \
              + [(p, 50) for p in range(60, 70)] \
              + [(p, 100) for p in range(70, 200)]
    peak = _StubPeak(heights)
    regions = [(0, 60), (70, 150)]
    out = merge_close_or_low_prominence(
        regions, peak, _LambdaLikeStrategy(10.0),
        min_pas_spacing=0, min_pas_prominence=999.0,  # would force merge if static
    )
    assert out == regions


# ---------------------------------------------------------------------------
# Chains
# ---------------------------------------------------------------------------


def test_three_way_chain_collapses_to_one():
    """A → B by Tier 1; (A∪B) → C by Tier 1 → single merged region."""
    peak = _StubPeak([(p, 100) for p in range(0, 200)])
    regions = [(0, 30), (31, 60), (61, 90)]
    out = merge_close_or_low_prominence(
        regions, peak, _StaticStrategy(),
        min_pas_spacing=50, min_pas_prominence=-1.0,  # Tier 2 disabled
    )
    assert out == [(0, 90)]


def test_unsorted_input_is_sorted_before_merging():
    """Input regions in reverse order still merge correctly."""
    peak = _StubPeak([(p, 100) for p in range(0, 200)])
    regions = [(100, 150), (0, 50)]
    out = merge_close_or_low_prominence(
        regions, peak, _StaticStrategy(),
        min_pas_spacing=100, min_pas_prominence=-1.0,
    )
    # After sort: [(0,50), (100,150)] — gap 50 < 100 → merge.
    assert out == [(0, 150)]
