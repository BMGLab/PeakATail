import numpy as np
from typing import List, Tuple

from ema.strategies.base import PeakFinderStrategy
from ema.strategies import register
from ema.strategies.utils import reconstruct_cb_dict, partition_peak_region
from ema.statistics.background import compute_lambda
from ema.statistics.significance import poisson_pvalue


def _moving_average(data: list, window: int) -> np.ndarray:
    """Smooth data with a centered moving average.

    Uses numpy convolution in 'same' mode so the returned array has the same
    length as the input.  When the window is larger than the data, or the
    window is 1, the raw data is returned unchanged.

    Args:
        data: Raw coverage heights (one value per position).
        window: Number of positions to average across.

    Returns:
        Smoothed array with the same length as ``data``.
    """
    if window <= 1 or len(data) <= window:
        return np.array(data, dtype=float)
    kernel = np.ones(window) / window
    return np.convolve(data, kernel, mode="same")


def _find_local_maxima(smoothed: np.ndarray) -> List[int]:
    """Return indices where the gradient crosses from positive to non-positive.

    A positive-to-non-positive crossing of the first derivative corresponds to
    a local maximum in the smoothed coverage profile and is a candidate PAS
    location.

    Args:
        smoothed: Smoothed coverage array.

    Returns:
        Sorted list of array indices that are local maxima.
    """
    if len(smoothed) < 3:
        return []
    gradient = np.diff(smoothed)
    maxima: List[int] = []
    for i in range(len(gradient) - 1):
        if gradient[i] > 0 and gradient[i + 1] <= 0:
            maxima.append(i + 1)
    return maxima


def _compute_prominence(heights: list, idx: int) -> float:
    """Compute topographic prominence for a local maximum.

    Prominence is defined as the height at ``idx`` minus the higher of the
    minimum values found while walking left and right until a taller neighbour
    is encountered (or the array boundary is reached).

    Args:
        heights: Full list of raw coverage heights.
        idx: Index of the candidate local maximum.

    Returns:
        Prominence value (>= 0).
    """
    h = heights[idx]

    left_min = h
    for i in range(idx - 1, -1, -1):
        if heights[i] < left_min:
            left_min = heights[i]
        if heights[i] > h:
            break

    right_min = h
    for i in range(idx + 1, len(heights)):
        if heights[i] < right_min:
            right_min = heights[i]
        if heights[i] > h:
            break

    # The base is the higher of the two valley floors
    base = max(left_min, right_min)
    return h - base


def _second_derivative_magnitude(smoothed: np.ndarray, idx: int) -> float:
    """Return the absolute value of the second derivative at ``idx``.

    A large negative second derivative at a local maximum means the profile is
    sharply peaked there — a high-confidence PAS location.  The result is
    returned as a positive magnitude so that higher values mean more confidence.

    Args:
        smoothed: Smoothed coverage array.
        idx: Index at which to evaluate the curvature.

    Returns:
        Absolute second derivative magnitude (>= 0).
    """
    if idx == 0 or idx >= len(smoothed) - 1:
        return 0.0
    d2 = smoothed[idx - 1] - 2.0 * smoothed[idx] + smoothed[idx + 1]
    # At a maximum d2 is negative; we want a positive confidence contribution
    return float(abs(d2))


@register("lambda_gradient")
class LambdaGradientStrategy(PeakFinderStrategy):
    """Hybrid strategy: lambda validates the region, gradient finds multiple PAS.

    This is PeakATail's novel approach combining two independent signals:

    Phase 1 — Region significance (lambda gate):
        Background lambda is estimated from floor positions within the peak
        region (positions below ``floor_fraction`` of the maximum height).
        A Poisson test on the regional maximum rejects noise regions before
        any gradient work is done.

    Phase 2 — Multi-PAS localisation (gradient zero-crossings):
        The coverage profile is smoothed with a moving average, then the first
        derivative is computed.  Positions where the gradient transitions from
        positive to non-positive are local maxima and become PAS candidates.
        Candidates are filtered by:
          - minimum raw height (``min_height``)
          - height above local background (> ``lambda_local``)
          - topographic prominence (``min_prominence``)

    Confidence scoring:
        Each surviving candidate is scored as::

            confidence = prominence * (height / lambda_local)
                       + second_derivative_magnitude

        This rewards candidates that are both prominent and sharply peaked.
        The top ``max_pas`` candidates are returned, sorted by genomic position.
    """

    def __init__(
        self,
        smoothing_window: int = 50,
        min_prominence: float = 5.0,
        lambda_method: str = "median",
        max_pas: int = 5,
        min_height: int = 20,
        alpha: float = 0.05,
        floor_fraction: float = 0.1,
    ) -> None:
        """Initialise the strategy with tunable hyper-parameters.

        Args:
            smoothing_window: Width of the moving-average kernel in base pairs.
                Larger values suppress noise but may merge closely spaced PAS.
            min_prominence: Minimum topographic prominence required for a
                candidate to pass the filter (in coverage units).
            lambda_method: Method forwarded to ``compute_lambda``.
                One of ``"median"``, ``"percentile20"``, or ``"global"``.
            max_pas: Maximum number of PAS to return per peak region.
            min_height: Absolute minimum raw coverage height for a candidate.
            alpha: Poisson significance threshold for the region gate.
                Regions with p-value >= alpha are discarded entirely.
            floor_fraction: Fraction of the regional maximum below which a
                position is considered background (forwarded to
                ``compute_lambda``).
        """
        self.smoothing_window = smoothing_window
        self.min_prominence = min_prominence
        self.lambda_method = lambda_method
        self.max_pas = max_pas
        self.min_height = min_height
        self.alpha = alpha
        self.floor_fraction = floor_fraction

    # ------------------------------------------------------------------
    # PeakFinderStrategy interface
    # ------------------------------------------------------------------

    def find_pas(self, peak) -> List[Tuple[int, int]]:
        """Find PAS positions within a peak region.

        Args:
            peak: Peak object exposing ``peak_list`` (list of [position,
                height]) and ``peak_strand`` (bool).

        Returns:
            List of ``(pas_start, pas_end)`` tuples, one per detected PAS,
            sorted by ascending genomic position.  Returns an empty list when
            the region fails the significance gate or no candidates survive
            filtering.
        """
        if not peak.peak_list:
            return []

        positions = [p for p, _ in peak.peak_list]
        heights = [h for _, h in peak.peak_list]

        if not heights:
            return []

        max_h = max(heights)
        if max_h < self.min_height:
            return []

        # ------------------------------------------------------------------
        # Phase 1: region significance gate
        # ------------------------------------------------------------------
        lambda_local = compute_lambda(
            heights,
            method=self.lambda_method,
            floor_fraction=self.floor_fraction,
        )

        p_val = poisson_pvalue(max_h, lambda_local)
        if p_val >= self.alpha:
            return []

        # ------------------------------------------------------------------
        # Phase 2: gradient zero-crossing localisation
        # ------------------------------------------------------------------
        smoothed = _moving_average(heights, self.smoothing_window)
        maxima_indices = _find_local_maxima(smoothed)

        if not maxima_indices:
            # Fallback: global maximum with full peak region
            return [(min(positions), max(positions))]

        # ------------------------------------------------------------------
        # Candidate filtering and scoring
        # ------------------------------------------------------------------
        candidates: List[Tuple[int, int, float]] = []  # (array_idx, genomic_pos, score)

        for idx in maxima_indices:
            if idx >= len(heights):
                continue

            h = heights[idx]

            if h < self.min_height:
                continue

            if h <= lambda_local:
                # Must exceed local background, not merely equal it
                continue

            # Per-PAS Poisson significance test
            pas_pval = poisson_pvalue(h, lambda_local)
            if pas_pval >= self.alpha:
                continue

            prominence = _compute_prominence(heights, idx)
            if prominence < self.min_prominence:
                continue

            # Composite confidence: prominence * enrichment + curvature
            enrichment = h / max(lambda_local, 1.0)
            curvature = _second_derivative_magnitude(smoothed, idx)
            confidence = prominence * enrichment + curvature

            candidates.append((idx, positions[idx], confidence))

        if not candidates:
            # Fallback: return the global maximum with full peak region
            return [(min(positions), max(positions))]

        # Top N by confidence, then re-sort by genomic position for output
        candidates.sort(key=lambda x: x[2], reverse=True)
        top = candidates[: self.max_pas]
        top.sort(key=lambda x: x[1])

        # Partition the peak region among found PAS by midpoint splitting
        summits = [pos for _, pos, _ in top]
        regions = partition_peak_region(summits, positions)
        return regions

    def get_cb_dict_for_pas(self, peak, pas_1: int, pas_2: int) -> dict:
        """Return per-barcode counts from reads near a specific PAS.

        Delegates to ``reconstruct_cb_dict`` which applies a window around
        ``[pas_1, pas_2]`` to pull only locally supporting reads.

        Args:
            peak: Peak object with ``cb_positions``.
            pas_1: PAS start position.
            pas_2: PAS end position.

        Returns:
            Dict of ``{barcode: count}`` for reads in the PAS neighbourhood.
        """
        return reconstruct_cb_dict(peak.cb_positions, pas_1, pas_2)

    def valley_threshold(self, peak, user_min_pas_prominence: float) -> float:
        """Return ``compute_lambda(heights)`` — the strategy's own background.

        Overrides the static default in :class:`PeakFinderStrategy` so the
        post-detection merger uses the same dynamic background estimate this
        strategy uses for the region significance gate and per-PAS Poisson
        tests.  The user-supplied ``min_pas_prominence`` is ignored for
        lambda-based strategies — the valley test becomes "is the dip at or
        below the local lambda?".
        """
        if not peak.peak_list:
            return float(user_min_pas_prominence)
        heights = [h for _, h in peak.peak_list]
        return float(compute_lambda(
            heights,
            method=self.lambda_method,
            floor_fraction=self.floor_fraction,
        ))

    def get_params(self) -> dict:
        """Return the current parameter set for logging and benchmarking.

        Returns:
            Dict with strategy name and all hyper-parameter values.
        """
        return {
            "strategy": "lambda_gradient",
            "smoothing_window": self.smoothing_window,
            "min_prominence": self.min_prominence,
            "lambda_method": self.lambda_method,
            "max_pas": self.max_pas,
            "min_height": self.min_height,
            "alpha": self.alpha,
            "floor_fraction": self.floor_fraction,
        }
