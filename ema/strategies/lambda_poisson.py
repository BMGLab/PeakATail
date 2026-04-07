from typing import List, Tuple

from ema.strategies import register
from ema.strategies.base import PeakFinderStrategy
from ema.strategies.utils import reconstruct_cb_dict
from ema.statistics.background import compute_lambda
from ema.statistics.significance import poisson_pvalue


@register("lambda_poisson")
class LambdaPoissonStrategy(PeakFinderStrategy):
    """MACS2-style: Poisson p-value with local lambda from floor positions.

    Finds the single most significant peak in a region using:
    1. Compute local lambda from floor positions (< 10% of max height)
    2. Test summit height against Poisson(lambda)
    3. Return summit as PAS if significant (p < alpha)

    Args:
        lambda_method: Background estimation method passed to compute_lambda.
            One of "median", "percentile20", or "global". Default "median".
        min_height: Minimum absolute read depth required at the summit before
            any statistical test is applied. Default 20.
        alpha: Significance threshold for the one-sided Poisson p-value.
            Peaks with p >= alpha are discarded. Default 0.05.
        floor_fraction: Fraction of max height below which a position is
            treated as background (floor). Passed to compute_lambda.
            Default 0.1 (i.e. positions at < 10 % of max height).
    """

    def __init__(
        self,
        lambda_method: str = "median",
        min_height: int = 20,
        alpha: float = 0.05,
        floor_fraction: float = 0.1,
    ) -> None:
        self.lambda_method = lambda_method
        self.min_height = min_height
        self.alpha = alpha
        self.floor_fraction = floor_fraction

    def find_pas(self, peak) -> List[Tuple[int, int]]:
        """Find the single most significant PAS within a peak region.

        Algorithm:
            1. Guard against empty input and insufficient depth.
            2. Estimate local background lambda from floor positions.
            3. Locate the summit (position of maximum height).
            4. Test summit height against Poisson(lambda); discard if not
               significant.
            5. Walk outward from the summit to determine peak boundaries
               (stop where height drops to or below lambda).
            6. Return the boundary span as a single (pas_start, pas_end) tuple.

        Args:
            peak: Peak object with attributes:
                - peak_list: list of [position, height] pairs
                - cb_dict: {barcode: count} total per-barcode counts
                - cb_positions: {position: {barcode: count}} per-position counts
                - peak_strand: bool, strand direction

        Returns:
            A list containing at most one (pas_start, pas_end) tuple, or an
            empty list when no significant PAS is detected.
        """
        if not peak.peak_list:
            return []

        positions: List[int] = [p for p, _ in peak.peak_list]
        heights: List[int] = [h for _, h in peak.peak_list]

        max_h = max(heights)
        if max_h < self.min_height:
            return []

        # Estimate local background from floor positions.
        lambda_local = compute_lambda(
            heights,
            method=self.lambda_method,
            floor_fraction=self.floor_fraction,
        )

        # Locate summit.
        max_idx = heights.index(max_h)
        summit_pos = positions[max_idx]  # noqa: F841 — kept for readability

        # Significance test: reject if p-value is at or above alpha.
        p_val = poisson_pvalue(max_h, lambda_local)
        if p_val >= self.alpha:
            return []

        # Walk left from summit until height drops to lambda or the edge.
        left = max_idx
        while left > 0 and heights[left - 1] > lambda_local:
            left -= 1

        # Walk right from summit until height drops to lambda or the edge.
        right = max_idx
        while right < len(heights) - 1 and heights[right + 1] > lambda_local:
            right += 1

        pas_start = positions[left]
        pas_end = positions[right]

        return [(pas_start, pas_end)]

    def get_cb_dict_for_pas(self, peak, pas_1: int, pas_2: int) -> dict:
        """Return per-barcode counts for reads within the PAS region.

        Delegates to reconstruct_cb_dict, which aggregates counts from
        peak.cb_positions for positions falling within the PAS window.

        Args:
            peak: Peak object with a cb_positions attribute.
            pas_1: PAS start position (genomic coordinate).
            pas_2: PAS end position (genomic coordinate).

        Returns:
            Dict of {barcode: count} for this PAS.
        """
        return reconstruct_cb_dict(peak.cb_positions, pas_1, pas_2)

    def get_params(self) -> dict:
        """Return current strategy parameters for logging and benchmarking.

        Returns:
            Dict with strategy name and all constructor parameters.
        """
        return {
            "strategy": "lambda_poisson",
            "lambda_method": self.lambda_method,
            "min_height": self.min_height,
            "alpha": self.alpha,
            "floor_fraction": self.floor_fraction,
        }
