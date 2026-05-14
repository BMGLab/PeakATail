from abc import ABC, abstractmethod
from typing import List, Tuple


class PeakFinderStrategy(ABC):
    """Base class for PAS finding within a detected peak region.

    Each strategy receives a Peak object (with peak_list, cb_dict, cb_positions)
    and returns a list of (pas_start, pas_end) tuples — one per detected PAS.

    The 'original' strategy wraps peak.pasfind() and returns at most 1 PAS.
    New strategies can return multiple PAS per peak region.
    """

    @abstractmethod
    def find_pas(self, peak) -> List[Tuple[int, int]]:
        """Find PAS positions within a peak region.

        Args:
            peak: Peak object with peak_list, cb_dict, cb_positions, peak_strand

        Returns:
            List of (pas_start, pas_end) tuples. Empty list = no valid PAS found.
        """
        pass

    @abstractmethod
    def get_cb_dict_for_pas(self, peak, pas_1: int, pas_2: int) -> dict:
        """Get the cell barcode count dict for a specific PAS.

        For original strategy: returns peak.cb_dict (full dict).
        For new strategies: reconstructs from cb_positions for the PAS region only.

        Args:
            peak: Peak object
            pas_1, pas_2: PAS start and end positions

        Returns:
            Dict of {barcode: count} for this PAS.
        """
        pass

    @abstractmethod
    def get_params(self) -> dict:
        """Return current strategy parameters for logging/benchmarking."""
        pass

    def valley_threshold(self, peak, user_min_pas_prominence: float) -> float:
        """Return the valley height below which two adjacent PAS should be merged.

        Used by the post-detection PAS merger (Tier 2) to decide whether the
        dip between two emitted summits is "background" — if the valley is
        below this threshold the pair is merged.

        Default implementation returns ``user_min_pas_prominence`` unchanged
        — a static threshold supplied by the user.  Lambda-based strategies
        (``lambda_poisson``, ``lambda_gradient``) override this to return
        their own ``compute_lambda(heights)`` so the merger uses the same
        dynamic background estimate the strategy itself uses for peak
        significance.

        Args:
            peak: ``Peak`` instance providing ``peak_list``.
            user_min_pas_prominence: Static fallback supplied from CLI/YAML.
                Returned as-is by this base implementation.

        Returns:
            Valley height threshold in coverage units.
        """
        return float(user_min_pas_prominence)
