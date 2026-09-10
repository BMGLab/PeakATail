import numpy as np
from typing import List, Tuple
from ema.strategies.base import PeakFinderStrategy
from ema.strategies import register
from ema.strategies.utils import reconstruct_cb_dict, partition_peak_region


@register("sierra_iterative")
class SierraIterativeStrategy(PeakFinderStrategy):
    """Sierra-style iterative peak subtraction.

    Based on Sierra (Patrick et al., Genome Biology 2020):
    1. Find the dominant peak (maximum coverage)
    2. Record it as a PAS
    3. Zero out coverage within ±subtraction_window
    4. Repeat for next-highest peak
    5. Stop when remaining max < min_height or max_pas reached

    This naturally handles multiple PAS per UTR — once the dominant
    peak is subtracted, smaller peaks become visible.
    """

    def __init__(
        self,
        subtraction_window: int = 300,
        min_height: int = 20,
        max_pas: int = 5,
        min_fraction: float = 0.1,
    ) -> None:
        """Initialise the Sierra iterative strategy.

        Args:
            subtraction_window: Nucleotides zeroed out on each side of a found
                summit (±bp).  Sierra default is 300.
            min_height: Absolute minimum coverage for a position to be called a
                PAS.  Iterations stop when the working maximum drops below this.
            max_pas: Hard cap on the number of PAS returned per peak region.
            min_fraction: Relative floor — stop if the current maximum is below
                this fraction of the original maximum coverage.  Prevents
                calling tiny residual bumps after the dominant peak is removed.
        """
        self.subtraction_window = subtraction_window
        self.min_height = min_height
        self.max_pas = max_pas
        self.min_fraction = min_fraction

    def find_pas(self, peak) -> List[Tuple[int, int]]:
        """Find PAS positions using iterative peak subtraction.

        Args:
            peak: Peak object with peak_list (list of [position, height]),
                cb_positions, and peak_strand.

        Returns:
            List of (pas_start, pas_end) tuples sorted by genomic position.
            Returns an empty list when no valid PAS are found.
        """
        if not peak.peak_list:
            return []

        positions = [p for p, _ in peak.peak_list]
        heights = [h for _, h in peak.peak_list]

        if not heights:
            return []

        original_max = max(heights)
        if original_max < self.min_height:
            return []

        # Work on a mutable copy — never modify the original peak_list data.
        working: List[int] = list(heights)
        summits: List[int] = []

        for _ in range(self.max_pas):
            max_val = max(working)

            # Absolute floor check.
            if max_val < self.min_height:
                break

            # Relative floor check — remaining peaks are negligible.
            if max_val < original_max * self.min_fraction:
                break

            max_idx = working.index(max_val)
            summit_pos = positions[max_idx]
            summits.append(summit_pos)

            # Zero out all positions within ±subtraction_window of the summit.
            for i, pos in enumerate(positions):
                if abs(pos - summit_pos) <= self.subtraction_window:
                    working[i] = 0

        if not summits:
            return []

        # Partition the peak region among found PAS by midpoint splitting.
        # Every read is assigned to exactly one PAS — no reads lost.
        summits.sort()
        regions = partition_peak_region(summits, positions)
        return regions

    def get_cb_dict_for_pas(self, peak, pas_1: int, pas_2: int) -> dict:
        """Return the cell-barcode count dict for a specific PAS.

        Delegates to reconstruct_cb_dict which aggregates read counts from
        cb_positions within a window around the PAS coordinates.

        Args:
            peak: Peak object with cb_positions attribute.
            pas_1: PAS start position.
            pas_2: PAS end position.

        Returns:
            Dict of {barcode: count} for reads near this PAS.
        """
        return reconstruct_cb_dict(peak.cb_positions, pas_1, pas_2)

    def get_params(self) -> dict:
        """Return current strategy parameters for logging and benchmarking.

        Returns:
            Dict with strategy name and all tunable parameters.
        """
        return {
            "strategy": "sierra_iterative",
            "subtraction_window": self.subtraction_window,
            "min_height": self.min_height,
            "max_pas": self.max_pas,
            "min_fraction": self.min_fraction,
        }
