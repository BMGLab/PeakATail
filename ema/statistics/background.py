import numpy as np
from typing import List, Optional

def compute_lambda(heights: List[float], method: str = "median",
                   floor_fraction: float = 0.1) -> float:
    """Compute local background lambda from coverage heights.

    Lambda is estimated from FLOOR positions only — positions where coverage
    is below a fraction of the maximum. This prevents peaks from inflating
    the background estimate.

    Args:
        heights: List of coverage heights at each position in the peak region
        method: "median", "percentile20", or "global"
        floor_fraction: Fraction of max_height below which positions are considered floor

    Returns:
        Estimated background lambda (minimum 1.0 to avoid division by zero)
    """
    if not heights:
        return 1.0

    max_h = max(heights)
    if max_h == 0:
        return 1.0

    # Extract floor positions (below floor_fraction of max)
    floor_heights = [h for h in heights if h < max_h * floor_fraction]

    # If no floor positions (entire region is peak), use all heights
    if not floor_heights:
        floor_heights = heights

    if method == "median":
        result = float(np.median(floor_heights))
    elif method == "percentile20":
        result = float(np.percentile(floor_heights, 20))
    elif method == "global":
        result = float(np.mean(heights))
    else:
        raise ValueError(f"Unknown lambda method: {method}. Use 'median', 'percentile20', or 'global'")

    return max(result, 1.0)  # minimum 1.0


def compute_lambda_flanking(heights: List[float], peak_start_idx: int,
                            peak_end_idx: int, window: int = 100) -> float:
    """Compute lambda from flanking regions only (excluding the peak itself).

    Args:
        heights: Full coverage heights array
        peak_start_idx: Index where the peak starts
        peak_end_idx: Index where the peak ends
        window: Number of positions to use on each flank

    Returns:
        Estimated background lambda from flanking regions
    """
    left_flank = heights[max(0, peak_start_idx - window):peak_start_idx]
    right_flank = heights[peak_end_idx:min(len(heights), peak_end_idx + window)]

    flank = left_flank + right_flank
    if not flank:
        return 1.0

    return max(float(np.median(flank)), 1.0)
