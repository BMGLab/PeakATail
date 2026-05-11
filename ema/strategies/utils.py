from typing import List, Tuple, Dict


def partition_peak_region(summits: List[int], peak_positions: List[int]) -> List[Tuple[int, int]]:
    """Partition a peak region among multiple PAS by midpoint splitting.

    Every read position in the peak is assigned to exactly one PAS.
    Adjacent PAS are divided at the midpoint between their summits.
    The leftmost PAS owns everything from peak start to the first midpoint.
    The rightmost PAS owns everything from the last midpoint to peak end.

    Args:
        summits: Sorted list of PAS summit positions (genomic coordinates).
        peak_positions: All positions in the peak region (from peak_list).

    Returns:
        List of (region_start, region_end) tuples, one per summit,
        covering the entire peak region with no gaps.
    """
    if not summits or not peak_positions:
        return []

    peak_start = min(peak_positions)
    peak_end = max(peak_positions)

    if len(summits) == 1:
        return [(peak_start, peak_end)]

    # Compute midpoints between adjacent summits
    midpoints = []
    for i in range(len(summits) - 1):
        mid = (summits[i] + summits[i + 1]) // 2
        midpoints.append(mid)

    # Build regions
    regions = []
    for i, summit in enumerate(summits):
        if i == 0:
            # Leftmost PAS: peak_start to first midpoint
            region_start = peak_start
            region_end = midpoints[0]
        elif i == len(summits) - 1:
            # Rightmost PAS: last midpoint to peak_end
            region_start = midpoints[-1] + 1
            region_end = peak_end
        else:
            # Middle PAS: between two midpoints
            region_start = midpoints[i - 1] + 1
            region_end = midpoints[i]

        regions.append((region_start, region_end))

    return regions


def reconstruct_cb_dict(cb_positions: dict, region_start: int, region_end: int) -> dict:
    """Extract cell barcode counts for reads within a specific PAS region.

    Uses the actual region boundaries (from partition_peak_region) instead
    of a hardcoded window. Every read is assigned to exactly one PAS.

    Args:
        cb_positions: {position: {barcode: count}} from Peak object.
        region_start: Start of this PAS's assigned region.
        region_end: End of this PAS's assigned region.

    Returns:
        Dict of {barcode: count} for reads in this region only.
    """
    result = {}
    for pos, cb_counts in cb_positions.items():
        if region_start <= pos <= region_end:
            for cb, count in cb_counts.items():
                result[cb] = result.get(cb, 0) + count
    return result
