from typing import List, Tuple, Dict


def merge_close_or_low_prominence(
    regions: List[Tuple[int, int]],
    peak,
    strategy,
    min_pas_spacing: int,
    min_pas_prominence: float,
) -> List[Tuple[int, int]]:
    """Merge adjacent PAS regions whose summits are too close or share a shallow valley.

    Strategy-agnostic post-detection merger.  Applied to ``find_pas`` output at
    the caller layer (``peackcalling.py``) so every strategy benefits, including
    future strategies whose contract is ``find_pas -> [(start, end), ...]``.

    Two tiers, applied in order to each adjacent pair:

    Tier 1 — Distance gate (hard floor):
        If ``next.start - prev.end < min_pas_spacing``, the pair is merged
        unconditionally.  Two summits within ~one read length are
        biologically indistinguishable for short-read scRNA — one read can
        straddle both, so they cannot be resolved as distinct PAS.

    Tier 2 — Valley gate (per-strategy threshold):
        If the gap passes Tier 1, the valley between the two regions is
        compared against a threshold supplied by the strategy itself via
        :meth:`~ema.strategies.base.PeakFinderStrategy.valley_threshold`:

            - **Lambda-based strategies** (``lambda_poisson``,
              ``lambda_gradient``) override the method to return their own
              ``compute_lambda(heights)`` — fully dynamic, in lock-step with
              whatever background the strategy uses for peak significance.
              The user's ``min_pas_prominence`` is ignored for them.
            - **Other strategies** (``original``, ``sierra_iterative``) fall
              through to the base implementation which returns
              ``min_pas_prominence`` unchanged — a static user-set threshold.

        The pair is merged when ``valley < threshold``: the dip between the
        candidates is not significantly above background.  Set
        ``min_pas_prominence < 0`` to disable Tier 2 entirely.

    Reads inside the gap are absorbed automatically: the next pipeline call
    ``strategy.get_cb_dict_for_pas(peak, merged_start, merged_end)`` routes
    through ``reconstruct_cb_dict``, which slices ``cb_positions`` by
    inclusive range and picks up every position in the widened span.

    Args:
        regions: List of ``(start, end)`` tuples from ``find_pas``.
        peak: ``Peak`` instance providing ``peak_list`` (list of
            ``(position, height)``).
        strategy: The strategy instance that produced ``regions``.  Its
            :meth:`~ema.strategies.base.PeakFinderStrategy.valley_threshold`
            decides whether Tier 2 uses a dynamic (lambda) or static
            (user-set) valley threshold.
        min_pas_spacing: Distance floor in bp.  Pairs with gap below this
            value are always merged.  ``0`` disables the distance tier.
        min_pas_prominence: Static fallback valley threshold passed to
            ``strategy.valley_threshold`` for non-lambda strategies.
            Lambda strategies ignore this value.  Negative disables Tier 2.

    Returns:
        List of merged ``(start, end)`` tuples sorted by genomic position.
        Returns ``regions`` unchanged when it has 0 or 1 entries, or when
        both tiers are disabled.
    """
    if len(regions) <= 1:
        return regions

    tier2_active = min_pas_prominence >= 0
    tier1_active = min_pas_spacing > 0
    if not tier1_active and not tier2_active:
        return regions

    regions = sorted(regions)
    pos = [p for p, _ in peak.peak_list]
    h = [hi for _, hi in peak.peak_list]

    # Ask the strategy for the right threshold.  Lambda strategies return
    # compute_lambda(heights); the base implementation returns the user's
    # min_pas_prominence unchanged.
    if tier2_active and strategy is not None:
        threshold = strategy.valley_threshold(peak, min_pas_prominence)
    else:
        threshold = float(min_pas_prominence)

    def _min_h_between(re: int, rs: int) -> int:
        vals = [h[i] for i, p in enumerate(pos) if re < p < rs]
        # No measured positions inside the gap → treat as zero-height valley
        # (i.e. maximally deep) so an empty gap never forces a Tier-2 merge.
        return min(vals) if vals else 0

    merged: List[Tuple[int, int]] = [regions[0]]
    for s, e in regions[1:]:
        ps, pe = merged[-1]
        gap = s - pe

        # Tier 1: distance floor
        if tier1_active and gap < min_pas_spacing:
            merged[-1] = (ps, e)
            continue

        # Tier 2: valley vs threshold
        if tier2_active:
            valley = _min_h_between(pe, s)
            if valley < threshold:
                merged[-1] = (ps, e)
                continue

        merged.append((s, e))

    return merged


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
