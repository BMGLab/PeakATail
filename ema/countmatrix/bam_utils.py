"""Small BAM-inspection helpers used by the peak-calling driver."""
from __future__ import annotations

import pysam


def infer_median_read_length(bam_path: str, n: int = 1000) -> int:
    """Return the median read length detected from the first ``n`` reads of a BAM.

    Used by the post-detection PAS merger to set its Tier-1 distance floor.
    The default ``-1`` value of ``variable_config.min_pas_spacing`` triggers
    this auto-detection once per BAM in ``ema/main.py``; the result is cached
    in ``variable_config.dataset_read_lengths``.

    Args:
        bam_path: Absolute path to an indexed BAM file.
        n: Maximum number of reads to sample.  Iterating the first ``n``
            reads is effectively free even for very large BAMs.

    Returns:
        Median read length in base pairs.  Returns ``100`` as a safe
        fallback when no read has a populated ``query_length`` (extremely
        unusual but possible for header-only or hard-clipped BAMs).
    """
    bam = pysam.AlignmentFile(bam_path, "rb", threads=1)
    lengths: list[int] = []
    try:
        for read in bam:
            ql = read.query_length
            if ql:
                lengths.append(ql)
                if len(lengths) >= n:
                    break
    finally:
        bam.close()

    if not lengths:
        return 100

    lengths.sort()
    return lengths[len(lengths) // 2]
