"""Runner function for BAM merging.

Provides the library-level ``run_merge`` entry point consumed by
``ema/cli/merge.py``.  The old argparse ``cli()`` has been removed;
use ``ema merge`` instead.
"""
from __future__ import annotations

import logging

from ema.merge_bam.merge import merge

log = logging.getLogger(__name__)


def run_merge(bam_files: list[str], output: str, threads: int = 4) -> None:
    """Library-level entry point for BAM merging.

    Merges one or more BAM files into a single sorted, indexed BAM.  Body is
    the post-argparse logic of the legacy ``cli()`` function, promoted so
    ``ema/cli/merge.py`` can call it without going through argparse.

    Args:
        bam_files: List of absolute paths to input BAM files.
        output: Output BAM path.
        threads: Number of samtools threads (default: 4).
    """
    log.info("run_merge: merging %d BAM(s) -> %s (threads=%d)", len(bam_files), output, threads)
    merge(bam_files=bam_files, threads=threads)
    log.info("run_merge: done.")
