"""CLI entry point for pre-warming the global GTF cache.

Usage::

    ema_parse_gtf --gtf /path/to/annotation.gtf
    ema_parse_gtf --gtf /path/to/annotation.gtf --cache-dir /tmp/my_cache

Running this before the main pipeline ensures the first pipeline run skips
the GTF parse entirely and loads from the global cache instead.
"""

import argparse
import logging
import os
import sys
import tempfile
import time
from pathlib import Path

log = logging.getLogger(__name__)


def cli() -> int:
    """Entry point for the ``ema_parse_gtf`` console script.

    Returns:
        int: Exit code — 0 on success, 1 on error.
    """
    parser = argparse.ArgumentParser(
        prog="ema_parse_gtf",
        description="Pre-warm the PeakATail global GTF cache.",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        metavar="PATH",
        help="Path to the GTF annotation file.",
    )
    parser.add_argument(
        "--cache-dir",
        metavar="PATH",
        default=None,
        help=(
            "Override the global cache directory "
            "(default: $XDG_CACHE_HOME/peakatail/gtf or ~/.cache/peakatail/gtf)."
        ),
    )
    args = parser.parse_args()

    gtf_path = Path(args.gtf)
    if not gtf_path.exists():
        log.error("GTF file not found: %s", gtf_path)
        return 1

    # Import here so startup is fast even when module tree is large
    from ema.annotate.gtf_cache import (
        global_cache_dir,
        gtf_global_fingerprint,
        lookup_global_cache,
        populate_global_cache,
    )
    from ema.annotate.gtftobed import gtf_bed

    # Allow caller to override cache root (useful for tests / CI)
    if args.cache_dir is not None:
        import ema.annotate.gtf_cache as _mod
        _orig_global_cache_dir = _mod.global_cache_dir

        override = Path(args.cache_dir)

        def _patched_global_cache_dir() -> Path:  # noqa: ANN202
            return override

        _mod.global_cache_dir = _patched_global_cache_dir  # type: ignore[method-assign]
        # Re-import after monkey-patch so the rest of this function sees it
        from ema.annotate import gtf_cache as _gc
        _gc.global_cache_dir = _patched_global_cache_dir  # type: ignore[method-assign]

    fp = gtf_global_fingerprint(gtf_path)
    cache_root = global_cache_dir()
    entry_dir = cache_root / fp

    # --- Probe for existing entry ---
    hit = lookup_global_cache(gtf_path)
    if hit is not None:
        log.info("cache HIT — %s", entry_dir)
        log.info("  endbed    : %s", hit['endbed_path'])
        log.info("  features  : %s", hit['features_path'])
        log.info("  utr_lengths: %s", hit['isoform_utrs_path'])
        log.info("  genes     : %d", len(hit['utr_lengths']))
        return 0

    log.info("cache MISS — parsing %s ...", gtf_path)
    t0 = time.perf_counter()

    # Parse into a temp output directory (pipeline artifacts written there,
    # then promoted to the global cache).
    with tempfile.TemporaryDirectory(prefix="ema_parse_gtf_") as tmp_out:
        endbed = os.path.join(tmp_out, "gene_ends.bed")
        features = os.path.join(tmp_out, "raw_features.tsv")
        utr_tsv = os.path.join(tmp_out, "utr_lengths.tsv")

        try:
            utr_lengths = gtf_bed(
                endbeddir=endbed,
                gtfdir=str(gtf_path),
                featuresdir=features,
                utr_lengths_dir=utr_tsv,
            )
        except Exception as exc:
            log.error("ERROR during GTF parsing: %s", exc)
            return 1

        populate_global_cache(
            gtf_path=str(gtf_path),
            utr_lengths=utr_lengths,
            endbed_path=endbed,
            features_path=features,
            utr_lengths_path=utr_tsv,
        )

    elapsed = time.perf_counter() - t0
    log.info("done in %.1fs — cache written to %s", elapsed, entry_dir)
    log.info("  genes parsed: %d", len(utr_lengths))
    return 0


if __name__ == "__main__":
    sys.exit(cli())
