"""CLI entry point for pre-warming the global GTF cache.

Usage::

    ema_parse_gtf --gtf /path/to/annotation.gtf
    ema_parse_gtf --gtf /path/to/annotation.gtf --cache-dir /tmp/my_cache

Running this before the main pipeline ensures the first pipeline run skips
the GTF parse entirely and loads from the global cache instead.
"""

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path


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
        print(f"ERROR: GTF file not found: {gtf_path}", file=sys.stderr)
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
        print(f"[ema_parse_gtf] cache HIT — {entry_dir}")
        print(f"  endbed    : {hit['endbed_path']}")
        print(f"  features  : {hit['features_path']}")
        print(f"  utr_lengths: {hit['isoform_utrs_path']}")
        print(f"  genes     : {len(hit['utr_lengths'])}")
        return 0

    print(f"[ema_parse_gtf] cache MISS — parsing {gtf_path} …")
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
            print(f"ERROR during GTF parsing: {exc}", file=sys.stderr)
            return 1

        populate_global_cache(
            gtf_path=str(gtf_path),
            utr_lengths=utr_lengths,
            endbed_path=endbed,
            features_path=features,
            utr_lengths_path=utr_tsv,
        )

    elapsed = time.perf_counter() - t0
    print(f"[ema_parse_gtf] done in {elapsed:.1f}s — cache written to {entry_dir}")
    print(f"  genes parsed: {len(utr_lengths)}")
    return 0


if __name__ == "__main__":
    sys.exit(cli())
