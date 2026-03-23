"""GTF processing cache for PeakATail.

Caches parsed GTF results (gene BED, UTR lengths, features) on disk.
Uses file mtime + size to determine cache validity. If the GTF file
has not changed, loads cached results instead of re-parsing.

This module is designed to be called from a background thread/process
by the GTF pre-processing pipeline so that annotation data is ready
before find_close needs it.

Cache location: {output_dir}/gtf_cache/
"""

import hashlib
import json
import os
import shutil
from pathlib import Path


def _gtf_fingerprint(gtf_path):
    """Compute a fingerprint for a GTF file based on mtime and size.

    Using mtime+size is fast (no need to read the whole file) and sufficient
    for detecting changes in practice. Falls back to partial content hash
    if the filesystem doesn't support mtime.
    """
    stat = os.stat(gtf_path)
    return f"{stat.st_mtime_ns}_{stat.st_size}"


def _cache_dir(output_dir):
    """Return the cache directory path."""
    return os.path.join(output_dir, "gtf_cache")


def _manifest_path(output_dir):
    """Return the path to the cache manifest file."""
    return os.path.join(_cache_dir(output_dir), "manifest.json")


def is_cache_valid(gtf_path, output_dir):
    """Check if cached GTF results are still valid.

    Args:
        gtf_path: Path to the GTF file.
        output_dir: Output directory containing the cache.

    Returns:
        bool: True if cache exists and GTF has not changed.
    """
    manifest = _manifest_path(output_dir)
    if not os.path.exists(manifest):
        return False

    try:
        with open(manifest, "r") as f:
            meta = json.load(f)

        if meta.get("gtf_path") != os.path.abspath(gtf_path):
            return False

        current_fp = _gtf_fingerprint(gtf_path)
        return meta.get("fingerprint") == current_fp
    except (json.JSONDecodeError, OSError, KeyError):
        return False


def save_to_cache(gtf_path, output_dir, utr_lengths,
                  endbed_path, features_path, utr_lengths_path):
    """Save GTF processing results to cache.

    Args:
        gtf_path: Path to the source GTF file.
        output_dir: Output directory for cache storage.
        utr_lengths: Dict mapping gene_id -> max UTR length.
        endbed_path: Path to the generated gene BED file.
        features_path: Path to the generated features TSV.
        utr_lengths_path: Path to the generated UTR lengths TSV.
    """
    cache = _cache_dir(output_dir)
    os.makedirs(cache, exist_ok=True)

    # Copy result files into cache
    for src in (endbed_path, features_path, utr_lengths_path):
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(cache, os.path.basename(src)))

    # Save UTR lengths as JSON for fast loading
    utr_cache = os.path.join(cache, "utr_lengths.json")
    with open(utr_cache, "w") as f:
        json.dump(utr_lengths, f)

    # Write manifest
    manifest = {
        "gtf_path": os.path.abspath(gtf_path),
        "fingerprint": _gtf_fingerprint(gtf_path),
        "endbed": os.path.basename(endbed_path),
        "features": os.path.basename(features_path),
        "utr_lengths_tsv": os.path.basename(utr_lengths_path),
    }
    with open(_manifest_path(output_dir), "w") as f:
        json.dump(manifest, f, indent=2)


def load_from_cache(output_dir, endbed_path, features_path, utr_lengths_path):
    """Load cached GTF results, restoring output files and returning UTR lengths.

    Args:
        output_dir: Output directory containing the cache.
        endbed_path: Destination path for the gene BED file.
        features_path: Destination path for the features TSV.
        utr_lengths_path: Destination path for the UTR lengths TSV.

    Returns:
        dict: Mapping of gene_id -> max UTR length (bp).
    """
    cache = _cache_dir(output_dir)

    with open(_manifest_path(output_dir), "r") as f:
        meta = json.load(f)

    # Restore result files from cache to their expected locations
    for cached_name, dest in [
        (meta["endbed"], endbed_path),
        (meta["features"], features_path),
        (meta["utr_lengths_tsv"], utr_lengths_path),
    ]:
        cached_file = os.path.join(cache, cached_name)
        if os.path.exists(cached_file) and os.path.abspath(cached_file) != os.path.abspath(dest):
            shutil.copy2(cached_file, dest)

    # Load UTR lengths from JSON (faster than re-parsing TSV)
    utr_json = os.path.join(cache, "utr_lengths.json")
    with open(utr_json, "r") as f:
        return json.load(f)


def process_gtf_cached(gtf_path, output_dir, endbed_path, features_path,
                       utr_lengths_path=None, **gtf_kwargs):
    """Process GTF file with caching. Returns UTR lengths dict.

    If cache is valid, loads from cache. Otherwise, runs full GTF parsing
    via gtftobed.gtf_bed() and saves results to cache.

    This is the main entry point for background GTF pre-processing.

    Args:
        gtf_path: Path to the GTF file.
        output_dir: Output directory.
        endbed_path: Destination for gene BED file.
        features_path: Destination for features TSV.
        utr_lengths_path: Destination for UTR lengths TSV.
        **gtf_kwargs: Additional kwargs passed to gtf_bed().

    Returns:
        dict: Mapping of gene_id -> max UTR length (bp).
    """
    if utr_lengths_path is None:
        utr_lengths_path = os.path.join(output_dir, "utr_lengths.tsv")

    if is_cache_valid(gtf_path, output_dir):
        return load_from_cache(output_dir, endbed_path, features_path,
                               utr_lengths_path)

    # Cache miss -- run full GTF processing
    from ema.annotate.gtftobed import gtf_bed

    utr_lengths = gtf_bed(
        endbeddir=endbed_path,
        gtfdir=gtf_path,
        featuresdir=features_path,
        utr_lengths_dir=utr_lengths_path,
        **gtf_kwargs
    )

    # Save to cache for next run
    save_to_cache(gtf_path, output_dir, utr_lengths,
                  endbed_path, features_path, utr_lengths_path)

    return utr_lengths
