"""GTF processing cache for PeakATail.

Caches parsed GTF results (gene BED, UTR lengths, features) on disk.
Uses file mtime + size to determine per-output-dir cache validity. Also
maintains a *global* cache keyed by SHA256(first 4 KB) + mtime_ns + size so
that parsed results survive across different output directories and runs.

Cache hierarchy (checked in order):
  1. Global cache  — ~/.cache/peakatail/gtf/<fingerprint>/
  2. Per-output-dir cache — {output_dir}/gtf_cache/

This module is designed to be called from a background thread/process
by the GTF pre-processing pipeline so that annotation data is ready
before find_close needs it.
"""

import hashlib
import json
import logging
import os
import shutil
import tempfile
from pathlib import Path

log = logging.getLogger(__name__)



# ---------------------------------------------------------------------------
# Low-level helpers shared by both cache tiers
# ---------------------------------------------------------------------------

def _gtf_fingerprint(gtf_path: str | os.PathLike) -> str:
    """Compute a per-output-dir fingerprint (mtime_ns + size only, cheap)."""
    stat = os.stat(gtf_path)
    return f"{stat.st_mtime_ns}_{stat.st_size}"


def _cache_dir(output_dir: str | os.PathLike) -> str:
    """Return the per-output-dir cache directory path."""
    return os.path.join(output_dir, "gtf_cache")


def _manifest_path(output_dir: str | os.PathLike) -> str:
    """Return the path to the per-output-dir cache manifest file."""
    return os.path.join(_cache_dir(output_dir), "manifest.json")


# ---------------------------------------------------------------------------
# Global cache helpers
# ---------------------------------------------------------------------------

def global_cache_dir() -> Path:
    """Return the global GTF cache directory, honouring XDG_CACHE_HOME.

    Returns:
        Path: ``$XDG_CACHE_HOME/peakatail/gtf`` if ``$XDG_CACHE_HOME`` is set,
        otherwise ``~/.cache/peakatail/gtf``.
    """
    xdg = os.environ.get("XDG_CACHE_HOME", "").strip()
    base = Path(xdg) if xdg else Path.home() / ".cache"
    return base / "peakatail" / "gtf"


def gtf_global_fingerprint(gtf_path: str | os.PathLike) -> str:
    """Compute a robust global fingerprint for a GTF file.

    Combines SHA256 of the first 4 KB of the file with mtime_ns and size so
    the key is both cheap to compute and resistant to content changes that
    preserve mtime (e.g. same-second overwrites, NFS with coarse timestamps).

    Args:
        gtf_path: Path to the GTF file.

    Returns:
        str: ``<sha256_prefix>_<mtime_ns>_<size>``
    """
    stat = os.stat(gtf_path)
    with open(gtf_path, "rb") as fh:
        head = fh.read(4096)
    sha = hashlib.sha256(head).hexdigest()[:16]
    return f"{sha}_{stat.st_mtime_ns}_{stat.st_size}"


def lookup_global_cache(gtf_path: str | os.PathLike) -> dict | None:
    """Probe the global cache for a pre-parsed GTF.

    Args:
        gtf_path: Path to the GTF file.

    Returns:
        dict with keys ``utr_lengths``, ``endbed_path``, ``features_path``,
        ``isoform_utrs_path`` if a valid cache entry exists; ``None`` otherwise.
    """
    fp = gtf_global_fingerprint(gtf_path)
    entry_dir = global_cache_dir() / fp
    manifest_file = entry_dir / "manifest.json"

    if not manifest_file.exists():
        return None

    try:
        with open(manifest_file, "r") as fh:
            meta = json.load(fh)

        # Verify fingerprint still matches (guards against hash collisions)
        if meta.get("fingerprint") != fp:
            return None

        utr_json = entry_dir / "utr_lengths.json"
        if not utr_json.exists():
            return None

        with open(utr_json, "r") as fh:
            utr_lengths = json.load(fh)

        return {
            "utr_lengths": utr_lengths,
            "endbed_path": str(entry_dir / meta["endbed"]),
            "features_path": str(entry_dir / meta["features"]),
            "isoform_utrs_path": str(entry_dir / meta.get("utr_lengths_tsv", "utr_lengths.tsv")),
        }
    except (json.JSONDecodeError, OSError, KeyError):
        return None


def populate_global_cache(
    gtf_path: str | os.PathLike,
    utr_lengths: dict,
    endbed_path: str | os.PathLike,
    features_path: str | os.PathLike,
    utr_lengths_path: str | os.PathLike,
) -> dict:
    """Write parsed GTF artifacts to the global cache atomically.

    Uses a temporary sibling directory then renames it so concurrent writers
    or a mid-write crash cannot produce a partial cache entry.

    Args:
        gtf_path: Path to the source GTF file.
        utr_lengths: Dict mapping gene_id -> max UTR length.
        endbed_path: Path to the generated gene BED file.
        features_path: Path to the generated features TSV.
        utr_lengths_path: Path to the generated UTR lengths TSV.

    Returns:
        dict: Same shape as ``lookup_global_cache`` return value.
    """
    fp = gtf_global_fingerprint(gtf_path)
    entry_dir = global_cache_dir() / fp
    parent = entry_dir.parent
    parent.mkdir(parents=True, exist_ok=True)

    if entry_dir.exists():
        # Already populated by a concurrent process — nothing to do.
        hit = lookup_global_cache(gtf_path)
        if hit is not None:
            return hit

    # Write to a temp dir then rename for atomicity
    tmp_dir = Path(tempfile.mkdtemp(dir=parent, prefix=".tmp_gtf_"))
    try:
        for src in (endbed_path, features_path, utr_lengths_path):
            src = Path(src)
            if src.exists():
                shutil.copy2(src, tmp_dir / src.name)

        utr_json = tmp_dir / "utr_lengths.json"
        with open(utr_json, "w") as fh:
            json.dump(utr_lengths, fh)

        manifest = {
            "fingerprint": fp,
            "gtf_path": os.path.abspath(gtf_path),
            "endbed": Path(endbed_path).name,
            "features": Path(features_path).name,
            "utr_lengths_tsv": Path(utr_lengths_path).name,
        }
        with open(tmp_dir / "manifest.json", "w") as fh:
            json.dump(manifest, fh, indent=2)

        # Atomic rename (same filesystem guaranteed — both under global_cache_dir)
        try:
            tmp_dir.rename(entry_dir)
        except OSError:
            # Another process beat us to it — clean up and use theirs
            shutil.rmtree(tmp_dir, ignore_errors=True)
    except Exception:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise

    return {
        "utr_lengths": utr_lengths,
        "endbed_path": str(entry_dir / Path(endbed_path).name),
        "features_path": str(entry_dir / Path(features_path).name),
        "isoform_utrs_path": str(entry_dir / Path(utr_lengths_path).name),
    }


# ---------------------------------------------------------------------------
# Per-output-dir cache (unchanged API, kept for backward compat)
# ---------------------------------------------------------------------------

def is_cache_valid(gtf_path: str | os.PathLike, output_dir: str | os.PathLike) -> bool:
    """Check if per-output-dir cached GTF results are still valid.

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
        with open(manifest, "r") as fh:
            meta = json.load(fh)

        if meta.get("gtf_path") != os.path.abspath(gtf_path):
            return False

        current_fp = _gtf_fingerprint(gtf_path)
        return meta.get("fingerprint") == current_fp
    except (json.JSONDecodeError, OSError, KeyError):
        return False


def save_to_cache(
    gtf_path: str | os.PathLike,
    output_dir: str | os.PathLike,
    utr_lengths: dict,
    endbed_path: str | os.PathLike,
    features_path: str | os.PathLike,
    utr_lengths_path: str | os.PathLike,
) -> None:
    """Save GTF processing results to the per-output-dir cache.

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

    for src in (endbed_path, features_path, utr_lengths_path):
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(cache, os.path.basename(src)))

    utr_cache = os.path.join(cache, "utr_lengths.json")
    with open(utr_cache, "w") as fh:
        json.dump(utr_lengths, fh)

    manifest = {
        "gtf_path": os.path.abspath(gtf_path),
        "fingerprint": _gtf_fingerprint(gtf_path),
        "endbed": os.path.basename(endbed_path),
        "features": os.path.basename(features_path),
        "utr_lengths_tsv": os.path.basename(utr_lengths_path),
    }
    with open(_manifest_path(output_dir), "w") as fh:
        json.dump(manifest, fh, indent=2)


def load_from_cache(
    output_dir: str | os.PathLike,
    endbed_path: str | os.PathLike,
    features_path: str | os.PathLike,
    utr_lengths_path: str | os.PathLike,
) -> dict:
    """Load cached GTF results from the per-output-dir cache.

    Restores output files to their expected locations and returns UTR lengths.

    Args:
        output_dir: Output directory containing the cache.
        endbed_path: Destination path for the gene BED file.
        features_path: Destination path for the features TSV.
        utr_lengths_path: Destination path for the UTR lengths TSV.

    Returns:
        dict: Mapping of gene_id -> max UTR length (bp).
    """
    cache = _cache_dir(output_dir)

    with open(_manifest_path(output_dir), "r") as fh:
        meta = json.load(fh)

    for cached_name, dest in [
        (meta["endbed"], endbed_path),
        (meta["features"], features_path),
        (meta["utr_lengths_tsv"], utr_lengths_path),
    ]:
        cached_file = os.path.join(cache, cached_name)
        if os.path.exists(cached_file) and os.path.abspath(cached_file) != os.path.abspath(dest):
            shutil.copy2(cached_file, dest)

    utr_json = os.path.join(cache, "utr_lengths.json")
    with open(utr_json, "r") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Main entry point — tries global cache, then per-output-dir, then parses
# ---------------------------------------------------------------------------

def process_gtf_cached(
    gtf_path: str | os.PathLike,
    output_dir: str | os.PathLike,
    endbed_path: str | os.PathLike,
    features_path: str | os.PathLike,
    utr_lengths_path: str | os.PathLike | None = None,
    **gtf_kwargs,
) -> dict:
    """Process GTF file with a two-tier cache. Returns UTR lengths dict.

    Cache lookup order:
      1. Global cache  (``~/.cache/peakatail/gtf/<sha_mtime_size>/``)
      2. Per-output-dir cache  (``{output_dir}/gtf_cache/``)
      3. Full parse via ``gtftobed.gtf_bed()`` → populate both caches.

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

    # --- Tier 1: global cache ---
    global_hit = lookup_global_cache(gtf_path)
    if global_hit is not None:
        fp_short = gtf_global_fingerprint(gtf_path)[:16]
        log.info("global hit at %s", global_cache_dir() / fp_short[:8])
        # Restore files to output_dir locations the rest of the pipeline expects
        for src_key, dest in [
            ("endbed_path", endbed_path),
            ("features_path", features_path),
            ("isoform_utrs_path", utr_lengths_path),
        ]:
            src = global_hit[src_key]
            if os.path.exists(src) and os.path.abspath(src) != os.path.abspath(dest):
                os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)
                shutil.copy2(src, dest)
        # Also refresh per-output-dir cache from global data so next run skips
        # even this copy step when output_dir is the same.
        save_to_cache(gtf_path, output_dir, global_hit["utr_lengths"],
                      endbed_path, features_path, utr_lengths_path)
        return global_hit["utr_lengths"]

    # --- Tier 2: per-output-dir cache ---
    if is_cache_valid(gtf_path, output_dir):
        log.info("local hit")
        return load_from_cache(output_dir, endbed_path, features_path, utr_lengths_path)

    # --- Tier 3: full parse ---
    log.info("miss — parsing GTF ...")
    # Import lazily to avoid circular / argparse side-effects at module load
    # time.  The module-level name _gtf_bed_fn is set here so tests can patch
    # ema.annotate.gtf_cache._gtf_bed_fn before the function runs.
    _fn = globals().get("_gtf_bed_fn")
    if _fn is None:
        from ema.annotate.gtftobed import gtf_bed as _fn  # noqa: PLC0415
    utr_lengths = _fn(
        endbeddir=endbed_path,
        gtfdir=gtf_path,
        featuresdir=features_path,
        utr_lengths_dir=utr_lengths_path,
        **gtf_kwargs,
    )

    # Populate both caches so future runs in any output dir benefit
    populate_global_cache(gtf_path, utr_lengths, endbed_path, features_path,
                          utr_lengths_path)
    save_to_cache(gtf_path, output_dir, utr_lengths, endbed_path,
                  features_path, utr_lengths_path)

    return utr_lengths
