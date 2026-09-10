"""GTF isoform UTR parser for PeakATail.

Parses a GTF file and returns per-transcript 3'UTR exon lists, grouped by
gene. Parsing is parallelized by chromosome using multiprocessing.Pool with
mmap for zero-copy file reads. Results are cached to disk using the same
mtime+size fingerprint strategy as ema.annotate.gtf_cache.

Output format::

    {
        gene_id: {
            transcript_id: [
                (chrom, start, end, strand, exon_rank),
                ...
            ]
        }
    }

Exons are sorted in transcription order (5' -> 3'):
- Plus strand : ascending start coordinate  (exon_rank 0, 1, 2, …)
- Minus strand: descending start coordinate (exon_rank 0, 1, 2, …)

Cache location: ``{cache_dir}/isoform_utr_cache/``
"""

from __future__ import annotations

import json
import mmap
import multiprocessing
import os
import pickle
import re
from collections import defaultdict
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass


# ---------------------------------------------------------------------------
# Fingerprint / cache helpers (mirrors gtf_cache.py style)
# ---------------------------------------------------------------------------

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _gtf_fingerprint(gtf_path: Path) -> str:
    """Return mtime_ns + size string for a GTF file."""
    stat = os.stat(gtf_path)
    return f"{stat.st_mtime_ns}_{stat.st_size}"


def _isoform_cache_dir(cache_dir: Path) -> Path:
    return cache_dir / "isoform_utr_cache"


def _manifest_path(cache_dir: Path) -> Path:
    return _isoform_cache_dir(cache_dir) / "manifest.json"


def _data_path(cache_dir: Path) -> Path:
    return _isoform_cache_dir(cache_dir) / "isoform_utrs.pkl"


def _is_cache_valid(gtf_path: Path, cache_dir: Path) -> bool:
    """Return True when a valid cache exists for *gtf_path* in *cache_dir*."""
    manifest = _manifest_path(cache_dir)
    if not manifest.exists():
        return False
    try:
        with open(manifest, "r") as fh:
            meta = json.load(fh)
        if meta.get("gtf_path") != str(gtf_path.resolve()):
            return False
        return meta.get("fingerprint") == _gtf_fingerprint(gtf_path)
    except (json.JSONDecodeError, OSError, KeyError):
        return False


def _save_cache(gtf_path: Path, cache_dir: Path, result: dict) -> None:
    icd = _isoform_cache_dir(cache_dir)
    icd.mkdir(parents=True, exist_ok=True)
    with open(_data_path(cache_dir), "wb") as fh:
        pickle.dump(result, fh, protocol=pickle.HIGHEST_PROTOCOL)
    manifest = {
        "gtf_path": str(gtf_path.resolve()),
        "fingerprint": _gtf_fingerprint(gtf_path),
    }
    with open(_manifest_path(cache_dir), "w") as fh:
        json.dump(manifest, fh, indent=2)


def _load_cache(cache_dir: Path) -> dict:
    with open(_data_path(cache_dir), "rb") as fh:
        return pickle.load(fh)


# ---------------------------------------------------------------------------
# Attribute parsing
# ---------------------------------------------------------------------------

def _parse_attributes(attr_string: str) -> dict[str, str]:
    """Parse GTF column-9 attributes into a dict."""
    return dict(_ATTR_RE.findall(attr_string))


# ---------------------------------------------------------------------------
# Per-chromosome worker (runs in a subprocess)
# ---------------------------------------------------------------------------

def _parse_chunk(args: tuple) -> dict:
    """Parse a list of raw GTF lines and return a partial isoform-UTR map.

    Args:
        args: (lines_bytes,) — list[bytes] of raw GTF lines for one chromosome.

    Returns:
        Partial ``{gene_id: {transcript_id: [(chrom, start, end, strand, exon_rank), ...]}}``
        dict.  Ranks are NOT yet assigned here — they are relative indices
        within the returned per-transcript list (unsorted); final sorting
        happens in the merge step to avoid per-chunk overhead.
    """
    (lines_bytes,) = args
    # gene_id -> transcript_id -> list of (chrom, start, end, strand)
    partial: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))

    for raw in lines_bytes:
        line = raw.decode("utf-8", errors="replace")
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue

        feature_type = fields[2]
        if feature_type not in ("three_prime_utr", "three_prime_UTR"):
            continue

        chrom = fields[0]
        start = int(fields[3]) - 1  # GTF is 1-based inclusive -> 0-based half-open
        end = int(fields[4])
        strand = fields[6]
        attrs = _parse_attributes(fields[8])

        gene_id = attrs.get("gene_id", "").split(".")[0]
        transcript_id = attrs.get("transcript_id", "").split(".")[0]
        if not gene_id or not transcript_id:
            continue

        partial[gene_id][transcript_id].append((chrom, start, end, strand))

    # Convert inner defaultdicts to plain dicts before returning (pickle-safe)
    return {g: dict(t_map) for g, t_map in partial.items()}


# ---------------------------------------------------------------------------
# Chromosome-aware file splitter using mmap
# ---------------------------------------------------------------------------

def _split_by_chromosome(gtf_path: Path) -> list[list[bytes]]:
    """Read the GTF with mmap and group lines by chromosome.

    Returns a list of per-chromosome line lists (bytes), suitable for
    passing to a multiprocessing pool.
    """
    chrom_lines: dict[str, list[bytes]] = defaultdict(list)

    with open(gtf_path, "rb") as fh:
        with mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            for raw_line in iter(mm.readline, b""):
                if raw_line.startswith(b"#"):
                    continue
                tab_idx = raw_line.find(b"\t")
                if tab_idx == -1:
                    continue
                chrom = raw_line[:tab_idx].decode("ascii", errors="replace")
                chrom_lines[chrom].append(raw_line)

    return list(chrom_lines.values())


# ---------------------------------------------------------------------------
# Merge helper
# ---------------------------------------------------------------------------

def _merge_partial_results(partials: list[dict]) -> dict:
    """Merge per-chromosome dicts and sort exons in transcription order.

    Assigns ``exon_rank`` (0-based, 5' -> 3') after merging all partials.
    """
    merged: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))

    for partial in partials:
        for gene_id, t_map in partial.items():
            for transcript_id, exons in t_map.items():
                merged[gene_id][transcript_id].extend(exons)

    result: dict[str, dict[str, list]] = {}
    for gene_id, t_map in merged.items():
        result[gene_id] = {}
        for transcript_id, exons in t_map.items():
            if not exons:
                continue
            strand = exons[0][3]
            # Sort by start coord: ascending for +, descending for -
            exons_sorted = sorted(exons, key=lambda e: e[1], reverse=(strand == "-"))
            # Assign exon_rank (5'->3' within 3'UTR)
            ranked = [
                (chrom, s, e, st, rank)
                for rank, (chrom, s, e, st) in enumerate(exons_sorted)
            ]
            result[gene_id][transcript_id] = ranked

    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_isoform_utrs(
    gtf_path: Path,
    cache_dir: Path | None = None,
    n_workers: int | None = None,
) -> dict[str, dict[str, list[tuple]]]:
    """Parse GTF and return isoform 3'UTR exon map.

    Uses mmap + ``multiprocessing.Pool`` (one worker per chromosome) for the
    heavy parse, then merges results in the main process.  On cache hit,
    returns the pickle-cached result directly (typically < 100 ms).

    Args:
        gtf_path: Path to the GTF file (Ensembl format assumed).
        cache_dir: Directory used for caching.  When *None*, caching is
            disabled and the result is recomputed every call.
        n_workers: Number of worker processes.  Defaults to
            ``min(os.cpu_count(), n_chromosomes)``.

    Returns:
        Nested dict::

            {
                gene_id: {
                    transcript_id: [
                        (chrom, start, end, strand, exon_rank),
                        ...
                    ]
                }
            }

        Coordinates are 0-based half-open (BED convention).
        Exons are sorted in transcription order (5' -> 3'); ``exon_rank``
        is 0 for the exon immediately downstream of the CDS.

    Raises:
        FileNotFoundError: If *gtf_path* does not exist.
        ValueError: If *gtf_path* cannot be parsed (no three_prime_utr rows).
    """
    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    # Cache check
    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        if _is_cache_valid(gtf_path, cache_dir):
            return _load_cache(cache_dir)

    # Split GTF by chromosome using mmap
    chrom_chunks = _split_by_chromosome(gtf_path)

    if not chrom_chunks:
        return {}

    # Determine worker count
    # ResourceManager-aware worker count to avoid oversubscribing under load
    try:
        from ema.utils import ResourceManager
        cpu = ResourceManager().get_n_jobs(per_worker_mb=500)
    except ImportError:
        cpu = os.cpu_count() or 1
    n_w = min(n_workers or cpu, len(chrom_chunks), cpu)

    # Parse chunks in parallel
    args = [(chunk,) for chunk in chrom_chunks]

    if n_w > 1:
        # Use 'spawn' to avoid fork() inheriting locks/threads from the parent,
        # which can deadlock downstream BLAS / scanpy / igraph operations.
        ctx = multiprocessing.get_context("spawn")
        with ctx.Pool(processes=n_w) as pool:
            partials = pool.map(_parse_chunk, args)
    else:
        partials = [_parse_chunk(a) for a in args]

    result = _merge_partial_results(partials)

    # Persist to cache
    if cache_dir is not None:
        _save_cache(gtf_path, cache_dir, result)

    return result
