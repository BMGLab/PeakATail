"""Internal priming filter for PeakATail.

Filters peaks near genomic A-rich stretches that cause false positive
poly(A) signals from internal priming during reverse transcription.

Reference implementations:
- SAPAS: >=6 A's genomic check
- scraps: A-content analysis
- scAPAtrap: 6 continuous A's within -140 to +10 bp
- scPAISO: 6-mer AAAAAA within -5/+20 bp

Strand convention
-----------------
Internal priming is caused by a genomic A-rich stretch DOWNSTREAM of the
cleavage site in TRANSCRIPT orientation (the oligo-dT primer anneals to the
genome-encoded A's that the reverse transcriptase then reads as a tail). The
window is therefore defined relative to the transcript: ``window_left`` nt
upstream and ``window_right`` nt downstream of the cleavage site, on BOTH
strands. In genomic (forward) coordinates that means (``pos`` is the BED
``end`` on '+' and the BED ``start`` on '-'; see :func:`ip_window`)::

    '+'  [pos - left,  pos + right)   scanned as-is      (A-run / A-fraction)
    '-'  [pos - right, pos + left )   reverse-complemented, then the same scan
                                      (== T-run / T-fraction on the forward strand)

Up to and including the 4efeb12 line the '-' strand used the SAME forward
window as '+', i.e. it tested 30 nt upstream / 10 nt downstream in transcript
orientation -- mostly the wrong side. Fixed in fix/ip-filter-strand.

This filter is OPTIONAL and configurable via CLI:
  --internal-priming-filter     Enable the filter
  --genome-fasta PATH           Path to genome FASTA (required if filter enabled)
  --ip-window-left INT          Left window from PAS (default: -10)
  --ip-window-right INT         Right window from PAS (default: +30)
  --ip-a-stretch INT            Minimum consecutive A's to flag (default: 6)
  --ip-a-fraction FLOAT         Alternative: fraction of A's in window (default: 0.7)
"""

import logging
from typing import Tuple

logger = logging.getLogger("peakatail.filters.internal_priming")

_RC_TABLE = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def reverse_complement(seq: str) -> str:
    """Reverse complement of *seq* (IUPAC letters other than ACGTUN pass through)."""
    return seq.translate(_RC_TABLE)[::-1]


def ip_window(pas_pos: int, strand: str,
              window_left: int = 10, window_right: int = 30) -> Tuple[int, int]:
    """Genomic half-open ``[start, end)`` window to fetch for one PAS.

    *pas_pos* is the BED ``end`` on '+' and the BED ``start`` on '-' (the
    convention :func:`filter_internal_priming` has always used). The window
    covers *window_left* nt upstream and *window_right* nt downstream of the
    cleavage site **in transcript orientation**, so on '-' it is the mirror
    image of the '+' window in genomic coordinates. ``start`` is clamped at
    0; the caller's slice clamps ``end`` at the contig length.
    """
    if strand == "+":
        return max(0, pas_pos - window_left), pas_pos + window_right
    return max(0, pas_pos - window_right), pas_pos + window_left


def call_internal_priming(window_seq: str, strand: str,
                          a_stretch: int = 6, a_fraction: float = 0.7) -> Tuple[bool, str]:
    """Decide internal priming from the FORWARD-strand genomic sequence of
    the :func:`ip_window` of one PAS.

    Returns ``(flag, tested_seq)`` where *tested_seq* is the sequence in
    transcript orientation (upper-cased; reverse-complemented on '-'). The
    flag is set when *tested_seq* contains ``'A' * a_stretch`` or its
    A-fraction is ``>= a_fraction``. Pure function -- no I/O.
    """
    seq = window_seq.upper()
    tested = seq if strand == "+" else reverse_complement(seq)
    if not tested:
        return False, tested
    if "A" * a_stretch in tested:
        return True, tested
    return tested.count("A") / len(tested) >= a_fraction, tested


def check_internal_priming(sequence, pas_pos: int, strand: str,
                           window_left: int = 10, window_right: int = 30,
                           a_stretch: int = 6, a_fraction: float = 0.7) -> Tuple[bool, str]:
    """:func:`ip_window` + :func:`call_internal_priming` on a sliceable contig.

    *sequence* may be a plain ``str`` or a ``pyfaidx`` record (anything whose
    ``[start:end]`` slice stringifies to the forward-strand bases; slices past
    the contig end are expected to truncate, as both do).
    """
    seq_start, seq_end = ip_window(pas_pos, strand, window_left, window_right)
    return call_internal_priming(str(sequence[seq_start:seq_end]), strand,
                                 a_stretch, a_fraction)


class _NullWriter:
    """A write sink for the features-only scan (``output_path=None``)."""

    def write(self, _line: str) -> None:  # pragma: no cover - trivial
        return

    def __enter__(self):
        return self

    def __exit__(self, *exc) -> bool:
        return False


def _open_or_null(path):
    return _NullWriter() if path is None else open(path, "w")


def _bed_score(parts: list) -> int:
    """BED column 5 as an int (0 when absent or non-numeric)."""
    try:
        return int(parts[4])
    except (IndexError, ValueError):
        return 0


def filter_internal_priming(bed_path: str, genome_fasta: str,
                            output_path: str | None,
                            window_left: int = 10,
                            window_right: int = 30,
                            a_stretch: int = 6,
                            a_fraction: float = 0.7,
                            mode: str = "annotate",
                            features=None) -> dict:
    """Check peaks for genomic A-rich stretches indicating internal priming.

    Reads a BED file, checks the genomic sequence around each PAS for
    A-rich stretches that indicate internal priming artifacts.

    Args:
        bed_path: Input BED file with PAS peaks.
        genome_fasta: Path to genome FASTA file (must be indexed with .fai).
        output_path: Path to write output BED.
        window_left: bp upstream of PAS to check (default 10).
        window_right: bp downstream of PAS to check (default 30).
        a_stretch: Minimum consecutive A's to flag as internal priming (default 6).
        a_fraction: Alternative: flag if A-fraction in window exceeds this (default 0.7).
        output_path: Where to write the (possibly filtered) BED.  ``None``
            scans without writing a BED at all -- the peakAtail-prime
            features-only pass, which needs the genome open over the same
            rows but changes nothing.  ``mode`` must be ``"annotate"`` then.
        features: Optional
            :class:`ema.countmatrix.pas_features.FeatureCollector`.  When
            given, each row also contributes its per-site scoring features,
            computed INSIDE this loop from the already-open genome -- there
            is no second pass over the FASTA and no second pass over the BED.
            Collecting features never changes which rows are written or
            flagged; ``ip_tool_flag`` is emitted **in addition to** the veto,
            never instead of it.
        mode: ``"annotate"`` (default) or ``"filter"``.
            * ``"annotate"`` -- an internally-primed peak is likely real
              alternative-PAS signal, not noise, so EVERY peak is written to
              ``output_path`` unchanged; the internal-priming flag is only
              returned (see ``stats["flags"]``), never used to drop a row.
            * ``"filter"`` -- today's pre-D6+ behaviour: flagged peaks are
              dropped from ``output_path``.

    Returns:
        Dict with filtering statistics:
        - total: total peaks processed
        - passed: peaks written to output_path
        - filtered: peaks removed (always 0 in "annotate" mode)
        - filtered_fraction: fraction removed
        - flagged: peaks flagged as likely internal priming (both modes)
        - flagged_fraction: fraction flagged
        - mode: the mode this call ran in
        - flags: ``{pas_id: bool}`` (BED column 4 -> internal_priming),
          covering every well-formed (>=6 column) row processed.
    """
    try:
        from pyfaidx import Fasta
    except ImportError:
        logger.error("pyfaidx not installed. Run: pip install pyfaidx")
        logger.error("Skipping internal priming filter.")
        # Copy input to output unchanged
        if output_path is not None:
            import shutil
            shutil.copy2(bed_path, output_path)
        return {"total": 0, "passed": 0, "filtered": 0, "filtered_fraction": 0,
                "flagged": 0, "flagged_fraction": 0, "mode": mode, "flags": {},
                "error": "pyfaidx not installed"}

    if mode not in ("annotate", "filter"):
        raise ValueError(f"filter_internal_priming: mode must be 'annotate' or 'filter', got {mode!r}")
    if output_path is None and mode != "annotate":
        raise ValueError(
            "filter_internal_priming: output_path=None only makes sense with "
            f"mode='annotate' (a scan that drops nothing), got mode={mode!r}"
        )

    if features is not None:
        from ema.countmatrix.pas_features import (
            bed_cleavage, feature_window, ip_covariates, na_features,
            orient_window, sequence_features,
        )

    genome = Fasta(genome_fasta)

    total = 0
    passed = 0
    filtered = 0
    flagged = 0
    flags: dict[str, bool] = {}

    with open(bed_path) as infile, _open_or_null(output_path) as outfile:
        for line in infile:
            total += 1
            parts = line.strip().split('\t')
            if len(parts) < 6:
                outfile.write(line)
                passed += 1
                continue

            pas_id = parts[3]
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = parts[5]

            # PAS position depends on strand
            if strand == '+':
                pas_pos = end  # 3' end for positive strand
            else:
                pas_pos = start  # 3' end for negative strand

            # Extract the genomic window (transcript-oriented; mirrored on '-')
            # and test it -- see ip_window() / call_internal_priming().
            try:
                is_internal_priming, seq = check_internal_priming(
                    genome[chrom], pas_pos, strand,
                    window_left, window_right, a_stretch, a_fraction,
                )
            except (KeyError, ValueError):
                # Chromosome not in FASTA or out of range — keep the peak,
                # flag unknown/not-flagged (we couldn't check it).
                outfile.write(line)
                passed += 1
                flags[pas_id] = False
                if features is not None:
                    _f = na_features()
                    features.features[str(pas_id)] = _f
                    features.rows.append((chrom, strand,
                                          bed_cleavage(start, end, strand),
                                          _bed_score(parts), str(pas_id)))
                continue

            flags[pas_id] = is_internal_priming
            if features is not None:
                # Sequence features from the SAME open contig record, one
                # extra slice of 71 nt.  `seq` is the string the veto just
                # tested, so ip_tool_afrac/arun describe exactly that window
                # under whatever --ip-window-left/right the run used.
                _afrac, _arun = ip_covariates(seq)
                _fs, _fe = feature_window(pas_pos, strand)
                _tseq, _r0 = orient_window(
                    str(genome[chrom][_fs:_fe]), pas_pos, strand, _fs)
                _feat = sequence_features(_tseq, _r0)
                _feat["ip_tool_flag"] = int(is_internal_priming)
                _feat["ip_tool_afrac"] = round(_afrac, 4)
                _feat["ip_tool_arun"] = _arun
                features.add(pas_id, chrom, strand,
                             bed_cleavage(start, end, strand),
                             _bed_score(parts), _feat)
            if is_internal_priming:
                flagged += 1
                logger.debug(f"Internal priming: {chrom}:{pas_pos} strand={strand} seq={seq[:20]}...")

            if mode == "filter" and is_internal_priming:
                filtered += 1
            else:
                outfile.write(line)
                passed += 1

    if features is not None:
        stats_features = features.stats()
    filtered_fraction = filtered / total if total > 0 else 0
    flagged_fraction = flagged / total if total > 0 else 0

    stats = {
        "total": total,
        "passed": passed,
        "filtered": filtered,
        "filtered_fraction": round(filtered_fraction, 4),
        "flagged": flagged,
        "flagged_fraction": round(flagged_fraction, 4),
        "mode": mode,
        "flags": flags,
    }
    if features is not None:
        stats["features"] = stats_features

    logger.info(f"Internal priming filter (mode={mode}): {total} total, {passed} passed, "
                f"{filtered} filtered ({filtered_fraction:.1%}), {flagged} flagged ({flagged_fraction:.1%})")

    return stats
