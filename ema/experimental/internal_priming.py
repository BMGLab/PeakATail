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


def filter_internal_priming(bed_path: str, genome_fasta: str,
                            output_path: str,
                            window_left: int = 10,
                            window_right: int = 30,
                            a_stretch: int = 6,
                            a_fraction: float = 0.7,
                            mode: str = "annotate") -> dict:
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
        import shutil
        shutil.copy2(bed_path, output_path)
        return {"total": 0, "passed": 0, "filtered": 0, "filtered_fraction": 0,
                "flagged": 0, "flagged_fraction": 0, "mode": mode, "flags": {},
                "error": "pyfaidx not installed"}

    if mode not in ("annotate", "filter"):
        raise ValueError(f"filter_internal_priming: mode must be 'annotate' or 'filter', got {mode!r}")

    genome = Fasta(genome_fasta)

    total = 0
    passed = 0
    filtered = 0
    flagged = 0
    flags: dict[str, bool] = {}

    with open(bed_path) as infile, open(output_path, 'w') as outfile:
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
                continue

            flags[pas_id] = is_internal_priming
            if is_internal_priming:
                flagged += 1
                logger.debug(f"Internal priming: {chrom}:{pas_pos} strand={strand} seq={seq[:20]}...")

            if mode == "filter" and is_internal_priming:
                filtered += 1
            else:
                outfile.write(line)
                passed += 1

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

    logger.info(f"Internal priming filter (mode={mode}): {total} total, {passed} passed, "
                f"{filtered} filtered ({filtered_fraction:.1%}), {flagged} flagged ({flagged_fraction:.1%})")

    return stats
