"""Data-driven 3' cleavage-site offset correction for called PAS.

Motivation
----------
Atlas-independent motif analysis of PeakATail's called peaks (Laughney
cohort, 22,629 PAS; see GitHub issue #72) showed that the reported peak
**3' end systematically stops ~90-105 nt short of the true cleavage site**:

* AATAAA positional density peaks at **+75 nt** downstream of the reported
  peak end (canonical AATAAA -> cleavage spacing is 15-30 nt, implying the
  cleavage site sits at roughly +90..+105).
* Genomic A-fraction crests at 43% at **+98 nt** then cliffs to background --
  the poly(A)-junction signature -- exactly where 10x R2 coverage runs out.

The offset is a chemistry artefact (R2 read length runs out before the
poly(A) junction), so downstream benchmarks scored at tight cutoffs punish
the *offset*, not the calls.  Shifting the reported 3' end downstream by the
observed offset recovers a large chunk of benchmark precision for free, for
every strategy at once.

This module provides the *pure* coordinate transform plus a (currently
stub) data-driven estimator.  The correction is **opt-in** and **defaults
to legacy behaviour** (offset 0 = no shift); it is gated behind the
``--cleavage-offset`` CLI flag / ``cleavage_offset`` YAML key.

Strand convention
-----------------
BED strand is encoded as a character in column 6.  For a ``+`` strand PAS
the transcript 3' end is the **larger** genomic coordinate (``bed_end``);
for a ``-`` strand PAS it is the **smaller** genomic coordinate
(``bed_start``).  "Downstream" always means *in the direction of
transcription* -- increasing coordinate on ``+``, decreasing on ``-``.
"""

from __future__ import annotations

from typing import Iterable

# Sane fallback constant (bp) used when a per-run data-driven estimate is
# not available.  Midpoint of the observed +90..+105 nt window (issue #72).
DEFAULT_CLEAVAGE_OFFSET: int = 95


def shift_3prime_end(
    bed_start: int,
    bed_end: int,
    strand: str,
    offset: int,
) -> tuple[int, int]:
    """Extend the 3' end of a BED interval downstream by ``offset`` bp.

    The 5' end of the peak (where R2 coverage is real) is preserved; only
    the 3' end -- the inferred cleavage site -- moves downstream.

    Args:
        bed_start: 0-based BED start coordinate (inclusive).
        bed_end: 0-based BED end coordinate (exclusive).
        strand: ``"+"`` or ``"-"``.  Any other value is treated as ``"+"``
            (forward), matching the permissive behaviour of the rest of the
            pipeline.
        offset: Non-negative shift in bp.  ``0`` (or negative) is a no-op and
            returns the interval unchanged -- the legacy default.

    Returns:
        ``(new_start, new_end)`` with the 3' end shifted downstream.  Never
        returns a negative coordinate (``-`` strand shifts clamp at 0) and
        never returns a degenerate interval (``new_end > new_start`` is
        preserved).
    """
    if offset <= 0:
        return bed_start, bed_end

    if strand == "-":
        # 3' end is the smaller coordinate; downstream = decreasing coord.
        new_start = max(0, bed_start - offset)
        return new_start, bed_end

    # "+" (and the permissive fallback): 3' end is the larger coordinate;
    # downstream = increasing coord.
    return bed_start, bed_end + offset


def inferred_cleavage_site(
    bed_start: int,
    bed_end: int,
    strand: str,
    offset: int,
) -> int:
    """Return the single inferred cleavage-site coordinate (point mode).

    This is the point most benchmarks score against ("strand-matched
    point-mode scoring").  For ``+`` it is ``bed_end + offset``; for ``-``
    it is ``max(0, bed_start - offset)``.
    """
    new_start, new_end = shift_3prime_end(bed_start, bed_end, strand, offset)
    return new_start if strand == "-" else new_end


def estimate_cleavage_offset(
    *,
    aataaa_profile: "Iterable[float] | None" = None,
    a_fraction_profile: "Iterable[float] | None" = None,
    fallback: int = DEFAULT_CLEAVAGE_OFFSET,
) -> int:
    """Estimate the per-run cleavage offset from downstream signal profiles.

    TODO(issue #72): implement the real data-driven estimator.  The intended
    algorithm (both signals are computed over a downstream window anchored at
    the reported peak 3' end):

      1. AATAAA route -- take the mode/argmax of the AATAAA positional
         density profile downstream of the peak end and add the canonical
         AATAAA->cleavage spacing (15-30 nt, use ~22 nt midpoint).  The
         Laughney data put the AATAAA mode at +75 nt => ~+97 nt cleavage.
      2. A-fraction route -- find the position where the genomic A-fraction
         profile crests before cliffing to background (the poly(A) junction).
         The Laughney data put this at +98 nt.
      3. Reconcile the two estimates (e.g. average, or prefer AATAAA when the
         profile has a clear mode) and clamp to a plausible [60, 130] window.

    Until that lands, this function returns ``fallback`` (the sane constant)
    whenever the profiles are not supplied, so callers get correct,
    well-defined behaviour today.  When profiles *are* supplied we still fall
    back (the estimator is a stub) but the signature is stable so wiring the
    real implementation later is a drop-in change.

    Args:
        aataaa_profile: Optional AATAAA positional-density profile, one value
            per bp downstream of the peak end.  Currently unused (stub).
        a_fraction_profile: Optional genomic A-fraction profile, one value per
            bp downstream of the peak end.  Currently unused (stub).
        fallback: Constant returned when no estimate can be made.

    Returns:
        Estimated offset in bp (currently always ``fallback``).
    """
    # --- STUB -------------------------------------------------------------
    # The real estimator (see docstring TODO) will consume the profiles.
    # For now we return the sane fallback constant so the opt-in flag has
    # correct, deterministic behaviour without the (larger) profiling code.
    _ = (aataaa_profile, a_fraction_profile)
    return int(fallback)


def rewrite_bed_3prime_offset(path, offset: int) -> int:
    """Rewrite a 6-column PAS BED in place, shifting each 3' end downstream.

    Reads a strand-aware BED (strand in column 6), applies
    :func:`shift_3prime_end` to every record, and rewrites the file.  Blank
    lines and malformed rows (fewer than 6 columns, non-integer coordinates)
    are passed through unchanged so the function is safe to run on any BED.

    Args:
        path: Path to the BED file (str or :class:`pathlib.Path`).
        offset: Non-negative shift in bp.  ``<= 0`` is a no-op (returns 0)
            and leaves the file byte-identical -- the legacy default.

    Returns:
        Number of records whose coordinates were shifted.
    """
    import os

    if offset <= 0:
        return 0

    shifted = 0
    tmp = f"{os.fspath(path)}.offset.tmp"
    with open(path) as src, open(tmp, "w") as dst:
        for line in src:
            if not line.strip():
                dst.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                dst.write(line)
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                dst.write(line)
                continue
            strand = parts[5]
            new_start, new_end = shift_3prime_end(start, end, strand, offset)
            if (new_start, new_end) != (start, end):
                shifted += 1
            parts[1] = str(new_start)
            parts[2] = str(new_end)
            dst.write("\t".join(parts) + "\n")
    os.replace(tmp, path)
    return shifted
