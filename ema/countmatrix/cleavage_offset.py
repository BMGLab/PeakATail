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

This module provides the *pure* coordinate transform plus a data-driven
estimator (:func:`estimate_cleavage_offset_diag` /
:func:`build_downstream_profiles`).  The correction is **opt-in** and
**defaults to legacy behaviour** (offset 0 = no shift): a manual constant is
set with ``--cleavage-offset N`` / ``cleavage_offset: N``, and the estimator
is engaged with ``--auto-cleavage-offset`` (which reads the called peaks +
``--genome-fasta`` and infers the offset per run).

Strand convention
-----------------
BED strand is encoded as a character in column 6.  For a ``+`` strand PAS
the transcript 3' end is the **larger** genomic coordinate (``bed_end``);
for a ``-`` strand PAS it is the **smaller** genomic coordinate
(``bed_start``).  "Downstream" always means *in the direction of
transcription* -- increasing coordinate on ``+``, decreasing on ``-``.
"""

from __future__ import annotations

import logging
import os
import random
from typing import Iterable, Sequence

log = logging.getLogger(__name__)

# Sane fallback constant (bp) used when a per-run data-driven estimate is
# not available.  Midpoint of the observed +90..+105 nt window (issue #72).
DEFAULT_CLEAVAGE_OFFSET: int = 95

# Default downstream window (bp) profiled from each peak 3' end.
DEFAULT_MAX_OFFSET: int = 150
# Plausible band (bp downstream) the crest is searched within.  The observed
# A-fraction crest sits at +98 nt (issue #72); 60..120 brackets it with room.
DEFAULT_BAND: tuple[int, int] = (60, 120)
# Canonical AATAAA -> cleavage spacing (nt) added to the AATAAA density peak
# for the secondary cross-check estimate (issue #72: AATAAA mode +75 -> ~+97).
AATAAA_CLEAVAGE_SPACING: int = 22

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


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


def _smooth(profile: Sequence[float], window: int) -> list[float]:
    """Centred moving-average smoother (pure Python, edge-truncated)."""
    n = len(profile)
    if window <= 1 or n == 0:
        return list(profile)
    half = window // 2
    out: list[float] = []
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        seg = profile[lo:hi]
        out.append(sum(seg) / len(seg))
    return out


def _median(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    m = len(s) // 2
    if len(s) % 2:
        return s[m]
    return 0.5 * (s[m - 1] + s[m])


def estimate_cleavage_offset(
    *,
    aataaa_profile: "Iterable[float] | None" = None,
    a_fraction_profile: "Iterable[float] | None" = None,
    fallback: int = DEFAULT_CLEAVAGE_OFFSET,
    band: tuple[int, int] = DEFAULT_BAND,
    smooth_window: int = 5,
    min_prominence: float = 0.05,
) -> int:
    """Estimate the per-run cleavage offset from downstream signal profiles.

    Thin ``int``-returning wrapper around :func:`estimate_cleavage_offset_diag`
    (kept for backwards compatibility with callers that only want the number).
    See that function for the full algorithm and the diagnostics dict.
    """
    offset, _ = estimate_cleavage_offset_diag(
        aataaa_profile=aataaa_profile,
        a_fraction_profile=a_fraction_profile,
        fallback=fallback,
        band=band,
        smooth_window=smooth_window,
        min_prominence=min_prominence,
    )
    return offset


def estimate_cleavage_offset_diag(
    *,
    aataaa_profile: "Iterable[float] | None" = None,
    a_fraction_profile: "Iterable[float] | None" = None,
    fallback: int = DEFAULT_CLEAVAGE_OFFSET,
    band: tuple[int, int] = DEFAULT_BAND,
    smooth_window: int = 5,
    min_prominence: float = 0.05,
    n_peaks_used: int = 0,
) -> tuple[int, dict]:
    """Estimate the cleavage offset from downstream profiles + diagnostics.

    Algorithm (both profiles are indexed by bp downstream of the reported
    peak 3' end; index ``i`` == offset ``i``):

      1. **A-fraction crest (primary).** Lightly smooth the genomic A-fraction
         profile and take its ``argmax`` restricted to the plausible band
         ``[band_lo, band_hi]``.  This is the poly(A)-junction crest that
         cliffs to background (issue #72 put it at +98 nt).  The crest is
         accepted only if it rises at least ``min_prominence`` above the
         profile background (its median) -- otherwise the profile is deemed
         flat / inconclusive.
      2. **AATAAA density peak (secondary cross-check).** ``argmax`` of the
         AATAAA positional-density profile plus the canonical AATAAA->cleavage
         spacing (~22 nt).  Reported in the diagnostics for cross-checking; it
         is used as the estimate only when the A-fraction route is
         inconclusive but the AATAAA route lands inside the band.
      3. **Fallback.** If neither route yields a confident in-band estimate,
         return ``fallback`` (the documented ~95 nt constant) and mark the
         method ``"fallback"`` so the caller can log that estimation was
         inconclusive.

    Args:
        aataaa_profile: AATAAA positional-density profile (one value per bp
            downstream of the peak end), or ``None``.
        a_fraction_profile: Genomic A-fraction profile (one value per bp
            downstream of the peak end), or ``None``.
        fallback: Constant returned when no confident estimate can be made.
        band: ``(lo, hi)`` plausible offset window (inclusive) to search.
        smooth_window: Moving-average width used to smooth the A-fraction
            profile before taking the argmax.
        min_prominence: Minimum A-fraction rise above background for the crest
            to be accepted.
        n_peaks_used: Number of peaks the profiles were aggregated over (passed
            straight into the diagnostics dict; not used in the maths).

    Returns:
        ``(offset, diagnostics)`` where ``diagnostics`` is a small dict with
        ``method`` in ``{"a_fraction_crest", "aataaa_spacing", "fallback"}``,
        ``crest_pos``, ``crest_value``, ``background``, ``aataaa_peak`` and
        ``aataaa_estimate`` (any of which may be ``None``).
    """
    band_lo, band_hi = int(band[0]), int(band[1])
    diag: dict = {
        "method": "fallback",
        "n_peaks_used": int(n_peaks_used),
        "band": [band_lo, band_hi],
        "crest_pos": None,
        "crest_value": None,
        "background": None,
        "aataaa_peak": None,
        "aataaa_estimate": None,
        "fallback": int(fallback),
    }

    af = list(a_fraction_profile) if a_fraction_profile is not None else []
    aa = list(aataaa_profile) if aataaa_profile is not None else []

    # --- Secondary: AATAAA density peak + canonical spacing ---------------
    aataaa_estimate = None
    if aa:
        aa_lo = max(0, band_lo - AATAAA_CLEAVAGE_SPACING)
        aa_hi = min(len(aa) - 1, band_hi - AATAAA_CLEAVAGE_SPACING)
        if aa_hi >= aa_lo:
            seg = aa[aa_lo:aa_hi + 1]
            if any(v > 0 for v in seg):
                aa_peak = aa_lo + max(range(len(seg)), key=seg.__getitem__)
                diag["aataaa_peak"] = int(aa_peak)
                aataaa_estimate = int(aa_peak + AATAAA_CLEAVAGE_SPACING)
                diag["aataaa_estimate"] = aataaa_estimate

    # --- Primary: A-fraction crest ----------------------------------------
    if af and len(af) > band_lo:
        smoothed = _smooth(af, smooth_window)
        background = _median(smoothed)
        hi = min(band_hi, len(smoothed) - 1)
        seg = smoothed[band_lo:hi + 1]
        if seg:
            crest_rel = max(range(len(seg)), key=seg.__getitem__)
            crest_pos = band_lo + crest_rel
            crest_value = smoothed[crest_pos]
            diag["crest_pos"] = int(crest_pos)
            diag["crest_value"] = float(crest_value)
            diag["background"] = float(background)
            if crest_value - background >= min_prominence:
                diag["method"] = "a_fraction_crest"
                return int(crest_pos), diag

    # --- A-fraction inconclusive: fall back to AATAAA route if in-band ----
    if aataaa_estimate is not None and band_lo <= aataaa_estimate <= band_hi:
        diag["method"] = "aataaa_spacing"
        return int(aataaa_estimate), diag

    return int(fallback), diag


def build_downstream_profiles(
    bed_paths: "str | os.PathLike | Sequence[str | os.PathLike]",
    genome_fasta: "str | os.PathLike",
    *,
    max_offset: int = DEFAULT_MAX_OFFSET,
    max_sample: int = 3000,
    max_n_fraction: float = 0.2,
    rng_seed: int = 0,
) -> tuple[list[float], list[float], int]:
    """Aggregate downstream A-fraction and AATAAA-density profiles from BEDs.

    For each confidently-stranded PAS record the genomic window ``+0..
    +max_offset`` nt *downstream* of the reported 3' end is extracted
    (strand-aware; ``-`` strand windows are reverse-complemented so
    "downstream" is consistent), upper-cased, and accumulated into two
    per-position profiles:

      * ``a_fraction[i]`` -- fraction of ``A`` at downstream offset ``i``.
      * ``aataaa[i]``     -- fraction of sampled peaks with ``AATAAA`` starting
        at downstream offset ``i`` (positional density).

    Windows that run off the chromosome or contain more than ``max_n_fraction``
    ``N`` bases are skipped.

    Args:
        bed_paths: One BED path or a sequence of them (6-column, strand in
            col 6).  Missing files are skipped.
        genome_fasta: Genome FASTA indexed for pyfaidx (``.fai``).
        max_offset: Length of the downstream window (bp).
        max_sample: Cap on the number of peaks profiled (reservoir-sampled
            across all inputs; all peaks are used when fewer than the cap).
        max_n_fraction: Skip a window if the ``N`` fraction exceeds this.
        rng_seed: Seed for the reservoir sampler (determinism).

    Returns:
        ``(a_fraction_profile, aataaa_profile, n_peaks_used)``.  Both profiles
        have length ``max_offset``; ``n_peaks_used`` is how many windows were
        actually aggregated.

    Raises:
        ImportError: if pyfaidx is not installed.
    """
    from pyfaidx import Fasta

    if isinstance(bed_paths, (str, os.PathLike)):
        bed_paths = [bed_paths]

    rng = random.Random(rng_seed)
    # Reservoir sample of (chrom, start, end, strand) across all BEDs.
    reservoir: list[tuple[str, int, int, str]] = []
    seen = 0
    for path in bed_paths:
        try:
            fh = open(path)
        except FileNotFoundError:
            continue
        with fh:
            for line in fh:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                strand = parts[5]
                if strand not in ("+", "-"):
                    continue  # only confidently-stranded peaks
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                rec = (parts[0], start, end, strand)
                seen += 1
                if len(reservoir) < max_sample:
                    reservoir.append(rec)
                else:
                    j = rng.randint(0, seen - 1)
                    if j < max_sample:
                        reservoir[j] = rec

    a_count = [0] * max_offset
    total = [0] * max_offset
    aataaa_count = [0] * max_offset
    n_used = 0

    genome = Fasta(genome_fasta)
    for chrom, start, end, strand in reservoir:
        if strand == "+":
            anchor = end
            lo, hi = anchor, anchor + max_offset
        else:
            anchor = start
            lo, hi = anchor - max_offset, anchor
        if lo < 0:
            continue
        try:
            seq = str(genome[chrom][lo:hi]).upper()
        except (KeyError, ValueError):
            continue
        if len(seq) < max_offset:
            continue  # ran off the chromosome end
        if strand == "-":
            seq = _revcomp(seq)
        if seq.count("N") / max_offset > max_n_fraction:
            continue
        n_used += 1
        for i, base in enumerate(seq):
            if base == "N":
                continue
            total[i] += 1
            if base == "A":
                a_count[i] += 1
        # AATAAA positional density.
        idx = seq.find("AATAAA")
        while idx != -1:
            if idx < max_offset:
                aataaa_count[idx] += 1
            idx = seq.find("AATAAA", idx + 1)

    a_fraction = [
        (a_count[i] / total[i]) if total[i] else 0.0 for i in range(max_offset)
    ]
    denom = n_used if n_used else 1
    aataaa_profile = [aataaa_count[i] / denom for i in range(max_offset)]
    return a_fraction, aataaa_profile, n_used


def resolve_cleavage_offset(
    explicit_offset: int,
    auto: bool,
    bed_paths: "Sequence[str | os.PathLike] | None" = None,
    genome_fasta: "str | os.PathLike | None" = None,
    *,
    max_offset: int = DEFAULT_MAX_OFFSET,
    max_sample: int = 3000,
    fallback: int = DEFAULT_CLEAVAGE_OFFSET,
    band: tuple[int, int] = DEFAULT_BAND,
) -> tuple[int, dict | None]:
    """Resolve the effective cleavage offset for a run.

    This is the single decision point the pipeline calls:

      * ``auto`` is ``False``  -> return ``explicit_offset`` unchanged
        (``0`` = legacy no-op; a positive int = manual constant).  No FASTA or
        BED access; ``diagnostics`` is ``None``.
      * ``auto`` is ``True``   -> build downstream profiles from ``bed_paths``
        + ``genome_fasta``, run :func:`estimate_cleavage_offset_diag`, and
        return the data-driven estimate together with its diagnostics dict.
        If the FASTA is missing or pyfaidx is unavailable, fall back to
        ``fallback`` with a ``method="fallback"`` diagnostics dict.

    Returns:
        ``(offset, diagnostics)``.  ``diagnostics`` is ``None`` in manual mode.
    """
    if not auto:
        return int(explicit_offset), None

    if not genome_fasta or not os.path.exists(str(genome_fasta)):
        log.warning(
            "auto cleavage-offset requested but genome FASTA is missing (%r); "
            "falling back to the default constant %d bp.",
            genome_fasta, fallback,
        )
        return int(fallback), {"method": "fallback", "reason": "no_fasta",
                               "fallback": int(fallback), "n_peaks_used": 0}

    try:
        a_fraction, aataaa, n_used = build_downstream_profiles(
            list(bed_paths or []), genome_fasta,
            max_offset=max_offset, max_sample=max_sample,
        )
    except ImportError:
        log.warning(
            "auto cleavage-offset requested but pyfaidx is not installed; "
            "falling back to the default constant %d bp.", fallback,
        )
        return int(fallback), {"method": "fallback", "reason": "no_pyfaidx",
                               "fallback": int(fallback), "n_peaks_used": 0}

    if n_used == 0:
        log.warning(
            "auto cleavage-offset: no usable downstream windows; falling back "
            "to the default constant %d bp.", fallback,
        )
        return int(fallback), {"method": "fallback", "reason": "no_windows",
                               "fallback": int(fallback), "n_peaks_used": 0}

    offset, diag = estimate_cleavage_offset_diag(
        aataaa_profile=aataaa, a_fraction_profile=a_fraction,
        fallback=fallback, band=band, n_peaks_used=n_used,
    )
    if diag.get("method") == "fallback":
        log.warning(
            "auto cleavage-offset: A-fraction/AATAAA profiles inconclusive "
            "(n=%d peaks); falling back to the default constant %d bp.",
            n_used, fallback,
        )
    else:
        log.info(
            "auto cleavage-offset: estimated %d bp from %d peaks (method=%s, "
            "crest_value=%s).", offset, n_used, diag.get("method"),
            diag.get("crest_value"),
        )
    return int(offset), diag


def rewrite_bed_3prime_offset(path, offset: int, *, skip_supported: bool = False) -> int:
    """Rewrite a 6-column PAS BED in place, shifting each 3' end downstream.

    Reads a strand-aware BED (strand in column 6), applies
    :func:`shift_3prime_end` to every record, and rewrites the file.  Blank
    lines and malformed rows (fewer than 6 columns, non-integer coordinates)
    are passed through unchanged so the function is safe to run on any BED.

    Args:
        path: Path to the BED file (str or :class:`pathlib.Path`).
        offset: Non-negative shift in bp.  ``<= 0`` is a no-op (returns 0)
            and leaves the file byte-identical -- the legacy default.
        skip_supported: When True, rows whose BED score (column 5) is > 0
            are passed through UNSHIFTED.  The ``clip_seeded`` strategy
            writes the poly(A) clip-read count there, and those PAS are
            already placed at the observed cleavage site -- shifting them
            would push them *past* it.  Only coverage-only rows (score 0)
            carry the R2-read-length offset this correction exists for.

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
            if skip_supported:
                try:
                    if float(parts[4]) > 0:
                        dst.write(line)
                        continue
                except ValueError:
                    pass
            new_start, new_end = shift_3prime_end(start, end, strand, offset)
            if (new_start, new_end) != (start, end):
                shifted += 1
            parts[1] = str(new_start)
            parts[2] = str(new_end)
            dst.write("\t".join(parts) + "\n")
    os.replace(tmp, path)
    return shifted


# ===========================================================================
# peakAtail-prime (TASK E): the CLIP-ANCHORED offset, and a signed rewrite
# ===========================================================================
#
# Everything above this line is the COVERAGE caller's estimator (issue #72):
# it looks for a genomic A-fraction crest 60-120 bp downstream of a peak's
# reported 3' end, because a coverage peak stops where R2 coverage runs out.
# It cannot return a small or a negative offset -- the search band forbids it
# -- and the +95 bp it returns is measurably destructive for ``clip_seeded``,
# whose tier-1 PAS are already ON the cleavage base (P@10 0.5209 -> 0.0551,
# results/algo_headroom/A4_resolution table T7).
#
# What follows is the different quantity TASK E asked for: the offset between
# the base the caller REPORTS and the bases the poly(A) tails actually pinned,
# estimated from the run's own clip-anchored subset -- no genome, no long
# reads, no atlas.  The caller writes one number per tier-1 PAS
# (``clip_offset_mean`` in ``pas_support.tsv``, ``--pas-features on``); this
# module aggregates them into the per-library estimate.
#
# MEASURED, so nobody has to re-derive it: on the PBMC chr19+21 slice the
# read-weighted estimate is -0.334 bp over 267,520 clip reads in 16,338
# clusters (76.10 % of clip reads sit EXACTLY on the reported base); on the
# GSE104556 mouse 1 chr18+19 slice it is +0.214 bp over 124,308 reads in
# 8,722 clusters (56.90 % exactly on it).  Both round to ZERO, and both are
# an average of two opposite per-strand values (+0.808/-0.748 PBMC,
# +0.999/-1.102 mouse) produced by ``cluster_clip_sites``'s tie-break, which
# resolves ties toward the LOWEST COORDINATE on both strands and is therefore
# not strand-symmetric.
#
# The offset external truths prefer is NOT this number: on the same slice the
# base-pair-exact optimum is -1 bp against the atlas and -2 bp against Kinnex
# long reads.  The clip channel cannot see it -- which is exactly why
# ``--cleavage-offset`` defaults to ``none`` and why ``auto`` is honest about
# estimating ~0.

#: Accepted spellings of ``--cleavage-offset`` beyond a signed integer.
CLEAVAGE_OFFSET_MODES = ("none", "auto")

#: The v2-compatibility value of ``--cleavage-offset``.
V2_CLEAVAGE_OFFSET = "none"


def parse_cleavage_offset(spec) -> tuple[str, int]:
    """Resolve ``--cleavage-offset`` into ``(mode, constant)``.

    Accepts, in order of precedence:

    * ``"none"`` / ``None`` / ``""``      -> ``("none", 0)`` -- v2, no shift.
    * ``"auto"``                          -> ``("auto", 0)`` -- use this run's
      clip-anchored estimate (resolved later, once the clusters exist).
    * an ``int`` or an integer string     -> ``("const", value)``.

    ``int`` is accepted because ``variable_config.cleavage_offset`` was an int
    in v2 and a library caller may still assign one; ``0`` is spelled
    ``("none", 0)`` so the two v2 spellings collapse to one no-op.

    Raises:
        ValueError: on anything else, so a typo cannot be read as "no shift".
    """
    if spec is None:
        return ("none", 0)
    if isinstance(spec, bool):                       # guard: True is not 1 here
        raise ValueError("cleavage_offset must be none/auto/<int>, got %r" % spec)
    if isinstance(spec, int):
        return ("none", 0) if spec == 0 else ("const", int(spec))
    text = str(spec).strip().lower()
    if text in ("", "none", "off", "0"):
        return ("none", 0)
    if text == "auto":
        return ("auto", 0)
    try:
        value = int(text)
    except ValueError:
        raise ValueError(
            "cleavage_offset must be 'none', 'auto' or a signed integer "
            "number of bp, got %r" % (spec,)
        ) from None
    return ("none", 0) if value == 0 else ("const", value)


def shift_cleavage_point(bed_start: int, bed_end: int, strand: str,
                         offset: int) -> tuple[int, int]:
    """Move a PAS interval's 3' end by a SIGNED *offset*, transcript-oriented.

    Positive is downstream and reproduces :func:`shift_3prime_end` exactly
    (the v2 behaviour: the 5' end stays put and the interval grows).  Negative
    is upstream, which v2 could not express -- it returned the interval
    unchanged -- and is what a base-pair-resolution correction needs.

    A 1-bp interval (every ``clip_seeded`` tier-1 PAS) is carried whole, so
    the reported base really moves; a wider interval keeps its 5' end unless
    that would make it degenerate.  ``offset == 0`` returns the input tuple
    unchanged, byte-for-byte.

    Args:
        bed_start: 0-based BED start (inclusive).
        bed_end: 0-based BED end (exclusive).
        strand: ``"+"`` or ``"-"``; anything else is treated as ``"+"``.
        offset: Signed bp in TRANSCRIPT orientation (+ = downstream).

    Returns:
        ``(new_start, new_end)``, never negative and never degenerate.
    """
    if offset == 0:
        return bed_start, bed_end
    if strand == "-":
        new_start = max(0, bed_start - offset)
        new_end = max(bed_end, new_start + 1)
        return new_start, new_end
    new_end = max(1, bed_end + offset)
    new_start = min(bed_start, new_end - 1)
    return max(0, new_start), new_end


def cleavage_point(bed_start: int, bed_end: int, strand: str) -> int:
    """The single reported cleavage base of a PAS interval (0-based)."""
    return bed_start if strand == "-" else bed_end - 1


#: Which rows a signed offset is applied to under ``clip_seeded``.  The sign
#: decides, and the measurement is why: a POSITIVE offset is the coverage
#: correction (issue #72, ~+95 bp) and must not touch the clip-anchored tier
#: -- that is v2's ``skip_supported`` rule, and applying +95 to tier 1 costs
#: 46 points of P@10 -- while a NEGATIVE offset is the base-pair resolution
#: correction, which is a tier-1 quantity: tier 2's own base-pair-exact
#: precision is 0.0016 and moves by 0.0004 over the whole -8..+8 sweep.
def rows_for_offset(offset: int, strategy: str) -> str:
    """``"all"``, ``"tier1"`` or ``"tier2"`` -- which BED rows an offset moves."""
    if str(strategy) != "clip_seeded":
        return "all"
    if offset > 0:
        return "tier2"
    if offset < 0:
        return "tier1"
    return "all"


def rewrite_bed_cleavage_offset(path, offset: int, *,
                                rows: str = "all") -> int:
    """Rewrite a 6-column PAS BED in place, shifting 3' ends by a SIGNED offset.

    The signed generalisation of :func:`rewrite_bed_3prime_offset`.  For
    ``offset > 0`` with ``rows="tier2"`` it is byte-for-byte equivalent to
    ``rewrite_bed_3prime_offset(path, offset, skip_supported=True)``; that
    equivalence is pinned by
    ``tests/test_cleavage_offset_prime.py::test_the_positive_path_is_v2``.

    Args:
        path: Path to the BED (str or :class:`pathlib.Path`).
        offset: Signed bp, transcript orientation.  ``0`` is a no-op that
            leaves the file byte-identical and returns 0.
        rows: ``"all"``, ``"tier1"`` (BED score > 0) or ``"tier2"`` (score 0).
            Malformed rows are always passed through unchanged.

    Returns:
        Number of records whose coordinates changed.
    """
    import os as _os

    if offset == 0:
        return 0
    if rows not in ("all", "tier1", "tier2"):
        raise ValueError("rows must be all/tier1/tier2, got %r" % (rows,))

    shifted = 0
    tmp = f"{_os.fspath(path)}.offset.tmp"
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
            if rows != "all":
                try:
                    is_tier1 = float(parts[4]) > 0
                except ValueError:
                    is_tier1 = False
                if (rows == "tier1") != is_tier1:
                    dst.write(line)
                    continue
            new_start, new_end = shift_cleavage_point(start, end, parts[5], offset)
            if (new_start, new_end) != (start, end):
                shifted += 1
            parts[1] = str(new_start)
            parts[2] = str(new_end)
            dst.write("\t".join(parts) + "\n")
    _os.replace(tmp, path)
    return shifted


def estimate_clip_anchored_offset(support_paths) -> tuple[float, dict]:
    """Aggregate the caller's per-site ``clip_offset_mean`` into one number.

    The per-site column is the read-weighted mean of ``member - call`` over a
    tier-1 cluster's own clip positions, in transcript orientation; this
    function weights each site by its ``clip_reads`` so the run-level estimate
    is the read-weighted mean over every clip read in the run -- the same
    quantity a single-pass estimator would compute, recovered from a file the
    caller already writes.

    Tier-2 rows carry ``NA`` and are skipped, as are rows written without
    ``--pas-features on`` (no such column at all).

    Args:
        support_paths: Iterable of ``(path, strand)`` pairs -- the per-strand
            ``*.support.tsv`` sidecars.  ``strand`` is used only to report the
            per-strand split, which is not symmetric (see the module note).

    Returns:
        ``(offset, diagnostics)``.  ``offset`` is a float in bp (positive =
        the tails pinned bases DOWNSTREAM of the reported one) and is ``0.0``
        when nothing could be aggregated; ``diagnostics`` always carries
        ``n_sites`` / ``n_clip_reads`` / ``available`` so a caller can tell
        "measured zero" from "could not measure".
    """
    num = 0.0
    den = 0
    n_sites = 0
    per_strand: dict[str, list] = {}
    missing_column = False
    for path, strand in support_paths:
        try:
            fh = open(path)
        except (FileNotFoundError, TypeError):
            continue
        with fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                i_off = header.index("clip_offset_mean")
                i_reads = header.index("clip_reads")
                i_tier = header.index("tier")
            except ValueError:
                missing_column = True
                continue
            s_num, s_den, s_n = 0.0, 0, 0
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= i_off:
                    continue
                if parts[i_tier] != "1":
                    continue
                raw = parts[i_off]
                if raw == "NA":
                    continue
                try:
                    d = float(raw)
                    w = int(parts[i_reads])
                except ValueError:
                    continue
                if w <= 0:
                    continue
                s_num += d * w
                s_den += w
                s_n += 1
            num += s_num
            den += s_den
            n_sites += s_n
            if s_den:
                acc = per_strand.setdefault(strand, [0.0, 0, 0])
                acc[0] += s_num
                acc[1] += s_den
                acc[2] += s_n
    diag = {
        "available": bool(den),
        "n_sites": n_sites,
        "n_clip_reads": den,
        "missing_column": missing_column,
        "by_strand": {
            s: {"offset_bp": round(v[0] / v[1], 4), "n_sites": v[2],
                "n_clip_reads": v[1]}
            for s, v in sorted(per_strand.items()) if v[1]
        },
    }
    if not den:
        return 0.0, diag
    return num / den, diag


#: Sidecar column appended by ``--emit-inferred-cleavage on``.
INFERRED_CLEAVAGE_COLUMN = "inferred_cleavage"


def append_inferred_cleavage(bed_path, support_path,
                             offset_by_tier: dict) -> int:
    """Append ``inferred_cleavage`` to one per-strand ``*.support.tsv``.

    ``inferred_cleavage`` is the 0-based genomic coordinate this run's
    cleavage offset implies for the PAS: the reported cleavage base moved by
    ``offset_by_tier[tier]`` in transcript orientation.  It is REPORTED, never
    substituted -- ``pasbed.bed`` keeps whatever coordinate the run's
    ``--cleavage-offset`` produced, so with the default (``none``) the column
    equals the reported base and says so explicitly rather than by omission.

    The join is on the PAS id (BED column 4 == sidecar column 1), not on line
    order, so it survives any future reordering of either file.  Idempotent:
    a sidecar that already carries the column is left alone.

    Args:
        bed_path: The per-strand PAS BED whose coordinates are the source.
        support_path: Its sidecar (``paswrite.support_path_for(bed_path)``).
        offset_by_tier: ``{1: int, 2: int}`` signed bp, transcript oriented.

    Returns:
        Number of sidecar rows given a coordinate (0 if nothing was written).
    """
    import os as _os

    if not _os.path.exists(support_path) or not _os.path.exists(bed_path):
        return 0
    coords: dict[str, tuple[int, int, str]] = {}
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                coords[parts[3]] = (int(parts[1]), int(parts[2]), parts[5])
            except ValueError:
                continue

    tmp = f"{_os.fspath(support_path)}.infcleav.tmp"
    written = 0
    with open(support_path) as src, open(tmp, "w") as dst:
        header = src.readline().rstrip("\n")
        cols = header.split("\t")
        if INFERRED_CLEAVAGE_COLUMN in cols:
            _os.unlink(tmp)
            return 0
        try:
            i_tier = cols.index("tier")
        except ValueError:
            _os.unlink(tmp)
            return 0
        dst.write(header + "\t" + INFERRED_CLEAVAGE_COLUMN + "\n")
        for line in src:
            row = line.rstrip("\n")
            if not row:
                continue
            parts = row.split("\t")
            rec = coords.get(parts[0])
            if rec is None or len(parts) <= i_tier:
                dst.write(row + "\tNA\n")
                continue
            start, end, strand = rec
            try:
                off = int(offset_by_tier.get(int(parts[i_tier]), 0))
            except ValueError:
                off = 0
            base = cleavage_point(start, end, strand)
            dst.write("%s\t%d\n" % (row, max(0, base - off if strand == "-"
                                             else base + off)))
            written += 1
    _os.replace(tmp, support_path)
    return written
