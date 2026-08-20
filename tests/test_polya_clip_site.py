"""Unit tests for ``ema.countmatrix.polya.clip_site`` (Stage 1 Phase 1).

``clip_site`` is the strand-aware terminal poly(A) soft-clip detector that
gives PeakATail its first read-level PAS evidence.  Its exact semantics are
load-bearing: the measured 1.152% clip rate / 92x wrong-end specificity /
73.9%-within-100bp statistics on the PBMC BAM (manuscript/10_caller_fix_plan.md
section 1) were produced by this rule set, and Phase 2's F1 depends on the
returned coordinate being the cleavage site, not the clip's outer edge.

Rules under test:
  * forward read -> clip must be the LAST cigar op, bases must be A-rich,
    site = reference_end - 1
  * reverse read -> clip must be the FIRST cigar op, bases must be T-rich
    (revcomp of the poly(A) tail), site = reference_start
  * the homopolymer run ADJACENT to the alignment boundary must itself be
    >= min_clip (this is what rejects templated A-rich clips)
  * purity is a fraction over the WHOLE clipped segment
  * no soft clip, wrong-end clip, unmapped/CIGAR-less reads -> None
"""
from __future__ import annotations

import pytest

from ema.countmatrix.polya import clip_site

CSOFT = 4   # BAM_CSOFT_CLIP
CMATCH = 0  # BAM_CMATCH
CHARD = 5   # BAM_CHARD_CLIP
CREF_SKIP = 3  # BAM_CREF_SKIP (N)


class _Read:
    """Minimal ``pysam.AlignedSegment`` stand-in."""

    def __init__(self, cigartuples, seq, start=1000, end=None,
                 is_reverse=False, is_unmapped=False, umi=None):
        self.cigartuples = cigartuples
        self.query_sequence = seq
        self.reference_start = start
        self.is_reverse = is_reverse
        self.is_unmapped = is_unmapped
        self._umi = umi
        if end is not None:
            self.reference_end = end
        else:
            span = sum(
                n for op, n in (cigartuples or [])
                if op in (CMATCH, 2, CREF_SKIP, 7, 8)
            )
            self.reference_end = start + span

    def get_tag(self, tag):
        if tag == "UB" and self._umi is not None:
            return self._umi
        raise KeyError(tag)


def _fwd(clip_seq, aligned_len=60, start=1000):
    """Forward read: aligned_len M then a terminal soft clip."""
    return _Read(
        [(CMATCH, aligned_len), (CSOFT, len(clip_seq))],
        "C" * aligned_len + clip_seq,
        start=start,
    )


def _rev(clip_seq, aligned_len=60, start=1000):
    """Reverse read: leading soft clip then aligned_len M."""
    return _Read(
        [(CSOFT, len(clip_seq)), (CMATCH, aligned_len)],
        clip_seq + "C" * aligned_len,
        start=start,
        is_reverse=True,
    )


# ---------------------------------------------------------------------------
# Plus strand
# ---------------------------------------------------------------------------

def test_plus_strand_pure_clip_returns_reference_end_minus_one():
    read = _fwd("A" * 12, aligned_len=60, start=1000)
    # reference_end = 1060 -> cleavage site is the last ALIGNED base, 1059
    assert read.reference_end == 1060
    assert clip_site(read) == 1059


def test_plus_strand_site_is_independent_of_clip_length():
    """The whole point: the clip is non-templated, so lengthening it must
    not move the inferred cleavage coordinate downstream."""
    short = clip_site(_fwd("A" * 6))
    long = clip_site(_fwd("A" * 40))
    assert short == long == 1059


def test_plus_strand_impure_but_above_threshold_accepted():
    # 9 A + 1 G = 0.9 purity, and the A-run adjacent to the boundary is 9
    assert clip_site(_fwd("A" * 9 + "G")) == 1059


def test_plus_strand_below_purity_rejected():
    # 6 A + 4 G = 0.6 purity < 0.8
    assert clip_site(_fwd("A" * 6 + "G" * 4)) is None


def test_plus_strand_wrong_end_clip_rejected():
    """A poly(A) run on the 5' end of a forward read is the wrong-end control
    that bounds alignment artefacts (0.0125% -> 92x specificity)."""
    read = _Read([(CSOFT, 12), (CMATCH, 60)], "A" * 12 + "C" * 60, start=1000)
    assert clip_site(read) is None


def test_plus_strand_t_rich_clip_rejected():
    assert clip_site(_fwd("T" * 12)) is None


# ---------------------------------------------------------------------------
# Minus strand
# ---------------------------------------------------------------------------

def test_minus_strand_pure_clip_returns_reference_start():
    read = _rev("T" * 12, aligned_len=60, start=1000)
    assert clip_site(read) == 1000


def test_minus_strand_a_rich_clip_rejected():
    """On a reverse read the poly(A) tail appears as T; an A-run there is the
    wrong-end control."""
    assert clip_site(_rev("A" * 12)) is None


def test_minus_strand_wrong_end_clip_rejected():
    read = _Read([(CMATCH, 60), (CSOFT, 12)], "C" * 60 + "T" * 12,
                 start=1000, is_reverse=True)
    assert clip_site(read) is None


def test_minus_strand_run_must_be_adjacent_to_alignment():
    """For a reverse read the run adjacent to the alignment is the clip's
    TRAILING T-run: 'TTTTTTGGGG' has 6 T but 0 adjacent to the boundary."""
    assert clip_site(_rev("T" * 6 + "G" * 4)) is None
    # ...whereas the mirrored sequence has its run flush against the alignment
    assert clip_site(_rev("G" * 2 + "T" * 8)) == 1000


def test_plus_strand_run_must_be_adjacent_to_alignment():
    """For a forward read the run adjacent to the alignment is the clip's
    LEADING A-run."""
    assert clip_site(_fwd("G" * 4 + "A" * 6)) is None
    assert clip_site(_fwd("A" * 8 + "G" * 2)) == 1059


# ---------------------------------------------------------------------------
# Boundaries
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n,expected", [(5, None), (6, 1059), (7, 1059)])
def test_min_clip_boundary(n, expected):
    assert clip_site(_fwd("A" * n)) == expected


@pytest.mark.parametrize("min_clip,expected", [(6, 1059), (8, 1059), (9, None)])
def test_min_clip_parameter_is_honoured(min_clip, expected):
    assert clip_site(_fwd("A" * 8), min_clip=min_clip) == expected


def test_purity_boundary_exactly_at_threshold_accepted():
    # 8 A + 2 G over 10 bases = exactly 0.8 -> accepted (>= threshold)
    assert clip_site(_fwd("A" * 8 + "G" * 2), min_purity=0.8) == 1059


def test_purity_boundary_just_below_threshold_rejected():
    # 7 A + 3 G = 0.7 < 0.8
    assert clip_site(_fwd("A" * 7 + "G" * 3), min_purity=0.8) is None


def test_min_purity_parameter_is_honoured():
    read = _fwd("A" * 7 + "G" * 3)  # 0.7 purity, adjacent run 7
    assert clip_site(read, min_purity=0.8) is None
    assert clip_site(read, min_purity=0.7) == 1059


# ---------------------------------------------------------------------------
# CIGAR / read-state guards
# ---------------------------------------------------------------------------

def test_cigar_without_softclip_returns_none():
    assert clip_site(_Read([(CMATCH, 91)], "C" * 91)) is None


def test_spliced_read_with_terminal_clip_is_accepted_and_uses_reference_end():
    """A junction (N) read must still work — reference_end already accounts
    for the skipped intron, so no manual CIGAR walk is needed."""
    read = _Read(
        [(CMATCH, 30), (CREF_SKIP, 500), (CMATCH, 30), (CSOFT, 10)],
        "C" * 60 + "A" * 10,
        start=1000,
    )
    assert read.reference_end == 1560
    assert clip_site(read) == 1559


def test_hard_clip_instead_of_soft_clip_returns_none():
    """Hard clips carry no SEQ, so the evidence is gone — must not be
    mistaken for a soft clip."""
    read = _Read([(CMATCH, 60), (CHARD, 12)], "C" * 60, start=1000)
    assert clip_site(read) is None


def test_empty_cigar_returns_none():
    assert clip_site(_Read(None, "AAAAAAAAAA")) is None
    assert clip_site(_Read([], "AAAAAAAAAA")) is None


def test_missing_sequence_returns_none():
    read = _Read([(CMATCH, 60), (CSOFT, 12)], None, start=1000)
    assert clip_site(read) is None


def test_unmapped_read_with_cigar_and_none_reference_end_returns_none():
    """read_check's unmapped guard runs first in the caller, but clip_site
    must be safe standalone: a CIGAR-less unmapped read yields None, and a
    read whose reference_end is None must not raise."""
    unmapped = _Read(None, "A" * 20, is_unmapped=True)
    assert clip_site(unmapped) is None

    weird = _Read([(CMATCH, 60), (CSOFT, 12)], "C" * 60 + "A" * 12)
    weird.reference_end = None
    assert clip_site(weird) is None


def test_reverse_unmapped_guard_uses_reference_start_only():
    """A reverse read only needs reference_start, so it must not be affected
    by a missing reference_end."""
    read = _rev("T" * 12)
    read.reference_end = None
    assert clip_site(read) == 1000
