"""``read_check`` was reordered for speed (strand/mapped/length checks before
the CB tag lookup) and now interns the ``"<RG>_<CB>"`` composite.

Both are supposed to be pure optimisations, so this module pins the result
against a verbatim copy of the pre-change implementation over an exhaustive
grid of read shapes: wrong strand, unmapped, CIGAR-less, ignored contig,
too-long / too-short spans, missing or GEM-suffixed barcodes, missing RG,
and underscore-bearing RGs.
"""
from __future__ import annotations

import itertools

import pytest

from ema.countmatrix import read as read_mod
from ema.countmatrix.read import read_check, reset_composite_cache

BARCODE = "AAAACCCCGGGGTTTT"


def _reference_read_check(read, direction, barcode, barcode_len, seq_len,
                          ignore_chro, *, sample_id=None):
    """Verbatim pre-optimisation body (order: CB, mapped, RG, strand, chro, len)."""
    try:
        cb = read.get_tag(barcode)
        if len(cb) != barcode_len:
            dash = cb.rfind("-")
            if dash == barcode_len and cb[dash + 1:].isdigit():
                cb = cb[:dash]
        if len(cb) != barcode_len:
            return 0, 0, 0, 0, 0
    except KeyError:
        return 0, 0, 0, 0, 0

    if read.is_unmapped or read.reference_end is None or read.reference_name is None:
        return 0, 0, 0, 0, 0

    read_chro, read_start, read_end, read_strand = (
        read.reference_name, read.reference_start, read.reference_end, read.is_reverse,
    )
    try:
        rg = read.get_tag('RG')
    except (KeyError, ValueError):
        rg = sample_id if sample_id is not None else read_mod._default_sample_id

    if direction != read_strand:
        return 0, 0, 0, 0, 0
    if read_chro in ignore_chro:
        return 0, 0, 0, 0, 0
    if read_end - read_start > seq_len:
        return 0, 0, 0, 0, 0
    elif read_end - read_start < seq_len:
        read_end = read_start + seq_len
    return read_chro, read_start, read_end, read_strand, f"{rg}_{cb}"


class _Read:
    def __init__(self, cb, chrom, start, end, is_reverse, rg, is_unmapped):
        self._cb, self._rg = cb, rg
        self.reference_name = chrom
        self.reference_start = start
        self.reference_end = end
        self.is_reverse = is_reverse
        self.is_unmapped = is_unmapped

    def get_tag(self, tag):
        if tag == "CB":
            if self._cb is None:
                raise KeyError(tag)
            return self._cb
        if tag == "RG":
            if self._rg is None:
                raise KeyError(tag)
            return self._rg
        raise KeyError(tag)


SEQ_LEN = 150
IGNORE = ("MT", "mt")
_GRID = list(itertools.product(
    [BARCODE, BARCODE + "-1", BARCODE + "-x", "AAAA", None],   # CB variants
    ["chr1", "MT", None],                                      # contig
    [(100, 250), (100, 200), (100, 400), (100, None)],         # span
    [False, True],                                             # is_reverse
    ["s1", "sampleA_rep1", None],                              # RG
    [False, True],                                             # is_unmapped
))


@pytest.fixture(autouse=True)
def _clean_cache():
    reset_composite_cache()
    yield
    reset_composite_cache()


@pytest.mark.parametrize("direction", [False, True])
def test_matches_reference_implementation_on_every_shape(direction):
    seen_accepted = 0
    for cb, chrom, (start, end), is_reverse, rg, is_unmapped in _GRID:
        read = _Read(cb, chrom, start, end, is_reverse, rg, is_unmapped)
        kwargs = dict(barcode="CB", barcode_len=16, seq_len=SEQ_LEN, ignore_chro=IGNORE)
        got = read_check(read, direction=direction, sample_id="ds", **kwargs)
        want = _reference_read_check(read, direction, sample_id="ds", **kwargs)
        assert got == want, (cb, chrom, start, end, is_reverse, rg, is_unmapped)
        seen_accepted += got[0] != 0
    assert seen_accepted, "the grid never produced an accepted read"


def test_gem_suffix_and_short_span_padding_survive():
    read = _Read(BARCODE + "-1", "chr1", 100, 200, False, "s1", False)
    chro, start, end, strand, cb = read_check(
        read, direction=False, barcode="CB", barcode_len=16,
        seq_len=SEQ_LEN, ignore_chro=(),
    )
    assert (chro, start, end, strand) == ("chr1", 100, 250, False)
    assert cb == f"s1_{BARCODE}"


def test_composite_is_interned_per_rg_cb_pair():
    """One string object per cell, not one per read (the memory point)."""
    reads = [_Read(BARCODE, "chr1", 100 + i, 250 + i, False, "s1", False) for i in range(5)]
    cbs = [read_check(r, direction=False, barcode="CB", barcode_len=16,
                      seq_len=SEQ_LEN, ignore_chro=())[4] for r in reads]
    assert all(c == cbs[0] for c in cbs)
    assert all(c is cbs[0] for c in cbs)
    other = read_check(_Read(BARCODE, "chr1", 100, 250, False, "s2", False),
                       direction=False, barcode="CB", barcode_len=16,
                       seq_len=SEQ_LEN, ignore_chro=())[4]
    assert other == f"s2_{BARCODE}" and other is not cbs[0]
    reset_composite_cache()
    again = read_check(reads[0], direction=False, barcode="CB", barcode_len=16,
                       seq_len=SEQ_LEN, ignore_chro=())[4]
    assert again == cbs[0] and again is not cbs[0]


def test_sample_id_fallback_still_applies_without_rg(monkeypatch):
    read = _Read(BARCODE, "chr1", 100, 250, False, None, False)
    kwargs = dict(barcode="CB", barcode_len=16, seq_len=SEQ_LEN, ignore_chro=())
    assert read_check(read, direction=False, sample_id="explicit", **kwargs)[4] == \
        f"explicit_{BARCODE}"
    reset_composite_cache()
    monkeypatch.setattr(read_mod, "_default_sample_id", "singleton")
    assert read_check(read, direction=False, **kwargs)[4] == f"singleton_{BARCODE}"
