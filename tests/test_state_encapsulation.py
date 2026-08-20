"""Verify encapsulated state instances don't share mutation.

Phase 1 regression tests for:
- BarcodeIndex instance isolation
- PeakCallingState instance isolation
- Backward-compat module-level shim functions
- matrix_write() backward-compat and explicit-index paths
- read_check() sample_id parameter
"""
from __future__ import annotations

import io

import pytest

from ema.countmatrix.indexing import BarcodeIndex
from ema.countmatrix.peak_state import PeakCallingState


# ---------------------------------------------------------------------------
# BarcodeIndex isolation
# ---------------------------------------------------------------------------


def test_two_barcode_indexes_independent() -> None:
    """Two BarcodeIndex instances must not share any state."""
    a = BarcodeIndex()
    b = BarcodeIndex()
    a.get_index("sample1_AAAACCCCGGGGTTTT")
    a.get_index("sample1_TTTTGGGGCCCCAAAA")
    assert len(a.mapping) == 2
    assert len(b.mapping) == 0


def test_barcode_index_reset_isolated() -> None:
    """Resetting one BarcodeIndex must not affect another."""
    a = BarcodeIndex()
    b = BarcodeIndex()
    a.get_index("sample1_AAAACCCCGGGGTTTT")
    b.get_index("sample1_AAAACCCCGGGGTTTT")
    a.reset()
    assert len(a.mapping) == 0
    assert len(b.mapping) == 1  # b unaffected


def test_barcode_index_columns_independent() -> None:
    """Column assignments must be independent across instances."""
    a = BarcodeIndex()
    b = BarcodeIndex()
    # Add one entry to a first so its counter is ahead of b's
    a.get_index("sample1_AAAACCCCGGGGTTTT")
    # Both b and a should assign column 1 to the *first* CB they see
    col_a2 = a.get_index("sample1_TTTTGGGGCCCCAAAA")
    col_b1 = b.get_index("sample1_TTTTGGGGCCCCAAAA")
    assert col_a2 == 2  # a already had one entry
    assert col_b1 == 1  # b starts fresh


# ---------------------------------------------------------------------------
# PeakCallingState isolation
# ---------------------------------------------------------------------------


def test_two_peak_states_independent() -> None:
    """Two PeakCallingState instances must have independent counters."""
    s1 = PeakCallingState()
    s2 = PeakCallingState()
    s1.bump_pasnumber()
    s1.bump_pasnumber()
    s1.bump_pasnumber()
    assert s1.pasnumber == 3
    assert s2.pasnumber == 0


def test_peak_state_bump_returns_new_value() -> None:
    """bump_pasnumber() must return the post-increment value."""
    s = PeakCallingState()
    assert s.bump_pasnumber() == 1
    assert s.bump_pasnumber() == 2
    assert s.bump_pasnumber() == 3
    assert s.pasnumber == 3


def test_reset_clears_state() -> None:
    """reset() must bring the counter back to 0."""
    s = PeakCallingState()
    s.bump_pasnumber()
    s.reset()
    assert s.pasnumber == 0


def test_reset_then_bump_restarts_from_one() -> None:
    """After reset(), bump_pasnumber() must return 1."""
    s = PeakCallingState()
    for _ in range(10):
        s.bump_pasnumber()
    s.reset()
    assert s.bump_pasnumber() == 1


def test_peak_state_initial_value_configurable() -> None:
    """PeakCallingState can be seeded with a non-zero starting value."""
    s = PeakCallingState(pasnumber=42)
    assert s.bump_pasnumber() == 43


# ---------------------------------------------------------------------------
# Module-level shim functions (backward compat)
# ---------------------------------------------------------------------------


def test_module_shims_still_work_for_backward_compat() -> None:
    """Module-level indexing / reset_index / get_mapping must still work."""
    from ema.countmatrix.indexing import get_mapping, indexing, reset_index

    reset_index()
    col = indexing("sample1_AAAACCCCGGGGTTTT")
    assert col == 1
    assert len(get_mapping()) == 1


def test_module_shim_reset_clears_mapping() -> None:
    """reset_index() must clear the module-level singleton."""
    from ema.countmatrix.indexing import get_mapping, indexing, reset_index

    reset_index()
    indexing("sample1_AAAACCCCGGGGTTTT")
    assert len(get_mapping()) == 1
    reset_index()
    assert len(get_mapping()) == 0


def test_get_default_index_returns_singleton() -> None:
    """get_default_index() must return the same object as the module singleton."""
    from ema.countmatrix.indexing import _index, get_default_index

    assert get_default_index() is _index


# ---------------------------------------------------------------------------
# matrix_write — backward compat (no index arg)
# ---------------------------------------------------------------------------


def test_matrix_write_backward_compat_no_index_arg() -> None:
    """matrix_write without explicit index should use the module singleton."""
    from ema.countmatrix.indexing import reset_index
    from ema.countmatrix.paswrite import matrix_write

    reset_index()
    output = io.StringIO()
    matrix_write(
        {"sample1_AAAACCCCGGGGTTTT": 5, "sample1_TTTTGGGGCCCCAAAA": 3},
        pasnumber=42,
        output=output,
    )
    lines = [l for l in output.getvalue().strip().split("\n") if l]
    assert len(lines) == 2


def test_matrix_write_empty_cb_dict_no_output() -> None:
    """matrix_write with an empty cb_dict should write nothing."""
    from ema.countmatrix.paswrite import matrix_write

    output = io.StringIO()
    matrix_write({}, pasnumber=1, output=output)
    assert output.getvalue() == ""


# ---------------------------------------------------------------------------
# matrix_write — explicit index does NOT touch the module singleton
# ---------------------------------------------------------------------------


def test_matrix_write_with_explicit_index() -> None:
    """Explicit index argument must leave the module singleton untouched."""
    from ema.countmatrix.indexing import get_mapping, reset_index
    from ema.countmatrix.paswrite import matrix_write

    reset_index()
    explicit = BarcodeIndex()
    output = io.StringIO()
    matrix_write(
        {"sample1_AAAACCCCGGGGTTTT": 5},
        pasnumber=42,
        output=output,
        index=explicit,
    )
    # Module singleton must be untouched
    assert len(get_mapping()) == 0
    # Explicit index must have received the assignment
    assert len(explicit.mapping) == 1


def test_matrix_write_explicit_index_output_format() -> None:
    """matrix_write output lines must follow 'row col count' MMX format."""
    from ema.countmatrix.paswrite import matrix_write

    explicit = BarcodeIndex()
    output = io.StringIO()
    matrix_write(
        {"s1_AAAACCCCGGGGTTTT": 7},
        pasnumber=3,
        output=output,
        index=explicit,
    )
    line = output.getvalue().strip()
    parts = line.split()
    assert len(parts) == 3
    assert parts[0] == "3"   # pasnumber (row)
    assert parts[2] == "7"   # count
    assert int(parts[1]) >= 1  # column index


# ---------------------------------------------------------------------------
# read_check — sample_id keyword argument
# ---------------------------------------------------------------------------
# ema.countmatrix.read imports ema.config which calls argparse.parse_args()
# at module-load time.  Under pytest sys.argv contains test filenames that
# argparse rejects.  We patch sys.argv to ["ema"] before importing the
# module so the import succeeds without side-effects.


class _MockRead:
    """Minimal pysam.AlignedSegment mock for testing read_check."""

    def __init__(
        self,
        cb: str,
        chrom: str,
        start: int,
        end: int,
        is_reverse: bool,
        rg: str | None = None,
        cb_len: int = 16,
        is_unmapped: bool = False,
    ) -> None:
        self._cb = cb
        self.reference_name = chrom
        self.reference_start = start
        self.reference_end = end
        self.is_reverse = is_reverse
        self._rg = rg
        self._cb_len = cb_len
        # pysam.AlignedSegment always exposes is_unmapped; read_check reads it
        # to skip unmapped/CIGAR-less records (reference_end is None there).
        self.is_unmapped = is_unmapped

    def get_tag(self, tag: str) -> str:
        if tag == "CB":
            return self._cb
        if tag == "RG":
            if self._rg is None:
                raise KeyError(tag)
            return self._rg
        raise KeyError(tag)


def _make_read_check():
    """Import read_check with a patched sys.argv so argparse does not crash."""
    import sys
    from unittest.mock import patch

    with patch.object(sys, "argv", ["ema"]):
        # Force a fresh import if the module was not yet loaded
        import importlib

        import ema.countmatrix.read as _read_mod

        importlib.reload(_read_mod)
        return _read_mod.read_check, _read_mod._default_sample_id


# Shared keyword args for read_check that override the config-derived defaults
# (variable_config values may be None when pytest runs outside the CLI context)
_READ_CHECK_KWARGS: dict = dict(
    barcode="CB",
    barcode_len=16,
    seq_len=150,
    ignore_chro=set(),
)


def test_read_check_uses_rg_tag_when_present() -> None:
    """read_check must use the RG tag from the read when it exists."""
    read_check, _ = _make_read_check()

    read = _MockRead(
        cb="AAAACCCCGGGGTTTT",
        chrom="chr1",
        start=100,
        end=250,  # 150 bp exactly — matches seq_len
        is_reverse=False,
        rg="dataset_A",
    )
    result = read_check(read, direction=False, **_READ_CHECK_KWARGS)
    _, _, _, _, cb = result
    assert cb.startswith("dataset_A_"), f"Expected dataset_A_ prefix, got {cb!r}"


def test_read_check_uses_explicit_sample_id_when_no_rg() -> None:
    """read_check must use sample_id kwarg as RG fallback (not module global)."""
    read_check, default_sample_id = _make_read_check()

    # Ensure the module global is different from our explicit value
    assert default_sample_id != "explicit_sample"

    read = _MockRead(
        cb="AAAACCCCGGGGTTTT",
        chrom="chr1",
        start=100,
        end=250,
        is_reverse=False,
        rg=None,  # no RG tag — triggers fallback
    )
    result = read_check(
        read,
        direction=False,
        **_READ_CHECK_KWARGS,
        sample_id="explicit_sample",
    )
    _, _, _, _, cb = result
    assert cb.startswith("explicit_sample_"), f"Expected explicit_sample_ prefix, got {cb!r}"


def test_read_check_falls_back_to_module_global_when_no_sample_id() -> None:
    """Without sample_id kwarg, read_check must use _default_sample_id."""
    read_check, default_sample_id = _make_read_check()

    read = _MockRead(
        cb="AAAACCCCGGGGTTTT",
        chrom="chr1",
        start=100,
        end=250,
        is_reverse=False,
        rg=None,
    )
    result = read_check(read, direction=False, **_READ_CHECK_KWARGS)
    _, _, _, _, cb = result
    assert cb.startswith(f"{default_sample_id}_"), (
        f"Expected {default_sample_id}_ prefix, got {cb!r}"
    )


def test_read_check_invalid_cb_returns_zeros() -> None:
    """read_check must return (0,0,0,0,0) when CB tag is missing."""
    read_check, _ = _make_read_check()

    class _NoCBRead(_MockRead):
        def get_tag(self, tag: str) -> str:
            if tag == "CB":
                raise KeyError(tag)
            return super().get_tag(tag)

    read = _NoCBRead(
        cb="AAAACCCCGGGGTTTT",
        chrom="chr1",
        start=100,
        end=250,
        is_reverse=False,
    )
    result = read_check(read, direction=False, **_READ_CHECK_KWARGS)
    assert result == (0, 0, 0, 0, 0)


def test_read_check_wrong_direction_returns_zeros() -> None:
    """read_check must skip reads on the wrong strand."""
    read_check, _ = _make_read_check()

    read = _MockRead(
        cb="AAAACCCCGGGGTTTT",
        chrom="chr1",
        start=100,
        end=250,
        is_reverse=True,  # reverse
        rg="s1",
    )
    result = read_check(
        read,
        direction=False,  # asking for forward — should return zeros
        **_READ_CHECK_KWARGS,
    )
    assert result == (0, 0, 0, 0, 0)
