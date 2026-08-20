"""Issue #72: data-driven 3' cleavage-site offset correction.

Covers the pure coordinate transform, the point-mode cleavage-site helper,
the (stub) data-driven estimator, and the in-place BED rewriter -- including
that offset 0 is byte-identical legacy behaviour.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from ema.countmatrix.cleavage_offset import (
    DEFAULT_CLEAVAGE_OFFSET,
    estimate_cleavage_offset,
    inferred_cleavage_site,
    rewrite_bed_3prime_offset,
    shift_3prime_end,
)


# --------------------------------------------------------------------------
# shift_3prime_end
# --------------------------------------------------------------------------

def test_plus_strand_extends_downstream_end():
    # + strand: 3' end is the larger coord -> end moves downstream (+).
    assert shift_3prime_end(1000, 1100, "+", 95) == (1000, 1195)


def test_minus_strand_extends_downstream_start():
    # - strand: 3' end is the smaller coord -> start moves downstream (-).
    assert shift_3prime_end(1000, 1100, "-", 95) == (905, 1100)


def test_minus_strand_clamps_at_zero():
    assert shift_3prime_end(40, 200, "-", 95) == (0, 200)


def test_zero_offset_is_noop_both_strands():
    assert shift_3prime_end(1000, 1100, "+", 0) == (1000, 1100)
    assert shift_3prime_end(1000, 1100, "-", 0) == (1000, 1100)


def test_negative_offset_is_noop():
    assert shift_3prime_end(1000, 1100, "+", -10) == (1000, 1100)


def test_unknown_strand_treated_as_forward():
    # Permissive fallback matches the rest of the pipeline.
    assert shift_3prime_end(1000, 1100, ".", 50) == (1000, 1150)


def test_5prime_end_preserved():
    # The peak 5' end (real R2 coverage) must not move.
    s_plus, _ = shift_3prime_end(1000, 1100, "+", 95)
    assert s_plus == 1000
    _, e_minus = shift_3prime_end(1000, 1100, "-", 95)
    assert e_minus == 1100


# --------------------------------------------------------------------------
# inferred_cleavage_site (point mode)
# --------------------------------------------------------------------------

def test_inferred_cleavage_point_plus():
    assert inferred_cleavage_site(1000, 1100, "+", 95) == 1195


def test_inferred_cleavage_point_minus():
    assert inferred_cleavage_site(1000, 1100, "-", 95) == 905


def test_inferred_cleavage_point_minus_clamped():
    assert inferred_cleavage_site(40, 200, "-", 95) == 0


def test_inferred_cleavage_point_zero_offset():
    assert inferred_cleavage_site(1000, 1100, "+", 0) == 1100
    assert inferred_cleavage_site(1000, 1100, "-", 0) == 1000


# --------------------------------------------------------------------------
# estimate_cleavage_offset (stub)
# --------------------------------------------------------------------------

def test_estimator_returns_default_constant():
    assert estimate_cleavage_offset() == DEFAULT_CLEAVAGE_OFFSET


def test_estimator_default_in_observed_window():
    # The fallback must sit inside the empirically observed +90..+105 window.
    assert 90 <= DEFAULT_CLEAVAGE_OFFSET <= 105


def test_estimator_respects_custom_fallback():
    assert estimate_cleavage_offset(fallback=88) == 88


def test_estimator_stub_ignores_profiles_for_now():
    # Signature is stable; profiles are accepted but currently unused.
    assert estimate_cleavage_offset(
        aataaa_profile=[0, 0, 1, 5, 1],
        a_fraction_profile=[0.1, 0.2, 0.43, 0.1],
    ) == DEFAULT_CLEAVAGE_OFFSET


# --------------------------------------------------------------------------
# rewrite_bed_3prime_offset
# --------------------------------------------------------------------------

def _write(path: Path, rows) -> None:
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def test_rewrite_shifts_per_strand(tmp_path):
    bed = tmp_path / "pas.bed"
    _write(bed, [
        ("chr1", 1000, 1100, "1", 0, "+"),
        ("chr1", 2000, 2100, "2", 0, "-"),
    ])
    n = rewrite_bed_3prime_offset(bed, 95)
    assert n == 2
    lines = [ln.split("\t") for ln in bed.read_text().splitlines()]
    # + strand: end extended
    assert lines[0][1:3] == ["1000", "1195"]
    # - strand: start extended downstream
    assert lines[1][1:3] == ["1905", "2100"]
    # name/score/strand columns untouched
    assert lines[0][3:] == ["1", "0", "+"]
    assert lines[1][3:] == ["2", "0", "-"]


def test_rewrite_zero_offset_is_byte_identical(tmp_path):
    bed = tmp_path / "pas.bed"
    original = "chr1\t1000\t1100\t1\t0\t+\nchr1\t2000\t2100\t2\t0\t-\n"
    bed.write_text(original)
    n = rewrite_bed_3prime_offset(bed, 0)
    assert n == 0
    assert bed.read_text() == original


def test_rewrite_passes_through_malformed_and_blank_lines(tmp_path):
    bed = tmp_path / "pas.bed"
    with open(bed, "w") as f:
        f.write("chr1\t1000\t1100\t1\t0\t+\n")
        f.write("\n")                      # blank
        f.write("chr1\tnope\t1100\t2\t0\t+\n")  # bad coord
        f.write("chr1\t500\t600\n")        # too few columns
    n = rewrite_bed_3prime_offset(bed, 50)
    assert n == 1
    out = bed.read_text().splitlines()
    assert out[0].split("\t")[1:3] == ["1000", "1150"]
    assert out[1] == ""
    assert out[2] == "chr1\tnope\t1100\t2\t0\t+"
    assert out[3] == "chr1\t500\t600"


def test_rewrite_minus_strand_clamps(tmp_path):
    bed = tmp_path / "pas.bed"
    _write(bed, [("chr1", 40, 200, "1", 0, "-")])
    rewrite_bed_3prime_offset(bed, 95)
    assert bed.read_text().split("\t")[1:3] == ["0", "200"]


def test_rewrite_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        rewrite_bed_3prime_offset(tmp_path / "nope.bed", 95)
