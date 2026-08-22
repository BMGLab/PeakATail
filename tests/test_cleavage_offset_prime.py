"""peakAtail-prime TASK E item 1 — the clip-anchored cleavage offset.

Three separable things are pinned here:

1. ``--cleavage-offset`` parses ``none`` / ``auto`` / a signed int, and every
   spelling v2 accepted still means what it meant.
2. a SIGNED offset moves the right rows: positive is the coverage correction
   and, under ``clip_seeded``, must not touch the clip-anchored tier -- the
   rule v2 spelled ``skip_supported`` -- while negative is the base-pair
   resolution correction and must move only that tier.
3. the per-library estimate is aggregated from the caller's own
   ``clip_offset_mean`` column, and ``inferred_cleavage`` REPORTS a corrected
   coordinate without moving ``pasbed.bed``.
"""
from __future__ import annotations

import pytest

from ema.countmatrix.cleavage_offset import (
    append_inferred_cleavage,
    cleavage_point,
    estimate_clip_anchored_offset,
    parse_cleavage_offset,
    rewrite_bed_3prime_offset,
    rewrite_bed_cleavage_offset,
    rows_for_offset,
    shift_3prime_end,
    shift_cleavage_point,
)

ROWS = [
    ("1", 1000, 1001, "t1_plus", 7, "+"),    # clip-anchored tier-1 point
    ("1", 2000, 2001, "t1_minus", 3, "-"),
    ("1", 3000, 3300, "t2_plus", 0, "+"),    # coverage-only interval
    ("1", 4000, 4300, "t2_minus", 0, "-"),
]


def _write(tmp_path, rows=ROWS, name="pas.bed"):
    p = tmp_path / name
    p.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return p


def _read(p):
    return [l.split("\t") for l in p.read_text().splitlines()]


# --------------------------------------------------------------------------
# 1. the spec parser
# --------------------------------------------------------------------------
@pytest.mark.parametrize("spec,expected", [
    (None, ("none", 0)),
    ("", ("none", 0)),
    ("none", ("none", 0)),
    ("NONE", ("none", 0)),
    ("off", ("none", 0)),
    (0, ("none", 0)),              # v2 spelled "no shift" as the integer 0
    ("0", ("none", 0)),
    ("auto", ("auto", 0)),
    (95, ("const", 95)),           # v2's documented constant, still an int
    ("95", ("const", 95)),
    ("-2", ("const", -2)),
    (-1, ("const", -1)),
])
def test_the_spec_parses_every_v2_spelling_and_the_two_new_ones(spec, expected):
    assert parse_cleavage_offset(spec) == expected


@pytest.mark.parametrize("bad", ["ninety-five", "auto2", "1.5", True, "yes"])
def test_a_typo_is_an_error_not_silently_no_shift(bad):
    with pytest.raises(ValueError):
        parse_cleavage_offset(bad)


# --------------------------------------------------------------------------
# 2. the signed point shift
# --------------------------------------------------------------------------
def test_zero_returns_the_interval_untouched():
    assert shift_cleavage_point(1000, 1100, "+", 0) == (1000, 1100)
    assert shift_cleavage_point(1000, 1100, "-", 0) == (1000, 1100)


def test_a_positive_offset_is_exactly_v2s_shift():
    for start, end, strand in ((1000, 1100, "+"), (1000, 1100, "-"),
                               (5, 6, "+"), (5, 6, "-")):
        for off in (1, 25, 95, 150):
            assert (shift_cleavage_point(start, end, strand, off)
                    == shift_3prime_end(start, end, strand, off)), (
                        (start, end, strand, off))


def test_a_negative_offset_moves_a_one_bp_pas_upstream():
    # '+': the reported base is end-1; upstream is a smaller coordinate.
    assert shift_cleavage_point(1000, 1001, "+", -2) == (998, 999)
    assert cleavage_point(998, 999, "+") == 998
    # '-': the reported base is start; upstream is a LARGER coordinate.
    assert shift_cleavage_point(2000, 2001, "-", -2) == (2002, 2003)
    assert cleavage_point(2002, 2003, "-") == 2002


def test_a_negative_offset_never_inverts_a_wide_interval():
    assert shift_cleavage_point(3000, 3300, "+", -2) == (3000, 3298)
    assert shift_cleavage_point(3000, 3300, "-", -2) == (3002, 3300)
    # ...even when it is narrower than the shift
    s, e = shift_cleavage_point(10, 11, "+", -50)
    assert e > s >= 0


def test_coordinates_never_go_negative():
    s, e = shift_cleavage_point(3, 4, "+", -20)
    assert s >= 0 and e >= 1
    s, e = shift_cleavage_point(3, 4, "-", 100)
    assert s == 0


# --------------------------------------------------------------------------
# 3. which rows a signed offset moves
# --------------------------------------------------------------------------
def test_the_sign_decides_the_tier_under_clip_seeded():
    assert rows_for_offset(95, "clip_seeded") == "tier2"
    assert rows_for_offset(-2, "clip_seeded") == "tier1"
    assert rows_for_offset(0, "clip_seeded") == "all"


def test_a_coverage_strategy_has_only_one_tier_so_it_moves_everything():
    for strategy in ("lambda_gradient", "original", "sierra_iterative"):
        assert rows_for_offset(95, strategy) == "all"
        assert rows_for_offset(-2, strategy) == "all"


def test_the_positive_path_is_v2(tmp_path):
    """``rewrite_bed_cleavage_offset(+N, rows="tier2")`` must be byte-for-byte
    what v2's ``rewrite_bed_3prime_offset(+N, skip_supported=True)`` produced."""
    a = _write(tmp_path, name="a.bed")
    b = _write(tmp_path, name="b.bed")
    n_new = rewrite_bed_cleavage_offset(a, 90, rows="tier2")
    n_v2 = rewrite_bed_3prime_offset(b, 90, skip_supported=True)
    assert n_new == n_v2 == 2
    assert a.read_bytes() == b.read_bytes()


def test_the_positive_path_with_rows_all_is_v2s_unskipped_rewrite(tmp_path):
    a = _write(tmp_path, name="a.bed")
    b = _write(tmp_path, name="b.bed")
    assert rewrite_bed_cleavage_offset(a, 90, rows="all") == 4
    assert rewrite_bed_3prime_offset(b, 90) == 4
    assert a.read_bytes() == b.read_bytes()


def test_a_negative_offset_moves_the_clip_tier_only(tmp_path):
    p = _write(tmp_path)
    before = _read(p)
    n = rewrite_bed_cleavage_offset(p, -2, rows="tier1")
    out = _read(p)
    assert n == 2
    assert int(out[0][1]) == 998 and int(out[0][2]) == 999      # + tier1
    assert int(out[1][1]) == 2002 and int(out[1][2]) == 2003    # - tier1
    assert out[2] == before[2] and out[3] == before[3]          # tier2 frozen


def test_zero_leaves_the_file_byte_identical(tmp_path):
    p = _write(tmp_path)
    before = p.read_bytes()
    assert rewrite_bed_cleavage_offset(p, 0, rows="all") == 0
    assert p.read_bytes() == before


def test_malformed_rows_pass_through(tmp_path):
    p = tmp_path / "odd.bed"
    p.write_text("\n1\tnotanint\t10\tx\t0\t+\nshort\n"
                 "1\t100\t101\tok\t5\t+\n")
    before = p.read_text().splitlines()
    rewrite_bed_cleavage_offset(p, -3, rows="tier1")
    after = p.read_text().splitlines()
    assert after[:3] == before[:3]
    assert after[3].split("\t")[1] == "97"


def test_an_unknown_row_selector_is_an_error(tmp_path):
    with pytest.raises(ValueError):
        rewrite_bed_cleavage_offset(_write(tmp_path), -2, rows="tier9")


# --------------------------------------------------------------------------
# 4. the per-library estimate
# --------------------------------------------------------------------------
HEADER = ("pas_id\tclip_reads\tclip_umis\tclip_reads_f3844\tclip_umis_f3844"
          "\twindow_reads\ttier\tclip_positions\tclip_span\tclip_offset_mean")


def _support(tmp_path, name, rows, header=HEADER):
    p = tmp_path / name
    p.write_text(header + "\n" + "".join(r + "\n" for r in rows))
    return p


def test_the_estimate_is_read_weighted_over_tier_one_rows(tmp_path):
    pos = _support(tmp_path, "pos.support.tsv", [
        "1\t10\t5\t10\t5\t100\t1\t3\t9\t1.00",     # weight 10
        "2\t0\t0\t0\t0\t50\t2\t0\t0\tNA",          # tier 2: ignored
    ])
    neg = _support(tmp_path, "neg.support.tsv", [
        "3\t30\t9\t30\t9\t200\t1\t4\t12\t-3.00",   # weight 30
    ])
    off, diag = estimate_clip_anchored_offset([(pos, "+"), (neg, "-")])
    # (10*1 + 30*-3) / 40 = -2.0 -- NOT the unweighted mean (-1.0)
    assert off == pytest.approx(-2.0)
    assert diag["available"] is True
    assert diag["n_sites"] == 2 and diag["n_clip_reads"] == 40
    assert diag["by_strand"]["+"]["offset_bp"] == pytest.approx(1.0)
    assert diag["by_strand"]["-"]["offset_bp"] == pytest.approx(-3.0)


def test_a_v2_sidecar_has_no_such_column_and_says_so(tmp_path):
    v2_header = ("pas_id\tclip_reads\tclip_umis\tclip_reads_f3844"
                 "\tclip_umis_f3844\twindow_reads\ttier")
    p = _support(tmp_path, "v2.support.tsv",
                 ["1\t10\t5\t10\t5\t100\t1"], header=v2_header)
    off, diag = estimate_clip_anchored_offset([(p, "+")])
    assert off == 0.0
    assert diag["available"] is False and diag["missing_column"] is True


def test_measured_zero_is_distinguishable_from_could_not_measure(tmp_path):
    p = _support(tmp_path, "z.support.tsv",
                 ["1\t10\t5\t10\t5\t100\t1\t1\t0\t0.00"])
    off, diag = estimate_clip_anchored_offset([(p, "+")])
    assert off == 0.0 and diag["available"] is True and diag["n_sites"] == 1


def test_a_missing_sidecar_is_skipped_not_fatal(tmp_path):
    off, diag = estimate_clip_anchored_offset([(tmp_path / "nope.tsv", "+")])
    assert off == 0.0 and diag["available"] is False


# --------------------------------------------------------------------------
# 5. the inferred_cleavage column
# --------------------------------------------------------------------------
def test_inferred_cleavage_reports_without_moving_the_bed(tmp_path):
    bed = _write(tmp_path)
    sup = _support(tmp_path, "pas.support.tsv", [
        "t1_plus\t10\t7\t10\t7\t100\t1\t2\t5\t0.10",
        "t1_minus\t8\t3\t8\t3\t90\t1\t2\t4\t-0.20",
        "t2_plus\t0\t0\t0\t0\t60\t2\t0\t0\tNA",
        "t2_minus\t0\t0\t0\t0\t60\t2\t0\t0\tNA",
    ])
    before = bed.read_bytes()
    n = append_inferred_cleavage(bed, sup, {1: -2, 2: 0})
    assert n == 4
    assert bed.read_bytes() == before, "the column must not move pasbed.bed"
    lines = sup.read_text().splitlines()
    assert lines[0].split("\t")[-1] == "inferred_cleavage"
    got = {l.split("\t")[0]: l.split("\t")[-1] for l in lines[1:]}
    assert got["t1_plus"] == "998"       # + : reported base 1000, -2 upstream
    assert got["t1_minus"] == "2002"     # - : reported base 2000, -2 upstream
    assert got["t2_plus"] == "3299"      # tier-2 offset 0 == the reported base
    assert got["t2_minus"] == "4000"


def test_the_default_offset_makes_the_column_equal_the_reported_base(tmp_path):
    bed = _write(tmp_path)
    sup = _support(tmp_path, "pas.support.tsv", [
        "t1_plus\t10\t7\t10\t7\t100\t1\t2\t5\t0.10",
        "t1_minus\t8\t3\t8\t3\t90\t1\t2\t4\t-0.20",
        "t2_plus\t0\t0\t0\t0\t60\t2\t0\t0\tNA",
        "t2_minus\t0\t0\t0\t0\t60\t2\t0\t0\tNA",
    ])
    append_inferred_cleavage(bed, sup, {1: 0, 2: 0})
    for line in sup.read_text().splitlines()[1:]:
        f = line.split("\t")
        row = [r for r in _read(bed) if r[3] == f[0]][0]
        assert int(f[-1]) == cleavage_point(int(row[1]), int(row[2]), row[5])


def test_appending_twice_is_a_no_op(tmp_path):
    bed = _write(tmp_path)
    sup = _support(tmp_path, "pas.support.tsv",
                   ["t1_plus\t10\t7\t10\t7\t100\t1\t2\t5\t0.10"])
    append_inferred_cleavage(bed, sup, {1: -2, 2: 0})
    once = sup.read_bytes()
    assert append_inferred_cleavage(bed, sup, {1: -5, 2: 0}) == 0
    assert sup.read_bytes() == once


def test_a_sidecar_row_with_no_bed_row_gets_NA(tmp_path):
    bed = _write(tmp_path, rows=ROWS[:1])
    sup = _support(tmp_path, "pas.support.tsv", [
        "t1_plus\t10\t7\t10\t7\t100\t1\t2\t5\t0.10",
        "ghost\t1\t1\t1\t1\t10\t1\t1\t0\t0.00",
    ])
    append_inferred_cleavage(bed, sup, {1: -2, 2: 0})
    got = {l.split("\t")[0]: l.split("\t")[-1]
           for l in sup.read_text().splitlines()[1:]}
    assert got["ghost"] == "NA"


# --------------------------------------------------------------------------
# 6. the decision point — where a regression actually hid
# --------------------------------------------------------------------------
from ema.countmatrix.cleavage_offset import resolve_offset_for_run  # noqa: E402


def test_the_default_moves_nothing():
    assert resolve_offset_for_run("none", strategy="clip_seeded",
                                  auto_legacy=False) == ("none", 0, "all")
    assert resolve_offset_for_run(0, strategy="lambda_gradient",
                                  auto_legacy=False) == ("none", 0, "all")


def test_an_explicit_constant_wins_and_picks_its_tier():
    assert resolve_offset_for_run(95, strategy="clip_seeded",
                                  auto_legacy=False) == ("const", 95, "tier2")
    assert resolve_offset_for_run(-2, strategy="clip_seeded",
                                  auto_legacy=False) == ("const", -2, "tier1")
    assert resolve_offset_for_run(95, strategy="lambda_gradient",
                                  auto_legacy=False) == ("const", 95, "all")


def test_auto_uses_this_runs_clip_anchored_estimate():
    assert resolve_offset_for_run("auto", strategy="clip_seeded",
                                  auto_legacy=False,
                                  clip_estimate=-0.3339) == ("auto", 0, "all")
    assert resolve_offset_for_run("auto", strategy="clip_seeded",
                                  auto_legacy=False,
                                  clip_estimate=-2.7) == ("auto", -3, "tier1")


def test_auto_without_an_estimate_raises_rather_than_silently_using_zero():
    with pytest.raises(ValueError, match="clip_offset_mean"):
        resolve_offset_for_run("auto", strategy="clip_seeded",
                               auto_legacy=False, clip_estimate=None)


def test_the_legacy_coverage_estimator_is_still_wired_up():
    """The regression this function exists to prevent: a first cut of main.py's
    block dropped the issue-#72 path entirely, so --auto-cleavage-offset with a
    COVERAGE strategy silently did nothing -- and the whole suite stayed green,
    because the estimator has unit tests and the glue above it had none."""
    assert resolve_offset_for_run(
        "none", strategy="lambda_gradient", auto_legacy=True,
        legacy_estimate=95) == ("auto_legacy", 95, "all")
    # ...and it beats an explicit constant, exactly as v2 had it
    assert resolve_offset_for_run(
        50, strategy="lambda_gradient", auto_legacy=True,
        legacy_estimate=95) == ("auto_legacy", 95, "all")
    # ...but not --cleavage-offset auto, which is the clip-anchored estimator
    assert resolve_offset_for_run(
        "auto", strategy="lambda_gradient", auto_legacy=True,
        legacy_estimate=95, clip_estimate=-0.3) == ("auto", 0, "all")


def test_the_legacy_estimator_falling_over_leaves_the_explicit_constant():
    assert resolve_offset_for_run(
        7, strategy="lambda_gradient", auto_legacy=True,
        legacy_estimate=None) == ("const", 7, "all")
