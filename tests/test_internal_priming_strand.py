"""Strand-aware internal-priming window (fix/ip-filter-strand).

Internal priming happens when the oligo-dT primer anneals to a genome-encoded
A-stretch DOWNSTREAM of the cleavage site in TRANSCRIPT orientation. On the
'-' strand that stretch is genomically UPSTREAM of the site (lower
coordinates) and reads as a T-run on the forward strand, so the window must be
mirrored: ``[pos - right, pos + left)`` instead of ``[pos - left, pos + right)``.

REGRESSION (4efeb12 and earlier): ``filter_internal_priming`` used the '+'
window ``[pos-10, pos+30)`` on BOTH strands, so on '-' it tested 30 nt
*upstream* / 10 nt downstream of the cleavage site in transcript orientation
-- mostly the wrong side. ``_shipped_rule_4efeb12`` below re-implements that
rule verbatim so the two '-' cases can assert that the fixed helper and the
shipped behaviour disagree exactly where they should.

The synthetic genome is C/G-only so no A/T run (or A/T-rich window) can occur
by accident; every signal is implanted explicitly.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from ema.experimental.internal_priming import (
    call_internal_priming,
    check_internal_priming,
    filter_internal_priming,
    ip_window,
    reverse_complement,
)

L, R, K, F = 10, 30, 6, 0.7  # shipped defaults: window_left/right, a_stretch, a_fraction


def _cg(n: int) -> str:
    """C/G-only background of length *n* (no A/T anywhere)."""
    return ("CG" * (n // 2 + 1))[:n]


def _implant(genome: str, start: int, motif: str) -> str:
    assert start >= 0 and start + len(motif) <= len(genome)
    return genome[:start] + motif + genome[start + len(motif):]


def _shipped_rule_4efeb12(genome: str, pos: int, strand: str) -> bool:
    """The rule as shipped up to 4efeb12 (same forward window on both strands)."""
    seq = genome[max(0, pos - L):pos + R].upper()
    base = "A" if strand == "+" else "T"
    return (base * K in seq) or (len(seq) > 0 and seq.count(base) / len(seq) >= F)


# ---------------------------------------------------------------------------
# ip_window: the genomic window is mirrored on '-'
# ---------------------------------------------------------------------------
def test_ip_window_is_mirrored_on_minus_strand():
    assert ip_window(100, "+", L, R) == (90, 130)
    assert ip_window(100, "-", L, R) == (70, 110)
    # both windows have the same length and the same layout in transcript
    # orientation: `left` nt before the cleavage base, `right` nt after it
    assert (130 - 90) == (110 - 70) == L + R


def test_ip_window_clamps_start_at_zero():
    assert ip_window(5, "+", L, R) == (0, 35)
    assert ip_window(5, "-", L, R) == (0, 15)


# ---------------------------------------------------------------------------
# (a) '+' strand: A-run downstream flags, A-run upstream does not
# ---------------------------------------------------------------------------
def test_plus_a_run_downstream_is_flagged():
    g = _implant(_cg(300), 110, "A" * 6)  # BED end=100 -> window [90, 130)
    flag, tested = check_internal_priming(g, 100, "+", L, R, K, F)
    assert flag is True
    assert tested == g[90:130]
    assert "AAAAAA" in tested


def test_plus_a_run_upstream_is_not_flagged():
    g = _implant(_cg(300), 70, "A" * 6)  # upstream of BED end=100, outside [90, 130)
    flag, tested = check_internal_priming(g, 100, "+", L, R, K, F)
    assert flag is False
    assert "A" not in tested


def test_plus_behaviour_is_identical_to_shipped_rule():
    """'+' must stay byte-identical: the fix only touches '-'."""
    for start in range(60, 140, 3):
        g = _implant(_cg(300), start, "A" * 6)
        assert check_internal_priming(g, 100, "+", L, R, K, F)[0] == \
            _shipped_rule_4efeb12(g, 100, "+")


# ---------------------------------------------------------------------------
# (b) '-' strand: T-run genomically UPSTREAM (= downstream in transcript
#     orientation) flags; T-run genomically downstream does not.
# ---------------------------------------------------------------------------
def test_minus_t_run_genomically_upstream_is_flagged():
    # BED start=200 on '-': cleavage base 200, transcript runs toward lower
    # coordinates, so the A-stretch that primes oligo-dT is a T-run BELOW 200.
    # Place it at [175, 181): inside the mirrored window [170, 210) but
    # outside the shipped window [190, 230).
    g = _implant(_cg(400), 175, "T" * 6)
    flag, tested = check_internal_priming(g, 200, "-", L, R, K, F)
    assert flag is True
    assert tested == reverse_complement(g[170:210])
    assert "AAAAAA" in tested  # the T-run reads as an A-run in transcript orientation
    # REGRESSION: the shipped rule looked at [190, 230) and MISSED this site.
    assert _shipped_rule_4efeb12(g, 200, "-") is False


def test_minus_t_run_genomically_downstream_is_not_flagged():
    # T-run at [215, 221): genomically downstream of the site = UPSTREAM in
    # transcript orientation, i.e. inside the transcript body, not in the
    # primer-annealing region. Inside the shipped window [190, 230), outside
    # the mirrored window [170, 210).
    g = _implant(_cg(400), 215, "T" * 6)
    flag, tested = check_internal_priming(g, 200, "-", L, R, K, F)
    assert flag is False
    assert "A" not in tested and "T" not in tested
    # REGRESSION: the shipped rule flagged this site (wrong side).
    assert _shipped_rule_4efeb12(g, 200, "-") is True


def test_minus_equals_plus_on_reverse_complemented_genome():
    """Strand symmetry: a '-' site on g is the same call as the '+' site at the
    mirrored coordinate on revcomp(g) (BED start on '-' <-> BED end on '+')."""
    n = 400
    for start in (150, 175, 195, 205, 215):
        g = _implant(_cg(n), start, "T" * 6)
        minus = check_internal_priming(g, 200, "-", L, R, K, F)
        plus = check_internal_priming(reverse_complement(g), n - 200, "+", L, R, K, F)
        assert minus == plus


# ---------------------------------------------------------------------------
# (c) the A-fraction rule (>= 0.7 of the window) on both strands
# ---------------------------------------------------------------------------
# 40-nt window with 28 A's (exactly 0.70) and no run of 6: flagged by the
# fraction rule only; mutating one A -> 27/40 = 0.675 is below threshold.
_FRAC_WINDOW = "AAAAAC" * 5 + "AAAC" + "CCCCCC"
assert len(_FRAC_WINDOW) == L + R and _FRAC_WINDOW.count("A") == 28
assert "AAAAAA" not in _FRAC_WINDOW
_FRAC_WINDOW_BELOW = _FRAC_WINDOW.replace("AAAAAC", "AAAACC", 1)
assert _FRAC_WINDOW_BELOW.count("A") == 27


def test_fraction_rule_plus():
    g = _implant(_cg(300), 90, _FRAC_WINDOW)  # fills exactly [90, 130) for BED end=100
    assert check_internal_priming(g, 100, "+", L, R, K, F)[0] is True
    g = _implant(_cg(300), 90, _FRAC_WINDOW_BELOW)
    assert check_internal_priming(g, 100, "+", L, R, K, F)[0] is False


def test_fraction_rule_minus():
    # forward-strand bases of the mirrored window [170, 210) for BED start=200
    g = _implant(_cg(400), 170, reverse_complement(_FRAC_WINDOW))
    flag, tested = check_internal_priming(g, 200, "-", L, R, K, F)
    assert flag is True and tested == _FRAC_WINDOW
    g = _implant(_cg(400), 170, reverse_complement(_FRAC_WINDOW_BELOW))
    assert check_internal_priming(g, 200, "-", L, R, K, F)[0] is False


def test_call_internal_priming_is_pure_and_case_insensitive():
    assert call_internal_priming("cgcgaaaaaacg", "+", K, F) == (True, "CGCGAAAAAACG")
    assert call_internal_priming("cgttttttcgcg", "-", K, F) == (True, "CGCGAAAAAACG")
    assert call_internal_priming("", "+", K, F) == (False, "")
    assert call_internal_priming("", "-", K, F) == (False, "")


# ---------------------------------------------------------------------------
# (d) contig-edge handling
# ---------------------------------------------------------------------------
def test_contig_start_edge_minus_is_clamped_not_raised():
    g = _implant(_cg(60), 0, "T" * 6)  # T-run at the very start of the contig
    flag, tested = check_internal_priming(g, 5, "-", L, R, K, F)  # window [0, 15)
    assert flag is True
    assert tested == reverse_complement(g[0:15])


def test_contig_end_edge_plus_is_truncated_not_raised():
    g = _implant(_cg(60), 54, "A" * 6)  # run ends exactly at the contig end
    flag, tested = check_internal_priming(g, 50, "+", L, R, K, F)  # window [40, 80) -> [40, 60)
    assert flag is True
    assert tested == g[40:60]


def test_short_truncated_window_uses_its_own_length_for_the_fraction():
    g = _cg(60) + "AAA"  # 63 nt; '+' site at end=60 -> window [50, 63) = 10 C/G + 3 A
    assert check_internal_priming(g, 60, "+", L, R, K, F)[0] is False
    g = _cg(60) + "AAAAAAAAAAAAAAAAAAAAAAAAA"  # 25 A -> 25/35 = 0.714 >= 0.7
    assert check_internal_priming(g, 60, "+", L, R, K, 0.7)[0] is True


def test_position_beyond_contig_gives_empty_window_and_no_flag():
    g = _cg(60)
    assert check_internal_priming(g, 500, "+", L, R, K, F) == (False, "")
    assert check_internal_priming(g, 500, "-", L, R, K, F) == (False, "")


# ---------------------------------------------------------------------------
# (e) end-to-end through the real entry point on a tiny FASTA + BED
# ---------------------------------------------------------------------------
def _write_fasta(path: Path, contigs: dict[str, str]) -> None:
    with open(path, "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


@pytest.fixture
def tiny_genome(tmp_path: Path) -> Path:
    pytest.importorskip("pyfaidx")
    g = _cg(600)
    g = _implant(g, 110, "A" * 6)   # P1 '+', BED end=100, window [90,130)  -> flagged
    g = _implant(g, 175, "A" * 6)   # P2 '+', BED end=200, window [190,230) -> upstream, clean
    g = _implant(g, 275, "T" * 6)   # M1 '-', BED start=300: new [270,310) hit, old [290,330) miss
    g = _implant(g, 415, "T" * 6)   # M2 '-', BED start=400: new [370,410) miss, old [390,430) hit
    fa = tmp_path / "genome.fa"
    _write_fasta(fa, {"chrT": g, "chrE": _implant(_cg(20), 0, "T" * 6)})
    return fa


_SITES = [
    ("chrT", 90, 100, "P1_plus_downstream_run", 3, "+"),
    ("chrT", 190, 200, "P2_plus_upstream_run", 3, "+"),
    ("chrT", 300, 310, "M1_minus_genomic_upstream_run", 3, "-"),
    ("chrT", 400, 410, "M2_minus_genomic_downstream_run", 3, "-"),
]
_EXPECTED = {
    "P1_plus_downstream_run": True,
    "P2_plus_upstream_run": False,
    "M1_minus_genomic_upstream_run": True,    # shipped rule: False (missed)
    "M2_minus_genomic_downstream_run": False,  # shipped rule: True (wrong side)
}


def _write_bed(path: Path, rows) -> None:
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")


def test_filter_internal_priming_annotate_mode_flags_per_strand(tiny_genome, tmp_path):
    bed = tmp_path / "in.bed"
    _write_bed(bed, _SITES)
    out = tmp_path / "out.bed"
    stats = filter_internal_priming(str(bed), str(tiny_genome), str(out),
                                    window_left=L, window_right=R,
                                    a_stretch=K, a_fraction=F, mode="annotate")
    assert stats["flags"] == _EXPECTED
    assert stats["flagged"] == 2 and stats["filtered"] == 0
    assert out.read_text() == bed.read_text()  # annotate never drops a row


def test_filter_internal_priming_filter_mode_drops_the_right_sites(tiny_genome, tmp_path):
    bed = tmp_path / "in.bed"
    _write_bed(bed, _SITES)
    out = tmp_path / "out.bed"
    stats = filter_internal_priming(str(bed), str(tiny_genome), str(out),
                                    window_left=L, window_right=R,
                                    a_stretch=K, a_fraction=F, mode="filter")
    survivors = [ln.split("\t")[3] for ln in out.read_text().splitlines()]
    assert survivors == ["P2_plus_upstream_run", "M2_minus_genomic_downstream_run"]
    assert stats["filtered"] == 2 and stats["passed"] == 2 and stats["total"] == 4


def test_filter_internal_priming_edge_and_missing_contig_rows_are_kept(tiny_genome, tmp_path):
    bed = tmp_path / "in.bed"
    _write_bed(bed, [
        ("chrE", 5, 6, "E_minus_at_contig_start", 1, "-"),   # window [0,15): T-run at 0 -> flag
        ("chrE", 14, 15, "E_plus_at_contig_end", 1, "+"),    # window [5,45) truncated: 'TCGCG...' clean
        ("chrMissing", 10, 11, "missing_contig", 1, "+"),     # not in FASTA -> kept, False
    ])
    out = tmp_path / "out.bed"
    stats = filter_internal_priming(str(bed), str(tiny_genome), str(out), mode="filter")
    assert stats["flags"] == {
        "E_minus_at_contig_start": True,
        "E_plus_at_contig_end": False,
        "missing_contig": False,
    }
    survivors = [ln.split("\t")[3] for ln in out.read_text().splitlines()]
    assert survivors == ["E_plus_at_contig_end", "missing_contig"]


def test_apply_filters_seam_threads_strand_correct_flags(tiny_genome, tmp_path):
    """Through ema.experimental.peak_filters.apply_filters, the seam ema.main
    calls per strand BED -- the '-' BED must now be judged on the mirrored window."""
    from ema.experimental.peak_filters import apply_filters

    neg = tmp_path / "negbed.bed"
    _write_bed(neg, [s for s in _SITES if s[5] == "-"])
    out = tmp_path / "neg.filtered.bed"
    stats = apply_filters(str(neg), str(out), genome_fasta=str(tiny_genome),
                          enable_internal_priming=True, ip_mode="filter",
                          ip_window_left=L, ip_window_right=R,
                          ip_a_stretch=K, ip_a_fraction=F)
    assert stats["internal_priming_flags"] == {
        "M1_minus_genomic_upstream_run": True,
        "M2_minus_genomic_downstream_run": False,
    }
    survivors = [ln.split("\t")[3] for ln in out.read_text().splitlines()]
    assert survivors == ["M2_minus_genomic_downstream_run"]
