"""Adversarial strand verification for the peakAtail-prime feature emission.

WHY THIS FILE EXISTS
--------------------
This project has been burned twice by strand bugs in exactly this seam (the
internal-priming window was mirrored the wrong way until ``fix/ip-filter-strand``
/ #96).  ``tests/test_pas_features.py`` checks the same ground, but it does so
partly by comparing the tool against helper functions from the same module
(``feature_window`` / ``orient_window``), which cannot catch a convention that
is wrong *consistently*.

Everything here is derived from FIRST PRINCIPLES instead, from the one
definition the docs state:

    r = g - c   on '+'          c = BED end - 1 on '+'
    r = c - g   on '-'          c = BED start   on '-'
    r > 0 is DOWNSTREAM in transcript orientation.

:func:`implant` writes a transcript-oriented motif into a FORWARD-strand genome
using only that definition (and a complement table).  Nothing in this file
imports the window helpers it is testing, and every expected value is a
hand-written literal.

The decisive assertions are the ANTI-implants: a signal placed on the wrong
side in transcript space must produce a *zero* feature.  A strand flip passes
the "mirror" style of test and fails these.
"""
from __future__ import annotations

import pytest

from ema.countmatrix.pas_features import FeatureCollector
from ema.experimental.internal_priming import filter_internal_priming

_COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

CONTIG_LEN = 4000
CONTIG = "1"


def implant(genome: list, c: int, strand: str, r_lo: int, tseq: str) -> None:
    """Write *tseq* (transcript orientation) at transcript offsets
    ``r_lo .. r_lo + len(tseq) - 1`` around cleavage base *c*.

    Derived only from ``r = g - c`` ('+') / ``r = c - g`` ('-'): on '-' the
    forward genome gets the COMPLEMENT of each transcript base, at a
    coordinate that DECREASES as r increases.
    """
    for k, base in enumerate(tseq):
        r = r_lo + k
        if strand == "+":
            genome[c + r] = base
        else:
            genome[c - r] = _COMP[base]


def _fasta(tmp_path, genome: list):
    fa = tmp_path / "synthetic.fa"
    seq = "".join(genome)
    with open(fa, "w") as fh:
        fh.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")
    return str(fa)


def _bed(tmp_path, rows):
    """rows: (pas_id, cleavage_base c, strand, score)."""
    bed = tmp_path / "pas.bed"
    with open(bed, "w") as fh:
        for pas_id, c, strand, score in rows:
            # BED for a 1-bp PAS whose cleavage base is c:
            #   '+': end - 1 == c  ->  (c, c+1)
            #   '-': start   == c  ->  (c, c+1)
            fh.write(f"{CONTIG}\t{c}\t{c + 1}\t{pas_id}\t{score}\t{strand}\n")
    return str(bed)


def _collect(tmp_path, genome, rows):
    fa = _fasta(tmp_path, genome)
    bed = _bed(tmp_path, rows)
    coll = FeatureCollector()
    filter_internal_priming(
        bed, fa, None, window_left=10, window_right=30,
        a_stretch=6, a_fraction=0.7, mode="annotate", features=coll,
    )
    coll.finish()
    return {pid: coll.parsed(pid) for pid, *_ in rows}


# ---------------------------------------------------------------------------
# 1. the hexamer window, both strands, from first principles
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("strand", ["+", "-"])
def test_a_hexamer_is_found_at_the_transcript_offset_it_was_implanted_at(
        tmp_path, strand):
    """AATAAA with its LAST base at r = -10 must report hex_strong_off -10 on
    BOTH strands, and only the reverse-complemented forward bases can give it
    on '-'."""
    genome = ["C"] * CONTIG_LEN
    c = 2000
    # AATAAA occupies r -15..-10, so its last base sits at r = -10.
    implant(genome, c, strand, -15, "AATAAA")
    f = _collect(tmp_path, genome, [("H1", c, strand, 5)])["H1"]
    assert f["seq_ok"] == "1"
    assert f["hex_strong"] == "1"
    assert f["hex_any12"] == "1"
    assert f["hex_strong_off"] == "-10"
    assert f["hex_best_off"] == "-10"
    assert f["hex_n_types"] == "1"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_a_hexamer_on_the_wrong_side_is_not_found(tmp_path, strand):
    """The ANTI-implant: AATAAA placed DOWNSTREAM (r +10..+15) is outside the
    hexamer window -40..-5 and must not be reported.  A mirrored window would
    find it here and miss it above."""
    genome = ["C"] * CONTIG_LEN
    c = 2000
    implant(genome, c, strand, 10, "AATAAA")
    f = _collect(tmp_path, genome, [("H2", c, strand, 5)])["H2"]
    assert f["seq_ok"] == "1"
    assert f["hex_strong"] == "0"
    assert f["hex_any12"] == "0"
    assert f["hex_strong_off"] == "0"
    assert f["hex_best_off"] == "0"


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("r_last,expected", [(-5, "-5"), (-35, "-35")])
def test_the_hexamer_window_edges_are_inclusive(tmp_path, strand, r_last,
                                                expected):
    genome = ["C"] * CONTIG_LEN
    c = 2000
    implant(genome, c, strand, r_last - 5, "AATAAA")
    f = _collect(tmp_path, genome, [("H3", c, strand, 5)])["H3"]
    assert f["hex_strong"] == "1", f
    assert f["hex_strong_off"] == expected


@pytest.mark.parametrize("strand", ["+", "-"])
def test_a_hexamer_one_base_past_the_window_is_not_found(tmp_path, strand):
    genome = ["C"] * CONTIG_LEN
    c = 2000
    implant(genome, c, strand, -41, "AATAAA")   # last base at r = -36
    f = _collect(tmp_path, genome, [("H4", c, strand, 5)])["H4"]
    assert f["hex_strong"] == "0", f
    genome2 = ["C"] * CONTIG_LEN
    implant(genome2, c, strand, -9, "AATAAA")   # last base at r = -4
    f2 = _collect(tmp_path, genome2, [("H5", c, strand, 5)])["H5"]
    assert f2["hex_strong"] == "0", f2


# ---------------------------------------------------------------------------
# 2. the downstream A window and the caller's own IP window, both strands
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("strand", ["+", "-"])
def test_a_downstream_A_run_is_measured_on_both_strands(tmp_path, strand):
    """Seven A's at r +2..+8, transcript orientation.  On '-' that is seven
    T's at DECREASING genomic coordinates on the forward strand."""
    genome = ["C"] * CONTIG_LEN
    c = 2000
    implant(genome, c, strand, 2, "A" * 7)
    f = _collect(tmp_path, genome, [("A1", c, strand, 5)])["A1"]
    assert f["seq_ok"] == "1"
    assert f["a_count_d18"] == "7"
    assert f["a_run_d18"] == "7"
    assert f["a_frac_d18"] == str(round(7 / 18, 4))
    assert f["a_run_d30"] == "7"
    assert f["a_frac_d30"] == str(round(7 / 30, 4))
    assert f["kin_ip_flag"] == "0"          # 7 < 12
    # The caller's own veto: window is r -9..+30 (40 nt), the run is inside it.
    assert f["ip_tool_flag"] == "1"         # a_stretch 6 <= 7
    assert f["ip_tool_arun"] == "7"
    assert f["ip_tool_afrac"] == str(round(7 / 40, 4))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_an_upstream_A_run_is_not_a_downstream_A_run(tmp_path, strand):
    """The ANTI-implant that matters most: the same seven A's placed UPSTREAM
    (r -20..-14) must leave every downstream-A column at zero AND must not
    trip the internal-priming veto, whose window starts at r = -9.

    A mirrored window would report 7 here and 0 above."""
    genome = ["C"] * CONTIG_LEN
    c = 2000
    implant(genome, c, strand, -20, "A" * 7)
    f = _collect(tmp_path, genome, [("A2", c, strand, 5)])["A2"]
    assert f["seq_ok"] == "1"
    assert f["a_count_d18"] == "0"
    assert f["a_run_d18"] == "0"
    assert f["a_run_d30"] == "0"
    assert f["ip_tool_flag"] == "0"
    assert f["ip_tool_arun"] == "0"
    assert f["ip_tool_afrac"] == "0.0"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_the_kinnex_flag_needs_twelve_A_in_the_first_eighteen(tmp_path, strand):
    genome = ["C"] * CONTIG_LEN
    c = 2000
    implant(genome, c, strand, 1, "A" * 12)     # r +1..+12
    f = _collect(tmp_path, genome, [("K1", c, strand, 5)])["K1"]
    assert f["a_count_d18"] == "12"
    assert f["kin_ip_flag"] == "1"
    genome2 = ["C"] * CONTIG_LEN
    implant(genome2, c, strand, 1, "A" * 11)
    f2 = _collect(tmp_path, genome2, [("K2", c, strand, 5)])["K2"]
    assert f2["a_count_d18"] == "11"
    assert f2["kin_ip_flag"] == "0"
    # r +19 is outside the d18 window but inside d30.
    genome3 = ["C"] * CONTIG_LEN
    implant(genome3, c, strand, 19, "A" * 6)
    f3 = _collect(tmp_path, genome3, [("K3", c, strand, 5)])["K3"]
    assert f3["a_count_d18"] == "0"
    assert f3["a_run_d30"] == "6"


@pytest.mark.parametrize("strand", ["+", "-"])
def test_the_ip_veto_window_does_not_reach_r_minus_ten(tmp_path, strand):
    """--ip-window-left 10 means the tested string starts at r = -9.  Six A's
    at r -15..-10 must NOT flag; the same six at r -9..-4 must."""
    genome = ["C"] * CONTIG_LEN
    c = 2000
    implant(genome, c, strand, -15, "A" * 6)
    out = _collect(tmp_path, genome, [("W1", c, strand, 5)])["W1"]
    assert out["ip_tool_flag"] == "0", out
    genome2 = ["C"] * CONTIG_LEN
    implant(genome2, c, strand, -9, "A" * 6)
    out2 = _collect(tmp_path, genome2, [("W2", c, strand, 5)])["W2"]
    assert out2["ip_tool_flag"] == "1", out2
    assert out2["ip_tool_arun"] == "6"


def test_the_two_strands_agree_on_an_identical_transcript_neighbourhood(tmp_path):
    """Same motif set in transcript space on '+' and on '-' at two different
    loci: every sequence column must be identical."""
    genome = ["C"] * CONTIG_LEN
    implant(genome, 1000, "+", -15, "AATAAA")
    implant(genome, 1000, "+", 3, "AAAAAAAA")
    implant(genome, 3000, "-", -15, "AATAAA")
    implant(genome, 3000, "-", 3, "AAAAAAAA")
    got = _collect(tmp_path, genome,
                   [("P", 1000, "+", 5), ("M", 3000, "-", 5)])
    seq_cols = [k for k in got["P"] if not k.startswith(("d_", "n_cand", "mol_",
                                                         "is_local"))]
    for col in seq_cols:
        assert got["P"][col] == got["M"][col], (
            f"{col}: '+' {got['P'][col]} != '-' {got['M'][col]}")


# ---------------------------------------------------------------------------
# 3. clip_offset_mean: the sign is decided by the strand, at the PRODUCER
# ---------------------------------------------------------------------------
def _seeder_offset(direction: bool, positions_and_reads):
    """Drive a real ClipSeeder and return the tier-1 row's clip_offset_mean."""
    from ema.countmatrix.polya import ClipSeeder

    s = ClipSeeder(direction, seed_window=25, min_umis=1, window=100,
                   features=True)
    n = 0
    for pos, reads in positions_and_reads:
        for k in range(reads):
            n += 1
            s.add_read(pos - 50, pos + 41, f"CB{n}")
            s.add_clip(pos, f"CB{n}", f"UMI{n}", pos + 41, True)
    rows: list = []
    recs = s.flush(support_out=rows)
    tier1 = [r for r in rows if r["tier"] == 1]
    assert len(tier1) == 1, (recs, rows)
    return tier1[0]


def test_clip_offset_mean_is_positive_downstream_on_the_plus_strand():
    """Mode at 1000 (5 reads), one member at 1010.  On '+' the member is
    DOWNSTREAM of the call, so the read-weighted mean is +10/6."""
    row = _seeder_offset(False, [(1000, 5), (1010, 1)])
    assert row["clip_positions"] == 2
    assert row["clip_span"] == 10
    assert row["clip_offset_mean"] == "%.2f" % (10 / 6)


def test_clip_offset_mean_is_negative_for_the_same_geometry_on_minus():
    """Identical genomic geometry on '-': 1010 is now UPSTREAM of the call at
    1000, so the sign must flip.  A strand-blind implementation returns
    +1.67 here and passes the '+' test above."""
    row = _seeder_offset(True, [(1000, 5), (1010, 1)])
    assert row["clip_positions"] == 2
    assert row["clip_span"] == 10
    assert row["clip_offset_mean"] == "%.2f" % (-10 / 6)


def test_a_single_position_cluster_has_offset_exactly_zero():
    for direction in (False, True):
        row = _seeder_offset(direction, [(1000, 4)])
        assert row["clip_positions"] == 1
        assert row["clip_span"] == 0
        assert row["clip_offset_mean"] == "0.00"


def test_the_offset_is_read_weighted_not_position_weighted():
    """Mode at 1000 with 1 read, a member at 1005 with 9 reads would be a
    different number under position weighting (2.5) than under read weighting
    (4.5).  cluster_clip_sites puts the mode on the heaviest position, so use
    a three-position cluster to separate the two."""
    # positions 1000 (10 reads, the mode), 1004 (1 read), 1020 (9 reads)
    row = _seeder_offset(False, [(1000, 10), (1004, 1), (1020, 9)])
    expected = (0 * 10 + 4 * 1 + 20 * 9) / 20
    assert row["clip_offset_mean"] == "%.2f" % expected


# ---------------------------------------------------------------------------
# 4. the signed cleavage shift, both strands, from first principles
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("strand,offset,expect", [
    # '+': downstream is a HIGHER coordinate; the cleavage base is end - 1.
    ("+", +3, (1000, 1004)),      # v2 grows the interval downstream
    ("+", -3, (997, 998)),        # a 1-bp PAS moves upstream, still 1 bp
    # '-': downstream is a LOWER coordinate; the cleavage base is start.
    ("-", +3, (997, 1001)),
    ("-", -3, (1003, 1004)),
])
def test_the_shift_moves_the_cleavage_base_in_transcript_orientation(
        strand, offset, expect):
    from ema.countmatrix.cleavage_offset import (
        cleavage_point, shift_cleavage_point,
    )

    start, end = 1000, 1001               # a 1-bp PAS; c == 1000 on both strands
    assert cleavage_point(start, end, strand) == 1000
    got = shift_cleavage_point(start, end, strand, offset)
    assert got == expect
    # And the reported base really landed where transcript arithmetic says.
    moved_c = cleavage_point(got[0], got[1], strand)
    assert moved_c == (1000 + offset if strand == "+" else 1000 - offset)
