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
    build_downstream_profiles,
    estimate_cleavage_offset,
    estimate_cleavage_offset_diag,
    inferred_cleavage_site,
    resolve_cleavage_offset,
    rewrite_bed_3prime_offset,
    shift_3prime_end,
)

MAX_OFFSET = 150
PLANTED_ARUN = (92, 109)      # genomic A-run downstream offsets [lo, hi)
PLANTED_AATAAA = 76           # AATAAA starts at this downstream offset
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _downstream_template(max_offset: int = MAX_OFFSET) -> str:
    """A 'C' background with an AATAAA at +76 and an A-run at +92..108."""
    t = ["C"] * max_offset
    for i, ch in enumerate("AATAAA"):
        t[PLANTED_AATAAA + i] = ch
    for i in range(*PLANTED_ARUN):
        t[i] = "A"
    return "".join(t)


def _make_synthetic_genome(tmp_path, strand: str, n_peaks: int = 60):
    """Write a synthetic FASTA + strand BED that plants the cleavage signal.

    Each peak's downstream window (strand-aware) is the same template, so the
    aggregated A-fraction profile crests inside the planted A-run and the
    AATAAA density spikes at the planted offset.  Returns (fasta_path,
    bed_path).
    """
    template = _downstream_template()
    stride = 400
    genome = list("G" * (stride * (n_peaks + 2)))
    bed_lines = []
    for k in range(n_peaks):
        base = stride * (k + 1)
        if strand == "+":
            peak_start, peak_end = base + 100, base + 200
            # downstream (increasing coord) = template, verbatim.
            for i, ch in enumerate(template):
                genome[peak_end + i] = ch
        else:
            peak_start, peak_end = base + 200, base + 300
            # 3' end is peak_start; downstream = decreasing coord,
            # reverse-complemented -> place revcomp(template) upstream of it.
            rc = _revcomp(template)
            for i, ch in enumerate(rc):
                genome[peak_start - MAX_OFFSET + i] = ch
        bed_lines.append(f"chr1\t{peak_start}\t{peak_end}\tpas{k}\t0\t{strand}")

    fasta = tmp_path / f"genome_{strand}.fa"
    fasta.write_text(">chr1\n" + "".join(genome) + "\n")
    bed = tmp_path / f"pas_{strand}.bed"
    bed.write_text("\n".join(bed_lines) + "\n")
    return fasta, bed


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
# estimate_cleavage_offset — fallback behaviour / backwards compat
# --------------------------------------------------------------------------

def test_estimator_returns_default_constant_no_profiles():
    assert estimate_cleavage_offset() == DEFAULT_CLEAVAGE_OFFSET


def test_estimator_default_in_observed_window():
    # The fallback must sit inside the empirically observed +90..+105 window.
    assert 90 <= DEFAULT_CLEAVAGE_OFFSET <= 105


def test_estimator_respects_custom_fallback():
    assert estimate_cleavage_offset(fallback=88) == 88


def test_estimator_too_short_profiles_fall_back():
    # Profiles shorter than the search band cannot yield a crest -> fallback.
    assert estimate_cleavage_offset(
        aataaa_profile=[0, 0, 1, 5, 1],
        a_fraction_profile=[0.1, 0.2, 0.43, 0.1],
    ) == DEFAULT_CLEAVAGE_OFFSET


def test_estimator_flat_profile_is_inconclusive():
    # A dead-flat A-fraction profile has no crest above background -> fallback,
    # and the diagnostics record the inconclusive method.
    flat = [0.25] * MAX_OFFSET
    offset, diag = estimate_cleavage_offset_diag(a_fraction_profile=flat)
    assert offset == DEFAULT_CLEAVAGE_OFFSET
    assert diag["method"] == "fallback"


# --------------------------------------------------------------------------
# estimate_cleavage_offset — recovers a planted crest (pure profiles)
# --------------------------------------------------------------------------

def _planted_a_fraction(crest: int, max_offset: int = MAX_OFFSET):
    """Background 0.05 A-fraction with a sharp crest to ~0.9 at `crest`."""
    prof = [0.05] * max_offset
    for i in range(crest - 4, crest + 5):
        if 0 <= i < max_offset:
            prof[i] = 0.9
    return prof


def test_estimator_recovers_planted_crest_from_profile():
    prof = _planted_a_fraction(98)
    offset, diag = estimate_cleavage_offset_diag(a_fraction_profile=prof)
    assert diag["method"] == "a_fraction_crest"
    assert abs(offset - 98) <= 6
    assert diag["crest_value"] > diag["background"]


def test_estimator_crest_restricted_to_band():
    # A crest OUTSIDE the band must be ignored (poly-A tails, exons, etc.).
    prof = _planted_a_fraction(20)  # below band_lo=60
    offset, diag = estimate_cleavage_offset_diag(a_fraction_profile=prof)
    # No in-band crest -> fallback (band default 60..120).
    assert offset == DEFAULT_CLEAVAGE_OFFSET
    assert diag["method"] == "fallback"


def test_estimator_aataaa_secondary_crosscheck():
    # AATAAA density peaks at +76 -> +76 + spacing(22) = +98 secondary est.
    aa = [0.0] * MAX_OFFSET
    aa[76] = 1.0
    _, diag = estimate_cleavage_offset_diag(
        a_fraction_profile=_planted_a_fraction(98), aataaa_profile=aa,
    )
    assert diag["aataaa_peak"] == 76
    assert diag["aataaa_estimate"] == 98


def test_estimator_falls_back_to_aataaa_when_afraction_flat():
    # Flat A-fraction but a clear AATAAA peak -> use the AATAAA route.
    aa = [0.0] * MAX_OFFSET
    aa[80] = 1.0                      # 80 + 22 = 102, in band
    offset, diag = estimate_cleavage_offset_diag(
        a_fraction_profile=[0.2] * MAX_OFFSET, aataaa_profile=aa,
    )
    assert diag["method"] == "aataaa_spacing"
    assert offset == 102


# --------------------------------------------------------------------------
# build_downstream_profiles + end-to-end recovery on a synthetic genome
# --------------------------------------------------------------------------

def test_build_profiles_recovers_planted_offset_plus_strand(tmp_path):
    pytest.importorskip("pyfaidx")
    fasta, bed = _make_synthetic_genome(tmp_path, "+")
    a_frac, aataaa, n_used = build_downstream_profiles(
        bed, str(fasta), max_offset=MAX_OFFSET,
    )
    assert n_used == 60
    # A-fraction is ~1.0 across the planted A-run and ~0 in the 'C' background.
    assert a_frac[100] > 0.9
    assert a_frac[50] < 0.1
    # AATAAA density spikes at the planted +76 offset.
    assert aataaa[PLANTED_AATAAA] > 0.9

    offset, diag = estimate_cleavage_offset_diag(
        a_fraction_profile=a_frac, aataaa_profile=aataaa, n_peaks_used=n_used,
    )
    assert diag["method"] == "a_fraction_crest"
    assert PLANTED_ARUN[0] - 4 <= offset <= PLANTED_ARUN[1] + 4
    assert diag["aataaa_estimate"] == 98  # 76 + 22


def test_build_profiles_strand_symmetry(tmp_path):
    pytest.importorskip("pyfaidx")
    fa_p, bed_p = _make_synthetic_genome(tmp_path, "+")
    fa_m, bed_m = _make_synthetic_genome(tmp_path, "-")
    prof_p = build_downstream_profiles(bed_p, str(fa_p), max_offset=MAX_OFFSET)
    prof_m = build_downstream_profiles(bed_m, str(fa_m), max_offset=MAX_OFFSET)
    off_p = estimate_cleavage_offset(
        a_fraction_profile=prof_p[0], aataaa_profile=prof_p[1])
    off_m = estimate_cleavage_offset(
        a_fraction_profile=prof_m[0], aataaa_profile=prof_m[1])
    # Strand-aware reverse-complement handling => identical recovered offset.
    assert off_p == off_m
    assert PLANTED_ARUN[0] - 4 <= off_p <= PLANTED_ARUN[1] + 4


# --------------------------------------------------------------------------
# resolve_cleavage_offset — the wiring decision point (auto vs manual)
# --------------------------------------------------------------------------

def test_resolve_manual_mode_returns_explicit_unchanged():
    # auto=False: byte-for-byte legacy behaviour, no FASTA/BED access.
    assert resolve_cleavage_offset(0, auto=False) == (0, None)
    assert resolve_cleavage_offset(95, auto=False) == (95, None)


def test_resolve_auto_mode_estimates_from_data(tmp_path):
    pytest.importorskip("pyfaidx")
    fasta, bed = _make_synthetic_genome(tmp_path, "+")
    offset, diag = resolve_cleavage_offset(
        0, auto=True, bed_paths=[bed], genome_fasta=str(fasta),
    )
    assert diag is not None
    assert diag["method"] == "a_fraction_crest"
    assert PLANTED_ARUN[0] - 4 <= offset <= PLANTED_ARUN[1] + 4


def test_resolve_auto_mode_missing_fasta_falls_back(tmp_path):
    offset, diag = resolve_cleavage_offset(
        0, auto=True, bed_paths=[], genome_fasta=str(tmp_path / "nope.fa"),
    )
    assert offset == DEFAULT_CLEAVAGE_OFFSET
    assert diag["method"] == "fallback"


def test_resolve_auto_applied_shifts_peaks_but_unset_is_noop(tmp_path):
    pytest.importorskip("pyfaidx")
    fasta, bed = _make_synthetic_genome(tmp_path, "+")
    original = bed.read_text()

    # Manual unset (auto=False, offset 0) -> resolver no-op AND rewrite no-op.
    off_unset, _ = resolve_cleavage_offset(0, auto=False, bed_paths=[bed],
                                           genome_fasta=str(fasta))
    assert off_unset == 0
    assert rewrite_bed_3prime_offset(bed, off_unset) == 0
    assert bed.read_text() == original  # byte-identical legacy no-op

    # Auto mode -> resolver estimates a positive offset that shifts every peak.
    off_auto, diag = resolve_cleavage_offset(0, auto=True, bed_paths=[bed],
                                             genome_fasta=str(fasta))
    assert off_auto > 0
    n = rewrite_bed_3prime_offset(bed, off_auto)
    assert n == 60
    # + strand: every peak end moved downstream by the estimated offset.
    first = bed.read_text().splitlines()[0].split("\t")
    assert int(first[2]) == 600 + off_auto  # peak_end (400*1+200) + offset


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
