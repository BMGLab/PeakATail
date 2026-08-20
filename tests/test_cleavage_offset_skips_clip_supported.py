"""The 3' cleavage-offset correction must not move clip-supported PAS.

``clip_seeded`` places tier-1 PAS at the observed poly(A) clip site (the true
cleavage coordinate) and stores the supporting clip-read count in BED column 5.
Applying the R2-read-length offset to those rows would push them past the
cleavage site. Coverage-only rows (score 0) still need the shift.
"""
from ema.countmatrix.cleavage_offset import rewrite_bed_3prime_offset


def _write(tmp_path, rows):
    p = tmp_path / "pas.bed"
    p.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return p


def _read(p):
    return [l.split("\t") for l in p.read_text().splitlines()]


ROWS = [
    ("1", 1000, 1001, "tier1_plus", 7, "+"),   # clip-supported, + strand
    ("1", 2000, 2001, "tier1_minus", 3, "-"),  # clip-supported, - strand
    ("1", 3000, 3300, "tier2_plus", 0, "+"),   # coverage-only
    ("1", 4000, 4300, "tier2_minus", 0, "-"),  # coverage-only
]


def test_default_behaviour_shifts_every_row(tmp_path):
    p = _write(tmp_path, ROWS)
    n = rewrite_bed_3prime_offset(p, 90)
    assert n == 4
    out = _read(p)
    assert int(out[0][2]) == 1001 + 90          # + strand: end moves downstream
    assert int(out[1][1]) == 2000 - 90          # - strand: start moves downstream


def test_skip_supported_exempts_clip_tier_only(tmp_path):
    p = _write(tmp_path, ROWS)
    n = rewrite_bed_3prime_offset(p, 90, skip_supported=True)
    assert n == 2, "only the two coverage-only rows may move"
    out = _read(p)
    assert (int(out[0][1]), int(out[0][2])) == (1000, 1001)   # untouched
    assert (int(out[1][1]), int(out[1][2])) == (2000, 2001)   # untouched
    assert int(out[2][2]) == 3300 + 90
    assert int(out[3][1]) == 4000 - 90


def test_skip_supported_is_a_noop_at_zero_offset(tmp_path):
    p = _write(tmp_path, ROWS)
    before = p.read_text()
    assert rewrite_bed_3prime_offset(p, 0, skip_supported=True) == 0
    assert p.read_text() == before
