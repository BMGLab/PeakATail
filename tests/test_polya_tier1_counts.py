"""Tier-1 PAS must be counted from the reads that pile up at the cleavage
site, not from the ~1% of reads that carry the poly(A) clip (Stage 2 fix).

The Stage-2 benchmark re-run surfaced a ~10x quantification regression:
``_flush_seeded`` wrote the ClipSeeder's clip-only per-CB dict as the
count-matrix row of every tier-1 PAS (mouse1 testis: 4.1M counts vs 41.9M
shipped; 23% of the real STARsolo cells fell under ``min_read``).  These
tests pin the corrected accounting on a synthetic BAM built in-process:

* ``PEAK`` locus (+): a 200-read coverage peak with only 2 clipped reads.
  The tier-1 PAS at the cleavage site must carry the peak's ~200 counts —
  exactly the mass the coverage strategy assigns to that locus — not ~2.
* ``ISOLATED`` locus (+) and its minus-strand mirror: a 3-read clip cluster
  with 30 non-clip reads inside the count window and no coverage candidate
  (max stack height < lambda_gradient's ``min_height``).  Must count 33.
* ``BARE`` locus (+): a coverage-only peak with no clip evidence.  Its
  tier-2 row must be byte-for-byte what the coverage strategy writes.

Coordinates and tier tags are a counting-independent function of the clip
sites and coverage candidates, so the BED of this BAM must not move at all
between the clip-only and the corrected counting.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

from tests.test_polya_two_tier_integration import _barcode, _rows  # noqa: E402

CHROM = "1"
CHROM_LEN = 200_000
SEQ_LEN = 91

PEAK_SITE = 2_400        # + strand, 200-read coverage peak, 2 clipped reads
ISOLATED_SITE = 12_400   # + strand, 3 clipped + 30 in-window reads, no peak
BARE_SITE = 22_400       # + strand, coverage-only peak, no clips
ISOLATED_SITE_NEG = 32_400  # - strand mirror of ISOLATED

N_PEAK_RAMP = 150
N_PEAK_SUMMIT = 50
N_PEAK_CLIPPED = 2
N_ISOLATED_CLIP = 3
N_ISOLATED_FLANK = 30
N_BARE_RAMP = 60
N_BARE_SUMMIT = 40


def _specs():
    """``(start, aligned_len, clip_len, reverse)`` for every synthetic read."""
    specs = []
    # PEAK (+): ramp + summit, first N_PEAK_CLIPPED summit reads clipped.
    for i in range(N_PEAK_RAMP):
        specs.append((PEAK_SITE - 350 + 2 * i, 80, 0, False))
    for i in range(N_PEAK_SUMMIT):
        start = PEAK_SITE - 80 + (i % 10)
        specs.append((start, PEAK_SITE - start + 1,
                      12 if i < N_PEAK_CLIPPED else 0, False))
    # ISOLATED (+): 3 clipped reads at the site; 30 flank reads whose end1
    # (start + SEQ_LEN) lies inside [site - SEQ_LEN, site + 25] but split so
    # no more than 18 reads ever overlap one position (< min_height 20).
    for i in range(N_ISOLATED_CLIP):
        start = ISOLATED_SITE - 80 + i
        specs.append((start, ISOLATED_SITE - start + 1, 12, False))
    for i in range(15):   # end1 in [site-91, site-77]
        specs.append((ISOLATED_SITE - 182 + i, 60, 0, False))
    for i in range(15):   # end1 in [site+14, site+28] -> keep <= site+25
        specs.append((ISOLATED_SITE - 77 + (i % 12), 60, 0, False))
    # BARE (+): coverage-only peak, no clips anywhere near.
    for i in range(N_BARE_RAMP):
        specs.append((BARE_SITE - 350 + i * 4, 80, 0, False))
    for i in range(N_BARE_SUMMIT):
        start = BARE_SITE - 80 + (i % 10)
        specs.append((start, BARE_SITE - start + 1, 0, False))
    # ISOLATED (-): clip at the alignment START (reference_start == site).
    for i in range(N_ISOLATED_CLIP):
        specs.append((ISOLATED_SITE_NEG, 70 + i, 12, True))
    for i in range(15):   # transcript 3' ends (starts) in [site-25, site-11]
        specs.append((ISOLATED_SITE_NEG - 25 + i, 60, 0, True))
    for i in range(15):   # starts in [site+80, site+91]
        specs.append((ISOLATED_SITE_NEG + 80 + (i % 12), 60, 0, True))
    return specs


def _write_bam(path: Path) -> Path:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CHROM, "LN": CHROM_LEN}],
        "RG": [{"ID": "testsample", "SM": "testsample"}],
    }
    records = []
    for n, (start, aligned, clip, reverse) in enumerate(_specs()):
        a = pysam.AlignedSegment()
        a.query_name = f"r{n}"
        if reverse:
            a.query_sequence = "T" * clip + "C" * aligned
            a.cigartuples = ([(4, clip)] if clip else []) + [(0, aligned)]
            a.flag = 16
        else:
            a.query_sequence = "C" * aligned + "A" * clip
            a.cigartuples = [(0, aligned)] + ([(4, clip)] if clip else [])
            a.flag = 0
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = 60
        a.query_qualities = pysam.qualitystring_to_array("I" * len(a.query_sequence))
        a.set_tag("CB", _barcode(n % 37))
        a.set_tag("UB", f"UMI{n:06d}")
        a.set_tag("RG", "testsample")
        records.append(a)
    records.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for rec in records:
            out.write(rec)
    pysam.index(str(path))
    return path


@pytest.fixture(scope="module")
def bam(tmp_path_factory) -> Path:
    return _write_bam(tmp_path_factory.mktemp("tier1_counts") / "synthetic.bam")


@pytest.fixture(autouse=True)
def _config():
    from ema.config import variable_config

    saved = {k: getattr(variable_config, k)
             for k in ("seqlen", "cb_len", "barcode_tag", "ignore_chro",
                       "default_threshold", "merge_len", "read_geometry")}
    variable_config.seqlen = SEQ_LEN
    variable_config.cb_len = 16
    variable_config.barcode_tag = "CB"
    variable_config.ignore_chro = []
    variable_config.default_threshold = 5
    variable_config.merge_len = 100
    # peakAtail-prime: this module's arithmetic ("30 in-window reads + 3 clip
    # = 33") is written against v2's FABRICATED read ends -- every accepted
    # read normalised to start + seq_len.  Under the branch default
    # (--read-geometry true) a synthetic read's end moves upstream by its own
    # unaligned tail, so the same [site - seq_len, site + 25] window catches a
    # different set (18, not 33, at the + ISOLATED locus).  That is the
    # geometry change working, not a counting bug -- so this module keeps
    # pinning the v2 contract and the branch behaviour is asserted in
    # tests/test_read_geometry.py.
    variable_config.read_geometry = "fixed"
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


def _call(bam: Path, out: Path, tag: str, strategy_name: str,
          direction: bool = False, **kwargs):
    from ema.countmatrix.indexing import reset_index
    from ema.countmatrix.peackcalling import peak_calling
    from ema.countmatrix.peak import Peak
    from ema.strategies import get_strategy

    reset_index()
    Peak.reset_pasnumber()
    bed = out / f"{tag}.bed"
    mtx = out / f"{tag}.mtx"
    with patch.object(sys, "argv", ["ema"]):
        peak_calling(
            direction,
            bedfilepath=str(bed),
            matrixpath=str(mtx),
            bamfile_dir=str(bam),
            strategy=get_strategy(strategy_name),
            **kwargs,
        )
    return bed, mtx


def _counts(mtx: Path) -> dict[int, dict[int, int]]:
    out: dict[int, dict[int, int]] = defaultdict(dict)
    for line in mtx.read_text().splitlines():
        if line.strip():
            r, c, n = line.split()
            out[int(r)][int(c)] = int(n)
    return out


def _rows_with_counts(bed: Path, mtx: Path):
    counts = _counts(mtx)
    rows = _rows(bed)
    for r in rows:
        r["cells"] = counts.get(int(r["name"]), {})
        r["total"] = sum(r["cells"].values())
    return rows


def _near(rows, site, window=400):
    return [r for r in rows if r["start"] <= site + window and r["end"] >= site - window]


# ---------------------------------------------------------------------------
# 1. cluster inside a 200-read coverage peak -> ~200 counts, not ~2
# ---------------------------------------------------------------------------

def test_tier1_inside_a_coverage_peak_takes_the_peaks_counts(bam, tmp_path):
    bed, mtx = _call(bam, tmp_path, "seeded", "clip_seeded")
    rows = _rows_with_counts(bed, mtx)
    seeds = [r for r in rows if r["start"] == PEAK_SITE and r["end"] == PEAK_SITE + 1]
    assert seeds, "no tier-1 PAS at the PEAK cleavage site"
    seed = seeds[0]
    assert seed["score"] == N_PEAK_CLIPPED, "score stays the clip-read support"
    assert seed["total"] >= 150, (
        f"tier-1 PAS inside a {N_PEAK_RAMP + N_PEAK_SUMMIT}-read peak carries "
        f"{seed['total']} counts — the clip-only regression"
    )
    assert seed["total"] > 10 * N_PEAK_CLIPPED


def test_tier1_conserves_the_coverage_strategys_mass_at_the_locus(bam, tmp_path):
    """Everything the coverage strategy counts at the PEAK locus must land on
    the tier-1 PAS that suppressed it — per cell — plus the reads at the
    cleavage site that the coverage candidate never covered (the summit
    reads end past ``max(peak_list)``, the ~seq_len tail the streaming loop
    drops when a peak closes)."""
    bed_s, mtx_s = _call(bam, tmp_path, "cons_seeded", "clip_seeded")
    bed_c, mtx_c = _call(bam, tmp_path, "cons_cov", "lambda_gradient")
    seeded = _near(_rows_with_counts(bed_s, mtx_s), PEAK_SITE)
    cov = _near(_rows_with_counts(bed_c, mtx_c), PEAK_SITE)
    assert cov, "lambda_gradient must call the PEAK locus"
    assert all(r["score"] >= 1 for r in seeded), "PEAK locus must be tier-1 only"
    assert max(r["end"] for r in cov) < PEAK_SITE, (
        "fixture: the coverage candidate must stop short of the summit reads"
    )

    def per_cb(rows):
        agg: dict[int, int] = defaultdict(int)
        for r in rows:
            for c, n in r["cells"].items():
                agg[c] += n
        return dict(agg)

    s_cb, c_cb = per_cb(seeded), per_cb(cov)
    assert all(s_cb.get(cb, 0) >= n for cb, n in c_cb.items()), (
        "a cell's coverage count must never shrink on the tier-1 row"
    )
    s_total, c_total = sum(r["total"] for r in seeded), sum(r["total"] for r in cov)
    assert s_total >= c_total + N_PEAK_SUMMIT, (
        f"the {N_PEAK_SUMMIT} summit reads at the cleavage site must be counted "
        f"(seeded {s_total} vs coverage {c_total})"
    )
    assert s_total <= N_PEAK_RAMP + N_PEAK_SUMMIT, "no read counted twice"


# ---------------------------------------------------------------------------
# 2. cluster outside every coverage peak -> 30 in-window reads + 3 clip = 33
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("direction,site", [(False, ISOLATED_SITE), (True, ISOLATED_SITE_NEG)])
def test_tier1_outside_any_peak_counts_every_read_end_in_the_window(bam, tmp_path, direction, site):
    tag = "iso_neg" if direction else "iso_pos"
    bed_c, _ = _call(bam, tmp_path, f"{tag}_cov", "lambda_gradient", direction=direction)
    assert not _near(_rows(bed_c), site), (
        "fixture error: the ISOLATED locus must produce no coverage candidate"
    )
    bed, mtx = _call(bam, tmp_path, tag, "clip_seeded", direction=direction)
    rows = _rows_with_counts(bed, mtx)
    seeds = [r for r in rows if r["start"] == site and r["end"] == site + 1]
    assert seeds, f"no tier-1 PAS at the ISOLATED site on {'-' if direction else '+'}"
    assert seeds[0]["score"] == N_ISOLATED_CLIP
    assert seeds[0]["total"] == N_ISOLATED_CLIP + N_ISOLATED_FLANK, (
        f"expected {N_ISOLATED_CLIP + N_ISOLATED_FLANK} read ends in the count "
        f"window, got {seeds[0]['total']}"
    )


def test_count_window_flag_narrows_the_outside_count(bam, tmp_path):
    """With a 0,0 window only read ends exactly at the anchor count — on the
    synthetic + locus no read end sits at the site itself, so the row falls
    back to the 3 clip reads; a generous window recovers all 33."""
    bed, mtx = _call(bam, tmp_path, "win0", "clip_seeded", polya_count_window=(0, 0))
    rows = _rows_with_counts(bed, mtx)
    (seed,) = [r for r in rows if r["start"] == ISOLATED_SITE]
    assert seed["total"] == N_ISOLATED_CLIP
    bed, mtx = _call(bam, tmp_path, "winbig", "clip_seeded", polya_count_window="200,50")
    rows = _rows_with_counts(bed, mtx)
    (seed,) = [r for r in rows if r["start"] == ISOLATED_SITE]
    assert seed["total"] == N_ISOLATED_CLIP + N_ISOLATED_FLANK


# ---------------------------------------------------------------------------
# 3. tier-2 rows are byte-for-byte the coverage strategy's rows
# ---------------------------------------------------------------------------

def test_tier2_rows_are_unchanged_from_the_coverage_strategy(bam, tmp_path):
    bed_s, mtx_s = _call(bam, tmp_path, "t2_seeded", "clip_seeded")
    bed_c, mtx_c = _call(bam, tmp_path, "t2_cov", "lambda_gradient")
    seeded = _rows_with_counts(bed_s, mtx_s)
    cov = {(r["start"], r["end"]): r for r in _rows_with_counts(bed_c, mtx_c)}
    tier2 = [r for r in seeded if r["score"] == 0]
    assert tier2 and _near(tier2, BARE_SITE), "BARE locus must be tier 2"
    for r in tier2:
        c = cov.get((r["start"], r["end"]))
        assert c is not None, f"tier-2 row {r['start']}-{r['end']} not a coverage PAS"
        assert r["cells"] == c["cells"], "tier-2 counts must be the coverage counts"


# ---------------------------------------------------------------------------
# coordinates never move; totals are >= the coverage strategy's
# ---------------------------------------------------------------------------

def test_counting_change_leaves_coordinates_and_scores_untouched(bam, tmp_path):
    """Counting is orthogonal to calling: the BED must be identical whatever
    count window is in force."""
    bed_a, _ = _call(bam, tmp_path, "coord_a", "clip_seeded")
    bed_b, _ = _call(bam, tmp_path, "coord_b", "clip_seeded", polya_count_window=(0, 0))
    bed_c, _ = _call(bam, tmp_path, "coord_c", "clip_seeded", polya_count_window="300,300")
    assert bed_a.read_bytes() == bed_b.read_bytes() == bed_c.read_bytes()


def test_every_tier1_row_has_counts_and_mass_is_not_lost(bam, tmp_path):
    bed_s, mtx_s = _call(bam, tmp_path, "mass_seeded", "clip_seeded")
    bed_c, mtx_c = _call(bam, tmp_path, "mass_cov", "lambda_gradient")
    seeded = _rows_with_counts(bed_s, mtx_s)
    cov = _rows_with_counts(bed_c, mtx_c)
    assert all(r["total"] > 0 for r in seeded), "a tier-1 row must never be empty"
    assert sum(r["total"] for r in seeded) >= sum(r["total"] for r in cov)


# ---------------------------------------------------------------------------
# ClipSeeder unit level: partition + no double counting
# ---------------------------------------------------------------------------

class _FakePeak:
    def __init__(self, cb_positions):
        self.cb_positions = cb_positions
        self.cb_dict = {}
        for d in cb_positions.values():
            for cb, n in d.items():
                self.cb_dict[cb] = self.cb_dict.get(cb, 0) + n


def _reconstruct(peak, s, e):
    from ema.strategies.utils import reconstruct_cb_dict
    return reconstruct_cb_dict(peak.cb_positions, s, e)


def test_seeder_splits_a_candidate_between_two_clusters_at_the_midpoint():
    from ema.countmatrix.polya import ClipSeeder

    s = ClipSeeder(direction=False, window=100)
    for i in range(3):
        s.add_clip(1000, cb="A", umi=f"a{i}")
    for i in range(2):
        s.add_clip(1300, cb="B", umi=f"b{i}")
    peak = _FakePeak({
        1000: {"A": 5}, 1100: {"A": 4}, 1150: {"X": 1},   # <= midpoint 1150
        1151: {"Y": 1}, 1250: {"B": 6}, 1300: {"B": 7},   # > midpoint
    })
    s.add_coverage_pas(900, 1400, dict(peak.cb_dict), peak)
    recs = s.flush()
    assert [(r[0], r[1]) for r in recs] == [(1000, 1001), (1300, 1301)]
    assert recs[0][3] == {"A": 9, "X": 1}
    assert recs[1][3] == {"Y": 1, "B": 13}
    assert recs[0][2] == 3 and recs[1][2] == 2, (
        "score stays the clip support (in molecules)"
    )


def test_seeder_does_not_invent_reads_for_a_non_positional_candidate():
    """``original`` hands back the whole peak's cb_dict for any interval; the
    shares then do not sum to the candidate and the mass stays whole on the
    nearest cluster instead of being duplicated."""
    from ema.countmatrix.polya import ClipSeeder

    s = ClipSeeder(direction=False, window=100)
    s.add_clip(1000, cb="A", umi="a")
    s.add_clip(1300, cb="B", umi="b")
    peak = _FakePeak({1000: {"A": 5}, 1300: {"B": 7}, 1600: {"C": 9}})
    # 'original'-style cb_dict: the whole peak, wider than the candidate
    s.add_coverage_pas(900, 1400, dict(peak.cb_dict), peak)
    recs = s.flush()
    totals = sorted(sum(r[3].values()) for r in recs)
    # whole candidate (21) on the nearest cluster; the other cluster has no
    # read ends of its own (no add_read calls) so it falls back to its 1
    # clip read — never a duplicated share of the candidate
    assert totals == [1, 21], totals


def test_seeder_single_cluster_takes_the_whole_candidate_even_without_a_peak():
    from ema.countmatrix.polya import ClipSeeder

    s = ClipSeeder(direction=False, window=100)
    s.add_clip(1000, cb="CB1", umi="u")
    s.add_coverage_pas(900, 1100, {"CB1": 7, "CB2": 3})
    (rec,) = s.flush()
    assert rec[3] == {"CB1": 7, "CB2": 3}


def test_seeder_window_count_excludes_candidate_intervals_and_neighbours():
    """A read end inside a coverage candidate belongs to that candidate
    (tier 2 or partitioned) and must not also be counted by a nearby
    outside cluster; two outside clusters split the gap at the midpoint."""
    from ema.countmatrix.polya import ClipSeeder

    s = ClipSeeder(direction=False, window=100, count_window=(300, 25))
    # reads: end1 = start1 + 91
    for end, cb, n in ((1150, "IN", 4), (1300, "OUT", 5), (1420, "FAR", 2),
                       (1500, "MID", 3), (1600, "NEXT", 6)):
        for _ in range(n):
            s.add_read(end - 91, end, cb)
    s.add_clip(1400, cb="OUT", umi="u1")       # outside: 1400 > 1200 + 100
    s.add_clip(1700, cb="NEXT", umi="u2")      # second outside cluster
    s.add_coverage_pas(1000, 1200, {"IN": 4})  # tier-2 candidate
    recs = s.flush()
    by_start = {r[0]: r for r in recs}
    assert by_start[1000][3] == {"IN": 4}
    # cluster 1400: window [1100, 1425] minus [1000,1200] -> [1201,1425],
    # clipped at midpoint to 1700 (1550) -> [1201, 1425]
    assert by_start[1400][3] == {"OUT": 5, "FAR": 2}
    # cluster 1700: window [1400, 1725] clipped to [1551, 1725]
    assert by_start[1700][3] == {"NEXT": 6}
    assert sum(sum(r[3].values()) for r in recs) == 4 + 5 + 2 + 6
    assert s.stats["reads_tier2"] == 4
    assert s.stats["reads_from_window"] == 5 + 2 + 6
    assert s.stats["clip_fallback_rows"] == 0


def test_seeder_inside_cluster_also_takes_the_peaks_uncovered_tail():
    """A cluster that suppresses a candidate gets the candidate's partition
    share AND the read ends in its window that no candidate covers (the
    dropped tail), never the candidate's reads twice."""
    from ema.countmatrix.polya import ClipSeeder

    s = ClipSeeder(direction=False, window=100, count_window=(91, 25))
    # candidate [1000, 1200] built from these positions ...
    peak = _FakePeak({1000: {"A": 3}, 1100: {"A": 2}, 1200: {"B": 4},
                      1250: {"TAIL": 5}, 1310: {"TAIL": 1}})
    for pos, d in peak.cb_positions.items():
        for cb, n in d.items():
            for _ in range(n):
                s.add_read(pos - 91, pos, cb)
    s.add_clip(1290, cb="TAIL", umi="u")     # inside [1000-100, 1200+100]
    s.add_coverage_pas(1000, 1200, _reconstruct(peak, 1000, 1200), peak)
    (rec,) = s.flush()
    assert rec[:2] == (1290, 1291)
    # partition share {A:5, B:4} + window [1199, 1315] minus [1000,1200]
    # -> [1201, 1315] -> TAIL 5 + 1
    assert rec[3] == {"A": 5, "B": 4, "TAIL": 6}
    assert s.stats["reads_from_candidates"] == 9
    assert s.stats["reads_from_window"] == 6


def test_seeder_minus_strand_window_is_shifted_by_seq_len():
    from ema.countmatrix.polya import ClipSeeder

    s = ClipSeeder(direction=True, window=100, count_window=(-1, 25))
    seq_len = 91
    site = 5000
    # clipped read on '-': start == site, end1 == site + seq_len
    s.add_read(site, site + seq_len, "C")
    s.add_clip(site, cb="C", umi="u")
    # transcript-upstream read ends (starts > site) within seq_len
    s.add_read(site + 60, site + 60 + seq_len, "UP")
    # transcript-downstream within 25 bp
    s.add_read(site - 20, site - 20 + seq_len, "DOWN")
    # too far upstream (start > site + seq_len)
    s.add_read(site + 120, site + 120 + seq_len, "TOOFAR")
    (rec,) = s.flush()
    assert rec[:2] == (site, site + 1)
    assert rec[3] == {"C": 1, "UP": 1, "DOWN": 1}


def test_parse_count_window_forms():
    from ema.countmatrix.polya import parse_count_window

    assert parse_count_window("auto,25") == (-1, 25)
    assert parse_count_window("98,25") == (98, 25)
    assert parse_count_window("30") == (-1, 30)
    assert parse_count_window((40, 10)) == (40, 10)
    with pytest.raises(ValueError):
        parse_count_window("1,2,3")
    with pytest.raises(ValueError):
        parse_count_window("10,-5")
