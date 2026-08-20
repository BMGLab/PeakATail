"""Unit tests for the clip-cluster core of the ``clip_seeded`` strategy.

``cluster_clip_sites`` and ``ClipSeeder`` are the Phase-2 machinery ported
from the scoring prototype (``denovo.py``): single-linkage clustering at
``--polya-seed-window`` with a read-weighted modal position, plus the
two-tier emission that keeps coverage-only peaks as an explicit second tier.

These tests pin the properties the F1 measurement depends on:
  * the emitted position is the read-weighted MODE, not the cluster mean or
    edge (the prototype measured P@25 0.34 -> 0.73 from this alone);
  * linkage is single-linkage on the gap, so a chain of close sites is one
    cluster however long the chain is;
  * ties resolve deterministically (multiprocessing-safe, byte-identical
    reruns);
  * BED intervals are end-exclusive and strand-correct.
"""
from __future__ import annotations

from ema.countmatrix.polya import (
    ClipAccumulator,
    ClipSeeder,
    cluster_clip_sites,
)


def _sites(spec: dict) -> dict:
    """``{site: n_reads}`` -> accumulator layout, one UMI per read."""
    acc = ClipAccumulator()
    for site, n in spec.items():
        for i in range(n):
            acc.add(site, cb=f"CB{i}", umi=f"U{site}_{i}")
    return acc.sites


# ---------------------------------------------------------------------------
# cluster_clip_sites
# ---------------------------------------------------------------------------

def test_empty_input_returns_empty():
    assert cluster_clip_sites({}, 25) == []


def test_single_site_returns_itself():
    out = cluster_clip_sites(_sites({1000: 3}), 25)
    assert len(out) == 1
    assert out[0]["mode"] == 1000
    assert out[0]["nreads"] == 3
    assert out[0]["numis"] == 3


def test_sites_within_window_merge_and_report_the_read_weighted_mode():
    out = cluster_clip_sites(_sites({1000: 2, 1010: 9, 1020: 1}), 25)
    assert len(out) == 1
    assert out[0]["mode"] == 1010, "position must be the highest-read site"
    assert out[0]["nreads"] == 12
    assert out[0]["numis"] == 12


def test_gap_larger_than_window_splits_clusters():
    out = cluster_clip_sites(_sites({1000: 5, 1026: 3}), 25)
    assert [c["mode"] for c in out] == [1000, 1026]


def test_gap_exactly_at_window_still_merges():
    """The linkage test is ``gap > window`` -> split, so gap == window links."""
    out = cluster_clip_sites(_sites({1000: 5, 1025: 3}), 25)
    assert len(out) == 1
    assert out[0]["mode"] == 1000


def test_single_linkage_chains_beyond_the_window():
    """A chain of 20 bp steps is ONE cluster even though its span is 100 bp —
    that is what single linkage means, and it is what the prototype did."""
    out = cluster_clip_sites(
        _sites({1000: 1, 1020: 1, 1040: 9, 1060: 1, 1080: 1, 1100: 1}), 25
    )
    assert len(out) == 1
    assert out[0]["mode"] == 1040
    assert out[0]["nreads"] == 14


def test_seed_window_parameter_is_honoured():
    spec = {1000: 5, 1030: 3}
    assert len(cluster_clip_sites(_sites(spec), 25)) == 2
    assert len(cluster_clip_sites(_sites(spec), 50)) == 1


def test_mode_tie_breaks_to_the_lowest_coordinate_deterministically():
    for _ in range(5):
        out = cluster_clip_sites(_sites({1000: 4, 1010: 4}), 25)
        assert out[0]["mode"] == 1000


def test_clusters_come_out_sorted_by_position():
    out = cluster_clip_sites(_sites({5000: 2, 1000: 2, 3000: 2}), 25)
    assert [c["mode"] for c in out] == [1000, 3000, 5000]


def test_umi_dedup_collapses_a_pcr_stack_but_reads_do_not():
    """A duplicate stack is exactly the sharp single-coordinate spike the
    coverage caller mistakes for a PAS; molecule counts must not inflate."""
    acc = ClipAccumulator()
    for _ in range(50):
        acc.add(1000, cb="CB1", umi="SAME")
    out = cluster_clip_sites(acc.sites, 25)
    assert out[0]["nreads"] == 50
    assert out[0]["numis"] == 1


def test_same_umi_in_different_cells_counts_as_distinct_molecules():
    acc = ClipAccumulator()
    acc.add(1000, cb="CB1", umi="SAME")
    acc.add(1000, cb="CB2", umi="SAME")
    out = cluster_clip_sites(acc.sites, 25)
    assert out[0]["numis"] == 2


def test_reads_without_a_umi_each_count_as_one_molecule():
    acc = ClipAccumulator()
    for _ in range(4):
        acc.add(1000, cb="CB1", umi=None)
    out = cluster_clip_sites(acc.sites, 25)
    assert out[0]["nreads"] == 4
    assert out[0]["numis"] == 4, "no UB tag -> conservative: one molecule/read"


def test_cb_counts_are_summed_across_the_cluster():
    acc = ClipAccumulator()
    acc.add(1000, cb="CB1", umi="a")
    acc.add(1000, cb="CB1", umi="b")
    acc.add(1010, cb="CB2", umi="c")
    out = cluster_clip_sites(acc.sites, 25)
    assert out[0]["cb_counts"] == {"CB1": 2, "CB2": 1}


# ---------------------------------------------------------------------------
# ClipSeeder — two-tier emission
# ---------------------------------------------------------------------------

def test_seeder_emits_one_bp_end_exclusive_intervals():
    s = ClipSeeder(direction=False)
    for i in range(3):
        s.add_clip(1000, cb="CB1", umi=f"u{i}")
    (start, end, score, cbs), = s.flush()
    assert (start, end) == (1000, 1001)
    assert score == 3


def test_seeder_suppresses_a_coverage_peak_overlapping_a_cluster():
    s = ClipSeeder(direction=False, window=100)
    s.add_clip(1000, cb="CB1", umi="u")
    s.add_coverage_pas(900, 1100, {"CB1": 7})
    recs = s.flush()
    assert [(r[0], r[1]) for r in recs] == [(1000, 1001)], (
        "the coverage peak covering the cluster must not be emitted twice"
    )


def test_seeder_keeps_a_distant_coverage_peak_as_tier_two():
    s = ClipSeeder(direction=False, window=100)
    s.add_clip(1000, cb="CB1", umi="u")
    s.add_coverage_pas(50_000, 50_200, {"CB1": 7})
    recs = s.flush()
    assert len(recs) == 2
    tier2 = [r for r in recs if r[0] == 50_000][0]
    assert tier2[2] == 0, "coverage-only peaks must score 0 (the tier tag)"
    assert tier2[3] == {"CB1": 7}


def test_seeder_output_is_coordinate_sorted():
    s = ClipSeeder(direction=False, window=100)
    for pos in (9_000, 1_000, 5_000):
        s.add_clip(pos, cb="CB1", umi=f"u{pos}")
    s.add_coverage_pas(30_000, 30_100, {"CB1": 1})
    s.add_coverage_pas(20_000, 20_100, {"CB1": 1})
    starts = [r[0] for r in s.flush()]
    assert starts == sorted(starts)


def test_seeder_min_umis_gate_uses_molecules_not_reads():
    s = ClipSeeder(direction=False, min_reads=2)
    for _ in range(10):
        s.add_clip(1000, cb="CB1", umi="SAME")  # 10 reads, 1 molecule
    assert s.flush() == [], "10 PCR duplicates are not 2 molecules"

    s = ClipSeeder(direction=False, min_reads=2)
    s.add_clip(1000, cb="CB1", umi="a")
    s.add_clip(1000, cb="CB1", umi="b")
    assert len(s.flush()) == 1


def test_seeder_flush_resets_state():
    s = ClipSeeder(direction=False)
    s.add_clip(1000, cb="CB1", umi="u")
    s.add_coverage_pas(50_000, 50_100, {"CB1": 1})
    assert len(s.flush()) == 2
    assert s.flush() == [], "flush must clear the chromosome's state"


def test_minus_strand_tier_two_support_is_measured_at_the_bed_start():
    """On the minus strand the 3' base is bed_start, so a clip cluster just
    upstream of bed_end must NOT be credited as support."""
    near = ClipSeeder(direction=True, window=50, min_reads=99)
    near.add_clip(1000, cb="CB1", umi="u")      # sub-threshold: stays tier 2
    near.add_coverage_pas(1010, 1400, {"CB1": 1})
    (_, _, score_near, _), = near.flush()
    assert score_near == 1, "clip 10 bp from bed_start must count on '-'"

    far = ClipSeeder(direction=True, window=50, min_reads=99)
    far.add_clip(1390, cb="CB1", umi="u")
    far.add_coverage_pas(1010, 1400, {"CB1": 1})
    (_, _, score_far, _), = far.flush()
    assert score_far == 0, "clip near bed_end is the 5' side on '-'"


def test_plus_strand_tier_two_support_is_measured_at_the_bed_end():
    near = ClipSeeder(direction=False, window=50, min_reads=99)
    near.add_clip(1390, cb="CB1", umi="u")
    near.add_coverage_pas(1010, 1400, {"CB1": 1})
    (_, _, score_near, _), = near.flush()
    assert score_near == 1, "clip 10 bp from bed_end - 1 must count on '+'"

    far = ClipSeeder(direction=False, window=50, min_reads=99)
    far.add_clip(1010, cb="CB1", umi="u")
    far.add_coverage_pas(1010, 1400, {"CB1": 1})
    (_, _, score_far, _), = far.flush()
    assert score_far == 0, "clip at bed_start is the 5' side on '+'"


def test_at_default_min_reads_every_tier_two_peak_scores_exactly_zero():
    """The tier tag must be exact: with min_reads=1 no cluster is rejected,
    so a surviving coverage-only peak can have no clip support at all —
    otherwise ``score == 0`` would not mean coverage_only and
    ``--polya-mode filter`` would keep the wrong rows."""
    s = ClipSeeder(direction=False, window=100, min_reads=1)
    # a clip 120 bp past the coverage peak's 3' base: far enough that the
    # cluster mode escapes the expanded interval test on some geometries,
    # but close enough to be counted as support by the +/-window rule
    s.add_clip(1499, cb="CB1", umi="u")
    s.add_coverage_pas(1000, 1400, {"CB1": 3})
    recs = s.flush()
    tier2 = [r for r in recs if r[1] - r[0] > 1]
    assert all(r[2] == 0 for r in tier2), (
        f"tier-2 peaks must score 0 at min_reads=1, got {tier2}"
    )


def test_tier_two_may_score_above_zero_only_via_a_rejected_cluster():
    """Raising the gate is the one case where a coverage-only peak keeps a
    non-zero score — it does have clip evidence, just not enough for its own
    call, so --polya-mode filter is right to keep it."""
    s = ClipSeeder(direction=False, window=100, min_reads=5)
    s.add_clip(1390, cb="CB1", umi="u")      # 1 molecule -> cluster rejected
    s.add_coverage_pas(1000, 1400, {"CB1": 3})
    (start, end, score, _), = s.flush()
    assert (start, end) == (1000, 1400)
    assert score == 1


def test_support_window_counts_molecules_not_reads():
    """Stage 1c: the tier-2 score is in the same unit as the tier-1 score and
    the gate — distinct (CB, UMI) molecules.  Five reads from two molecules
    score 2, not 5 (the sidecar keeps the raw 5)."""
    s = ClipSeeder(direction=False, window=100, min_reads=99)
    for i in range(5):
        s.add_clip(1350, cb=f"CB{i % 2}", umi=f"u{i % 2}")
    s.add_coverage_pas(1300, 1400, {"CB1": 1})
    sup = []
    (_, _, score, _), = s.flush(support_out=sup)
    assert score == 2
    assert sup[0]["clip_reads"] == 5 and sup[0]["clip_umis"] == 2
