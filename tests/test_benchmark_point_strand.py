"""Tests for point-mode strand-matched benchmark scoring (ema.benchmark.metrics).

Coverage targets:
  - strand-aware 3'-most-base reduction; BED end is EXCLUSIVE
    ('+' -> end-1, '-' -> start); a 1-length interval at position p yields
    point p on BOTH strands
  - exact matching on '+' and '-' strands with synthetic mini-BEDs
  - strand-mismatch rejection: opposite-strand reference at distance 0
    must NOT match
  - cutoff boundary: distance exactly == cutoff matches (inclusive), one
    beyond does not
  - entries with strand not in {+, -} are excluded and counted
  - honest naming: the recall-like quantity is "reference_coverage";
    restricted-reference recall is "recall_restricted" with the restriction
    recorded
  - shuffled-null baseline: determinism (same seed -> identical result
    twice), per-seed + mean reporting, width/chrom/strand-preserving
    placement, gene-body constraint
  - compute_metrics: legacy keys stay byte-compatible when the new blocks
    are enabled (requires bedtools; skipped when unavailable)
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pytest

from ema.benchmark.metrics import (
    _intervals_to_points,
    _load_stranded_intervals,
    _shuffle_records,
    compute_null_baseline,
    compute_point_strand_metrics,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_bed(path: Path, rows: list[str]) -> str:
    path.write_text("".join(r + "\n" for r in rows))
    return str(path)


def _bed6(chrom: str, start: int, end: int, name: str, strand: str) -> str:
    return f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}"


# ---------------------------------------------------------------------------
# Point reduction: BED end-exclusive off-by-one
# ---------------------------------------------------------------------------


def test_one_length_interval_yields_point_p_on_both_strands():
    p = 1234
    records = [("chr1", p, p + 1, "+"), ("chr1", p, p + 1, "-")]
    points = _intervals_to_points(records)
    assert list(points[("chr1", "+")]) == [p]
    assert list(points[("chr1", "-")]) == [p]


def test_multi_base_interval_reduces_to_strand_aware_3prime_base():
    records = [("chr1", 100, 150, "+"), ("chr1", 100, 150, "-")]
    points = _intervals_to_points(records)
    # '+' strand: 3'-most covered base is end-1 (BED end is exclusive)
    assert list(points[("chr1", "+")]) == [149]
    # '-' strand: 3'-most covered base is start
    assert list(points[("chr1", "-")]) == [100]


# ---------------------------------------------------------------------------
# Strand handling and exclusion counting
# ---------------------------------------------------------------------------


def test_entries_without_valid_strand_are_excluded_and_counted(tmp_path):
    bed = _write_bed(tmp_path / "p.bed", [
        _bed6("chr1", 100, 101, "ok_plus", "+"),
        _bed6("chr1", 200, 201, "dot_strand", "."),
        "chr1\t300\t301",  # 3-column line: no strand field at all
        _bed6("chr1", 400, 401, "ok_minus", "-"),
    ])
    records, n_excluded = _load_stranded_intervals(bed)
    assert len(records) == 2
    assert n_excluded == 2


def test_exclusion_count_reported_in_block(tmp_path):
    pred = _write_bed(tmp_path / "p.bed", [
        _bed6("chr1", 100, 101, "a", "+"),
        _bed6("chr1", 200, 201, "b", "."),
    ])
    ref = _write_bed(tmp_path / "r.bed", [
        _bed6("chr1", 100, 101, "r1", "+"),
    ])
    block = compute_point_strand_metrics(pred, ref, distance_cutoffs=[10])
    assert block["n_predicted_points"] == 1
    assert block["n_predicted_excluded_no_strand"] == 1
    assert block["n_reference_excluded_no_strand"] == 0


# ---------------------------------------------------------------------------
# Matching: both strands, mismatch rejection, cutoff boundary
# ---------------------------------------------------------------------------


def test_exact_match_on_plus_and_minus_strands(tmp_path):
    # '+' prediction as a wide interval whose 3' end (end-1) hits the ref;
    # '-' prediction as a wide interval whose 3' end (start) hits the ref.
    pred = _write_bed(tmp_path / "p.bed", [
        _bed6("chr1", 950, 1001, "plus_peak", "+"),    # point = 1000
        _bed6("chr2", 2000, 2051, "minus_peak", "-"),  # point = 2000
    ])
    ref = _write_bed(tmp_path / "r.bed", [
        _bed6("chr1", 1000, 1001, "ref_plus", "+"),    # point = 1000
        _bed6("chr2", 2000, 2001, "ref_minus", "-"),   # point = 2000
    ])
    block = compute_point_strand_metrics(pred, ref, distance_cutoffs=[0])
    m = block["cutoffs"][0]
    assert m["precision"] == 1.0
    assert m["reference_coverage"] == 1.0
    assert m["matched_predicted"] == 2
    assert m["matched_reference"] == 2


def test_opposite_strand_reference_at_distance_zero_does_not_match(tmp_path):
    pred = _write_bed(tmp_path / "p.bed", [
        _bed6("chr1", 1000, 1001, "pred", "+"),   # point = 1000 on '+'
    ])
    ref = _write_bed(tmp_path / "r.bed", [
        _bed6("chr1", 1000, 1001, "ref", "-"),    # point = 1000 on '-'
    ])
    block = compute_point_strand_metrics(pred, ref, distance_cutoffs=[0, 5000])
    for cutoff in (0, 5000):
        m = block["cutoffs"][cutoff]
        assert m["precision"] == 0.0
        assert m["reference_coverage"] == 0.0
        assert m["matched_predicted"] == 0
        assert m["matched_reference"] == 0


def test_distance_exactly_at_cutoff_matches_but_one_beyond_does_not(tmp_path):
    cutoff = 50
    pred = _write_bed(tmp_path / "p.bed", [
        _bed6("chr1", 1000, 1001, "pred", "+"),            # point = 1000
    ])
    ref = _write_bed(tmp_path / "r.bed", [
        _bed6("chr1", 1000 + cutoff, 1001 + cutoff, "ref", "+"),  # dist = 50
    ])
    block = compute_point_strand_metrics(pred, ref,
                                         distance_cutoffs=[cutoff - 1, cutoff])
    assert block["cutoffs"][cutoff]["precision"] == 1.0        # dist == cutoff
    assert block["cutoffs"][cutoff - 1]["precision"] == 0.0    # dist > cutoff


def test_reference_coverage_key_replaces_recall_in_new_block(tmp_path):
    pred = _write_bed(tmp_path / "p.bed", [_bed6("chr1", 10, 11, "a", "+")])
    ref = _write_bed(tmp_path / "r.bed", [_bed6("chr1", 10, 11, "r", "+")])
    block = compute_point_strand_metrics(pred, ref, distance_cutoffs=[10])
    m = block["cutoffs"][10]
    assert "reference_coverage" in m
    assert "recall" not in m


# ---------------------------------------------------------------------------
# Restricted-reference recall
# ---------------------------------------------------------------------------


def test_recall_restricted_uses_restricted_reference_and_records_it(tmp_path):
    pred = _write_bed(tmp_path / "p.bed", [
        _bed6("chr1", 1000, 1001, "a", "+"),
    ])
    # Full reference: 2 sites, only one near the prediction -> coverage 0.5
    ref = _write_bed(tmp_path / "r.bed", [
        _bed6("chr1", 1000, 1001, "near", "+"),
        _bed6("chr1", 900000, 900001, "far", "+"),
    ])
    # Cohort-restricted reference: only the relevant site -> recall 1.0
    restricted = _write_bed(tmp_path / "restricted.bed", [
        _bed6("chr1", 1000, 1001, "near", "+"),
    ])
    block = compute_point_strand_metrics(pred, ref, distance_cutoffs=[10],
                                         restricted_reference_bed=restricted)
    m = block["cutoffs"][10]
    assert m["reference_coverage"] == 0.5
    assert m["recall_restricted"] == 1.0
    assert m["matched_reference_restricted"] == 1
    rr = block["restricted_reference"]
    assert rr["bed"] == restricted
    assert rr["n_reference_points"] == 1


# ---------------------------------------------------------------------------
# Shuffled-null baseline
# ---------------------------------------------------------------------------


def _null_beds(tmp_path):
    pred = _write_bed(tmp_path / "p.bed", [
        _bed6("chr1", 100, 150, f"p{i}", "+" if i % 2 else "-")
        for i in range(10)
    ] + [
        _bed6("chr2", 500 + 10 * i, 530 + 10 * i, f"q{i}", "+")
        for i in range(5)
    ])
    ref = _write_bed(tmp_path / "r.bed", [
        _bed6("chr1", 120, 121, "r1", "+"),
        _bed6("chr1", 5000, 5001, "r2", "-"),
        _bed6("chr2", 600, 601, "r3", "+"),
    ])
    return pred, ref


def test_null_shuffle_same_seed_is_deterministic(tmp_path):
    pred, ref = _null_beds(tmp_path)
    a = compute_null_baseline(pred, ref, distance_cutoffs=[50, 500],
                              n_seeds=3, base_seed=7)
    b = compute_null_baseline(pred, ref, distance_cutoffs=[50, 500],
                              n_seeds=3, base_seed=7)
    assert a == b
    assert a["seeds"] == [7, 8, 9]
    assert len(a["per_seed"]) == 3
    assert set(a["mean"]["cutoffs"]) == {50, 500}
    for entry in a["per_seed"]:
        for cutoff in (50, 500):
            assert "precision" in entry["cutoffs"][cutoff]
            assert "reference_coverage" in entry["cutoffs"][cutoff]


def test_null_mean_is_mean_of_per_seed(tmp_path):
    pred, ref = _null_beds(tmp_path)
    out = compute_null_baseline(pred, ref, distance_cutoffs=[500],
                                n_seeds=3, base_seed=0)
    per_seed = [s["cutoffs"][500]["precision"] for s in out["per_seed"]]
    assert out["mean"]["cutoffs"][500]["precision"] == pytest.approx(
        round(float(np.mean(per_seed)), 4))


def test_shuffle_preserves_width_chromosome_and_strand():
    records = [("chr1", 100, 150, "+"), ("chr2", 40, 45, "-")]
    rng = np.random.default_rng(0)
    shuffled, n_unplaced = _shuffle_records(
        records, rng, chrom_extent={"chr1": 100000, "chr2": 100000})
    assert n_unplaced == 0
    assert len(shuffled) == 2
    for (chrom, start, end, strand), (ochrom, ostart, oend, ostrand) in zip(
            shuffled, records):
        assert chrom == ochrom
        assert strand == ostrand
        assert end - start == oend - ostart
        assert start >= 0


def test_shuffle_gene_body_constraint_places_inside_gene_bodies(tmp_path):
    gene_body = {"chr1": [(1000, 2000), (5000, 5100)]}
    records = [("chr1", 0, 50, "+")] * 200
    rng = np.random.default_rng(1)
    shuffled, n_unplaced = _shuffle_records(records, rng, gene_body=gene_body)
    assert n_unplaced == 0
    assert len(shuffled) == 200
    for chrom, start, end, strand in shuffled:
        assert chrom == "chr1"
        assert end - start == 50
        assert any(gs <= start and end <= ge for gs, ge in gene_body["chr1"])


def test_shuffle_drops_and_counts_unplaceable_intervals():
    # Wider than every gene-body interval on its chromosome, and a
    # chromosome absent from the placement space entirely.
    gene_body = {"chr1": [(0, 30)]}
    records = [("chr1", 0, 50, "+"), ("chrX", 0, 10, "+")]
    rng = np.random.default_rng(0)
    shuffled, n_unplaced = _shuffle_records(records, rng, gene_body=gene_body)
    assert shuffled == []
    assert n_unplaced == 2


def test_null_records_gene_body_bed_and_placement_mode(tmp_path):
    pred, ref = _null_beds(tmp_path)
    gb = _write_bed(tmp_path / "genes.bed", [
        "chr1\t0\t100000", "chr2\t0\t100000",
    ])
    out = compute_null_baseline(pred, ref, distance_cutoffs=[50],
                                n_seeds=2, base_seed=0, gene_body_bed=gb)
    assert out["gene_body_bed"] == gb
    assert out["placement"] == "gene_body_constrained"
    unconstrained = compute_null_baseline(pred, ref, distance_cutoffs=[50],
                                          n_seeds=2, base_seed=0)
    assert unconstrained["gene_body_bed"] is None
    assert unconstrained["placement"] == "chromosome_extent"


# ---------------------------------------------------------------------------
# compute_metrics integration: legacy keys stay byte-compatible
# ---------------------------------------------------------------------------


needs_bedtools = pytest.mark.skipif(
    shutil.which("bedtools") is None, reason="bedtools binary not available")


@needs_bedtools
def test_compute_metrics_legacy_block_unchanged_when_new_blocks_enabled(tmp_path):
    pytest.importorskip("pybedtools")
    from ema.benchmark.metrics import compute_metrics

    pred = _write_bed(tmp_path / "p.bed", [
        _bed6("chr1", 950, 1001, "p1", "+"),
        _bed6("chr1", 3000, 3051, "p2", "-"),
    ])
    ref = _write_bed(tmp_path / "r.bed", [
        _bed6("chr1", 1000, 1001, "r1", "+"),
        _bed6("chr1", 9000, 9001, "r2", "-"),
    ])

    legacy = compute_metrics(pred, ref, distance_cutoffs=[50, 500])
    extended = compute_metrics(pred, ref, distance_cutoffs=[50, 500],
                               point_strand=True, null_shuffle=True,
                               null_n_seeds=2)

    # Every legacy key is present and byte-identical.
    for key in legacy:
        assert extended[key] == legacy[key]
    assert set(legacy) == {"n_predicted", "n_reference", "cutoffs"}

    # New blocks live alongside, under their own names.
    assert "point_strand" in extended
    assert "null" in extended
    assert extended["point_strand"]["mode"] == "point_strand_matched"
    assert extended["null"]["n_seeds"] == 2


@needs_bedtools
def test_compute_metrics_restricted_reference_implies_point_block(tmp_path):
    pytest.importorskip("pybedtools")
    from ema.benchmark.metrics import compute_metrics

    pred = _write_bed(tmp_path / "p.bed", [_bed6("chr1", 100, 101, "p", "+")])
    ref = _write_bed(tmp_path / "r.bed", [_bed6("chr1", 100, 101, "r", "+")])
    restricted = _write_bed(tmp_path / "rr.bed",
                            [_bed6("chr1", 100, 101, "r", "+")])

    out = compute_metrics(pred, ref, distance_cutoffs=[10],
                          restricted_reference_bed=restricted)
    assert "point_strand" in out
    assert out["point_strand"]["cutoffs"][10]["recall_restricted"] == 1.0
